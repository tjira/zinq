#include "viewer.h"
#include "viewer_math.hpp"
#include "viewer_mesh.hpp"
#include "viewer_ptable.hpp"
#include "viewer_shader.hpp"
#include "viewer_state.hpp"
#include "viewer_trajectory.hpp"

#include <glad/gl.h>

#define GLFW_INCLUDE_NONE

#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <stb_image_write.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

static const std::string vertex_shader_src = R"(
#version 330 core

layout(location = 0) in vec3 i_position;
layout(location = 1) in vec3   i_normal;

uniform mat4 u_model;
uniform mat4  u_view;
uniform mat4  u_proj;

out vec3 v_position;
out vec3   v_normal;

void main() {
    vec4 world_pos = u_model * vec4(i_position, 1.0);

    v_position = world_pos.xyz;

    v_normal = mat3(u_model) * i_normal;

    gl_Position = u_proj * u_view * world_pos;
}
)";

static const std::string fragment_shader_src = R"(
#version 330 core

in vec3 v_position;
in vec3   v_normal;

uniform vec3      u_color;
uniform vec3     u_camera;
uniform vec3  u_light_dir;
uniform float   u_ambient;
uniform float   u_diffuse;
uniform float  u_specular;
uniform float u_shininess;

out vec4 o_color;

void main() {
    vec3 n = normalize(v_normal             );
    vec3 l = normalize(u_light_dir          );
    vec3 v = normalize(u_camera - v_position);
    vec3 h = normalize(l + v                );

    float diff = max(dot(n, l), 0.0);

    float spec = (diff > 0.0 && u_shininess > 0.0) ? pow(max(dot(n, h), 0.0), u_shininess) : 0.0;

    vec3 ambient  = u_ambient * u_color;

    vec3 diffuse  = u_diffuse * diff * u_color;

    vec3 specular = u_specular * spec * vec3(1.0);

    o_color = vec4(ambient + diffuse + specular, 1.0);
}
)";

/// Manages GLFW window and context lifecycle with automated resource reclamation.
struct WindowContext {
    GLFWwindow* window = nullptr;

    /// Releases GLFW window resources and shuts down the GLFW subsystem.
    ~WindowContext() {
        if (window) {
            glfwSetWindowUserPointer(window, nullptr);
            glfwDestroyWindow(window);
        }

        glfwTerminate();
    }
};

/// Handles mouse button press and release events for camera interaction.
static void mouseButtonCallback(GLFWwindow* window, int button, int action, int) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (!state) return;

    if (action == GLFW_RELEASE) {
        if (button == GLFW_MOUSE_BUTTON_LEFT)  state->mouse.dragging_left  = false;
        if (button == GLFW_MOUSE_BUTTON_RIGHT) state->mouse.dragging_right = false;

        return;
    }

    if (state->gui.isVisible() && state->gui.wantsCaptureMouse()) return;

    double x = 0, y = 0;

    glfwGetCursorPos(window, &x, &y);

    state->mouse.last_x = static_cast<float>(x);
    state->mouse.last_y = static_cast<float>(y);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        state->mouse.dragging_left = true;
    }

    else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        state->mouse.dragging_right = true;
    }
}

/// Handles mouse motion events to orbit or pan the camera.
static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (!state) return;

    float fx = static_cast<float>(xpos);
    float fy = static_cast<float>(ypos);
    float dx = fx - state->mouse.last_x;
    float dy = fy - state->mouse.last_y;

    state->mouse.last_x = fx;
    state->mouse.last_y = fy;

    if (state->gui.isVisible() && state->gui.wantsCaptureMouse() && !state->mouse.dragging_left && !state->mouse.dragging_right) {
        return;
    }

    if (state->mouse.dragging_left) {
        state->camera.orbit(dx, dy);
    }

    else if (state->mouse.dragging_right) {
        state->camera.pan(dx, dy);
    }
}

/// Handles mouse wheel scroll events for interactive camera zoom.
static void scrollCallback(GLFWwindow* window, double, double dy) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (state && (!state->gui.isVisible() || !state->gui.wantsCaptureMouse())) {
        state->camera.zoom(static_cast<float>(dy));
    }
}

/// Adjusts the OpenGL viewport and projection matrix on window resize.
static void resizeCallback(GLFWwindow* window, int width, int height) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (state && width > 0 && height > 0) {
        state->width  =  width;
        state->height = height;

        state->proj = perspective(45 * static_cast<float>(M_PI) / 180, static_cast<float>(width) / height, 0.05, 2000);

        glViewport(0, 0, width, height);
    }
}

/// Handles keyboard input events for trajectory playback and rendering mode.
static void keyCallback(GLFWwindow* window, int key, int, int action, int) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (!state || action != GLFW_PRESS) return;

    if (key == GLFW_KEY_G) {
        state->gui.toggleVisible(); return;
    }

    if (key == GLFW_KEY_O) {
        if (state->gui.isFileDialogOpen() && state->gui.wantsTextInput()) return;

        state->gui.toggleFileDialog(); return;
    }

    if (key == GLFW_KEY_S) {
        if ((state->gui.isFileDialogOpen() || state->gui.isScreenshotDialogOpen()) && state->gui.wantsTextInput()) return;

        state->gui.toggleScreenshotDialog(); return;
    }

    if (state->gui.isVisible() && state->gui.wantsCaptureKeyboard()) return;

    switch (key) {
        case GLFW_KEY_SPACE:
            state->trajectory.togglePause(); break;
        case GLFW_KEY_LEFT:
            state->trajectory.stepFrame(-1); state->trajectory.setPaused(true);break;
        case GLFW_KEY_RIGHT:
            state->trajectory.stepFrame(1); state->trajectory.setPaused(true); break;
        case GLFW_KEY_UP:
            state->trajectory.adjustSpeed(0.8); break;
        case GLFW_KEY_DOWN:
            state->trajectory.adjustSpeed(1.25); break;
        case GLFW_KEY_R:
            state->camera.reset(); break;
        case GLFW_KEY_B:
            state->dark_mode = !state->dark_mode; break;
        case GLFW_KEY_TAB:
            state->style = static_cast<ViewStyle>((static_cast<int>(state->style) + 1) % 3); break;
        case GLFW_KEY_ESCAPE:
            if (state->gui.isFileDialogOpen()) {
                state->gui.closeFileDialog(); break;
            }

            if (state->gui.isScreenshotDialogOpen()) {
                state->gui.closeScreenshotDialog(); break;
            }

            glfwSetWindowShouldClose(window, GLFW_TRUE); break;
        default:
            break;
    }
}

/// Captures active OpenGL framebuffer raster pixels and encodes uncompressed RGBA image data to a PNG file.
static bool saveScreenshot(const std::string& filepath, int width, int height) {
    if (width <= 0 || height <= 0) return false;

    std::vector<unsigned char> pixels(static_cast<size_t>(width) * static_cast<size_t>(height) * 4);

    std::vector<float> depths(static_cast<size_t>(width) * static_cast<size_t>(height));

    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

    glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depths.data());

    if (glGetError() == GL_NO_ERROR) {
        for (size_t i = 0; i < depths.size(); ++i) {
            if (depths[i] >= 1) {
                pixels[i * 4 + 3] = 0;
            }
        }
    }

    for (size_t i = 0; i < pixels.size(); i += 4) {
        if (pixels[i + 3] == 0) {
            pixels[i] = 0; pixels[i + 1] = 0; pixels[i + 2] = 0;
        }
    }

    stbi_flip_vertically_on_write(1);

    int stride = width * 4;

    if (!stbi_write_png(filepath.c_str(), width, height, 4, pixels.data(), stride)) {
        std::cerr << "\nFAILED TO WRITE SCREENSHOT TO '" << filepath << "' IMAGE FILE" << std::endl;

        return false;
    }

    return true;
}

/// Uploads camera view, projection, and illumination uniforms to the shader.
static void setUniforms(const Shader& shader, const ViewerState& state) {
    shader.use();

    shader.set("u_view", state.camera.getViewMatrix());
    shader.set("u_proj",                   state.proj);
    shader.set("u_camera",      state.camera.getEye());

    vec3 f = normalize(state.camera.target - state.camera.getEye());
    vec3 r = normalize(cross(f, vec3(0, 1, 0))                    );

    vec3 u = cross(r, f);

    vec3 light_dir = normalize(-f + state.light_offset.x * r + state.light_offset.y * u + state.light_offset.z * (-f));

    shader.set("u_light_dir",  light_dir               );
    shader.set("u_ambient",    state.ambient_intensity );
    shader.set("u_diffuse",    state.diffuse_intensity );
    shader.set("u_specular",   state.specular_intensity);
    shader.set("u_shininess",  state.shininess         );
}

extern "C" {
    /// Opens an interactive OpenGL viewport rendering the molecular system or trajectory at the given file path.
    bool viewer_open_window(const char *filename) {
        if (!filename) return false;

        Trajectory traj = Trajectory::Load(filename);

        if (traj.empty()) {
            std::cerr << "NO MOLECULAR COORDINATES FOUND IN '" << filename << "' STRUCTURE FILE" << std::endl;

            return false;
        }

        if (!glfwInit()) return false;

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        glfwWindowHint(GLFW_SAMPLES, 4); glfwWindowHint(GLFW_ALPHA_BITS, 8);

#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

        WindowContext ctx;

        ctx.window = glfwCreateWindow(1280, 720, "Zinq Viewer", nullptr, nullptr);

        if (!ctx.window) return false;

        glfwMakeContextCurrent(ctx.window);

        if (!gladLoadGL(reinterpret_cast<GLADloadfunc>(glfwGetProcAddress))) return false;

        glEnable(GL_DEPTH_TEST); glEnable(GL_CULL_FACE); glfwSwapInterval(1);

        std::cout << "\nCONTROLS:\n"
                  << "  Left Drag             ORBIT ROTATION\n"
                  << "  Right Drag            PAN TRANSLATION\n"
                  << "  Scroll Wheel          SMOOTH ZOOM\n"
                  << "  Space                 PLAY / PAUSE TRAJECTORY\n"
                  << "  Left Arrow            PREVIOUS FRAME\n"
                  << "  Right Arrow           NEXT FRAME\n"
                  << "  Up Arrow              INCREASE SPEED\n"
                  << "  Down Arrow            DECREASE SPEED\n"
                  << "  R                     RESET CAMERA VIEW\n"
                  << "  B                     TOGGLE BACKGROUND\n"
                  << "  G                     TOGGLE GUI OVERLAY\n"
                  << "  O                     OPEN FILE DIALOG\n"
                  << "  S                     SAVE SCREENSHOT\n"
                  << "  Tab                   CYCLE STYLE\n"
                  << "  Esc                   QUIT VIEWER"
        << std::endl;

        ViewerState state; state.filename = filename;

        state.trajectory = std::move(traj);

        float camera_dist = std::max(6.0f, state.trajectory.getMaxRadius() * 2.5f);

        state.camera.distance         = camera_dist;
        state.camera.initial_distance = camera_dist;

        state.proj = perspective(45 * static_cast<float>(M_PI) / 180, 1280.0 / 720, 0.05, 2000);

        glfwSetWindowUserPointer(ctx.window, &state);

        glfwSetMouseButtonCallback(ctx.window, mouseButtonCallback);

        glfwSetCursorPosCallback(ctx.window, cursorPosCallback);

        glfwSetScrollCallback(ctx.window, scrollCallback);

        glfwSetWindowSizeCallback(ctx.window, resizeCallback);

        glfwSetKeyCallback(ctx.window, keyCallback);

        if (!state.gui.init(ctx.window)) return false;

        Shader shader(vertex_shader_src, fragment_shader_src);

        if (!shader.isValid()) return false;

        Mesh sphere_mesh   = Mesh::Icosphere(state.sphere_subdivisions, state.smooth_shading);
        Mesh cylinder_mesh = Mesh::Cylinder(state.cylinder_sectors,     state.smooth_shading);

        while (!glfwWindowShouldClose(ctx.window)) {
            if (state.mesh_dirty) {
                sphere_mesh   = Mesh::Icosphere(state.sphere_subdivisions, state.smooth_shading);
                cylinder_mesh = Mesh::Cylinder(state.cylinder_sectors,     state.smooth_shading);

                state.mesh_dirty = false;
            }

            vec3 bg = state.dark_mode ? vec3(0, 0, 0) : vec3(1, 1, 1);

            glClearColor(bg.x, bg.y, bg.z, 0);

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (state.style == ViewStyle::Wireframe) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            }

            if (state.style == ViewStyle::BallAndStick || state.style == ViewStyle::SpaceFilling) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }

            setUniforms(shader, state);

            float base_scale = (state.style == ViewStyle::SpaceFilling) ? 2.4f : 1.0f;

            float atom_scale = base_scale * state.atom_scale;

            bool show_bonds = (state.style != ViewStyle::SpaceFilling);

            vec3 bond_color = state.dark_mode ? vec3(0.85, 0.85, 0.85) : vec3(0.35, 0.35, 0.35);

            state.trajectory.render(shader, sphere_mesh, cylinder_mesh, atom_scale, show_bonds, bond_color);

            if (!state.pending_screenshot_path.empty()) {
                int fb_w = 0, fb_h = 0;

                glfwGetFramebufferSize(ctx.window, &fb_w, &fb_h);

                saveScreenshot(state.pending_screenshot_path, fb_w, fb_h);

                state.pending_screenshot_path.clear();
            }

            if (state.gui.isVisible()) {
                state.gui.newFrame(); state.gui.render(state); state.gui.draw();
            }

            state.updateTitle(ctx.window);

            glfwSwapBuffers(ctx.window); glfwPollEvents();
        }

        return true;
    }
}
