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
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

static const std::string vertex_shader_src = R"(
#version 330 core

layout(location = 0) in vec3 i_position;
layout(location = 1) in vec3   i_normal;
layout(location = 2) in mat4    i_model;
layout(location = 6) in vec3    i_color;

uniform mat4  u_view;
uniform mat4  u_proj;

out vec3 v_position;
out vec3   v_normal;
out vec3    v_color;

void main() {
    vec4 world_pos = i_model * vec4(i_position, 1.0);

    v_position = world_pos.xyz;

    v_normal = mat3(i_model) * i_normal;

    v_color = i_color;

    gl_Position = u_proj * u_view * world_pos;
}
)";

static const std::string fragment_shader_src = R"(
#version 330 core

in vec3 v_position;
in vec3   v_normal;
in vec3    v_color;

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

    vec3 ambient  = u_ambient * v_color;

    vec3 diffuse  = u_diffuse * diff * v_color;

    vec3 specular = u_specular * spec * vec3(1.0);

    o_color = vec4(ambient + diffuse + specular, 1.0);
}
)";

static const std::string outline_fragment_shader_src = R"(
#version 330 core

uniform vec3 u_outline_color;

out vec4 o_color;

void main() {
    o_color = vec4(u_outline_color, 1.0);
}
)";

/// Manages GLFW window and context lifecycle with automated resource reclamation.
struct WindowContext {
    GLFWwindow* window = nullptr;

    /// Releases GLFW window resources and shuts down the GLFW subsystem.
    ~WindowContext() {
        if (window) {
            glfwSetWindowUserPointer(window, nullptr); glfwDestroyWindow(window);
        }

        glfwTerminate();
    }
};

/// Stores camera-projected view depth and instance transform data for individual atom selection outlines.
struct SelectedOutline {
    float depth = 0.0f; InstanceData mask; InstanceData outline;
};

/// Handles mouse button press and release events for camera interaction and atom selection.
static void mouseButtonCallback(GLFWwindow* window, int button, int action, int) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (!state) return;

    if (action == GLFW_RELEASE) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            state->mouse.dragging_left = false;

            double rx = 0, ry = 0;

            glfwGetCursorPos(window, &rx, &ry);

            float dx = static_cast<float>(rx) - state->mouse.press_x;
            float dy = static_cast<float>(ry) - state->mouse.press_y;

            if (dx * dx + dy * dy < 25.0f && (!state->gui.isVisible(state) || !state->gui.wantsCaptureMouse())) {
                int picked = state->pickAtom(static_cast<float>(rx), static_cast<float>(ry));

                if (picked >= 0) {
                    state->toggleAtomSelection(picked);
                }
            }
        }

        if (button == GLFW_MOUSE_BUTTON_RIGHT) state->mouse.dragging_right = false;

        return;
    }

    if (state->gui.isVisible(state) && state->gui.wantsCaptureMouse()) return;

    double x = 0, y = 0;

    glfwGetCursorPos(window, &x, &y);

    state->mouse.last_x  = static_cast<float>(x);
    state->mouse.last_y  = static_cast<float>(y);
    state->mouse.press_x = static_cast<float>(x);
    state->mouse.press_y = static_cast<float>(y);

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

    if (state->gui.isVisible(state) && state->gui.wantsCaptureMouse() && !state->mouse.dragging_left && !state->mouse.dragging_right) {
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

    if (state && (!state->gui.isVisible(state) || !state->gui.wantsCaptureMouse())) {
        state->camera.zoom(static_cast<float>(dy));
    }
}

/// Adjusts the OpenGL viewport on window resize.
static void resizeCallback(GLFWwindow* window, int width, int height) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (state && width > 0 && height > 0) {
        state->width  =  width;
        state->height = height;

        glViewport(0, 0, width, height);
    }
}

/// Handles keyboard input events for trajectory playback and rendering mode.
static void keyCallback(GLFWwindow* window, int key, int, int action, int) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (!state || (action != GLFW_PRESS && action != GLFW_REPEAT)) return;

    if (action == GLFW_REPEAT && key != GLFW_KEY_LEFT && key != GLFW_KEY_RIGHT) return;

    if (key == GLFW_KEY_F1) {
        state->gui.toggleVisible(); return;
    }

    if (key == GLFW_KEY_ESCAPE) {
        if (state->gui.isFileDialogOpen()) {
            state->gui.closeFileDialog(); return;
        }

        if (state->gui.isScreenshotDialogOpen()) {
            state->gui.closeScreenshotDialog(); return;
        }

        glfwSetWindowShouldClose(window, GLFW_TRUE); return;
    }

    if (state->gui.wantsTextInput()) return;

    if (key == GLFW_KEY_O) {
        state->gui.toggleFileDialog(); return;
    }

    if (key == GLFW_KEY_S) {
        state->gui.toggleScreenshotDialog(); return;
    }

    if (key == GLFW_KEY_F12) {
        std::time_t raw_time = std::time(nullptr);

        std::tm* time_info = std::localtime(&raw_time);

        char time_buf[64];

        std::strftime(time_buf, sizeof(time_buf), "screenshot_%Y%m%d_%H%M%S.png", time_info);

        state->pending_screenshot_path = time_buf; return;
    }

    if (state->gui.wantsCaptureKeyboard()) return;

    switch (key) {
        case GLFW_KEY_SPACE:
            state->trajectory.togglePause(); break;
        case GLFW_KEY_LEFT:
            state->trajectory.stepFrame(-1); state->trajectory.setPaused(true); break;
        case GLFW_KEY_RIGHT:
            state->trajectory.stepFrame(1); state->trajectory.setPaused(true); break;
        case GLFW_KEY_HOME:
            state->trajectory.setFrame(0); state->trajectory.setPaused(true); break;
        case GLFW_KEY_END:
            state->trajectory.setFrame(state->trajectory.frameCount() - 1); state->trajectory.setPaused(true); break;
        case GLFW_KEY_UP:
            state->trajectory.adjustSpeed(0.8f); break;
        case GLFW_KEY_DOWN:
            state->trajectory.adjustSpeed(1.25f); break;
        case GLFW_KEY_R:
            state->camera.reset(); break;
        case GLFW_KEY_P:
            state->camera.toggleProjection(); break;
        case GLFW_KEY_B:
            state->toggleBackground(); break;
        case GLFW_KEY_TAB:
            state->style = static_cast<ViewStyle>((static_cast<int>(state->style) + 1) % 3); break;
        case GLFW_KEY_C:
            state->clearSelection(); break;
        default:
            break;
    }
}

/// Captures active OpenGL framebuffer raster pixels and encodes uncompressed RGBA image data to a PNG file.
static bool saveScreenshot(const std::string& filepath, int width, int height, bool transparent) {
    if (width <= 0 || height <= 0) return false;

    std::vector<unsigned char> pixels(static_cast<size_t>(width) * static_cast<size_t>(height) * 4);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

    for (size_t i = 0; i < pixels.size(); i += 4) {
        if (!transparent) {
            pixels[i + 3] = 255;
        }

        else if (pixels[i + 3] == 0) {
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
    shader.set("u_proj",        state.getProjMatrix());
    shader.set("u_camera",      state.camera.getEye());

    vec3 f = normalize(state.camera.target - state.camera.getEye());
    vec3 r = normalize(cross(f, vec3(0.0f, 1.0f, 0.0f))           );

    vec3 u = cross(r, f);

    vec3 light_dir = normalize(-f + state.light_offset.x * r + state.light_offset.y * u + state.light_offset.z * (-f));

    shader.set("u_light_dir",  light_dir               );
    shader.set("u_ambient",    state.ambient_intensity );
    shader.set("u_diffuse",    state.diffuse_intensity );
    shader.set("u_specular",   state.specular_intensity);
    shader.set("u_shininess",  state.shininess         );
}

/// Renders highlighted outline halos around selected atoms using depth-tested and back-to-front ordered stencil passes.
static void renderSelectionOutlines(const Mesh& sphere_mesh, const Shader& outline_shader, const ViewerState& state, float atom_scale) {
    if (state.selected_atoms.empty()) return;

    const Geometry* geom = state.trajectory.currentGeometry();

    if (!geom) return;

    const auto& atoms = geom->getAtoms();

    vec3 eye = state.camera.getEye(), f = normalize(state.camera.target - eye);

    std::vector<SelectedOutline> items;

    items.reserve(state.selected_atoms.size());

    for (int idx : state.selected_atoms) {
        if (idx >= 0 && idx < static_cast<int>(atoms.size())) {
            const auto& a = atoms[idx];

            mat4 m = translate(mat4(1.0f), a.position);

            float r = a.radius * atom_scale, outline_r = r * 1.07f + 0.02f;

            float depth = dot(a.position - eye, f);

            items.push_back({ depth, { scale(m, vec3(r)), vec3(1.0f) }, { scale(m, vec3(outline_r)), state.outline_color } });
        }
    }

    if (items.empty()) return;

    std::sort(items.begin(), items.end(), [](const SelectedOutline& a, const SelectedOutline& b) {
        return a.depth > b.depth;
    });

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_STENCIL_TEST);

    outline_shader.use();

    outline_shader.set("u_view", state.camera.getViewMatrix());
    outline_shader.set("u_proj",        state.getProjMatrix());
    outline_shader.set("u_outline_color", state.outline_color);

    std::vector<InstanceData> single_mask(1), single_outline(1);

    for (const auto& item : items) {
        glStencilMask(0xFF);

        glClear(GL_STENCIL_BUFFER_BIT);

        glStencilFunc(GL_ALWAYS, 1, 0xFF); glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

        glDepthMask(GL_FALSE); glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LEQUAL);

        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

        single_mask[0] = item.mask; sphere_mesh.renderInstanced(single_mask);

        glStencilMask(0x00); glStencilFunc(GL_NOTEQUAL, 1, 0xFF); glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

        single_outline[0] = item.outline; sphere_mesh.renderInstanced(single_outline);
    }

    glEnable(GL_DEPTH_TEST);

    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);

    glStencilMask(0xFF); glDisable(GL_STENCIL_TEST);

    if (state.style == ViewStyle::Wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
}

extern "C" {
    /// Opens an interactive OpenGL viewport rendering the molecular system or trajectory at the given file path.
    bool viewer_open_window(const char *filename) {
        Trajectory traj; std::string initial_path;

        if (filename && filename[0] != '\0') {
            traj = Trajectory::Load(filename);

            if (traj.empty()) {
                std::cerr << "NO MOLECULAR COORDINATES FOUND IN '" << filename << "' STRUCTURE FILE" << std::endl;

                return false;
            }

            initial_path = filename;
        }

        if (!glfwInit()) return false;

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        glfwWindowHint(GLFW_SAMPLES, 4); glfwWindowHint(GLFW_ALPHA_BITS, 8); glfwWindowHint(GLFW_STENCIL_BITS, 8);

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
                  << "  Left Click            PICK ATOM / TOGGLE SELECTION\n"
                  << "  Scroll Wheel          SMOOTH ZOOM\n"
                  << "  Space                 PLAY / PAUSE TRAJECTORY\n"
                  << "  Left Arrow            PREVIOUS FRAME\n"
                  << "  Right Arrow           NEXT FRAME\n"
                  << "  Home / End            FIRST / LAST FRAME\n"
                  << "  Up Arrow              INCREASE SPEED\n"
                  << "  Down Arrow            DECREASE SPEED\n"
                  << "  R                     RESET CAMERA VIEW\n"
                  << "  P                     TOGGLE PROJECTION\n"
                  << "  B                     TOGGLE BACKGROUND\n"
                  << "  C                     CLEAR ATOM SELECTION\n"
                  << "  F1                    TOGGLE GUI OVERLAY\n"
                  << "  F12                   QUICK SCREENSHOT\n"
                  << "  O                     OPEN FILE DIALOG\n"
                  << "  S                     SAVE SCREENSHOT DIALOG\n"
                  << "  Tab                   CYCLE STYLE\n"
                  << "  Esc                   QUIT VIEWER"
        << std::endl;

        ViewerState state; state.filename = initial_path;

        glfwGetFramebufferSize(ctx.window, &state.width, &state.height);

        glViewport(0, 0, state.width, state.height);

        state.trajectory = std::move(traj);

        if (!state.trajectory.empty()) {
            float camera_dist = std::max(6.0f, state.trajectory.getMaxRadius() * 2.5f);

            state.camera.distance         = camera_dist;
            state.camera.initial_distance = camera_dist;

            state.camera.reset();
        }

        glfwSetWindowUserPointer(ctx.window, &state);

        glfwSetMouseButtonCallback(ctx.window, mouseButtonCallback);

        glfwSetCursorPosCallback(ctx.window, cursorPosCallback);

        glfwSetScrollCallback(ctx.window, scrollCallback);

        glfwSetFramebufferSizeCallback(ctx.window, resizeCallback);

        glfwSetKeyCallback(ctx.window, keyCallback);

        if (!state.gui.init(ctx.window)) return false;

        if (state.trajectory.empty()) {
            state.gui.openFileDialog();
        }

        Shader shader(vertex_shader_src, fragment_shader_src);

        if (!shader.isValid()) return false;

        Shader outline_shader(vertex_shader_src, outline_fragment_shader_src);

        if (!outline_shader.isValid()) return false;

        Mesh sphere_mesh   = Mesh::Icosphere(state.sphere_subdivisions, state.smooth_shading);
        Mesh cylinder_mesh = Mesh::Cylinder(state.cylinder_sectors,     state.smooth_shading);

        while (!glfwWindowShouldClose(ctx.window)) {
            if (state.mesh_dirty) {
                sphere_mesh   = Mesh::Icosphere(state.sphere_subdivisions, state.smooth_shading);
                cylinder_mesh = Mesh::Cylinder(state.cylinder_sectors,     state.smooth_shading);

                state.mesh_dirty = false;
            }

            vec3 bg = state.bg_color;

            float alpha = state.transparent_bg ? 0.0f : 1.0f;

            glClearColor(bg.x, bg.y, bg.z, alpha);

            glStencilMask(0xFF);

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

            glDisable(GL_STENCIL_TEST);

            if (state.style == ViewStyle::Wireframe) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            }

            if (state.style == ViewStyle::BallAndStick || state.style == ViewStyle::SpaceFilling) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }

            setUniforms(shader, state);

            float atom_scale = state.getEffectiveAtomScale();

            bool show_bonds = (state.style != ViewStyle::SpaceFilling);

            state.trajectory.render(sphere_mesh, cylinder_mesh, atom_scale, show_bonds, state.bond_color);

            renderSelectionOutlines(sphere_mesh, outline_shader, state, atom_scale);

            if (!state.pending_screenshot_path.empty()) {
                int fb_w = 0, fb_h = 0;

                glfwGetFramebufferSize(ctx.window, &fb_w, &fb_h);

                saveScreenshot(state.pending_screenshot_path, fb_w, fb_h, state.transparent_bg);

                state.pending_screenshot_path.clear();
            }

            if (state.gui.isVisible(&state)) {
                state.gui.newFrame(); state.gui.render(state); state.gui.draw();
            }

            state.updateTitle(ctx.window);

            glfwSwapBuffers(ctx.window); glfwPollEvents();
        }

        return true;
    }
}
