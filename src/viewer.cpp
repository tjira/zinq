#include "viewer.h"
#include "viewer_math.hpp"
#include "viewer_mesh.hpp"
#include "viewer_ptable.hpp"
#include "viewer_shader.hpp"
#include "viewer_trajectory.hpp"

#include <glad/gl.h>

#define GLFW_INCLUDE_NONE

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

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

uniform vec3     u_color;
uniform vec3    u_camera;
uniform vec3 u_light_dir;

out vec4 o_color;

void main() {
    vec3 n = normalize(v_normal             );
    vec3 l = normalize(u_light_dir          );
    vec3 v = normalize(u_camera - v_position);
    vec3 h = normalize(l + v                );

    float diff = max(dot(n, l), 0.0);

    float spec = diff > 0.0 ? pow(max(dot(n, h), 0.0), 32.0) : 0.0;

    vec3 ambient  = 0.35 * u_color;

    vec3 diffuse  = 0.55 * diff * u_color;
    vec3 specular = 0.25 * spec * vec3(1);

    o_color = vec4(ambient + diffuse + specular, 1.0);
}
)";

enum class ViewStyle {
    BallAndStick = 0,
    SpaceFilling = 1,
    Wireframe    = 2
};

struct OrbitCamera {
    vec3 target = vec3(0);

    float distance = 10;
    float yaw      =  0;
    float pitch    =  0;

    vec3 initial_target = vec3(0);

    float initial_distance = 10;
    float initial_yaw      =  0;
    float initial_pitch    =  0;

    /// Resets camera position, orientation, and zoom distance to initial values.
    void reset() {
        target   =   initial_target;
        distance = initial_distance;
        yaw      =      initial_yaw;
        pitch    =    initial_pitch;
    }

    /// Computes the 3D Cartesian eye coordinates from spherical angles and target.
    vec3 getEye() const {
        float cp = std::cos(pitch);
        float sp = std::sin(pitch);

        float cy = std::cos(yaw);
        float sy = std::sin(yaw);

        return target + distance * vec3(cp * sy, sp, cp * cy);
    }

    /// Constructs the view matrix for camera transform looking toward target.
    mat4 getViewMatrix() const {
        return lookAt(getEye(), target, vec3(0, 1, 0));
    }

    /// Rotates the camera orbit angles around the target center.
    void orbit(float dx, float dy) {
        yaw   -= dx * 0.008;
        pitch += dy * 0.008;

        const float max_pitch = 1.55;

        pitch = std::clamp(pitch, -max_pitch, max_pitch);
    }

    /// Shifts the camera target along camera screen-aligned plane vectors.
    void pan(float dx, float dy) {
        vec3 eye = getEye();

        vec3 f = normalize(target - eye           );
        vec3 r = normalize(cross(f, vec3(0, 1, 0)));

        vec3 u = cross(r, f);

        float factor = distance * 0.0015;

        target += (-dx * r + dy * u) * factor;
    }

    /// Adjusts the camera view distance from the orbital center.
    void zoom(float dy) {
        distance *= std::pow(0.9, dy);

        if (distance <    0.2) distance =    0.2;
        if (distance > 2000.0) distance = 2000.0;
    }
};

struct MouseState {
    float last_x = 0;
    float last_y = 0;

    bool dragging_left  = false;
    bool dragging_right = false;
};

struct ViewerState {
    int width = 1280;
    int height = 720;

    std::string filename;

    OrbitCamera camera;
    MouseState   mouse;

    mat4 proj = mat4(1);

    ViewStyle style = ViewStyle::BallAndStick;

    bool dark_mode = true;

    Trajectory trajectory;

    int   last_reported_frame = -1;
    float last_reported_wait  =  0;

    bool last_reported_pause = false;

    /// Updates the window title with filename, frame progress, pause, and FPS status.
    void updateTitle(GLFWwindow* window) {
        int cur_frame = trajectory.currentFrame();

        int total = trajectory.frameCount();
        bool paused = trajectory.isPaused();
        float wait = trajectory.getWaitMs();

        if (cur_frame == last_reported_frame && paused == last_reported_pause && wait == last_reported_wait) {
            return;
        }

        last_reported_frame = cur_frame;
        last_reported_pause =    paused;
        last_reported_wait  =      wait;

        std::string title = "Zinq Viewer - " + filename;

        if (total > 1) {
            int fps = static_cast<int>(std::round(1000.0f / wait));

            title += " [Frame " + std::to_string(cur_frame + 1) + "/" + std::to_string(total) + "] ";

            title += (paused ? "[Paused]" : "[Playing]");

            title += " [" + std::to_string(fps) + " fps]";
        }

        glfwSetWindowTitle(window, title.c_str());
    }
};

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

    double x = 0, y = 0;

    glfwGetCursorPos(window, &x, &y);

    state->mouse.last_x = static_cast<float>(x);
    state->mouse.last_y = static_cast<float>(y);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        state->mouse.dragging_left = (action == GLFW_PRESS);
    }

    else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        state->mouse.dragging_right = (action == GLFW_PRESS);
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

    if (state) {
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
            glfwSetWindowShouldClose(window, GLFW_TRUE); break;
        default:
            break;
    }
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

    vec3 light_dir = normalize(-f + 0.4 * r + 0.7 * u);

    shader.set("u_light_dir", light_dir);
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

        glfwWindowHint(GLFW_SAMPLES, 4);

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

        Shader shader(vertex_shader_src, fragment_shader_src);

        if (!shader.isValid()) return false;

        Mesh sphere_mesh   = Mesh::Icosphere(2, true);
        Mesh cylinder_mesh = Mesh::Cylinder(16, true);

        while (!glfwWindowShouldClose(ctx.window)) {
            vec3 bg = state.dark_mode ? vec3(0, 0, 0) : vec3(1, 1, 1);

            glClearColor(bg.x, bg.y, bg.z, 1);

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (state.style == ViewStyle::Wireframe) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            }

            if (state.style == ViewStyle::BallAndStick || state.style == ViewStyle::SpaceFilling) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }

            setUniforms(shader, state);

            float atom_scale = (state.style == ViewStyle::SpaceFilling) ? 2.4 : 1;

            bool show_bonds = (state.style != ViewStyle::SpaceFilling);

            vec3 bond_color = state.dark_mode ? vec3(0.85, 0.85, 0.85) : vec3(0.35, 0.35, 0.35);

            state.trajectory.render(shader, sphere_mesh, cylinder_mesh, atom_scale, show_bonds, bond_color);

            state.updateTitle(ctx.window);

            glfwSwapBuffers(ctx.window); glfwPollEvents();
        }

        return true;
    }
}
