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
#include <iostream>
#include <string>
#include <unordered_map>

static const std::string vertex_shader_src = R"(
#version 330 core

layout(location = 0) in vec3 i_position;

layout(location = 1) in vec3 i_normal;

layout(location = 2) in vec3 i_color;

uniform mat4 u_model, u_view, u_proj;

out vec3 fragment, normal, color; out mat3 transform;

void main() {
    normal = normalize(mat3(transpose(inverse(u_model))) * i_normal);

    fragment = vec3(u_model * vec4(i_position, 1.0));

    color = i_color;

    gl_Position = u_proj * u_view * vec4(fragment, 1.0);

    transform = inverse(mat3(u_view));
}
)";

static const std::string fragment_shader_src = R"(
#version 330 core

struct Light { vec3 position; float ambient, diffuse, specular, shininess; };

uniform Light u_light;
uniform vec3 u_camera;

in vec3 fragment, normal, color; in mat3 transform;

out vec4 o_color;

void main() {
    vec3 lightPos = transform * u_light.position;

    vec3 reflection = reflect(-normalize(lightPos), normal);

    vec3 direction = normalize(u_camera - fragment);

    vec3 specular = vec3(pow(max(dot(direction, reflection), 0), u_light.shininess));

    vec3 diffuse = vec3(max(dot(normal, normalize(lightPos)), 0));

    vec3 phong = (vec3(u_light.ambient) + u_light.diffuse * diffuse + u_light.specular * specular);

    o_color = vec4(phong * color, 1.0);
}
)";

struct ViewerState {
    int width = 1280;
    int height = 720;

    vec2 mouse;

    mat4 view = mat4(1);
    mat4 proj = mat4(1);

    struct Light {
        float ambient = 0.45, diffuse = 0.55, specular = 0.35, shininess = 16; vec3 position = { 1, 1, 1 };
    } light;

    Trajectory trajectory;
};

static void positionCallback(GLFWwindow* window, double x, double y) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (state && glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        mat3 inv_view = inverse(to_mat3(state->view));

        vec3 xaxis = inv_view * vec3(0, 1, 0);
        vec3 yaxis = inv_view * vec3(1, 0, 0);

        state->view = rotate(state->view, 0.01 * (static_cast<float>(y) - state->mouse.y), yaxis);
        state->view = rotate(state->view, 0.01 * (static_cast<float>(x) - state->mouse.x), xaxis);
    }

    if (state) state->mouse = { static_cast<float>(x), static_cast<float>(y) };
}

static void scrollCallback(GLFWwindow* window, double, double dy) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (state) {
        float factor = 1 + 0.08 * static_cast<float>(dy);

        state->view = state->view * from_mat3(mat3(factor));
    }
}

static void resizeCallback(GLFWwindow* window, int width, int height) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (state && width > 0 && height > 0) {
        state->width = width; state->height = height;

        state->proj = perspective(45 * static_cast<float>(M_PI) / 180, static_cast<float>(width) / height, 0.01, 1000);

        glViewport(0, 0, width, height);
    }
}

static void keyCallback(GLFWwindow* window, int key, int, int action, int) {
    auto* state = static_cast<ViewerState*>(glfwGetWindowUserPointer(window));

    if (action == GLFW_PRESS && state) {
        if (key == GLFW_KEY_SPACE) {
            state->trajectory.togglePause();
        }

        else if (key == GLFW_KEY_Q || key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
    }
}

static void setUniforms(const Shader& shader, const ViewerState& state) {
    shader.use();

    mat3 inv_view = inverse(to_mat3(state.view));

    vec3 cam_pos = -(inv_view * state.view.col3(3));

    shader.set("u_view", state.view);
    shader.set("u_proj", state.proj);

    shader.set("u_camera", cam_pos);

    shader.set("u_light.position", state.light.position);

    shader.set("u_light.ambient", state.light.ambient);
    shader.set("u_light.diffuse", state.light.diffuse);

    shader.set("u_light.specular", state.light.specular);

    shader.set("u_light.shininess", state.light.shininess);
}

extern "C" {
    bool viewer_open_window(const char *filename) {
        if (!filename) return false;

        Trajectory traj = Trajectory::Load(filename);

        if (traj.empty()) {
            std::cerr << "Error: No molecular coordinates found in " << filename << std::endl;

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

        GLFWwindow* window = glfwCreateWindow(1280, 720, "Zinq Viewer", nullptr, nullptr);

        if (!window) {
            glfwTerminate(); return false;
        }

        glfwMakeContextCurrent(window);

        if (!gladLoadGL(reinterpret_cast<GLADloadfunc>(glfwGetProcAddress))) {
            glfwDestroyWindow(window); glfwTerminate(); return false;
        }

        glEnable(GL_DEPTH_TEST); glEnable(GL_CULL_FACE);

        glfwSwapInterval(1);

        ViewerState state;

        state.trajectory = std::move(traj);

        float camera_dist = std::max(6.0, state.trajectory.getMaxRadius() * 2.5);

        state.proj = perspective(45 * static_cast<float>(M_PI) / 180, 1280.0 / 720, 0.01, 1000);

        state.view = lookAt(vec3(0, 0, camera_dist), vec3(0, 0, 0), vec3(0, 1, 0));

        glfwSetWindowUserPointer(window, &state);

        glfwSetCursorPosCallback(window, positionCallback);

        glfwSetScrollCallback(window, scrollCallback);

        glfwSetWindowSizeCallback(window, resizeCallback);

        glfwSetKeyCallback(window, keyCallback);

        std::unordered_map<std::string, Mesh> meshes;

        for (const auto& pair : ptable) {
            Mesh sphere = Mesh::Icosphere(2, true);

            sphere.setColor(pair.second.color);

            meshes[pair.first] = std::move(sphere);
        }

        Mesh bond = Mesh::Cylinder(16, true);

        bond.setColor(vec3(0.85, 0.85, 0.85));

        meshes["bond"] = std::move(bond);

        Shader shader(vertex_shader_src, fragment_shader_src);

        while (!glfwWindowShouldClose(window)) {
            glClearColor(0, 0, 0, 1);

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            setUniforms(shader, state);

            state.trajectory.render(shader, meshes);

            glfwSwapBuffers(window); glfwPollEvents();
        }

        meshes.clear();

        glfwDestroyWindow(window); glfwTerminate();

        return true;
    }
}
