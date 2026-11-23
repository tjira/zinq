#include "viewer.h"

#include <glad/gl.h>

#define GLFW_INCLUDE_NONE

#include <GLFW/glfw3.h>

extern "C" {
    bool viewer_open_window(void) {
        if (!glfwInit()) return false;

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

        GLFWwindow* window = glfwCreateWindow(800, 600, "Zinq Viewer", nullptr, nullptr);

        if (!window) {
            glfwTerminate(); return false;
        }

        glfwMakeContextCurrent(window);

        if (!gladLoadGL(reinterpret_cast<GLADloadfunc>(glfwGetProcAddress))) {
            glfwDestroyWindow(window); glfwTerminate(); return false;
        }

        while (!glfwWindowShouldClose(window)) {
            glClearColor(0, 0, 0, 1);

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glfwSwapBuffers(window); glfwPollEvents();
        }

        glfwDestroyWindow(window); glfwTerminate(); return true;
    }
}
