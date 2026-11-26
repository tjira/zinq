#pragma once

#include "viewer_gui.hpp"
#include "viewer_math.hpp"
#include "viewer_trajectory.hpp"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

/// Specifies the geometric rendering topology for atomic spheres and molecular bond cylinders.
enum class ViewStyle {
    BallAndStick = 0,
    SpaceFilling = 1,
    Wireframe    = 2
};

/// Computes spherical coordinate camera orbits, target pans, and projection matrix transformations.
struct OrbitCamera {
    vec3 target = vec3(0);

    float distance = 10;
    float yaw      =  0;
    float pitch    =  0;

    vec3 initial_target = vec3(0);

    float initial_distance = 10;
    float initial_yaw      =  0;
    float initial_pitch    =  0;

    /// Restores the camera distance, target centroid, and spherical orbital angles to initial defaults.
    void reset() {
        target   =   initial_target;
        distance = initial_distance;
        yaw      =      initial_yaw;
        pitch    =    initial_pitch;
    }

    /// Evaluates the Cartesian 3D eye vector from spherical orbital angles and spatial target position.
    vec3 getEye() const {
        float cp = std::cos(pitch);
        float sp = std::sin(pitch);

        float cy = std::cos(yaw);
        float sy = std::sin(yaw);

        return target + distance * vec3(cp * sy, sp, cp * cy);
    }

    /// Evaluates the rigid-body view transformation matrix gazing at the target center.
    mat4 getViewMatrix() const {
        return lookAt(getEye(), target, vec3(0, 1, 0));
    }

    /// Increments spherical polar and azimuthal camera angles with clamped elevation boundaries.
    void orbit(float dx, float dy) {
        yaw   -= dx * 0.008f;
        pitch += dy * 0.008f;

        const float max_pitch = 1.55f;

        pitch = std::clamp(pitch, -max_pitch, max_pitch);
    }

    /// Displaces the focal target along screen-aligned camera plane basis vectors.
    void pan(float dx, float dy) {
        vec3 eye = getEye();

        vec3 f = normalize(target - eye           );
        vec3 r = normalize(cross(f, vec3(0, 1, 0)));

        vec3 u = cross(r, f);

        float factor = distance * 0.0015f;

        target += (-dx * r + dy * u) * factor;
    }

    /// Multiplies the radial orbit distance by an exponential zoom attenuation factor.
    void zoom(float dy) {
        distance *= std::pow(0.9f, dy);

        if (distance <    0.2f) distance =    0.2f;
        if (distance > 2000.0f) distance = 2000.0f;
    }
};

/// Records cursor Cartesian coordinates and active drag button states for viewport interaction.
struct MouseState {
    float last_x = 0;
    float last_y = 0;

    bool dragging_left  = false;
    bool dragging_right = false;
};

/// Maintains dynamic state for camera projection, molecular trajectory playback, and GUI overlay.
struct ViewerState {
    int width = 1280;
    int height = 720;

    std::string filename; std::string pending_screenshot_path;

    OrbitCamera camera;
    MouseState   mouse;

    mat4 proj = mat4(1);

    ViewStyle style = ViewStyle::BallAndStick;

    bool dark_mode = true;

    Trajectory trajectory;

    ViewerGui gui;

    int   last_reported_frame = -1;
    float last_reported_wait  =  0;

    bool last_reported_pause = false;

    /// Loads a new molecular trajectory from disk, recenters the viewport, and recalculates orbital bounding radius.
    bool loadTrajectory(const std::string& new_path) {
        Trajectory new_traj = Trajectory::Load(new_path);

        if (new_traj.empty()) {
            std::cerr << "NO MOLECULAR COORDINATES FOUND IN '" << new_path << "' STRUCTURE FILE" << std::endl;

            return false;
        }

        trajectory = std::move(new_traj);

        filename = new_path;

        float camera_dist = std::max(6.0f, trajectory.getMaxRadius() * 2.5f);

        camera.distance         = camera_dist;
        camera.initial_distance = camera_dist;

        camera.reset();

        last_reported_frame = -1;

        return true;
    }

    /// Formats GLFW window title string with active frame index, playback cadence, and timing metrics.
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
