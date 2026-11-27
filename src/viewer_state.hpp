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

/// Specifies the geometric projection mapping for 3D viewport rendering.
enum class ProjectionMode {
    Perspective  = 0,
    Orthographic = 1
};

/// Computes spherical coordinate camera orbits, target pans, and projection matrix transformations.
struct OrbitCamera {
    vec3 target = vec3(0.0f);

    float distance = 10.0f;
    float yaw      =  0.0f;
    float pitch    =  0.0f;

    vec3 initial_target = vec3(0.0f);

    float initial_distance = 10.0f;
    float initial_yaw      =  0.0f;
    float initial_pitch    =  0.0f;

    ProjectionMode projection = ProjectionMode::Perspective;

    float fov_deg =   45.00f;
    float z_near  =    0.05f;
    float z_far   = 2000.00f;

    /// Evaluates the perspective or orthographic projection matrix for the given viewport aspect ratio.
    mat4 getProjMatrix(float aspect) const {
        if (projection == ProjectionMode::Orthographic) {
            float half_h = distance * std::tan(fov_deg * static_cast<float>(M_PI) / 360.0f);

            float half_w = half_h * aspect;

            return ortho(-half_w, half_w, -half_h, half_h, z_near, z_far);
        }

        return perspective(fov_deg * static_cast<float>(M_PI) / 180.0f, aspect, z_near, z_far);
    }

    /// Restores the camera distance, target centroid, and spherical orbital angles to initial defaults.
    void reset() {
        target   =   initial_target;
        distance = initial_distance;
        yaw      =      initial_yaw;
        pitch    =    initial_pitch;
    }

    /// Toggles between perspective and orthographic camera projection modes.
    void toggleProjection() {
        projection = (projection == ProjectionMode::Perspective) ? ProjectionMode::Orthographic : ProjectionMode::Perspective;
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
        return lookAt(getEye(), target, vec3(0.0f, 1.0f, 0.0f));
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

        vec3 f = normalize(target - eye                    );
        vec3 r = normalize(cross(f, vec3(0.0f, 1.0f, 0.0f)));

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
    float last_x = 0.0f;
    float last_y = 0.0f;

    float press_x = 0.0f;
    float press_y = 0.0f;

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

    ViewStyle style = ViewStyle::BallAndStick;

    vec3 bg_color = vec3(0.0f);

    bool transparent_bg = true;

    vec3 bond_color = vec3(0.65f, 0.65f, 0.65f);

    vec3 outline_color = vec3(1.0f);

    float ambient_intensity  = 0.35f;
    float diffuse_intensity  = 0.55f;
    float specular_intensity = 0.25f;
    float shininess          = 32.0f;

    vec3 light_offset = vec3(0.4f, 0.7f, 0.0f);

    int sphere_subdivisions = 2;
    int cylinder_sectors    = 16;

    bool smooth_shading = true;
    bool mesh_dirty     = false;

    float binding_factor = 0.013f;
    float bond_radius    = 0.090f;
    float atom_scale     = 1.000f;

    Trajectory trajectory;

    ViewerGui gui;

    int   last_reported_frame =    -1;
    float last_reported_wait  =  0.0f;

    bool last_reported_pause = false;

    std::vector<int> selected_atoms;

    /// Evaluates the perspective projection matrix using camera settings and current viewport dimensions.
    mat4 getProjMatrix() const {
        float aspect = height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;

        return camera.getProjMatrix(aspect);
    }

    /// Loads a new molecular trajectory from disk, recenters the viewport, and recalculates orbital bounding radius.
    bool loadTrajectory(const std::string& new_path) {
        Trajectory new_traj = Trajectory::Load(new_path, binding_factor, bond_radius);

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

        clearSelection();

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

        std::string title = filename.empty() ? "Zinq Viewer" : ("Zinq Viewer - " + filename);

        if (total > 1) {
            int fps = static_cast<int>(std::round(1000.0f / wait));

            title += " [Frame " + std::to_string(cur_frame + 1) + "/" + std::to_string(total) + "] ";

            title += (paused ? "[Paused]" : "[Playing]");

            title += " [" + std::to_string(fps) + " fps]";
        }

        glfwSetWindowTitle(window, title.c_str());
    }

    /// Toggles the viewport background color between white and black presets and adjusts outline color.
    void toggleBackground() {
        bool was_dark = (bg_color.x + bg_color.y + bg_color.z) < 1.5f;

        bg_color = was_dark ? vec3(1.0f) : vec3(0.0f);

        outline_color = was_dark ? vec3(0.1f) : vec3(1.0f);
    }

    /// Projects a 3D Cartesian position to 2D screen viewport coordinates.
    bool projectToScreen(const vec3& pos, float& out_x, float& out_y) const {
        float aspect = height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;

        mat4 vp = camera.getProjMatrix(aspect) * camera.getViewMatrix();

        vec4 clip = vp * vec4(pos, 1.0f);

        if (clip.w <= 0.001f) return false;

        float ndc_x = clip.x / clip.w;
        float ndc_y = clip.y / clip.w;

        if (ndc_x < -1.5f || ndc_x > 1.5f || ndc_y < -1.5f || ndc_y > 1.5f) return false;

        out_x = (ndc_x * 0.5f + 0.5f) * static_cast<float>(width);

        out_y = (1.0f - (ndc_y * 0.5f + 0.5f)) * static_cast<float>(height);

        return true;
    }

    /// Casts an unprojected ray from viewport cursor coordinates into the 3D scene.
    Ray castRay(float mouse_x, float mouse_y) const {
        float aspect = height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;

        mat4 view_proj = camera.getProjMatrix(aspect) * camera.getViewMatrix();

        mat4 inv_vp = inverse(view_proj);

        float ndc_x = (2.0f * mouse_x) / static_cast<float>(width) - 1.0f;

        float ndc_y = 1.0f - (2.0f * mouse_y) / static_cast<float>(height);

        vec4 p_near = inv_vp * vec4(ndc_x, ndc_y, -1.0f, 1.0f);
        vec4 p_far  = inv_vp * vec4(ndc_x, ndc_y,  1.0f, 1.0f);

        vec3 world_near = p_near.w != 0.0f ? p_near.xyz() / p_near.w : p_near.xyz();
        vec3 world_far  = p_far.w  != 0.0f ? p_far.xyz()  / p_far.w  :  p_far.xyz();

        return Ray{ world_near, normalize(world_far - world_near) };
    }

    /// Finds the nearest atom intersecting an unprojected screen ray, or -1 if none intersect.
    int pickAtom(float mouse_x, float mouse_y) const {
        const Geometry* geom = trajectory.currentGeometry();

        if (!geom) return -1;

        Ray ray = castRay(mouse_x, mouse_y);

        const auto& atoms = geom->getAtoms();

        int best_idx = -1; float min_dist = 1e30f;

        for (size_t i = 0; i < atoms.size(); ++i) {
            const auto& atom = atoms[i];

            if (atom.is_dummy) continue;

            float t = 0.0f; float r = atom.radius * atom_scale;

            if (ray.intersectSphere(atom.position, r, t) && t < min_dist) {
                min_dist = t; best_idx = static_cast<int>(i);
            }
        }

        return best_idx;
    }

    /// Toggles the selection status of the atom at the specified index, retaining up to 4 atoms.
    void toggleAtomSelection(int idx) {
        if (idx < 0) return;

        auto it = std::find(selected_atoms.begin(), selected_atoms.end(), idx);

        if (it != selected_atoms.end()) {
            selected_atoms.erase(it);
        }

        else {
            if (selected_atoms.size() >= 4) {
                selected_atoms.erase(selected_atoms.begin());
            }

            selected_atoms.push_back(idx);
        }
    }

    /// Clears the collection of currently selected atom indices.
    void clearSelection() {
        selected_atoms.clear();
    }
};
