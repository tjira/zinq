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

    quat orientation = quat(1.0f, 0.0f, 0.0f, 0.0f);

    vec3 initial_target = vec3(0.0f);

    float initial_distance = 10.0f;

    quat initial_orientation = quat(1.0f, 0.0f, 0.0f, 0.0f);

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

    /// Restores camera distance, target centroid, and spatial orientation to initial defaults.
    void reset() {
        target      =      initial_target;
        distance    =    initial_distance;
        orientation = initial_orientation;
    }

    /// Toggles between perspective and orthographic camera projection modes.
    void toggleProjection() {
        projection = (projection == ProjectionMode::Perspective) ? ProjectionMode::Orthographic : ProjectionMode::Perspective;
    }

    /// Evaluates the normalized backward direction vector pointing from target to camera eye.
    vec3 getBack() const {
        return orientation.rotate(vec3(0.0f, 0.0f, 1.0f));
    }

    /// Evaluates the normalized forward line-of-sight vector pointing from camera eye to target.
    vec3 getForward() const {
        return -getBack();
    }

    /// Evaluates the screen-aligned vertical unit vector spanning the camera viewing plane.
    vec3 getUp() const {
        return orientation.rotate(vec3(0.0f, 1.0f, 0.0f));
    }

    /// Evaluates the screen-aligned horizontal unit vector spanning the camera viewing plane.
    vec3 getRight() const {
        return orientation.rotate(vec3(1.0f, 0.0f, 0.0f));
    }

    /// Evaluates the Cartesian 3D eye vector from camera target centroid, radial distance, and orientation.
    vec3 getEye() const {
        return target + getBack() * distance;
    }

    /// Evaluates the rigid-body view transformation matrix gazing at the target center.
    mat4 getViewMatrix() const {
        return lookAt(getEye(), target, getUp());
    }

    /// Rotates the camera orientation around screen-aligned vertical and horizontal viewport axes.
    void orbit(float dx, float dy) {
        quat rot_y = quat::fromAxisAngle(vec3(0.0f, 1.0f, 0.0f), -dx * 0.008f);
        quat rot_x = quat::fromAxisAngle(vec3(1.0f, 0.0f, 0.0f), -dy * 0.008f);

        orientation = normalize(orientation * (rot_x * rot_y));
    }

    /// Rotates the camera orientation around the line-of-sight view axis.
    void roll(float d) {
        quat rot_z = quat::fromAxisAngle(vec3(0.0f, 0.0f, 1.0f), -d * 0.008f);

        orientation = normalize(orientation * rot_z);
    }

    /// Displaces the focal target along screen-aligned camera plane basis vectors.
    void pan(float dx, float dy) {
        vec3 r = getRight();
        vec3 u =    getUp();

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

    bool dragging_left   = false;
    bool dragging_right  = false;
    bool dragging_middle = false;
    bool dragging_atoms  = false;
    bool drag_moved      = false;

    vec3 drag_plane_point  = vec3(0.0f);
    vec3 drag_last_hit     = vec3(0.0f);

    vec3 drag_plane_normal = vec3(0.0f, 0.0f, -1.0f);
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

    /// Evaluates the combined projection and view transformation matrix for the active viewport.
    mat4 getViewProjMatrix() const {
        return getProjMatrix() * camera.getViewMatrix();
    }

    /// Evaluates the effective atom rendering scale based on visual style and user scale factor.
    float getEffectiveAtomScale() const {
        float base = (style == ViewStyle::SpaceFilling) ? 2.4f : 1.0f;

        return base * atom_scale;
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
        vec4 clip = getViewProjMatrix() * vec4(pos, 1.0f);

        if (clip.w <= 0.001f) return false;

        float inv_w = 1.0f / clip.w;
        float ndc_x = clip.x * inv_w;
        float ndc_y = clip.y * inv_w;

        if (ndc_x < -1.5f || ndc_x > 1.5f || ndc_y < -1.5f || ndc_y > 1.5f) return false;

        out_x = (ndc_x * 0.5f + 0.5f) * static_cast<float>(width);

        out_y = (1.0f - (ndc_y * 0.5f + 0.5f)) * static_cast<float>(height);

        return true;
    }

    /// Casts an unprojected ray from viewport cursor coordinates into the 3D scene.
    Ray castRay(float mouse_x, float mouse_y) const {
        mat4 inv_vp = inverse(getViewProjMatrix());

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

        int best_idx = -1; float min_dist = 1e30f; float effective_scale = getEffectiveAtomScale();

        for (size_t i = 0; i < atoms.size(); ++i) {
            const auto& atom = atoms[i];

            if (atom.is_dummy) continue;

            float t = 0.0f; float r = atom.radius * effective_scale;

            if (ray.intersectSphere(atom.position, r, t) && t < min_dist) {
                min_dist = t; best_idx = static_cast<int>(i);
            }
        }

        return best_idx;
    }

    /// Checks whether the atom at the specified index is currently selected.
    bool isAtomSelected(int idx) const {
        return idx >= 0 && std::find(selected_atoms.begin(), selected_atoms.end(), idx) != selected_atoms.end();
    }

    /// Toggles the selection status of the atom at the specified index.
    void toggleAtomSelection(int idx) {
        if (idx < 0) return;

        auto it = std::find(selected_atoms.begin(), selected_atoms.end(), idx);

        if (it != selected_atoms.end()) {
            selected_atoms.erase(it);
        }

        else {
            selected_atoms.push_back(idx);
        }
    }

    /// Translates all selected atoms in the active frame and updates covalent bonds.
    void moveSelectedAtoms(const vec3& offset) {
        if (selected_atoms.empty()) return;

        Geometry* geom = trajectory.currentGeometry();

        if (!geom) return;

        geom->moveAtoms(selected_atoms, offset);

        geom->rebind(binding_factor, bond_radius);
    }

    /// Restores the active geometry frame to its loaded atom positions and updates covalent bonds.
    void resetGeometry() {
        Geometry* geom = trajectory.currentGeometry();

        if (!geom) return;

        geom->resetPositions();

        geom->rebind(binding_factor, bond_radius);

        mouse.dragging_atoms = false;
        mouse.drag_moved     = false;
    }

    /// Clears the collection of currently selected atom indices.
    void clearSelection() {
        selected_atoms.clear();

        mouse.dragging_atoms = false;
        mouse.drag_moved     = false;
    }

    /// Exports the current active molecular geometry frame to an XYZ file on disk.
    bool saveGeometry(const std::string& filepath) const {
        const Geometry* geom = trajectory.currentGeometry();

        if (!geom) return false;

        std::string path = filepath;

        if (path.size() < 4 || path.compare(path.size() - 4, 4, ".xyz") != 0) {
            path += ".xyz";
        }

        return geom->save(path);
    }

    /// Exports all trajectory frames to an XYZ file on disk.
    bool saveTrajectory(const std::string& filepath) const {
        if (trajectory.empty()) return false;

        std::string path = filepath;

        if (path.size() < 4 || path.compare(path.size() - 4, 4, ".xyz") != 0) {
            path += ".xyz";
        }

        return trajectory.save(path);
    }
};
