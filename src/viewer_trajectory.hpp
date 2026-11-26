#pragma once

#include "viewer_math.hpp"
#include "viewer_mesh.hpp"

#include <chrono>
#include <istream>
#include <string>
#include <vector>

/// Represents a single molecular geometry frame with spatial atom centers and covalent bonds.
class Geometry {
public:
    struct Atom {
        vec3 position; vec3 color = vec3(1.0f);

        float radius   =  0.5f;
        float covalent = 75.0f;

        bool is_dummy = false;
    };

    struct Bond {
        mat4 model = mat4(1.0f);
    };

    /// Constructs a default empty molecular geometry structure.
    Geometry() = default;

    /// Parses an XYZ format stream to construct a molecular geometry with connectivity bonds.
    static Geometry Load(std::istream& is, int natoms);

    /// Returns a constant reference to the collection of atoms in the geometry.
    const std::vector<Atom>& getAtoms() const { return atoms; }

    /// Returns a constant reference to the collection of bond models in the geometry.
    const std::vector<Bond>& getBonds() const { return bonds; }

    /// Calculates the geometric centroid of all atoms in the geometry.
    vec3 getCenter() const;

    /// Translates all atom positions and bond models by a spatial offset vector.
    void moveBy(const vec3& offset);

    /// Computes covalent bond connectivity models based on pairwise covalent radii.
    void rebind(float factor, float bond_size = 0.090f);

    /// Renders atoms as spheres and covalent bonds as cylinders using instanced draw calls.
    void render(const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color) const;

private:
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
};

/// Defines playback loop wrapping behavior for multi-frame molecular trajectory visualization.
enum class LoopMode {
    Loop     = 0,
    PingPong = 1,
    Once     = 2
};

/// Encapsulates a multi-frame molecular trajectory with frame timing, centering, and render dispatch.
class Trajectory {
public:
    /// Constructs a default empty molecular trajectory instance.
    Trajectory() = default;

    /// Reads a multi-frame molecular trajectory file from disk.
    static Trajectory Load(const std::string& filename);

    /// Checks whether the trajectory contains any geometry frames.
    bool empty() const { return geoms.empty(); }

    /// Returns the total number of frames in the trajectory.
    int frameCount() const { return static_cast<int>(geoms.size()); }

    /// Returns the zero-based index of the currently active frame.
    int currentFrame() const { return frame; }

    /// Checks whether trajectory playback is paused.
    bool isPaused() const { return paused; }

    /// Returns the active trajectory playback loop mode.
    LoopMode getLoopMode() const { return loop_mode; }

    /// Sets the trajectory playback loop mode.
    void setLoopMode(LoopMode m) { loop_mode = m; }

    /// Returns the playback wait interval between consecutive frames in milliseconds.
    float getWaitMs() const { return wait_ms; }

    /// Toggles between playback and paused states.
    void togglePause() { paused = !paused; }

    /// Explicitly sets the trajectory playback pause state.
    void setPaused(bool p) { paused = p; }

    /// Advances or rewinds the active frame by the specified frame count.
    void stepFrame(int delta);

    /// Clamps and sets the active trajectory frame index.
    void setFrame(int f);

    /// Multiplies playback frame delay by a scaling factor.
    void adjustSpeed(float factor);

    /// Centers and translates all geometry frames by a spatial offset vector.
    void moveBy(const vec3& offset);

    /// Calculates the maximum radial distance from origin across all atoms in all frames.
    float getMaxRadius() const;

    /// Recomputes covalent bond connectivity models across all trajectory frames using given bonding parameters.
    void rebind(float factor, float bond_size = 0.090f);

    /// Advances frame timing and dispatches geometry render calls for the current frame.
    void render(const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color);

private:
    std::chrono::high_resolution_clock::time_point timestamp;

    std::vector<Geometry> geoms;

    bool paused = false; float wait_ms = 16.0f; int frame = 0;

    LoopMode loop_mode = LoopMode::Loop;

    int direction = 1;
};
