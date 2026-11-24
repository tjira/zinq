#pragma once

#include "viewer_math.hpp"
#include "viewer_mesh.hpp"
#include "viewer_shader.hpp"

#include <chrono>
#include <istream>
#include <string>
#include <unordered_map>
#include <vector>

/// Represents a single molecular geometry frame with spatial atom centers and covalent bonds.
class Geometry {
public:
    struct Object {
        mat4 translate = mat4(1.0); mat4 rotate = mat4(1.0); mat4 scale = mat4(1.0); std::string name;

        mat4 getModel() const { return translate * rotate * scale; }

        vec3 getPosition() const { return translate.col3(3); }
    };

    Geometry() = default;

    static Geometry Load(std::istream& is, int natoms);

    const std::vector<Object>& getObjects() const { return objects; } vec3 getCenter() const;

    void moveBy(const vec3& offset);

    void rebind(float factor);

    void render(const Shader& shader, const std::unordered_map<std::string, Mesh>& meshes) const;

private:
    std::vector<Object> objects;
};

/// Encapsulates a multi-frame molecular trajectory with frame timing, centering, and render dispatch.
class Trajectory {
public:
    Trajectory() = default;

    static Trajectory Load(const std::string& filename);

    bool empty() const { return geoms.empty(); } void togglePause() { paused = !paused; }

    void moveBy(const vec3& offset);

    float getMaxRadius() const;

    void render(const Shader& shader, const std::unordered_map<std::string, Mesh>& meshes);

private:
    std::chrono::high_resolution_clock::time_point timestamp;

    std::vector<Geometry> geoms;

    bool paused = false; float wait = 16; int frame = 0;
};
