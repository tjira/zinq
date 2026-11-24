#pragma once

#include "viewer_math.hpp"
#include "viewer_mesh.hpp"
#include "viewer_shader.hpp"

#include <chrono>
#include <istream>
#include <string>
#include <vector>

/// Represents a single molecular geometry frame with spatial atom centers and covalent bonds.
class Geometry {
public:
    struct Atom {
        vec3 position; vec3 color = vec3(1);

        float radius   =  0.5;
        float covalent = 75.0;

        bool is_dummy = false;
    };

    struct Bond {
        mat4 model = mat4(1);
    };

    Geometry() = default;

    static Geometry Load(std::istream& is, int natoms);

    const std::vector<Atom>& getAtoms() const { return atoms; }
    const std::vector<Bond>& getBonds() const { return bonds; }

    vec3 getCenter() const;

    void moveBy(const vec3& offset);

    void rebind(float factor);

    void render(const Shader& shader, const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color) const;

private:
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
};

/// Encapsulates a multi-frame molecular trajectory with frame timing, centering, and render dispatch.
class Trajectory {
public:
    Trajectory() = default;

    static Trajectory Load(const std::string& filename);

    bool  empty()        const { return geoms.empty();                  }
    int   frameCount()   const { return static_cast<int>(geoms.size()); }
    int   currentFrame() const { return frame;                          }
    bool  isPaused()     const { return paused;                         }
    float getWaitMs()    const { return wait_ms;                        }

    void togglePause()     { paused = !paused; }
    void setPaused(bool p) { paused = p;       }

    void stepFrame(int delta); void setFrame(int f); void adjustSpeed(float factor);

    void moveBy(const vec3& offset);

    float getMaxRadius() const;

    void render(const Shader& shader, const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color);

private:
    std::chrono::high_resolution_clock::time_point timestamp;

    std::vector<Geometry> geoms;

    bool paused = false; float wait_ms = 16; int frame = 0;
};
