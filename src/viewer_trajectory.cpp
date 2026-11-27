#include "viewer_trajectory.hpp"
#include "viewer_ptable.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>

static constexpr float kBindingFactor  = 0.013f;
static constexpr float kBondSize       = 0.090f;
static constexpr float kAtomSizeFactor = 0.007f;

/// Parses an XYZ format stream to construct a molecular geometry with connectivity bonds.
Geometry Geometry::Load(std::istream& is, int natoms) {
    Geometry geom; geom.atoms.reserve(natoms);

    std::string line;

    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(is, line)) break;

        std::istringstream iss(line); std::string symbol;

        float x = 0.0f, y = 0.0f, z = 0.0f;

        if (!(iss >> symbol >> x >> y >> z)) continue;

        const AtomInfo& info = get_atom_info(symbol);

        bool dummy = (symbol == "El" || symbol == "el");

        geom.atoms.push_back({ vec3(x, y, z), info.color, kAtomSizeFactor * info.radius, info.covalent, dummy });
    }

    return geom;
}

/// Calculates the geometric centroid of all atoms in the geometry.
vec3 Geometry::getCenter() const {
    vec3 center(0.0f);

    float count = 0.0f;

    for (const auto& atom : atoms) {
        if (!atom.is_dummy) {
            center += atom.position; count += 1.0f;
        }
    }

    return count > 0.0f ? center / count : vec3(0.0f);
}

/// Translates all atom positions and bond models by a spatial offset vector.
void Geometry::moveBy(const vec3& offset) {
    for (auto& atom : atoms) {
        atom.position += offset;
    }

    for (auto& bond : bonds) {
        bond.model.at(0, 3) += offset.x;
        bond.model.at(1, 3) += offset.y;
        bond.model.at(2, 3) += offset.z;
    }

    instances_dirty = true;
}

/// Computes covalent bond connectivity models based on pairwise covalent radii.
void Geometry::rebind(float factor, float bond_size) {
    bonds.clear();

    instances_dirty = true;

    size_t n = atoms.size();

    if (n < 2) return;

    bonds.reserve(n * 2);

    for (size_t i = 0; i < n; ++i) {
        if (atoms[i].is_dummy) continue;

        for (size_t j = i + 1; j < n; ++j) {
            if (atoms[j].is_dummy) continue;

            vec3 pi = atoms[i].position;
            vec3 pj = atoms[j].position;

            vec3 vector = pj - pi;

            float dist = length(vector);

            if (dist < 1e-4f) continue;

            float thresh = factor * (atoms[i].covalent + atoms[j].covalent);

            if (dist < thresh) {
                vec3 mid = (pi + pj) * 0.5f; vec3 up(0.0f, 1.0f, 0.0f);

                vec3 v_norm = vector / dist;
                vec3 cr = cross(up, v_norm);

                float cr_len = length(cr);

                mat4 rot(1.0f);

                if (cr_len > 1e-6f) {
                    float angle = std::atan2(cr_len, dot(up, v_norm));

                    rot = rotate(mat4(1.0f), angle, cr / cr_len);
                }

                else if (v_norm.y < 0.0f) {
                    rot = rotate(mat4(1.0f), static_cast<float>(M_PI), vec3(1.0f, 0.0f, 0.0f));
                }

                mat4 sc = scale(mat4(1.0f), vec3(bond_size, dist * 0.5f, bond_size)); mat4 tr = translate(mat4(1.0f), mid);

                bonds.push_back({ tr * rot * sc });
            }
        }
    }
}

/// Renders atoms as spheres and covalent bonds as cylinders using instanced draw calls.
void Geometry::render(const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color) const {
    if (!atoms.empty()) {
        bool atom_dirty = instances_dirty || (cached_atom_scale != atom_scale) || (cached_atom_instances.size() != atoms.size());

        if (atom_dirty) {
            cached_atom_instances.clear();

            cached_atom_instances.reserve(atoms.size());

            for (const auto& atom : atoms) {
                mat4 model = translate(mat4(1.0f), atom.position);

                model = scale(model, vec3(atom.radius * atom_scale));

                cached_atom_instances.push_back({ model, atom.color });
            }

            cached_atom_scale = atom_scale;
        }

        sphere.renderInstanced(cached_atom_instances, atom_dirty);
    }

    if (show_bonds && !bonds.empty()) {
        bool bond_dirty = instances_dirty || (cached_bond_color != bond_color) || (cached_bond_instances.size() != bonds.size());

        if (bond_dirty) {
            cached_bond_instances.clear();

            cached_bond_instances.reserve(bonds.size());

            for (const auto& bond : bonds) {
                cached_bond_instances.push_back({ bond.model, bond_color });
            }

            cached_bond_color = bond_color;
        }

        cylinder.renderInstanced(cached_bond_instances, bond_dirty);
    }

    instances_dirty = false;
}

/// Reads a multi-frame molecular trajectory file from disk and calculates covalent connectivity.
Trajectory Trajectory::Load(const std::string& filename, float factor, float bond_size) {
    Trajectory traj;

    std::ifstream file(filename);

    if (!file) return traj;

    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        int natoms = 0;

        if (!(iss >> natoms) || natoms <= 0) continue;

        if (!std::getline(file, line)) break;

        traj.geoms.push_back(Geometry::Load(file, natoms));
    }

    if (!traj.geoms.empty()) {
        traj.moveBy(-traj.geoms[0].getCenter());

        traj.rebind(factor, bond_size);

        traj.timestamp = std::chrono::high_resolution_clock::now();
    }

    return traj;
}

/// Advances or rewinds the active frame by the specified frame count.
void Trajectory::stepFrame(int delta) {
    if (geoms.empty()) return;

    int n = static_cast<int>(geoms.size());

    frame = (frame + delta) % n;

    if (frame < 0) frame += n;

    timestamp = std::chrono::high_resolution_clock::now();

    time_accumulator_ms = 0.0;
}

/// Clamps and sets the active trajectory frame index.
void Trajectory::setFrame(int f) {
    if (geoms.empty()) return;

    int n = static_cast<int>(geoms.size());

    frame = std::clamp(f, 0, n - 1);

    timestamp = std::chrono::high_resolution_clock::now();

    time_accumulator_ms = 0.0;
}

/// Multiplies playback frame delay by a scaling factor.
void Trajectory::adjustSpeed(float factor) {
    wait_ms = std::clamp(wait_ms * factor, 4.0f, 500.0f);
}

/// Centers and translates all geometry frames by a spatial offset vector.
void Trajectory::moveBy(const vec3& offset) {
    for (auto& g : geoms) {
        g.moveBy(offset);
    }
}

/// Calculates the maximum radial distance from origin across all atoms in all frames.
float Trajectory::getMaxRadius() const {
    if (geoms.empty()) return 3.0f;

    float max_r = 0.0f;

    for (const auto& geom : geoms) {
        for (const auto& atom : geom.getAtoms()) {
            if (!atom.is_dummy) {
                max_r = std::max(max_r, length(atom.position));
            }
        }
    }

    return max_r > 0.0f ? max_r : 3.0f;
}

/// Recomputes covalent bond connectivity models across all trajectory frames using given bonding parameters.
void Trajectory::rebind(float factor, float bond_size) {
    for (auto& g : geoms) {
        g.rebind(factor, bond_size);
    }
}

/// Advances frame timing and dispatches geometry render calls for the current frame.
void Trajectory::render(const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color) {
    if (geoms.empty()) return;

    if (!paused && geoms.size() > 1) {
        auto now = std::chrono::high_resolution_clock::now();

        double elapsed = std::chrono::duration<double, std::milli>(now - timestamp).count();

        timestamp = now;

        time_accumulator_ms += elapsed;

        if (time_accumulator_ms >= wait_ms) {
            int advance = static_cast<int>(time_accumulator_ms / wait_ms);
            time_accumulator_ms -= static_cast<double>(advance) * wait_ms;

            int n = static_cast<int>(geoms.size());

            if (loop_mode == LoopMode::Loop) {
                frame = (frame + advance) % n;
            }

            else if (loop_mode == LoopMode::PingPong) {
                for (int i = 0; i < advance; ++i) {
                    if (direction > 0 && frame >= n - 1) {
                        direction = -1;
                    }

                    else if (direction < 0 && frame <= 0) {
                        direction = 1;
                    }

                    frame += direction;
                }

                frame = std::clamp(frame, 0, n - 1);
            }

            else if (loop_mode == LoopMode::Once) {
                if (frame + advance >= n - 1) {
                    frame = n - 1; paused = true; time_accumulator_ms = 0.0;
                }

                else {
                    frame += advance;
                }
            }
        }
    }

    else {
        timestamp = std::chrono::high_resolution_clock::now();

        time_accumulator_ms = 0.0;
    }

    geoms.at(frame).render(sphere, cylinder, atom_scale, show_bonds, bond_color);
}
