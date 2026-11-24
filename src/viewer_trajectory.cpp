#include "viewer_trajectory.hpp"
#include "viewer_ptable.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>

static constexpr float kBindingFactor  = 0.013f;
static constexpr float kBondSize       = 0.090f;
static constexpr float kAtomSizeFactor = 0.007f;

Geometry Geometry::Load(std::istream& is, int natoms) {
    Geometry geom; geom.atoms.reserve(natoms);

    std::string line;

    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(is, line)) break;

        std::istringstream iss(line); std::string symbol;

        float x = 0, y = 0, z = 0;

        if (!(iss >> symbol >> x >> y >> z)) continue;

        const AtomInfo& info = get_atom_info(symbol);

        bool dummy = (symbol == "El" || symbol == "el");

        geom.atoms.push_back({ vec3(x, y, z), info.color, kAtomSizeFactor * info.radius, info.covalent, dummy });
    }

    geom.rebind(kBindingFactor);

    return geom;
}

vec3 Geometry::getCenter() const {
    vec3 center(0);

    float count = 0;

    for (const auto& atom : atoms) {
        if (!atom.is_dummy) {
            center += atom.position; count += 1;
        }
    }

    return count > 0 ? center / count : vec3(0);
}

void Geometry::moveBy(const vec3& offset) {
    for (auto& atom : atoms) {
        atom.position += offset;
    }

    for (auto& bond : bonds) {
        bond.model.at(0, 3) += offset.x;
        bond.model.at(1, 3) += offset.y;
        bond.model.at(2, 3) += offset.z;
    }
}

void Geometry::rebind(float factor) {
    bonds.clear();

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

            if (dist < 1e-4) continue;

            float thresh = factor * (atoms[i].covalent + atoms[j].covalent);

            if (dist < thresh) {
                vec3 mid = (pi + pj) * 0.5; vec3 up(0, 1, 0);

                vec3 v_norm = vector / dist;
                vec3 cr = cross(up, v_norm);

                float cr_len = length(cr);

                mat4 rot(1);

                if (cr_len > 1e-6) {
                    float angle = std::atan2(cr_len, dot(up, v_norm));

                    rot = rotate(mat4(1), angle, cr / cr_len);
                }

                else if (v_norm.y < 0) {
                    rot = rotate(mat4(1), static_cast<float>(M_PI), vec3(1, 0, 0));
                }

                mat4 sc = scale(mat4(1), vec3(kBondSize, dist * 0.5, kBondSize)); mat4 tr = translate(mat4(1), mid);

                bonds.push_back({ tr * rot * sc });
            }
        }
    }
}

void Geometry::render(const Shader& shader, const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color) const {
    for (const auto& atom : atoms) {
        mat4 model = translate(mat4(1.0f), atom.position);

        model = scale(model, vec3(atom.radius * atom_scale));

        sphere.render(shader, model, atom.color);
    }

    if (show_bonds) {
        for (const auto& bond : bonds) {
            cylinder.render(shader, bond.model, bond_color);
        }
    }
}

Trajectory Trajectory::Load(const std::string& filename) {
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

        traj.timestamp = std::chrono::high_resolution_clock::now();
    }

    return traj;
}

void Trajectory::stepFrame(int delta) {
    if (geoms.empty()) return;

    int n = static_cast<int>(geoms.size());

    frame = (frame + delta) % n;

    if (frame < 0) frame += n;

    timestamp = std::chrono::high_resolution_clock::now();
}

void Trajectory::setFrame(int f) {
    if (geoms.empty()) return;

    int n = static_cast<int>(geoms.size());

    frame = std::clamp(f, 0, n - 1);

    timestamp = std::chrono::high_resolution_clock::now();
}

void Trajectory::adjustSpeed(float factor) {
    wait_ms = std::clamp(wait_ms * factor, 4.0f, 500.0f);
}

void Trajectory::moveBy(const vec3& offset) {
    for (auto& g : geoms) {
        g.moveBy(offset);
    }
}

float Trajectory::getMaxRadius() const {
    if (geoms.empty()) return 3;

    float max_r = 0;

    for (const auto& geom : geoms) {
        for (const auto& atom : geom.getAtoms()) {
            if (!atom.is_dummy) {
                max_r = std::max(max_r, length(atom.position));
            }
        }
    }

    return max_r > 0 ? max_r : 3;
}

void Trajectory::render(const Shader& shader, const Mesh& sphere, const Mesh& cylinder, float atom_scale, bool show_bonds, const vec3& bond_color) {
    if (geoms.empty()) return;

    if (!paused && geoms.size() > 1) {
        auto now = std::chrono::high_resolution_clock::now();

        float elapsed = std::chrono::duration<float, std::milli>(now - timestamp).count();

        if (elapsed >= wait_ms) {
            int advance = static_cast<int>(elapsed / wait_ms);

            frame = (frame + advance) % static_cast<int>(geoms.size());

            timestamp += std::chrono::milliseconds(static_cast<long long>(advance * wait_ms));
        }
    }

    else {
        timestamp = std::chrono::high_resolution_clock::now();
    }

    geoms.at(frame).render(shader, sphere, cylinder, atom_scale, show_bonds, bond_color);
}
