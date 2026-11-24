#include "viewer_trajectory.hpp"
#include "viewer_ptable.hpp"

#include <fstream>
#include <sstream>

#define BINDINGFACTOR 0.013
#define BONDSIZE 0.09
#define ATOMSIZEFACTOR 0.007

Geometry Geometry::Load(std::istream& is, int natoms) {
    Geometry geom; std::string line;

    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(is, line)) break;

        std::istringstream iss(line); std::string atom;

        float x = 0, y = 0, z = 0;

        if (!(iss >> atom >> x >> y >> z)) continue;

        const AtomInfo& info = get_atom_info(atom);

        mat4 sc = scale(mat4(1.0), vec3(ATOMSIZEFACTOR * info.radius)); mat4 tr = translate(mat4(1.0), vec3(x, y, z));

        geom.objects.push_back({ tr, mat4(1.0), sc, atom });
    }

    geom.rebind(BINDINGFACTOR);

    return geom;
}

vec3 Geometry::getCenter() const {
    vec3 center(0);

    float count = 0.0;

    for (const auto& obj : objects) {
        if (obj.name != "bond" && obj.name != "El") {
            center += obj.getPosition(); count += 1;
        }
    }

    return count > 0 ? center / count : vec3(0);
}

void Geometry::moveBy(const vec3& offset) {
    for (auto& obj : objects) {
        obj.translate = translate(obj.translate, offset);
    }
}

void Geometry::rebind(float factor) {
    std::vector<Object> atoms;

    for (const auto& obj : objects) {
        if (obj.name != "bond") {
            atoms.push_back(obj);
        }
    }

    std::vector<Object> all_objects = atoms;

    size_t n_atoms = atoms.size();

    for (size_t i = 0; i < n_atoms; ++i) {
        for (size_t j = i + 1; j < n_atoms; ++j) {
            if (atoms[i].name == "El" || atoms[j].name == "El") continue;

            vec3 pi = atoms[i].getPosition();
            vec3 pj = atoms[j].getPosition();

            vec3 vector = pj - pi;

            float dist = length(vector);

            if (dist < 1e-4) continue;

            const AtomInfo& info_i = get_atom_info(atoms[i].name);
            const AtomInfo& info_j = get_atom_info(atoms[j].name);

            float thresh = factor * (info_i.covalent + info_j.covalent);

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

                mat4 sc = scale(mat4(1), vec3(BONDSIZE, dist * 0.5, BONDSIZE));

                mat4 tr = translate(mat4(1), mid);

                all_objects.push_back({ tr, rot, sc, "bond" });
            }
        }
    }

    objects = std::move(all_objects);
}

void Geometry::render(const Shader& shader, const std::unordered_map<std::string, Mesh>& meshes) const {
    for (const auto& obj : objects) {
        auto it = meshes.find(obj.name);

        if (it != meshes.end()) {
            it->second.render(shader, obj.getModel());
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

void Trajectory::moveBy(const vec3& offset) {
    for (auto& g : geoms) {
        g.moveBy(offset);
    }
}

float Trajectory::getMaxRadius() const {
    if (geoms.empty()) return 3;

    float max_r = 0;

    for (const auto& obj : geoms[0].getObjects()) {
        if (obj.name != "bond") {
            max_r = std::max(max_r, length(obj.getPosition()));
        }
    }

    return max_r > 0 ? max_r : 3;
}

void Trajectory::render(const Shader& shader, const std::unordered_map<std::string, Mesh>& meshes) {
    if (geoms.empty()) return;

    auto now = std::chrono::high_resolution_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - timestamp).count();

    if (elapsed > wait) {
        if (!paused && geoms.size() > 1) {
            frame = (frame + (int)(elapsed / wait)) % (int)geoms.size();
        }

        timestamp = now;
    }

    geoms.at(frame).render(shader, meshes);
}
