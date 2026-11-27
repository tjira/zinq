#include "viewer_mesh.hpp"

#include <glad/gl.h>

#include <cmath>
#include <vector>

/// Initializes an empty OpenGL vertex array and buffer object.
Buffer::Buffer() {
    std::vector<Vertex> empty; generate(empty);
}

/// Allocates and uploads vertex coordinate and normal data to GPU buffers.
Buffer::Buffer(const std::vector<Vertex>& data) {
    generate(data);
}

/// Move constructor transferring GPU vertex array and buffer ownership.
Buffer::Buffer(Buffer&& other) noexcept : vao(other.vao), vbo(other.vbo), instance_vbo(other.instance_vbo), count(other.count), last_instances(other.last_instances), last_instance_count(other.last_instance_count) {
    other.vao          = 0;
    other.vbo          = 0;
    other.instance_vbo = 0;

    other.last_instances      = nullptr;
    other.last_instance_count =       0;

    other.count = 0;
}

/// Destroys OpenGL buffer objects and releases GPU memory.
Buffer::~Buffer() {
    destroy();
}

/// Move assignment transferring GPU buffer resource handles.
Buffer& Buffer::operator=(Buffer&& other) noexcept {
    if (this != &other) {
        destroy();

        vao = other.vao;
        vbo = other.vbo;

        instance_vbo = other.instance_vbo;

        count = other.count;

        last_instances      =      other.last_instances;
        last_instance_count = other.last_instance_count;

        other.vao          = 0;
        other.vbo          = 0;
        other.instance_vbo = 0;

        other.last_instances      = nullptr;
        other.last_instance_count =       0;

        other.count = 0;
    }

    return *this;
}

/// Generates GPU vertex attributes and uploads triangle vertex records.
void Buffer::generate(const std::vector<Vertex>& data) {
    count = data.size();

    last_instances      = nullptr;
    last_instance_count =       0;

    glGenVertexArrays(1, &vao         );
    glGenBuffers     (1, &vbo         );
    glGenBuffers     (1, &instance_vbo);

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal  ));

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    if (!data.empty()) {
        glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(Vertex), data.data(), GL_STATIC_DRAW);
    }

    glBindBuffer(GL_ARRAY_BUFFER, instance_vbo);

    for (int i = 0; i < 4; ++i) {
        unsigned int loc = 2 + i;

        glEnableVertexAttribArray(loc);

        glVertexAttribPointer(loc, 4, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)(i * 4 * sizeof(float)));

        glVertexAttribDivisor(loc, 1);
    }

    glEnableVertexAttribArray(6);

    glVertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, color));

    glVertexAttribDivisor(6, 1);
}

/// Deletes GPU vertex array and vertex buffer resources.
void Buffer::destroy() {
    if (instance_vbo) {
        glDeleteBuffers(1, &instance_vbo);

        instance_vbo = 0;
    }

    if (vbo) {
        glDeleteBuffers(1, &vbo);

        vbo = 0;
    }

    if (vao) {
        glDeleteVertexArrays(1, &vao);

        vao = 0;
    }

    count = 0;

    last_instances      = nullptr;
    last_instance_count =       0;
}

/// Dispatches instanced triangle draw calls using the uploaded instance buffer.
void Buffer::renderInstanced(const InstanceData* instances, size_t instance_count, bool dirty) const {
    if (count == 0 || instance_count == 0) return;

    if (dirty || last_instances != instances || last_instance_count != instance_count) {
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo);

        glBufferData(GL_ARRAY_BUFFER, instance_count * sizeof(InstanceData), instances, GL_STREAM_DRAW);

        last_instances      =      instances;
        last_instance_count = instance_count;
    }

    glBindVertexArray(vao);

    glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(count), static_cast<GLsizei>(instance_count));
}

/// Constructs a polygonal mesh initialized with given vertex geometry data.
Mesh::Mesh(const std::vector<Vertex>& data) : buffer(data) {}

/// Procedurally generates a cylindrical polygonal mesh for molecular bonds.
Mesh Mesh::Cylinder(int sectors, bool smooth) {
    std::vector<Vertex> data; data.reserve(sectors * 12);

    for (int j = 0; j < sectors; j++) {
        float theta0 = 2.0f * static_cast<float>(M_PI) / static_cast<float>(sectors) * static_cast<float>(j + 0);
        float theta1 = 2.0f * static_cast<float>(M_PI) / static_cast<float>(sectors) * static_cast<float>(j + 1);

        float c0 = std::cos(theta0), s0 = std::sin(theta0);
        float c1 = std::cos(theta1), s1 = std::sin(theta1);

        vec3 p0_top = { c0,  1.0f, s0 };
        vec3 p1_top = { c1,  1.0f, s1 };
        vec3 p1_bot = { c1, -1.0f, s1 };
        vec3 p0_bot = { c0, -1.0f, s0 };

        vec3 n0 = smooth ? vec3(c0, 0.0f, s0) : normalize(cross(p1_top - p0_top, p1_bot - p0_top));

        vec3 n1 = smooth ? vec3(c1, 0.0f, s1) : n0;

        data.push_back({ p0_top, n0 });
        data.push_back({ p1_top, n1 });
        data.push_back({ p1_bot, n1 });

        data.push_back({ p0_top, n0 });
        data.push_back({ p1_bot, n1 });
        data.push_back({ p0_bot, n0 });

        data.push_back({ { 0.0f, 1.0f, 0.0f }, { 0.0f, 1.0f, 0.0f } });

        data.push_back({ p0_top, { 0.0f, 1.0f, 0.0f } });
        data.push_back({ p1_top, { 0.0f, 1.0f, 0.0f } });

        data.push_back({ { 0.0f, -1.0f, 0.0f }, { 0.0f, -1.0f, 0.0f } });

        data.push_back({ p1_bot, { 0.0f, -1.0f, 0.0f } });
        data.push_back({ p0_bot, { 0.0f, -1.0f, 0.0f } });
    }

    return Mesh(data);
}

/// Procedurally generates an icosphere mesh for spherical atom centers.
Mesh Mesh::Icosphere(int subdivisions, bool smooth) {
    std::vector<vec3> positions;

    float k = (1.0f + std::sqrt(5.0f)) * 0.5f;

    auto push_tri = [&](vec3 a, vec3 b, vec3 c) {
        positions.push_back(normalize(a));
        positions.push_back(normalize(b));
        positions.push_back(normalize(c));
    };

    push_tri({ -1.0f,  k,  0.0f }, { -k,  0.0f,  1.0f }, {  0.0f,  1.0f,  k });
    push_tri({ -1.0f,  k,  0.0f }, {  0.0f,  1.0f,  k }, {  1.0f,  k,  0.0f });
    push_tri({ -1.0f,  k,  0.0f }, {  1.0f,  k,  0.0f }, {  0.0f,  1.0f, -k });
    push_tri({ -1.0f,  k,  0.0f }, {  0.0f,  1.0f, -k }, { -k,  0.0f, -1.0f });
    push_tri({ -1.0f,  k,  0.0f }, { -k,  0.0f, -1.0f }, { -k,  0.0f,  1.0f });

    push_tri({  1.0f,  k,  0.0f }, {  0.0f,  1.0f,  k }, {  k,  0.0f,  1.0f });
    push_tri({  0.0f,  1.0f,  k }, { -k,  0.0f,  1.0f }, {  0.0f, -1.0f,  k });
    push_tri({ -k,  0.0f,  1.0f }, { -k,  0.0f, -1.0f }, { -1.0f, -k,  0.0f });
    push_tri({ -k,  0.0f, -1.0f }, {  0.0f,  1.0f, -k }, {  0.0f, -1.0f, -k });
    push_tri({  0.0f,  1.0f, -k }, {  1.0f,  k,  0.0f }, {  k,  0.0f, -1.0f });

    push_tri({  0.0f, -1.0f,  k }, {  k,  0.0f,  1.0f }, {  0.0f,  1.0f,  k });
    push_tri({ -1.0f, -k,  0.0f }, {  0.0f, -1.0f,  k }, { -k,  0.0f,  1.0f });
    push_tri({  0.0f, -1.0f, -k }, { -1.0f, -k,  0.0f }, { -k,  0.0f, -1.0f });
    push_tri({  k,  0.0f, -1.0f }, {  0.0f, -1.0f, -k }, {  0.0f,  1.0f, -k });
    push_tri({  k,  0.0f,  1.0f }, {  k,  0.0f, -1.0f }, {  1.0f,  k,  0.0f });

    push_tri({  1.0f, -k,  0.0f }, {  0.0f, -1.0f,  k }, { -1.0f, -k,  0.0f });
    push_tri({  1.0f, -k,  0.0f }, { -1.0f, -k,  0.0f }, {  0.0f, -1.0f, -k });
    push_tri({  1.0f, -k,  0.0f }, {  0.0f, -1.0f, -k }, {  k,  0.0f, -1.0f });
    push_tri({  1.0f, -k,  0.0f }, {  k,  0.0f, -1.0f }, {  k,  0.0f,  1.0f });
    push_tri({  1.0f, -k,  0.0f }, {  k,  0.0f,  1.0f }, {  0.0f, -1.0f,  k });

    for (int s = 0; s < subdivisions; ++s) {
        std::vector<vec3> subdivided; subdivided.reserve(positions.size() * 4);

        for (size_t j = 0; j < positions.size(); j += 3) {
            vec3 p1 = positions[j + 0];
            vec3 p2 = positions[j + 1];
            vec3 p3 = positions[j + 2];

            vec3 p4 = normalize((p1 + p2) * 0.5f);
            vec3 p5 = normalize((p2 + p3) * 0.5f);
            vec3 p6 = normalize((p3 + p1) * 0.5f);

            subdivided.push_back(p1); subdivided.push_back(p4); subdivided.push_back(p6);
            subdivided.push_back(p4); subdivided.push_back(p2); subdivided.push_back(p5);
            subdivided.push_back(p6); subdivided.push_back(p5); subdivided.push_back(p3);
            subdivided.push_back(p4); subdivided.push_back(p5); subdivided.push_back(p6);
        }

        positions = std::move(subdivided);
    }

    std::vector<Vertex> data; data.reserve(positions.size());

    for (size_t i = 0; i < positions.size(); i += 3) {
        vec3 p1 = positions[i + 0];
        vec3 p2 = positions[i + 1];
        vec3 p3 = positions[i + 2];

        vec3 flat_n = normalize(cross(p2 - p1, p3 - p1));

        data.push_back({ p1, smooth ? p1 : flat_n });
        data.push_back({ p2, smooth ? p2 : flat_n });
        data.push_back({ p3, smooth ? p3 : flat_n });
    }

    return Mesh(data);
}

/// Dispatches instanced drawing commands for the mesh geometry.
void Mesh::renderInstanced(const std::vector<InstanceData>& instances, bool dirty) const {
    buffer.renderInstanced(instances.data(), instances.size(), dirty);
}
