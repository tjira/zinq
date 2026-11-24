#include "viewer_mesh.hpp"
#include "viewer_shader.hpp"

#include <glad/gl.h>

#include <algorithm>

Buffer::Buffer() {
    generate();
}

Buffer::Buffer(const std::vector<Vertex>& data) : data(data) {
    generate();
}

Buffer::Buffer(const Buffer& other) : data(other.data) {
    generate();
}

Buffer::Buffer(Buffer&& other) noexcept : data(std::move(other.data)), vao(other.vao), vbo(other.vbo) {
    other.vao = 0;
    other.vbo = 0;
}

Buffer::~Buffer() {
    destroy();
}

Buffer& Buffer::operator=(const Buffer& other) {
    if (this != &other) {
        destroy(); data = other.data; generate();
    }

    return *this;
}

Buffer& Buffer::operator=(Buffer&& other) noexcept {
    if (this != &other) {
        destroy();

        data = std::move(other.data);

        vao = other.vao;
        vbo = other.vbo;

        other.vao = 0;
        other.vbo = 0;
    }

    return *this;
}

void Buffer::bind() const {
    glBindVertexArray(vao);
}

void Buffer::generate() {
    glGenVertexArrays(1, &vao); glGenBuffers(1, &vbo); glBindVertexArray(vao); glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(Vertex), data.data(), GL_STATIC_DRAW);
}

void Buffer::destroy() {
    if (vbo) {
        glDeleteBuffers(1, &vbo); vbo = 0;
    }

    if (vao) {
        glDeleteVertexArrays(1, &vao); vao = 0;
    }
}

Mesh::Mesh(const std::vector<Vertex>& data) : buffer(data) {}

Mesh Mesh::Cylinder(int sectors, bool smooth) {
    std::vector<Vertex> data;

    for (int j = 0; j < sectors; j++) {
        float theta0 = 2 * static_cast<float>(M_PI) / static_cast<float>(sectors) * static_cast<float>(j + 0);
        float theta1 = 2 * static_cast<float>(M_PI) / static_cast<float>(sectors) * static_cast<float>(j + 1);

        data.push_back({ { std::cos(theta0),  1, std::sin(theta0) } });
        data.push_back({ { std::cos(theta1),  1, std::sin(theta1) } });
        data.push_back({ { std::cos(theta1), -1, std::sin(theta1) } });

        data.push_back({ { std::cos(theta0),  1, std::sin(theta0) } });
        data.push_back({ { std::cos(theta1), -1, std::sin(theta1) } });
        data.push_back({ { std::cos(theta0), -1, std::sin(theta0) } });
    }

    for (size_t i = 0; i < data.size(); i += 3) {
        vec3 v1 = data[i + 1].position - data[i].position;
        vec3 v2 = data[i + 2].position - data[i].position;

        vec3 n = smooth ? vec3(data[i].position.x, 0.0, data[i].position.z) : cross(v1, v2);

        data[i + 0].normal = normalize(smooth ? vec3(data[i + 0].position.x, 0.0, data[i + 0].position.z) : n);
        data[i + 1].normal = normalize(smooth ? vec3(data[i + 1].position.x, 0.0, data[i + 1].position.z) : n);
        data[i + 2].normal = normalize(smooth ? vec3(data[i + 2].position.x, 0.0, data[i + 2].position.z) : n);
    }

    return Mesh(data);
}

Mesh Mesh::Icosphere(int subdivisions, bool smooth) {
    std::vector<Vertex> data;

    float k = (1.0 + std::sqrt(5.0)) / 2.0;

    auto push_tri = [&](vec3 a, vec3 b, vec3 c) {
        data.push_back({ normalize(a) });
        data.push_back({ normalize(b) });
        data.push_back({ normalize(c) });
    };

    push_tri({ -1,  k,  0 }, { -k,  0,  1 }, {  0,  1,  k });
    push_tri({ -1,  k,  0 }, {  0,  1,  k }, {  1,  k,  0 });
    push_tri({ -1,  k,  0 }, {  1,  k,  0 }, {  0,  1, -k });
    push_tri({ -1,  k,  0 }, {  0,  1, -k }, { -k,  0, -1 });
    push_tri({ -1,  k,  0 }, { -k,  0, -1 }, { -k,  0,  1 });

    push_tri({  1,  k,  0 }, {  0,  1,  k }, {  k,  0,  1 });
    push_tri({  0,  1,  k }, { -k,  0,  1 }, {  0, -1,  k });
    push_tri({ -k,  0,  1 }, { -k,  0, -1 }, { -1, -k,  0 });
    push_tri({ -k,  0, -1 }, {  0,  1, -k }, {  0, -1, -k });
    push_tri({  0,  1, -k }, {  1,  k,  0 }, {  k,  0, -1 });

    push_tri({  0, -1,  k }, {  k,  0,  1 }, {  0,  1,  k });
    push_tri({ -1, -k,  0 }, {  0, -1,  k }, { -k,  0,  1 });
    push_tri({  0, -1, -k }, { -1, -k,  0 }, { -k,  0, -1 });
    push_tri({  k,  0, -1 }, {  0, -1, -k }, {  0,  1, -k });
    push_tri({  k,  0,  1 }, {  k,  0, -1 }, {  1,  k,  0 });

    push_tri({  1, -k,  0 }, {  0, -1,  k }, { -1, -k,  0 });
    push_tri({  1, -k,  0 }, { -1, -k,  0 }, {  0, -1, -k });
    push_tri({  1, -k,  0 }, {  0, -1, -k }, {  k,  0, -1 });
    push_tri({  1, -k,  0 }, {  k,  0, -1 }, {  k,  0,  1 });
    push_tri({  1, -k,  0 }, {  k,  0,  1 }, {  0, -1,  k });

    for (int s = 0; s < subdivisions; ++s) {
        std::vector<Vertex> subdivided;

        for (size_t j = 0; j < data.size(); j += 3) {
            vec3 p1 = data[j + 0].position;
            vec3 p2 = data[j + 1].position;
            vec3 p3 = data[j + 2].position;

            vec3 p4 = normalize((p1 + p2) * 0.5);
            vec3 p5 = normalize((p2 + p3) * 0.5);
            vec3 p6 = normalize((p3 + p1) * 0.5);

            subdivided.push_back({ p1 }); subdivided.push_back({ p4 }); subdivided.push_back({ p6 });
            subdivided.push_back({ p4 }); subdivided.push_back({ p2 }); subdivided.push_back({ p5 });
            subdivided.push_back({ p6 }); subdivided.push_back({ p5 }); subdivided.push_back({ p3 });
            subdivided.push_back({ p4 }); subdivided.push_back({ p5 }); subdivided.push_back({ p6 });
        }

        data = subdivided;
    }

    for (size_t i = 0; i < data.size(); i += 3) {
        vec3 v1 = data[i + 1].position - data[i].position;
        vec3 v2 = data[i + 2].position - data[i].position;

        vec3 n = normalize(cross(v1, v2));

        data[i + 0].normal = smooth ? data[i + 0].position : n;
        data[i + 1].normal = smooth ? data[i + 1].position : n;
        data[i + 2].normal = smooth ? data[i + 2].position : n;
    }

    return Mesh(data);
}

void Mesh::setColor(const vec3& color) {
    std::vector<Vertex> data = buffer.getData();

    for (auto& v : data) {
        v.color = color;
    }

    buffer = Buffer(data);
}

void Mesh::render(const Shader& shader, const mat4& transform) const {
    shader.use(); shader.set("u_model", transform);

    buffer.bind();

    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)buffer.getSize());
}
