#pragma once

#include "viewer_math.hpp"

#include <vector>

class Shader;

/// Represents a 3D vertex with spatial coordinates, surface normal, and RGB color.
struct Vertex {
    vec3 position;

    vec3 normal = vec3(0.0); vec3 color = vec3(1.0);
};

/// OpenGL vertex array and buffer object abstraction for managing GPU geometry data.
class Buffer {
public:
    Buffer(); ~Buffer();

    Buffer(const std::vector<Vertex>& data);

    Buffer(const Buffer& other); Buffer(Buffer&& other) noexcept;

    Buffer& operator=(const Buffer& other); Buffer& operator=(Buffer&& other) noexcept;

    std::vector<Vertex> getData() const { return data; }

    size_t getSize() const { return data.size(); }

    void bind() const;

private:
    std::vector<Vertex> data;

    unsigned int vao = 0; unsigned int vbo = 0;

    void generate(); void destroy();
};

/// 3D polygonal mesh container supporting procedural primitives and shader rendering.
class Mesh {
public:
    Mesh() = default;

    Mesh(const std::vector<Vertex>& data);

    static Mesh Cylinder(int sectors, bool smooth); static Mesh Icosphere(int subdivisions, bool smooth);

    void setColor(const vec3& color);

    void render(const Shader& shader, const mat4& transform = mat4(1.0)) const;

private:
    Buffer buffer;
};
