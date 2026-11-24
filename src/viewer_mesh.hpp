#pragma once

#include "viewer_math.hpp"

#include <vector>

class Shader;

/// Represents a 3D vertex with spatial coordinates and surface normal.
struct Vertex {
    vec3 position; vec3 normal = vec3(0.0f);
};

/// OpenGL vertex array and buffer object abstraction for managing GPU geometry data.
class Buffer {
public:
    Buffer(); ~Buffer();

    Buffer(const std::vector<Vertex>& data);

    Buffer(const Buffer&  other);
    Buffer(      Buffer&& other) noexcept;

    Buffer& operator=(const Buffer&  other);
    Buffer& operator=(      Buffer&& other) noexcept;

    size_t getSize() const { return count; }

    void bind() const;

private:
    unsigned int vao = 0;
    unsigned int vbo = 0;

    size_t count = 0;

    void generate(const std::vector<Vertex>& data);

    void destroy();
};

/// 3D polygonal mesh container supporting procedural primitives and shader rendering.
class Mesh {
public:
    Mesh() = default;

    Mesh(const std::vector<Vertex>& data);

    static Mesh Cylinder (int sectors,      bool smooth);
    static Mesh Icosphere(int subdivisions, bool smooth);

    void render(const Shader& shader, const mat4& transform, const vec3& color) const;

private:
    Buffer buffer;
};
