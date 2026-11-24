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
    /// Initializes an empty OpenGL vertex array and buffer object.
    Buffer();

    /// Destroys OpenGL buffer objects and releases GPU memory.
    ~Buffer();

    /// Allocates and uploads vertex coordinate and normal data to GPU buffers.
    Buffer(const std::vector<Vertex>& data);

    /// Copy constructor creating an independent vertex buffer instance.
    Buffer(const Buffer&  other);

    /// Move constructor transferring GPU vertex array and buffer ownership.
    Buffer(Buffer&& other) noexcept;

    /// Copy assignment reallocating vertex buffer state.
    Buffer& operator=(const Buffer&  other);

    /// Move assignment transferring GPU buffer resource handles.
    Buffer& operator=(Buffer&& other) noexcept;

    /// Returns the total number of vertices stored in the buffer.
    size_t getSize() const { return count; }

    /// Binds the OpenGL vertex array object for draw dispatch.
    void bind() const;

private:
    unsigned int vao = 0;
    unsigned int vbo = 0;

    size_t count = 0;

    /// Generates GPU vertex attributes and uploads triangle vertex records.
    void generate(const std::vector<Vertex>& data);

    /// Deletes GPU vertex array and vertex buffer resources.
    void destroy();
};

/// 3D polygonal mesh container supporting procedural primitives and shader rendering.
class Mesh {
public:
    /// Constructs a default empty polygonal mesh container.
    Mesh() = default;

    /// Constructs a polygonal mesh initialized with given vertex geometry data.
    Mesh(const std::vector<Vertex>& data);

    /// Procedurally generates a cylindrical polygonal mesh for molecular bonds.
    static Mesh Cylinder(int sectors, bool smooth);

    /// Procedurally generates an icosphere mesh for spherical atom centers.
    static Mesh Icosphere(int subdivisions, bool smooth);

    /// Dispatches drawing commands for the mesh using active shader uniforms.
    void render(const Shader& shader, const mat4& transform, const vec3& color) const;

private:
    Buffer buffer;
};
