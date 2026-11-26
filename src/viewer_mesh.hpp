#pragma once

#include "viewer_math.hpp"

#include <cstddef>
#include <vector>

/// Represents a 3D vertex with spatial coordinates and surface normal.
struct Vertex {
    vec3 position; vec3 normal = vec3(0.0f);
};

/// Per-instance model transformation matrix and RGB color for instanced mesh rendering.
struct InstanceData {
    mat4 model; vec3 color;
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

    /// Copy constructor and assignment are deleted to prevent accidental resource duplication.
    Buffer(const Buffer&) = delete;

    /// Copy assignment operator is deleted to prevent accidental resource duplication.
    Buffer& operator=(const Buffer&) = delete;

    /// Move constructor transferring GPU vertex array and buffer ownership.
    Buffer(Buffer&& other) noexcept;

    /// Move assignment transferring GPU buffer resource handles.
    Buffer& operator=(Buffer&& other) noexcept;

    /// Returns the total number of vertices stored in the buffer.
    size_t getSize() const { return count; }

    /// Dispatches instanced triangle draw calls using the uploaded instance buffer.
    void renderInstanced(const InstanceData* instances, size_t instance_count) const;

private:
    unsigned int vao          = 0;
    unsigned int vbo          = 0;
    unsigned int instance_vbo = 0;

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

    /// Copy constructor and assignment are deleted to prevent accidental resource duplication.
    Mesh(const Mesh&) = delete;

    /// Copy assignment operator is deleted to prevent accidental resource duplication.
    Mesh& operator=(const Mesh&) = delete;

    /// Move constructor transferring mesh buffer resources.
    Mesh(Mesh&&) noexcept = default;

    /// Move assignment transferring mesh buffer resources.
    Mesh& operator=(Mesh&&) noexcept = default;

    /// Procedurally generates a cylindrical polygonal mesh for molecular bonds.
    static Mesh Cylinder(int sectors, bool smooth);

    /// Procedurally generates an icosphere mesh for spherical atom centers.
    static Mesh Icosphere(int subdivisions, bool smooth);

    /// Dispatches instanced drawing commands for the mesh geometry.
    void renderInstanced(const std::vector<InstanceData>& instances) const;

private:
    Buffer buffer;
};
