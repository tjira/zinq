#pragma once

#include "viewer_math.hpp"

#include <string>

/// Compiles and encapsulates an OpenGL GLSL shader program and its active uniform bindings.
class Shader {
public:
    /// Compiles vertex and fragment shaders and links them into an executable shader program.
    Shader(const std::string& vertex_src, const std::string& fragment_src);

    /// Deletes the compiled OpenGL shader program handle.
    ~Shader();

    /// Activates the compiled shader program in the OpenGL rendering pipeline.
    void use() const;

    /// Checks whether the OpenGL shader program was successfully compiled and linked.
    bool isValid() const { return id != 0; }

    /// Sets a floating-point scalar uniform variable in the shader program.
    void set(const std::string& name, float value) const;

    /// Sets a 3D Cartesian vector uniform variable in the shader program.
    void set(const std::string& name, const vec3& value) const;

    /// Sets a 4x4 transformation matrix uniform variable in the shader program.
    void set(const std::string& name, const mat4& value) const;

private:
    unsigned int id = 0;

    /// Queries and verifies shader compilation or program link status.
    bool errorCheck(unsigned int shader, const std::string& stage) const;
};
