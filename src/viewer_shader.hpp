#pragma once

#include "viewer_math.hpp"

#include <string>

/// Compiles and encapsulates an OpenGL GLSL shader program and its active uniform bindings.
class Shader {
public:
    Shader(const std::string& vertex_src, const std::string& fragment_src); ~Shader();

    void use() const;

    void set(const std::string& name, float value) const;

    void set(const std::string& name, const vec3& value) const;
    void set(const std::string& name, const mat4& value) const;

private:
    unsigned int id = 0;

    void errorCheck(unsigned int shader, const std::string& stage) const;
};
