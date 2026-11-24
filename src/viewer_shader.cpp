#include "viewer_shader.hpp"

#include <glad/gl.h>

#include <vector>

Shader::Shader(const std::string& vertex_src, const std::string& fragment_src) : id(glCreateProgram()) {
    const char* vs_code = vertex_src.c_str(); const char* fs_code = fragment_src.c_str();

    unsigned int vs = glCreateShader(GL_VERTEX_SHADER); unsigned int fs = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vs, 1, &vs_code, nullptr); glCompileShader(vs); errorCheck(vs, "VERTEX");

    glShaderSource(fs, 1, &fs_code, nullptr); glCompileShader(fs); errorCheck(fs, "FRAGMENT");

    glAttachShader(id, vs);
    glAttachShader(id, fs);

    glLinkProgram(id);

    int success = 0;

    glGetProgramiv(id, GL_LINK_STATUS, &success);

    if (!success) {
        int length = 0;

        glGetProgramiv(id, GL_INFO_LOG_LENGTH, &length);

        std::vector<char> log(length);

        glGetProgramInfoLog(id, length, &length, log.data());

        glDeleteShader(vs);
        glDeleteShader(fs);

        throw std::runtime_error(std::string("Shader Program Link Error: ") + log.data());
    }

    glDetachShader(id, vs);
    glDetachShader(id, fs);

    glDeleteShader(vs);
    glDeleteShader(fs);
}

Shader::~Shader() {
    if (id) {
        glDeleteProgram(id);
    }
}

void Shader::errorCheck(unsigned int shader, const std::string& stage) const {
    int success = 0;

    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success) {
        int length = 0;

        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);

        std::vector<char> log(length);

        glGetShaderInfoLog(shader, length, &length, log.data());

        throw std::runtime_error(stage + " Shader Compile Error: " + log.data());
    }
}

void Shader::use() const {
    glUseProgram(id);
}

void Shader::set(const std::string& name, float value) const {
    glUniform1f(glGetUniformLocation(id, name.c_str()), value);
}

void Shader::set(const std::string& name, const vec3& value) const {
    glUniform3f(glGetUniformLocation(id, name.c_str()), value.x, value.y, value.z);
}

void Shader::set(const std::string& name, const mat4& value) const {
    glUniformMatrix4fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, value.ptr());
}
