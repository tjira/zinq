#pragma once

#include <cmath>

struct vec2 {
    float x = 0, y = 0;

    vec2() = default; vec2(float x_, float y_) : x(x_), y(y_) {}
};

struct vec3 {
    float x = 0, y = 0, z = 0;

    vec3() = default; vec3(float v) : x(v), y(v), z(v) {} vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    vec3 operator-() const { return vec3(-x, -y, -z); }

    vec3 operator+(const vec3& b) const { return vec3(x + b.x, y + b.y, z + b.z); }
    vec3 operator-(const vec3& b) const { return vec3(x - b.x, y - b.y, z - b.z); }

    vec3 operator*(float s) const { return vec3(x * s, y * s, z * s); }
    vec3 operator/(float s) const { return vec3(x / s, y / s, z / s); }

    vec3& operator+=(const vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
};

inline vec3 operator*(float s, const vec3& v) { return vec3(v.x * s, v.y * s, v.z * s); }

inline float dot(const vec3& a, const vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline float length(const vec3& v) { return std::sqrt(dot(v, v)); }

inline vec3 normalize(const vec3& v) {
    float len = length(v); return len > 1e-8 ? v / len : vec3(0);
}

inline vec3 cross(const vec3& a, const vec3& b) {
    return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

struct mat3 {
    float data[9] = {0};

    mat3() = default;

    mat3(float diag) {
        data[0] = diag; data[4] = diag; data[8] = diag;
    }

    float& at(int r, int c) { return data[c * 3 + r]; } float at(int r, int c) const { return data[c * 3 + r]; }

    vec3 operator*(const vec3& v) const {
        return vec3(at(0, 0) * v.x + at(0, 1) * v.y + at(0, 2) * v.z, at(1, 0) * v.x + at(1, 1) * v.y + at(1, 2) * v.z, at(2, 0) * v.x + at(2, 1) * v.y + at(2, 2) * v.z);
    }
};

inline mat3 inverse(const mat3& m) {
    float det1 = m.at(0, 0) * (m.at(1, 1) * m.at(2, 2) - m.at(1, 2) * m.at(2, 1));
    float det2 = m.at(0, 1) * (m.at(1, 0) * m.at(2, 2) - m.at(1, 2) * m.at(2, 0));
    float det3 = m.at(0, 2) * (m.at(1, 0) * m.at(2, 1) - m.at(1, 1) * m.at(2, 0));

    float det = det1 - det2 + det3;

    if (std::abs(det) < 1e-12) return mat3(1.0);

    float inv_det = 1.0 / det;

    mat3 r;

    r.at(0, 0) =  (m.at(1, 1) * m.at(2, 2) - m.at(1, 2) * m.at(2, 1)) * inv_det;
    r.at(0, 1) = -(m.at(0, 1) * m.at(2, 2) - m.at(0, 2) * m.at(2, 1)) * inv_det;
    r.at(0, 2) =  (m.at(0, 1) * m.at(1, 2) - m.at(0, 2) * m.at(1, 1)) * inv_det;
    r.at(1, 0) = -(m.at(1, 0) * m.at(2, 2) - m.at(1, 2) * m.at(2, 0)) * inv_det;
    r.at(1, 1) =  (m.at(0, 0) * m.at(2, 2) - m.at(0, 2) * m.at(2, 0)) * inv_det;
    r.at(1, 2) = -(m.at(0, 0) * m.at(1, 2) - m.at(0, 2) * m.at(1, 0)) * inv_det;
    r.at(2, 0) =  (m.at(1, 0) * m.at(2, 1) - m.at(1, 1) * m.at(2, 0)) * inv_det;
    r.at(2, 1) = -(m.at(0, 0) * m.at(2, 1) - m.at(0, 1) * m.at(2, 0)) * inv_det;
    r.at(2, 2) =  (m.at(0, 0) * m.at(1, 1) - m.at(0, 1) * m.at(1, 0)) * inv_det;

    return r;
}

struct mat4 {
    float data[16] = { 0 };

    mat4() = default;

    mat4(float diag) {
        data[0] = diag; data[5] = diag; data[10] = diag; data[15] = diag;
    }

    float& at(int r, int c) { return data[c * 4 + r]; } float at(int r, int c) const { return data[c * 4 + r]; }

    const float* ptr() const { return data; }

    vec3 col3(int c) const {
        return vec3(at(0, c), at(1, c), at(2, c));
    }

    mat4 operator*(const mat4& b) const {
        mat4 r;

        for (int c = 0; c < 4; ++c) {
            for (int r_idx = 0; r_idx < 4; ++r_idx) {
                r.at(r_idx, c) = at(r_idx, 0) * b.at(0, c) + at(r_idx, 1) * b.at(1, c) + at(r_idx, 2) * b.at(2, c) + at(r_idx, 3) * b.at(3, c);
            }
        }

        return r;
    }
};

inline mat3 to_mat3(const mat4& m) {
    mat3 r;

    for (int c = 0; c < 3; ++c) {
        for (int row = 0; row < 3; ++row) r.at(row, c) = m.at(row, c);
    }

    return r;
}

inline mat4 from_mat3(const mat3& m) {
    mat4 r(1);

    for (int c = 0; c < 3; ++c) {
        for (int row = 0; row < 3; ++row) r.at(row, c) = m.at(row, c);
    }

    return r;
}

inline mat4 scale(const mat4& m, const vec3& s) {
    mat4 r = m;

    for (int row = 0; row < 4; ++row) {
        r.at(row, 0) *= s.x;
        r.at(row, 1) *= s.y;
        r.at(row, 2) *= s.z;
    }

    return r;
}

inline mat4 translate(const mat4& m, const vec3& v) {
    mat4 r = m;

    for (int row = 0; row < 4; ++row) {
        r.at(row, 3) += m.at(row, 0) * v.x + m.at(row, 1) * v.y + m.at(row, 2) * v.z;
    }

    return r;
}

inline mat4 rotate(const mat4& m, float angle, const vec3& axis) {
    vec3 a = normalize(axis);

    float c = std::cos(angle);
    float s = std::sin(angle);

    float t = 1 - c;

    mat4 rot(1);

    rot.at(0, 0) = t * a.x * a.x + c;
    rot.at(1, 1) = t * a.y * a.y + c;
    rot.at(2, 2) = t * a.z * a.z + c;

    rot.at(1, 0) = t * a.x * a.y + s * a.z;
    rot.at(2, 0) = t * a.x * a.z - s * a.y;

    rot.at(0, 1) = t * a.x * a.y - s * a.z;
    rot.at(2, 1) = t * a.y * a.z + s * a.x;

    rot.at(0, 2) = t * a.x * a.z + s * a.y;
    rot.at(1, 2) = t * a.y * a.z - s * a.x;

    return m * rot;
}

inline mat4 perspective(float fovy_rad, float aspect, float zNear, float zFar) {
    float tan_half = std::tan(fovy_rad * 0.5);

    mat4 r;

    r.at(0, 0) = 1 / (aspect * tan_half); r.at(1, 1) = 1 / tan_half; r.at(2, 2) = -(zFar + zNear) / (zFar - zNear);

    r.at(3, 2) = -1; r.at(2, 3) = -(2 * zFar * zNear) / (zFar - zNear);

    return r;
}

inline mat4 lookAt(const vec3& eye, const vec3& center, const vec3& up) {
    vec3 f = normalize(center - eye);
    vec3 s = normalize(cross(f, up));

    vec3 u = cross(s, f);

    mat4 r(1.0);

    r.at(0, 0) = s.x; r.at(0, 1) = s.y; r.at(0, 2) = s.z;
    r.at(1, 0) = u.x; r.at(1, 1) = u.y; r.at(1, 2) = u.z;

    r.at(2, 0) = -f.x; r.at(2, 1) = -f.y; r.at(2, 2) = -f.z;

    r.at(0, 3) = -dot(s, eye);
    r.at(1, 3) = -dot(u, eye);

    r.at(2, 3) = dot(f, eye);

    return r;
}
