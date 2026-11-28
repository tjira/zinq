#pragma once

#include <algorithm>
#include <cmath>

/// 3D Cartesian floating-point vector for spatial coordinates, normals, and colors.
struct vec3 {
    float x = 0, y = 0, z = 0;

    /// Constructs a 3D Cartesian vector initialized to zero, uniform scalar, or specified spatial components.
    vec3() = default; vec3(float v) : x(v), y(v), z(v) {} vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    /// Compares two 3D vectors for equality of all spatial components.
    bool operator!=(const vec3& b) const { return x != b.x || y != b.y || z != b.z; }

    /// Computes the negation of the 3D Cartesian vector.
    vec3 operator-() const { return vec3(-x, -y, -z); }

    /// Computes the vector addition of two 3D vectors.
    vec3 operator+(const vec3& b) const { return vec3(x + b.x, y + b.y, z + b.z); }

    /// Computes the vector subtraction of two 3D vectors.
    vec3 operator-(const vec3& b) const { return vec3(x - b.x, y - b.y, z - b.z); }

    /// Multiplies the 3D vector by a scalar factor.
    vec3 operator*(float s) const { return vec3(x * s, y * s, z * s); }

    /// Divides the 3D vector by a scalar factor.
    vec3 operator/(float s) const { return vec3(x / s, y / s, z / s); }

    /// Adds another 3D vector in-place to this vector.
    vec3& operator+=(const vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
};

/// Scales a 3D vector by a scalar factor.
inline vec3 operator*(float s, const vec3& v) { return vec3(v.x * s, v.y * s, v.z * s); }

/// Computes the standard Euclidean dot product between two 3D vectors.
inline float dot(const vec3& a, const vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/// Computes the Euclidean norm or length of the 3D vector.
inline float length(const vec3& v) { return std::sqrt(dot(v, v)); }

/// Normalizes the 3D vector to unit Euclidean length.
inline vec3 normalize(const vec3& v) {
    float len = length(v); return len > 1e-8f ? v / len : vec3(0.0f);
}

/// Computes the right-handed Cartesian cross product of two 3D vectors.
inline vec3 cross(const vec3& a, const vec3& b) {
    return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/// 4x4 column-major floating-point matrix for 3D homogeneous affine transformations and projections.
struct mat4 {
    float data[16] = { 0 };

    /// Constructs a 4x4 column-major matrix initialized to zero or scalar diagonal values.
    mat4() = default;

    /// Constructs a 4x4 diagonal matrix with uniform diagonal elements.
    mat4(float diag) {
        data[0] = diag; data[5] = diag; data[10] = diag; data[15] = diag;
    }

    /// Returns a reference to the element at the specified row and column in the 4x4 matrix.
    float& at(int r, int c) { return data[c * 4 + r]; }

    /// Returns the element value at the specified row and column in the 4x4 matrix.
    float at(int r, int c) const { return data[c * 4 + r]; }

    /// Returns a raw pointer to the underlying column-major 4x4 floating-point buffer.
    const float* ptr() const { return data; }

    /// Extracts the first three spatial components of the specified column.
    vec3 col3(int c) const {
        return vec3(at(0, c), at(1, c), at(2, c));
    }

    /// Computes matrix multiplication of two 4x4 transformation matrices.
    mat4 operator*(const mat4& b) const {
        mat4 r;

        for (int c = 0; c < 4; ++c) {
            for (int r_idx = 0; r_idx < 4; ++r_idx) {
                float sum = 0.0f;

                for (int k = 0; k < 4; ++k) sum += at(r_idx, k) * b.at(k, c);

                r.at(r_idx, c) = sum;
            }
        }

        return r;
    }
};

/// Applies non-uniform spatial scaling along coordinate axes to a transformation matrix.
inline mat4 scale(const mat4& m, const vec3& s) {
    mat4 r = m;

    for (int row = 0; row < 4; ++row) {
        r.at(row, 0) *= s.x;
        r.at(row, 1) *= s.y;
        r.at(row, 2) *= s.z;
    }

    return r;
}

/// Applies spatial translation by a Cartesian vector to a transformation matrix.
inline mat4 translate(const mat4& m, const vec3& v) {
    mat4 r = m;

    for (int row = 0; row < 4; ++row) {
        r.at(row, 3) += m.at(row, 0) * v.x + m.at(row, 1) * v.y + m.at(row, 2) * v.z;
    }

    return r;
}

/// Applies arbitrary-axis spatial rotation by an angle in radians using Rodrigues rotation.
inline mat4 rotate(const mat4& m, float angle, const vec3& axis) {
    vec3 a = normalize(axis);

    float c = std::cos(angle);
    float s = std::sin(angle);

    float t = 1.0f - c;

    mat4 rot(1.0f);

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

/// Constructs a standard perspective projection matrix from field-of-view and clip planes.
inline mat4 perspective(float fovy_rad, float aspect, float zNear, float zFar) {
    float tan_half = std::tan(fovy_rad * 0.5f);

    mat4 r;

    r.at(0, 0) = 1.0f / (aspect * tan_half); r.at(1, 1) = 1.0f / tan_half; r.at(2, 2) = -(zFar + zNear) / (zFar - zNear);

    r.at(3, 2) = -1.0f; r.at(2, 3) = -(2.0f * zFar * zNear) / (zFar - zNear);

    return r;
}

/// Constructs an orthographic projection matrix defining a parallel-projection viewing box.
inline mat4 ortho(float left, float right, float bottom, float top, float zNear, float zFar) {
    mat4 r;

    r.at(0, 0) =  2.0f / (right - left);
    r.at(1, 1) =  2.0f / (top - bottom);
    r.at(2, 2) = -2.0f / (zFar - zNear);

    r.at(0, 3) = -(right + left) / (right - left);
    r.at(1, 3) = -(top + bottom) / (top - bottom);
    r.at(2, 3) = -(zFar + zNear) / (zFar - zNear);

    r.at(3, 3) =  1.0f;

    return r;
}

/// Constructs a camera view transformation matrix looking from eye toward center.
inline mat4 lookAt(const vec3& eye, const vec3& center, const vec3& up) {
    vec3 f = normalize(center - eye);
    vec3 s = normalize(cross(f, up));

    vec3 u = cross(s, f);

    mat4 r(1.0f);

    r.at(0, 0) = s.x; r.at(0, 1) = s.y; r.at(0, 2) = s.z;
    r.at(1, 0) = u.x; r.at(1, 1) = u.y; r.at(1, 2) = u.z;

    r.at(2, 0) = -f.x; r.at(2, 1) = -f.y; r.at(2, 2) = -f.z;

    r.at(0, 3) = -dot(s, eye);
    r.at(1, 3) = -dot(u, eye);

    r.at(2, 3) = dot(f, eye);

    return r;
}

/// 4D Cartesian floating-point vector for homogeneous spatial coordinates and projections.
struct vec4 {
    float x = 0, y = 0, z = 0, w = 0;

    /// Constructs a 4D Cartesian vector initialized to zero, uniform scalar, or specified spatial components.
    vec4() = default; vec4(float v) : x(v), y(v), z(v), w(v) {}

    /// Constructs a 4D Cartesian vector from specified spatial components and homogeneous coordinate.
    vec4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}

    /// Constructs a 4D Cartesian vector from a 3D vector and homogeneous coordinate.
    vec4(const vec3& v, float w_) : x(v.x), y(v.y), z(v.z), w(w_) {}

    /// Extracts the 3D spatial components of the 4D homogeneous vector.
    vec3 xyz() const { return vec3(x, y, z); }
};

/// Multiplies a 4x4 matrix by a 4D homogeneous vector.
inline vec4 operator*(const mat4& m, const vec4& v) {
    return vec4(
        m.at(0, 0) * v.x + m.at(0, 1) * v.y + m.at(0, 2) * v.z + m.at(0, 3) * v.w,
        m.at(1, 0) * v.x + m.at(1, 1) * v.y + m.at(1, 2) * v.z + m.at(1, 3) * v.w,
        m.at(2, 0) * v.x + m.at(2, 1) * v.y + m.at(2, 2) * v.z + m.at(2, 3) * v.w,
        m.at(3, 0) * v.x + m.at(3, 1) * v.y + m.at(3, 2) * v.z + m.at(3, 3) * v.w
    );
}

/// Computes the inverse of a 4x4 transformation matrix using cofactor expansion.
inline mat4 inverse(const mat4& m) {
    float a0 = m.at(0, 0) * m.at(1, 1) - m.at(0, 1) * m.at(1, 0);
    float a1 = m.at(0, 0) * m.at(1, 2) - m.at(0, 2) * m.at(1, 0);
    float a2 = m.at(0, 0) * m.at(1, 3) - m.at(0, 3) * m.at(1, 0);
    float a3 = m.at(0, 1) * m.at(1, 2) - m.at(0, 2) * m.at(1, 1);
    float a4 = m.at(0, 1) * m.at(1, 3) - m.at(0, 3) * m.at(1, 1);
    float a5 = m.at(0, 2) * m.at(1, 3) - m.at(0, 3) * m.at(1, 2);
    float b0 = m.at(2, 0) * m.at(3, 1) - m.at(2, 1) * m.at(3, 0);
    float b1 = m.at(2, 0) * m.at(3, 2) - m.at(2, 2) * m.at(3, 0);
    float b2 = m.at(2, 0) * m.at(3, 3) - m.at(2, 3) * m.at(3, 0);
    float b3 = m.at(2, 1) * m.at(3, 2) - m.at(2, 2) * m.at(3, 1);
    float b4 = m.at(2, 1) * m.at(3, 3) - m.at(2, 3) * m.at(3, 1);
    float b5 = m.at(2, 2) * m.at(3, 3) - m.at(2, 3) * m.at(3, 2);

    float det = a0 * b5 - a1 * b4 + a2 * b3 + a3 * b2 - a4 * b1 + a5 * b0;

    if (std::abs(det) < 1e-12f) return mat4(1.0f);

    float inv_det = 1.0f / det;

    mat4 r;

    r.at(0, 0) = (+m.at(1, 1) * b5 - m.at(1, 2) * b4 + m.at(1, 3) * b3) * inv_det;
    r.at(0, 1) = (-m.at(0, 1) * b5 + m.at(0, 2) * b4 - m.at(0, 3) * b3) * inv_det;
    r.at(0, 2) = (+m.at(3, 1) * a5 - m.at(3, 2) * a4 + m.at(3, 3) * a3) * inv_det;
    r.at(0, 3) = (-m.at(2, 1) * a5 + m.at(2, 2) * a4 - m.at(2, 3) * a3) * inv_det;

    r.at(1, 0) = (-m.at(1, 0) * b5 + m.at(1, 2) * b2 - m.at(1, 3) * b1) * inv_det;
    r.at(1, 1) = (+m.at(0, 0) * b5 - m.at(0, 2) * b2 + m.at(0, 3) * b1) * inv_det;
    r.at(1, 2) = (-m.at(3, 0) * a5 + m.at(3, 2) * a2 - m.at(3, 3) * a1) * inv_det;
    r.at(1, 3) = (+m.at(2, 0) * a5 - m.at(2, 2) * a2 + m.at(2, 3) * a1) * inv_det;

    r.at(2, 0) = (+m.at(1, 0) * b4 - m.at(1, 1) * b2 + m.at(1, 3) * b0) * inv_det;
    r.at(2, 1) = (-m.at(0, 0) * b4 + m.at(0, 1) * b2 - m.at(0, 3) * b0) * inv_det;
    r.at(2, 2) = (+m.at(3, 0) * a4 - m.at(3, 1) * a2 + m.at(3, 3) * a0) * inv_det;
    r.at(2, 3) = (-m.at(2, 0) * a4 + m.at(2, 1) * a2 - m.at(2, 3) * a0) * inv_det;

    r.at(3, 0) = (-m.at(1, 0) * b3 + m.at(1, 1) * b1 - m.at(1, 2) * b0) * inv_det;
    r.at(3, 1) = (+m.at(0, 0) * b3 - m.at(0, 1) * b1 + m.at(0, 2) * b0) * inv_det;
    r.at(3, 2) = (-m.at(3, 0) * a3 + m.at(3, 1) * a1 - m.at(3, 2) * a0) * inv_det;
    r.at(3, 3) = (+m.at(2, 0) * a3 - m.at(2, 1) * a1 + m.at(2, 2) * a0) * inv_det;

    return r;
}

/// Unit quaternion representing 3D spatial orientation and rotation in SO(3).
struct quat {
    float w = 1.0f, x = 0.0f, y = 0.0f, z = 0.0f;

    /// Constructs a unit quaternion initialized to identity or specified scalar components.
    quat() = default; quat(float w_, float x_, float y_, float z_) : w(w_), x(x_), y(y_), z(z_) {}

    /// Evaluates the Hamilton product of two quaternions representing compound spatial rotation.
    quat operator*(const quat& b) const {
        return quat(
            w * b.w - x * b.x - y * b.y - z * b.z,
            w * b.x + x * b.w + y * b.z - z * b.y,
            w * b.y - x * b.z + y * b.w + z * b.x,
            w * b.z + x * b.y - y * b.x + z * b.w
        );
    }

    /// Evaluates the quaternion conjugate representing the inverse spatial rotation.
    quat conjugate() const { return quat(w, -x, -y, -z); }

    /// Evaluates the Euclidean norm of the four-dimensional quaternion.
    float length() const { return std::sqrt(w * w + x * x + y * y + z * z); }

    /// Applies the 3D spatial rotation represented by the unit quaternion to a Cartesian vector.
    vec3 rotate(const vec3& v) const {
        vec3 qv(x, y, z), t = 2.0f * cross(qv, v);

        return v + w * t + cross(qv, t);
    }

    /// Constructs a unit quaternion representing spatial rotation around an arbitrary 3D Cartesian axis.
    static quat fromAxisAngle(const vec3& axis, float angle_rad) {
        vec3 a = normalize(axis);

        float half = angle_rad * 0.5f;

        float s = std::sin(half);

        return quat(std::cos(half), a.x * s, a.y * s, a.z * s);
    }
};

/// Normalizes a quaternion to unit length to preserve rigid-body rotational isometry.
inline quat normalize(const quat& q) {
    float len = q.length();

    return len > 1e-8f ? quat(q.w / len, q.x / len, q.y / len, q.z / len) : quat(1.0f, 0.0f, 0.0f, 0.0f);
}

/// Semi-infinite 3D ray parameterized by origin and unit direction.
struct Ray {
    vec3 origin, direction;

    /// Tests ray intersection with a sphere and returns the nearest positive distance.
    bool intersectSphere(const vec3& center, float radius, float& t_out) const {
        vec3 m = origin - center;

        float b = dot(m, direction); float c = dot(m, m) - radius * radius;

        if (c > 0.0f && b > 0.0f) return false;

        float disc = b * b - c;

        if (disc < 0.0f) return false;

        float t = -b - std::sqrt(disc);

        if (t < 0.0f) t = -b + std::sqrt(disc);

        if (t < 0.0f) return false;

        t_out = t;

        return true;
    }

    /// Tests ray intersection with a plane and outputs the 3D intersection point.
    bool intersectPlane(const vec3& plane_point, const vec3& plane_normal, vec3& hit_out) const {
        float denom = dot(direction, plane_normal);

        if (std::abs(denom) < 1e-6f) return false;

        float t = dot(plane_point - origin, plane_normal) / denom;

        if (t < 0.0f) return false;

        hit_out = origin + direction * t;

        return true;
    }
};

/// Computes the bond angle in radians formed by three spatial atom positions.
inline float computeAngle(const vec3& p1, const vec3& p2, const vec3& p3) {
    vec3 v1 = normalize(p1 - p2);
    vec3 v2 = normalize(p3 - p2);

    float cos_theta = std::clamp(dot(v1, v2), -1.0f, 1.0f);

    return std::acos(cos_theta);
}

/// Computes the dihedral angle in radians formed by four sequential spatial atom positions.
inline float computeDihedral(const vec3& p1, const vec3& p2, const vec3& p3, const vec3& p4) {
    vec3 b1 = p2 - p1;
    vec3 b2 = p3 - p2;
    vec3 b3 = p4 - p3;

    vec3 n1 = cross(b1, b2);
    vec3 n2 = cross(b2, b3);

    vec3 m1 = cross(n1, normalize(b2));

    return std::atan2(dot(m1, n2), dot(n1, n2));
}
