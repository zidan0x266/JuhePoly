// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Small vectors: integer ones for lattice positions, steps and bond vectors, and real ones for
// positions in continuous space.
#pragma once

namespace poly::geometry {

struct Vec3i {
  int x{0};
  int y{0};
  int z{0};

  constexpr bool operator==(const Vec3i&) const noexcept = default;
  constexpr Vec3i operator+(Vec3i o) const noexcept { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3i operator-(Vec3i o) const noexcept { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3i operator-() const noexcept { return {-x, -y, -z}; }
  constexpr Vec3i& operator+=(Vec3i o) noexcept {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  [[nodiscard]] constexpr int dot(Vec3i o) const noexcept { return x * o.x + y * o.y + z * o.z; }
  [[nodiscard]] constexpr int norm2() const noexcept { return dot(*this); }
  [[nodiscard]] constexpr int operator[](int axis) const noexcept { return axis == 0 ? x : (axis == 1 ? y : z); }
};

struct Vec3d {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  constexpr bool operator==(const Vec3d&) const noexcept = default;
  constexpr Vec3d operator+(Vec3d o) const noexcept { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3d operator-(Vec3d o) const noexcept { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3d operator-() const noexcept { return {-x, -y, -z}; }
  constexpr Vec3d operator*(double s) const noexcept { return {x * s, y * s, z * s}; }
  constexpr Vec3d& operator+=(Vec3d o) noexcept {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  [[nodiscard]] constexpr double dot(Vec3d o) const noexcept { return x * o.x + y * o.y + z * o.z; }
  [[nodiscard]] constexpr double norm2() const noexcept { return dot(*this); }
  [[nodiscard]] constexpr double operator[](int axis) const noexcept { return axis == 0 ? x : (axis == 1 ? y : z); }
};

}  // namespace poly::geometry
