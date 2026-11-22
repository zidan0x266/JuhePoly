// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "gusev/sim/packing.hpp"

#include <cmath>

namespace gusev {

namespace {

[[nodiscard]] Vec3d uniform_point(double box, poly::Rng& rng) noexcept {
  const double x = box * rng.uniform01();
  const double y = box * rng.uniform01();
  return {x, y, box * rng.uniform01()};
}

[[nodiscard]] Vec3d normal_step(double sigma, poly::Rng& rng) noexcept {
  const double x = sigma * rng.normal();
  const double y = sigma * rng.normal();
  return {x, y, sigma * rng.normal()};
}

}  // namespace

System place_crosslinkers(std::span<const std::size_t> lengths, std::size_t n_crosslinkers, int functionality, double box,
                          poly::Rng& rng) {
  System system = poly::endlink::make_system<Vec3d>(lengths, n_crosslinkers, functionality, box, 1.0);
  for (Vec3d& p : system.crosslinkers) p = uniform_point(box, rng);
  return system;
}

void brownian_bridge(const Vec3d& from, const Vec3d& to, double bond_length, poly::Rng& rng, std::span<Vec3d> out) {
  const double sigma = bond_length / std::sqrt(3.0);
  const std::size_t steps = out.size() + 1;
  Vec3d x = from;
  for (std::size_t k = 0; k < out.size(); ++k) {
    const auto m = static_cast<double>(steps - k);  // steps left from x
    x = x + (to - x) * (1.0 / m) + normal_step(sigma * std::sqrt((m - 1) / m), rng);
    out[k] = x;
  }
}

void place_beads(System& system, double bond_length, poly::Rng& rng) {
  const double sigma = bond_length / std::sqrt(3.0);
  system.beads.assign(system.chain_beads(), Vec3d{});
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    const std::span<Vec3d> beads(system.beads.data() + system.chain_start[c], system.chain_length(c));
    const Link& end0 = system.link(c, 0);
    const Link& end1 = system.link(c, 1);
    if (end0.linked() && end1.linked()) {
      const Vec3d from = system.image_of(end0.crosslinker, end0.image);
      brownian_bridge(from, from + system.span(c), bond_length, rng, beads);
    } else if (end0.linked()) {
      Vec3d x = system.image_of(end0.crosslinker, end0.image);
      for (Vec3d& bead : beads) bead = x = x + normal_step(sigma, rng);
    } else if (end1.linked()) {
      Vec3d x = system.image_of(end1.crosslinker, end1.image);
      for (std::size_t i = beads.size(); i-- > 0;) beads[i] = x = x + normal_step(sigma, rng);
    } else {
      Vec3d x = uniform_point(system.period, rng);
      beads[0] = x;
      for (std::size_t i = 1; i < beads.size(); ++i) beads[i] = x = x + normal_step(sigma, rng);
    }
  }
}

}  // namespace gusev
