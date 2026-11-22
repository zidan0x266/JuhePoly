// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include <format>
#include <utility>

#include "poly/endlink/system.hpp"

namespace poly::endlink {

template <class Point>
System<Point> make_system(std::span<const std::size_t> lengths, std::size_t n_crosslinkers, int functionality,
                          typename System<Point>::Coordinate period, double scale) {
  System<Point> s;
  s.period = period;
  s.scale = scale;
  s.functionality = functionality;
  s.chain_start.reserve(lengths.size() + 1);
  for (const std::size_t n : lengths) s.chain_start.push_back(s.chain_start.back() + n);
  s.crosslinkers.assign(n_crosslinkers, Point{});
  s.links.assign(2 * lengths.size(), Link{});
  s.used_sites.assign(n_crosslinkers, 0);
  return s;
}

template <class Point>
Statistics statistics(const System<Point>& system, double bond_length) {
  Statistics s;
  double ratio = 0.0;
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    const int ends = system.linked_ends(c);
    s.linked_ends += static_cast<std::size_t>(ends);
    if (ends == 0) {
      ++s.free;
    } else if (ends == 1) {
      ++s.dangling;
    } else {
      ++s.connected;
      const Link& a = system.link(c, 0);
      const Link& b = system.link(c, 1);
      if (a.crosslinker == b.crosslinker) ++(a.image == b.image ? s.primary_loops : s.image_loops);
      const Point r = system.span(c);
      const double r2 = system.scale * system.scale * static_cast<double>(r.norm2());
      ratio += r2 / (static_cast<double>(system.chain_length(c) + 1) * bond_length * bond_length);
    }
  }
  for (std::size_t k = 0; k < system.n_crosslinkers(); ++k) s.full_crosslinkers += system.used_sites[k] == system.functionality ? 1 : 0;
  const std::size_t minority = system.minority_groups();
  s.conversion = minority ? static_cast<double>(s.linked_ends) / static_cast<double>(minority) : 0.0;
  s.span_ratio = s.connected ? ratio / static_cast<double>(s.connected) : 0.0;
  return s;
}

template <class Point>
std::expected<void, std::string> validate(const System<Point>& system) {
  const auto fail = [](std::string message) { return std::unexpected(std::move(message)); };
  if (system.chain_start.empty() || system.chain_start.front() != 0) return fail("chain offsets must start at 0");
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    if (system.chain_start[c + 1] < system.chain_start[c] + 2) return fail(std::format("chain {} has fewer than two beads", c + 1));
  }
  if (system.links.size() != 2 * system.n_chains()) return fail("two links per chain are needed");
  if (system.used_sites.size() != system.n_crosslinkers()) return fail("one used-site count per cross-linker is needed");
  for (std::size_t k = 0; k < system.n_crosslinkers(); ++k) {
    const Point& p = system.crosslinkers[k];
    if (p.x < 0 || p.y < 0 || p.z < 0 || p.x >= system.period || p.y >= system.period || p.z >= system.period) {
      return fail(std::format("cross-linker {} lies outside the box", k + 1));
    }
  }
  std::vector<int> used(system.n_crosslinkers(), 0);
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    for (int end = 0; end < 2; ++end) {
      const Link& l = system.link(c, end);
      if (!l.linked()) {
        if (l.image != Vec3i{}) return fail(std::format("chain {} has an image at an unlinked end", c + 1));
        continue;
      }
      if (l.crosslinker >= system.n_crosslinkers()) return fail(std::format("chain {} is linked to a cross-linker that does not exist", c + 1));
      ++used[l.crosslinker];
    }
  }
  for (std::size_t k = 0; k < system.n_crosslinkers(); ++k) {
    if (used[k] != system.used_sites[k] || used[k] > system.functionality) {
      return fail(std::format("the sites of cross-linker {} do not match the chains linked to it", k + 1));
    }
  }
  if (!system.beads.empty() && !system.has_beads()) return fail("bead count does not match the chains");
  return {};
}

template System<Vec3i> make_system<Vec3i>(std::span<const std::size_t>, std::size_t, int, int, double);
template System<Vec3d> make_system<Vec3d>(std::span<const std::size_t>, std::size_t, int, double, double);
template Statistics statistics(const System<Vec3i>&, double);
template Statistics statistics(const System<Vec3d>&, double);
template std::expected<void, std::string> validate(const System<Vec3i>&);
template std::expected<void, std::string> validate(const System<Vec3d>&);

}  // namespace poly::endlink
