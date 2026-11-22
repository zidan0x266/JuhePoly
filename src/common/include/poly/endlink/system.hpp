// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// An end-linking system: precursor chains with a reactive bead at each end and f-functional
// cross-linkers in a periodic cubic box, with the links formed so far between chain ends and
// cross-linker sites. The packing modules build it (randomwalk packs unlinked chains, gusev links
// them by Gusev's procedure), and the writers of src/common turn it into LAMMPS data files.
//
// Positions are of type Point: lattice points (Vec3i) or points in continuous space (Vec3d), in
// units of `scale`. Chain beads are unwrapped, so every chain is continuous; cross-linkers lie in
// the box. A linked end records the periodic image of its cross-linker that it is bonded to, as
// seen from the unwrapped beads of its chain.
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "poly/geometry/vec3.hpp"

namespace poly::endlink {

using geometry::Vec3d;
using geometry::Vec3i;

/// Cross-linker of a chain end that is not linked.
inline constexpr std::uint32_t kFree = std::numeric_limits<std::uint32_t>::max();

/// Atom types of the data files, labeled as the end-linking inputs of LAMMPS expect them: chain
/// beads (linked ends included), unlinked chain ends, cross-linkers, and cross-linkers with all
/// sites linked.
enum class AtomType : int { Bead = 1, End = 2, Crosslinker = 3, FullCrosslinker = 4 };
inline constexpr std::array<std::string_view, 4> kTypeLabels{"bead", "end", "xlink", "xlink_full"};

struct Link {
  std::uint32_t crosslinker{kFree};
  Vec3i image{};  ///< periodic image of the cross-linker, in box edges, seen from the chain
  [[nodiscard]] bool linked() const noexcept { return crosslinker != kFree; }
};

template <class Point>
struct System {
  using Coordinate = std::remove_cvref_t<decltype(Point::x)>;

  Coordinate period{};  ///< box edge, in coordinate units
  double scale{1.0};    ///< length of one coordinate unit
  int functionality{4};
  /// The first bead of every chain, then the number of chain beads: chain c holds the beads
  /// chain_start[c] to chain_start[c + 1] - 1, from end 0 to end 1.
  std::vector<std::size_t> chain_start{0};
  std::vector<Point> beads;  ///< unwrapped; empty until placed
  std::vector<Point> crosslinkers;
  std::vector<Link> links;                ///< end e of chain c at 2c + e
  std::vector<std::uint8_t> used_sites;  ///< per cross-linker

  // ---- chains and cross-linkers ----------------------------------------------------------------
  [[nodiscard]] std::size_t n_chains() const noexcept { return chain_start.size() - 1; }
  [[nodiscard]] std::size_t chain_length(std::size_t c) const noexcept { return chain_start[c + 1] - chain_start[c]; }
  [[nodiscard]] std::size_t chain_beads() const noexcept { return chain_start.back(); }
  [[nodiscard]] std::size_t n_crosslinkers() const noexcept { return crosslinkers.size(); }
  [[nodiscard]] bool has_beads() const noexcept { return beads.size() == chain_beads(); }
  [[nodiscard]] double box_length() const noexcept { return scale * static_cast<double>(period); }

  [[nodiscard]] const Link& link(std::size_t c, int end) const noexcept { return links[2 * c + static_cast<std::size_t>(end)]; }
  [[nodiscard]] Link& link(std::size_t c, int end) noexcept { return links[2 * c + static_cast<std::size_t>(end)]; }
  [[nodiscard]] bool linked(std::size_t c, int end) const noexcept { return link(c, end).linked(); }
  [[nodiscard]] int linked_ends(std::size_t c) const noexcept { return int{linked(c, 0)} + int{linked(c, 1)}; }
  [[nodiscard]] int free_sites(std::uint32_t k) const noexcept { return functionality - used_sites[k]; }
  /// Reactive groups of the minority component: min(2 chains, f cross-linkers).
  [[nodiscard]] std::size_t minority_groups() const noexcept {
    const std::size_t sites = static_cast<std::size_t>(functionality) * n_crosslinkers();
    return 2 * n_chains() < sites ? 2 * n_chains() : sites;
  }
  /// Link end `end` of chain `c` to cross-linker `k`, through its periodic image `image`.
  void attach(std::size_t c, int end, std::uint32_t k, Vec3i image = {}) noexcept {
    link(c, end) = {k, image};
    ++used_sites[k];
  }

  /// Position of cross-linker `k` shifted by `image` box edges.
  [[nodiscard]] Point image_of(std::uint32_t k, Vec3i image) const noexcept {
    const Point& p = crosslinkers[k];
    return {p.x + static_cast<Coordinate>(image.x) * period, p.y + static_cast<Coordinate>(image.y) * period,
            p.z + static_cast<Coordinate>(image.z) * period};
  }
  /// From the cross-linker of end 0 to that of end 1 of a chain with both ends linked, through
  /// the images they are bonded to.
  [[nodiscard]] Point span(std::size_t c) const noexcept {
    return image_of(link(c, 1).crosslinker, link(c, 1).image) - image_of(link(c, 0).crosslinker, link(c, 0).image);
  }

  // ---- atoms: the chain beads, chain after chain, then the cross-linkers -----------------------
  [[nodiscard]] std::size_t atoms() const noexcept { return chain_beads() + n_crosslinkers(); }
  /// Bonds along the chains and from every linked end to its cross-linker.
  [[nodiscard]] std::size_t bonds() const noexcept {
    std::size_t n = chain_beads() - n_chains();
    for (const Link& l : links) n += l.linked() ? 1 : 0;
    return n;
  }
  [[nodiscard]] bool is_chain_bead(std::size_t i) const noexcept { return i < chain_beads(); }
  /// The chain that the chain bead `i` belongs to.
  [[nodiscard]] std::size_t chain_of(std::size_t i) const noexcept {
    return static_cast<std::size_t>(std::ranges::upper_bound(chain_start, i) - chain_start.begin()) - 1;
  }
  /// One-based molecule ID: chain c is c + 1, cross-linker k is n_chains + k + 1.
  [[nodiscard]] std::size_t molecule(std::size_t i) const noexcept {
    return is_chain_bead(i) ? chain_of(i) + 1 : i - chain_beads() + n_chains() + 1;
  }
  [[nodiscard]] AtomType crosslinker_type(std::size_t k) const noexcept {
    return used_sites[k] == functionality ? AtomType::FullCrosslinker : AtomType::Crosslinker;
  }
  /// Type of the bead `i` of chain `c`.
  [[nodiscard]] AtomType bead_type(std::size_t c, std::size_t i) const noexcept {
    const bool open_end = (i == chain_start[c] && !linked(c, 0)) || (i + 1 == chain_start[c + 1] && !linked(c, 1));
    return open_end ? AtomType::End : AtomType::Bead;
  }
  [[nodiscard]] AtomType type(std::size_t i) const noexcept {
    return is_chain_bead(i) ? bead_type(chain_of(i), i) : crosslinker_type(i - chain_beads());
  }
  [[nodiscard]] const Point& position(std::size_t i) const noexcept {
    return is_chain_bead(i) ? beads[i] : crosslinkers[i - chain_beads()];
  }
};

/// A system of unlinked chains of `lengths` beads and `n_crosslinkers` cross-linkers, with the
/// cross-linkers at the origin and no beads placed yet.
template <class Point>
[[nodiscard]] System<Point> make_system(std::span<const std::size_t> lengths, std::size_t n_crosslinkers, int functionality,
                                        typename System<Point>::Coordinate period, double scale);

struct Statistics {
  std::size_t free{0};           ///< chains with no end linked
  std::size_t dangling{0};       ///< chains with one end linked
  std::size_t connected{0};      ///< chains with both ends linked
  std::size_t primary_loops{0};  ///< both ends on one cross-linker, through the same image
  std::size_t image_loops{0};    ///< both ends on one cross-linker, through different images
  std::size_t linked_ends{0};
  std::size_t full_crosslinkers{0};
  double conversion{0.0};  ///< linked ends over the minority groups
  /// Mean of R^2 / (N b^2) over the connected chains, R the span and N the bonds from one
  /// cross-linker to the other (chain length + 1).
  double span_ratio{0.0};
};

template <class Point>
[[nodiscard]] Statistics statistics(const System<Point>& system, double bond_length);

/// Check the layout of the chains, the links and the sites, the cross-linker positions and, if
/// placed, the number of beads. Returns the first violation.
template <class Point>
[[nodiscard]] std::expected<void, std::string> validate(const System<Point>& system);

}  // namespace poly::endlink
