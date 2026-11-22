// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Packing of a network built by Gusev's procedure, into the end-linking system of src/common
// (poly/endlink/system.hpp) with positions in continuous space. The cross-linkers are placed at
// random in the box before the chains are linked (sim/linking.hpp); the beads of the chains are
// placed afterwards, as Gaussian random walks with steps of mean-square length b^2 (each
// component normal with variance b^2/3), without excluded volume. A chain with both ends linked
// is a Brownian bridge between its two cross-linkers (Wang, Ramkrishna and Narsimhan, J. Chem.
// Phys. 2020, 153, 034901), sampled exactly step by step: with m steps left from x to the target
// y, the next point is normal with mean x + (y - x)/m and variance (b^2/3)(m - 1)/m per component.
// A dangling chain is a free walk from its cross-linker, and a free chain one from a random point
// in the box.
#pragma once

#include <cstddef>
#include <span>

#include "poly/endlink/system.hpp"
#include "poly/random.hpp"

namespace gusev {

using poly::geometry::Vec3d;
using poly::geometry::Vec3i;
/// The network: chains and cross-linkers in continuous space, the chains unwrapped.
using System = poly::endlink::System<Vec3d>;
using poly::endlink::AtomType;
using poly::endlink::Link;

/// Unlinked chains of `lengths` beads and `n_crosslinkers` cross-linkers with `functionality`
/// sites, placed uniformly at random in a box of edge `box`; the beads are not placed yet.
[[nodiscard]] System place_crosslinkers(std::span<const std::size_t> lengths, std::size_t n_crosslinkers, int functionality,
                                        double box, poly::Rng& rng);

/// Write to `out` the out.size() interior points of a Gaussian random walk of out.size() + 1
/// steps from `from` to `to`.
void brownian_bridge(const Vec3d& from, const Vec3d& to, double bond_length, poly::Rng& rng, std::span<Vec3d> out);

/// Place the beads of every chain of `system`, unwrapped, chain after chain.
void place_beads(System& system, double bond_length, poly::Rng& rng);

}  // namespace gusev
