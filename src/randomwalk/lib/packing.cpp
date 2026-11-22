// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/sim/packing.hpp"

#include <algorithm>
#include <array>
#include <format>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace rw {

namespace {

using NeighborOrder = std::array<int, FccLattice::kCoordination>;

/// Order in which the original program visited the twelve neighbors of a head: by cubic cell
/// offset (x, then y, then z, each -1, 0, +1) and then by basis site (0,0,0), (1/2,1/2,0),
/// (1/2,0,1/2), (0,1/2,1/2). It depends only on the parities of the head's coordinates.
constexpr std::array<NeighborOrder, 4> make_legacy_orders() {
  std::array<NeighborOrder, 4> orders{};
  for (int parity = 0; parity < 4; ++parity) {
    const int px = parity & 1;
    const int py = parity >> 1;
    const int pz = (px + py) & 1;
    std::array<int, FccLattice::kCoordination> key{};
    for (int k = 0; k < FccLattice::kCoordination; ++k) {
      const Vec3i d = FccLattice::kNeighbors[k];
      const Vec3i c{px + d.x, py + d.y, pz + d.z};  // neighbor relative to the head's cell, doubled
      const int cell = (((c.x + 2) >> 1) * 3 + ((c.y + 2) >> 1)) * 3 + ((c.z + 2) >> 1);
      const int basis = (c.x & 1) ? ((c.y & 1) ? 1 : 2) : ((c.y & 1) ? 3 : 0);
      key[k] = cell * 4 + basis;
      orders[parity][k] = k;
    }
    std::ranges::sort(orders[parity], [&key](int a, int b) { return key[a] < key[b]; });
  }
  return orders;
}

constexpr std::array<NeighborOrder, 4> kLegacyOrders = make_legacy_orders();

constexpr NeighborOrder kPlainOrder{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

}  // namespace

const NeighborOrder& legacy_scan_order(Vec3i head) noexcept { return kLegacyOrders[(head.x & 1) | ((head.y & 1) << 1)]; }

std::string_view to_string(SiteWeighting w) noexcept {
  switch (w) {
    case SiteWeighting::Legacy: return "legacy";
    case SiteWeighting::FreeNeighbors: return "free_neighbors";
    case SiteWeighting::Uniform: return "uniform";
  }
  return "?";
}

std::optional<SiteWeighting> parse_site_weighting(std::string_view s) noexcept {
  if (s == "legacy") return SiteWeighting::Legacy;
  if (s == "free_neighbors") return SiteWeighting::FreeNeighbors;
  if (s == "uniform") return SiteWeighting::Uniform;
  return std::nullopt;
}

int cells_for_occupancy(std::size_t chain_beads, double occupancy) {
  if (!(occupancy > 0)) throw std::invalid_argument("occupancy must be positive");
  int n = 1;
  while (4.0 * n * n * n * occupancy <= static_cast<double>(chain_beads)) ++n;
  return n;
}

Packer::Packer(int cells, WalkRules rules, poly::Rng& rng) : lattice_(cells), rules_(rules), rng_(rng) {
  if (rules_.angles.empty()) throw std::invalid_argument("at least one bond angle must be allowed");
  if (rules_.max_attempts < 1) throw std::invalid_argument("max_attempts must be at least 1");
  if (rules_.weighting == SiteWeighting::Legacy) head_count_.assign(lattice_.num_sites(), 0);
}

void Packer::count_head(Site s, int delta) noexcept {
  const Vec3i p = lattice_.position(s);
  for (int k = 0; k < FccLattice::kCoordination; ++k) {
    std::uint8_t& count = head_count_[lattice_.neighbor(p, k)];
    count = static_cast<std::uint8_t>(count + delta);
  }
}

Site Packer::random_free_site() {
  // Rejection sampling; with the lattice at most half full it rarely takes more than two draws.
  constexpr int kMaxDraws = 1'000'000;
  const auto n = static_cast<std::uint32_t>(lattice_.num_sites());
  for (int draw = 0; draw < kMaxDraws; ++draw) {
    const Site s = rng_.below(n);
    if (!lattice_.occupied(s)) return s;
  }
  throw std::runtime_error("no empty lattice site found; the lattice is almost full");
}

bool Packer::grow_chain(std::span<Vec3i> beads) {
  struct Candidate {
    Vec3i step;
    Site site;
    std::uint32_t weight;
  };
  std::array<Candidate, FccLattice::kCoordination> candidates;

  path_.clear();
  path_.push_back(random_free_site());
  lattice_.occupy(path_.back());
  Vec3i head = lattice_.position(path_.back());  // wrapped
  beads[0] = head;
  Vec3i previous{};

  const bool legacy = rules_.weighting == SiteWeighting::Legacy;
  for (std::size_t v = 1; v < beads.size(); ++v) {
    if (legacy) count_head(path_.back(), +1);
    int n = 0;
    int n_free = 0;
    std::array<std::uint32_t, FccLattice::kCoordination> free_weights;
    const NeighborOrder& order = legacy ? legacy_scan_order(head) : kPlainOrder;
    for (const int k : order) {
      const Vec3i d = FccLattice::kNeighbors[k];
      const Site s = lattice_.site(lattice_.wrap(head + d));
      if (lattice_.occupied(s)) continue;
      std::uint32_t weight = 1;
      if (legacy) {
        weight = static_cast<std::uint32_t>(FccLattice::kCoordination - head_count_[s]);
      } else if (rules_.weighting == SiteWeighting::FreeNeighbors) {
        weight = static_cast<std::uint32_t>(FccLattice::kCoordination - lattice_.occupied_neighbors(s));
      }
      free_weights[n_free++] = weight;
      if (v > 1 && !rules_.angles.allows_steps(previous, d)) continue;
      candidates[n++] = {d, s, weight};
    }
    if (n == 0) {  // trapped: remove the partial chain
      for (const Site s : path_) {
        if (legacy) count_head(s, -1);
        lattice_.vacate(s);
      }
      return false;
    }
    std::uint32_t total = 0;
    for (int i = 0; i < n; ++i) {
      if (legacy) candidates[i].weight = free_weights[i];  // the original program's index mix-up
      total += candidates[i].weight;
    }

    int pick = 0;
    if (total == 0) {  // every candidate is enclosed; choose among them uniformly
      pick = static_cast<int>(rng_.below(static_cast<std::uint32_t>(n)));
    } else {
      std::uint32_t r = rng_.below(total);
      while (r >= candidates[pick].weight) r -= candidates[pick++].weight;
    }
    const Candidate& next = candidates[pick];
    lattice_.occupy(next.site);
    path_.push_back(next.site);
    beads[v] = beads[v - 1] + next.step;
    head = lattice_.wrap(head + next.step);
    previous = next.step;
  }
  return true;
}

System Packer::pack(std::size_t n_chains, std::size_t chain_length, std::size_t n_crosslinkers, double bond_length) {
  const std::vector<std::size_t> lengths(n_chains, chain_length);
  return pack(lengths, n_crosslinkers, bond_length);
}

System Packer::pack(std::span<const std::size_t> lengths, std::size_t n_crosslinkers, double bond_length, int functionality) {
  std::size_t longest = 0;
  for (const std::size_t n : lengths) {
    if (n < 2) throw std::invalid_argument("chains need at least two beads");
    longest = std::max(longest, n);
  }
  System out =
      poly::endlink::make_system<Vec3i>(lengths, n_crosslinkers, functionality, lattice_.extent(), bond_length / std::numbers::sqrt2);
  if (out.chain_beads() + n_crosslinkers >= lattice_.num_sites() - lattice_.occupied_count()) {
    throw std::invalid_argument("more beads than empty lattice sites");
  }
  out.beads.resize(out.chain_beads());
  path_.reserve(longest);

  const std::span<Vec3i> all(out.beads);
  for (std::size_t c = 0; c < lengths.size(); ++c) {
    const auto chain = all.subspan(out.chain_start[c], lengths[c]);
    for (int attempt = 1; !grow_chain(chain); ++attempt) {
      ++dead_ends_;
      if (attempt == rules_.max_attempts) {
        throw std::runtime_error(std::format(
            "chain {} was trapped in all {} attempts; lower the occupancy or allow more bond angles", c + 1, attempt));
      }
    }
  }
  for (std::size_t k = 0; k < n_crosslinkers; ++k) {
    const Site s = random_free_site();
    lattice_.occupy(s);
    out.crosslinkers[k] = lattice_.position(s);
  }
  return out;
}

std::expected<void, std::string> validate(const System& system, AngleSet angles) {
  const auto fail = [](std::string message) { return std::unexpected(std::move(message)); };
  if (auto ok = poly::endlink::validate(system); !ok) return ok;
  if (!system.has_beads()) return fail("the chain beads are not placed");

  FccLattice lattice(system.period / 2);
  for (std::size_t i = 0; i < system.atoms(); ++i) {
    const Vec3i p = system.position(i);
    if (((p.x + p.y + p.z) & 1) != 0) return fail(std::format("bead {} is not on a lattice site", i + 1));
    const Site s = lattice.site(lattice.wrap_any(p));
    if (lattice.occupied(s)) return fail(std::format("bead {} sits on an occupied site", i + 1));
    lattice.occupy(s);
  }

  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    const std::size_t first = system.chain_start[c];
    for (std::size_t v = 1; v < system.chain_length(c); ++v) {
      const Vec3i step = system.beads[first + v] - system.beads[first + v - 1];
      if (step.norm2() != 2) return fail(std::format("bond {}-{} is not a nearest-neighbor step", first + v, first + v + 1));
      if (v >= 2 && !angles.allows_steps(system.beads[first + v - 1] - system.beads[first + v - 2], step)) {
        return fail(std::format("bond angle at bead {} is not allowed", first + v));
      }
    }
  }
  return {};
}

}  // namespace rw
