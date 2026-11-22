// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Linking of the chains by Gusev's Monte Carlo procedure for end-linked networks, after the
// cross-linkers have been placed (sim/packing.hpp) and ignoring excluded volume. An end whose
// chain has no linked end yet reacts with an active site chosen at random in the whole box. An
// end whose other end is linked at X reacts with an active site within the contour length N b of
// X, chosen with probability proportional to exp[-3 d^2 / (2 N b^2)], the Gaussian distribution of
// the end-to-end vector of N segments of length b; d is measured to any periodic image within N b,
// and the link records the image. Two variants:
//
//   ends     (Gusev & Schwarz, Macromolecules 2019, 52, 9445): each step reacts an end chosen at
//            random from all remaining active ends, until the target conversion is reached;
//   strands  (Gusev, Macromolecules 2019, 52, 3244): the chains are taken one after another, and
//            each of their ends reacts with the probability xi.
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

#include "gusev/sim/packing.hpp"
#include "poly/random.hpp"

namespace gusev {

enum class Procedure { Ends, Strands };
[[nodiscard]] std::string_view to_string(Procedure p) noexcept;
[[nodiscard]] std::optional<Procedure> parse_procedure(std::string_view s) noexcept;

struct LinkRules {
  Procedure procedure{Procedure::Ends};
  /// ends: the target conversion of the minority groups; strands: the probability xi that an
  /// end reacts.
  double conversion{0.95};
};

struct LinkReport {
  std::size_t reactions{0};    ///< ends linked
  std::size_t unreachable{0};  ///< ends that found no active site within reach
};

/// The cross-linkers with active sites near a point, through all periodic images, from a cell list.
class SiteFinder {
 public:
  struct Candidate {
    std::uint32_t crosslinker;
    Vec3d displacement;  ///< from the point to this image of the cross-linker
    Vec3i image;         ///< the image, in box edges
  };

  explicit SiteFinder(const System& system);
  /// Append to `out` every image within `radius` of `point` (in the box) of a cross-linker with
  /// an active site.
  void within(const Vec3d& point, double radius, std::vector<Candidate>& out) const;

 private:
  const System* system_;
  int cells_;          ///< per box edge
  double cell_size_;
  std::vector<std::uint32_t> cell_start_;
  std::vector<std::uint32_t> members_;
};

/// Link the chain ends of `system` by Gusev's procedure, with segments of length `bond_length`.
/// The first end linked of a chain sees its cross-linker in the box, and a chain with both ends
/// linked is seen from end 0, so its beads start at that cross-linker.
LinkReport link_strands(System& system, const LinkRules& rules, double bond_length, poly::Rng& rng);

}  // namespace gusev
