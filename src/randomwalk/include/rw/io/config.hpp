// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration: an INI-style file (`key = value`, optional [sections], `#` comments) plus
// `key=value` command-line overrides. The keys of the chains, cross-linkers, seed and output are
// those of every packing module (poly/endlink/config.hpp); this module adds the lattice box and
// the rules of the walks. The defaults are the No = 20 system of the fracture paper
// (Macromolecules 2026, 59, 99).
#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "poly/chains/lengths.hpp"
#include "poly/endlink/config.hpp"
#include "rw/sim/packing.hpp"

namespace rw {

using poly::chains::ChainLengthDistribution;
using poly::chains::LengthDistribution;
using poly::chains::LengthSpec;
using poly::chains::LengthStats;
using poly::chains::length_stats;

struct Config : poly::endlink::ConfigBase {
  // [system]
  /// Sizes the box: the smallest number of cells in which the chain beads fill less than this
  /// fraction of the lattice sites (cross-linkers are not counted, as in the original program).
  double occupancy{0.05};
  /// Cells per box edge; 0 derives it from `occupancy`.
  int cells{0};

  // [walk]
  WalkRules walk{};

  /// Load from an INI file; unknown keys are an error so typos cannot silently change runs.
  static Config from_file(const std::filesystem::path& path);
  /// Apply `key=value` overrides (keys as in the file, optionally section-qualified `run.seed`).
  void apply_overrides(const std::vector<std::string>& overrides);
  /// Sanity checks; throws std::invalid_argument with a readable message.
  void validate() const;
  /// The configuration in the same INI syntax (written next to the results).
  [[nodiscard]] std::string to_ini() const;

  /// Cells per box edge for `chain_beads` chain beads: `cells` if set, else derived from
  /// `occupancy`.
  [[nodiscard]] int box_cells(std::size_t chain_beads) const;
  /// The same for the expected number of chain beads, n_chains times the number average.
  [[nodiscard]] int box_cells() const { return box_cells(expected_chain_beads()); }

  /// Set one key from its textual value; errors are returned, not thrown.
  [[nodiscard]] std::expected<void, std::string> set(std::string_view key, std::string_view value);
};

}  // namespace rw
