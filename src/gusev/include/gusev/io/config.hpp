// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration: an INI-style file (`key = value`, optional [sections], `#` comments) plus
// `key=value` command-line overrides. The keys of the chains, cross-linkers, seed and output are
// those of every packing module (poly/endlink/config.hpp); this module adds the box and the
// linking procedure. The defaults are the composition of the No = 20 network of the fracture
// paper (Macromolecules 2026, 59, 99) at the bead density of a Kremer-Grest melt.
#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "gusev/sim/linking.hpp"
#include "poly/chains/lengths.hpp"
#include "poly/endlink/config.hpp"

namespace gusev {

using poly::chains::ChainLengthDistribution;
using poly::chains::LengthDistribution;
using poly::chains::LengthSpec;
using poly::chains::LengthStats;
using poly::chains::length_stats;

struct Config : poly::endlink::ConfigBase {
  // [system]
  /// Number density of all beads, chains and cross-linkers; sets the box unless `box` is set.
  double density{0.85};
  double box{0.0};  ///< box edge; 0 derives it from `density`

  // [network]
  LinkRules linking{};

  /// Load from an INI file; unknown keys are an error so typos cannot silently change runs.
  static Config from_file(const std::filesystem::path& path);
  /// Apply `key=value` overrides (keys as in the file, optionally section-qualified `run.seed`).
  void apply_overrides(const std::vector<std::string>& overrides);
  /// Sanity checks; throws std::invalid_argument with a readable message.
  void validate() const;
  /// The configuration in the same INI syntax (written next to the results).
  [[nodiscard]] std::string to_ini() const;

  /// Box edge for `beads` beads: `box` if set, else from `density`.
  [[nodiscard]] double box_edge(std::size_t beads) const;

  /// Set one key from its textual value; errors are returned, not thrown.
  [[nodiscard]] std::expected<void, std::string> set(std::string_view key, std::string_view value);
};

}  // namespace gusev
