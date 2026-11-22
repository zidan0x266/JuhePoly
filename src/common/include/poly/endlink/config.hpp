// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The configuration keys shared by the packing modules: the chains and cross-linkers of
// [system], the seed of [run] and the files of [output]. Each module derives its Config from
// ConfigBase and adds the keys of its method. The defaults are the No = 20 system of the fracture
// paper (Macromolecules 2026, 59, 99).
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

#include "poly/chains/lengths.hpp"

namespace poly::endlink {

struct ConfigBase {
  // [system]
  /// Chain lengths: `chain_length` (the number average), `length_distribution`, `dispersity`,
  /// `min_chain_length` and `max_chain_length`.
  chains::LengthSpec lengths{};
  std::size_t n_chains{5680};
  std::size_t n_crosslinkers{2840};
  int functionality{4};  ///< sites per cross-linker
  double bond_length{1.0};
  double mass{1.0};

  // [run]
  std::uint64_t seed{0};  ///< 0 -> entropy seed, written to run_config.ini

  // [output]
  std::filesystem::path output_dir{"."};
  std::string data_file{"config.lmp"};
  bool write_xyz{true};

  /// Set one of the shared keys, `key` already normalized; std::nullopt if it is not one of them.
  [[nodiscard]] std::optional<std::expected<void, std::string>> set_shared(std::string_view key, std::string_view value);
  /// Check the shared keys; throws std::invalid_argument with a readable message.
  void validate_shared() const;
  /// The chain lengths to draw, with no chain longer than all chain beads together when
  /// max_chain_length is 0.
  [[nodiscard]] chains::LengthSpec length_spec() const;
  /// n_chains times the number-average chain length.
  [[nodiscard]] std::size_t expected_chain_beads() const;
  /// The shared keys of [system] in INI syntax, one line each.
  [[nodiscard]] std::string system_ini() const;
  /// The [run] and [output] sections in INI syntax.
  [[nodiscard]] std::string run_output_ini() const;
};

}  // namespace poly::endlink
