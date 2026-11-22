// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>

#include "poly/endlink/config.hpp"
#include "poly/ini.hpp"

namespace poly::endlink {

std::optional<std::expected<void, std::string>> ConfigBase::set_shared(std::string_view key, std::string_view value) {
  using ini::assign;
  using ini::parse_number;
  if (auto done = chains::set_length_key(lengths, key, value)) return done;
  if (key == "n_chains") return assign(n_chains, parse_number<std::size_t>(key, value));
  if (key == "n_crosslinkers") return assign(n_crosslinkers, parse_number<std::size_t>(key, value));
  if (key == "functionality") return assign(functionality, parse_number<int>(key, value));
  if (key == "bond_length") return assign(bond_length, parse_number<double>(key, value));
  if (key == "mass") return assign(mass, parse_number<double>(key, value));
  if (key == "seed") return assign(seed, parse_number<std::uint64_t>(key, value));
  if (key == "directory" || key == "output_dir") {
    output_dir = ini::trim(value);
    return std::expected<void, std::string>{};
  }
  if (key == "data_file") {
    data_file = ini::trim(value);
    return std::expected<void, std::string>{};
  }
  if (key == "write_xyz") return assign(write_xyz, ini::parse_bool(key, value));
  return std::nullopt;
}

void ConfigBase::validate_shared() const {
  const auto fail = [](std::string_view msg) { throw std::invalid_argument(std::format("invalid configuration: {}", msg)); };
  if (!(lengths.mean >= 2)) fail("chain_length must be at least 2");
  if (n_chains < 1) fail("n_chains must be at least 1");
  if (functionality < 1 || functionality > 255) fail("functionality must be from 1 to 255");
  if (static_cast<double>(n_crosslinkers) * functionality > std::numeric_limits<std::uint32_t>::max() ||
      2.0 * static_cast<double>(n_chains) > std::numeric_limits<std::uint32_t>::max()) {
    fail("too many cross-linker sites or chain ends");
  }
  [[maybe_unused]] const chains::ChainLengthDistribution check(length_spec());  // throws if the lengths cannot be drawn
  if (!(bond_length > 0)) fail("bond_length must be positive");
  if (!(mass > 0)) fail("mass must be positive");
  if (data_file.empty()) fail("data_file must not be empty");
}

chains::LengthSpec ConfigBase::length_spec() const { return chains::bounded_by_system(lengths, n_chains); }

std::size_t ConfigBase::expected_chain_beads() const {
  return static_cast<std::size_t>(std::llround(static_cast<double>(n_chains) * lengths.mean));
}

std::string ConfigBase::system_ini() const {
  return std::format("{}n_chains = {}\nn_crosslinkers = {}\nfunctionality = {}\nbond_length = {}\nmass = {}\n",
                     chains::length_keys_ini(lengths), n_chains, n_crosslinkers, functionality, bond_length, mass);
}

std::string ConfigBase::run_output_ini() const {
  return std::format("[run]\nseed = {}\n\n[output]\ndirectory = {}\ndata_file = {}\nwrite_xyz = {}\n", seed, output_dir.string(), data_file,
                     write_xyz ? "true" : "false");
}

}  // namespace poly::endlink
