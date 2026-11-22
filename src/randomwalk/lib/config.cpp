// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/io/config.hpp"

#include <format>
#include <sstream>
#include <stdexcept>

#include "poly/ini.hpp"

namespace rw {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_number;
using poly::ini::trim;

namespace {

/// Angles in degrees, separated by spaces or commas, e.g. "90, 120".
poly::ini::Parsed<AngleSet> parse_angles(std::string_view key, std::string_view value) {
  std::string text(value);
  for (char& ch : text) {
    if (ch == ',') ch = ' ';
  }
  AngleSet angles;
  std::istringstream in{text};
  std::string token;
  while (in >> token) {
    const auto degrees = parse_number<int>(key, token);
    if (!degrees || !angles.allow(*degrees)) {
      return std::unexpected(std::format("config key '{}': '{}' is not one of 60, 90, 120, 180", key, token));
    }
  }
  if (angles.empty()) return std::unexpected(std::format("config key '{}': no angle given", key));
  return angles;
}

std::string to_string(AngleSet angles) {
  std::string out;
  for (int degrees : AngleSet::kValues) {
    if (angles.allows(degrees)) out += std::format("{}{}", out.empty() ? "" : ", ", degrees);
  }
  return out;
}

}  // namespace

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (auto done = set_shared(key, value)) return *done;
  if (key == "occupancy") return assign(occupancy, parse_number<double>(key, value));
  if (key == "cells") return assign(cells, parse_number<int>(key, value));
  if (key == "bond_angles") return assign(walk.angles, parse_angles(key, value));
  if (key == "site_weighting") {
    const auto w = parse_site_weighting(lower(trim(value)));
    if (!w) return std::unexpected("config key 'site_weighting': expected 'legacy', 'free_neighbors' or 'uniform'");
    walk.weighting = *w;
    return {};
  }
  if (key == "max_attempts") return assign(walk.max_attempts, parse_number<int>(key, value));
  return std::unexpected(std::format("unknown configuration key '{}'", raw_key));
}

Config Config::from_file(const std::filesystem::path& path) {
  Config cfg;
  poly::ini::read_file(path, [&cfg](std::string_view key, std::string_view value) { return cfg.set(key, value); });
  return cfg;
}

void Config::apply_overrides(const std::vector<std::string>& overrides) {
  poly::ini::apply_overrides(overrides, [this](std::string_view key, std::string_view value) { return set(key, value); });
}

int Config::box_cells(std::size_t chain_beads) const { return cells > 0 ? cells : cells_for_occupancy(chain_beads, occupancy); }

void Config::validate() const {
  const auto fail = [](std::string_view msg) { throw std::invalid_argument(std::format("invalid configuration: {}", msg)); };
  validate_shared();
  if (cells < 0) fail("cells must be positive, or 0 to derive it from occupancy");
  if (cells == 0 && !(occupancy > 0 && occupancy < 1)) fail("occupancy must be in (0,1)");
  if (walk.max_attempts < 1) fail("max_attempts must be at least 1");

  const int n = box_cells();
  if (n < FccLattice::kMinCells || n > FccLattice::kMaxCells) {
    fail(std::format("the box must have {} to {} cells per edge, got {}", FccLattice::kMinCells, FccLattice::kMaxCells, n));
  }
  const double sites = 4.0 * n * n * n;
  if (static_cast<double>(expected_chain_beads() + n_crosslinkers) >= sites) fail("more beads than lattice sites");
}

std::string Config::to_ini() const {
  return std::format(
      "[system]\n"
      "{}"
      "occupancy = {}\n"
      "cells = {}\n\n"
      "[walk]\n"
      "bond_angles = {}\n"
      "site_weighting = {}\n"
      "max_attempts = {}\n\n"
      "{}",
      system_ini(), occupancy, cells, to_string(walk.angles), rw::to_string(walk.weighting), walk.max_attempts, run_output_ini());
}

}  // namespace rw
