// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "gusev/io/config.hpp"

#include <cmath>
#include <format>
#include <stdexcept>

#include "poly/ini.hpp"

namespace gusev {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_number;
using poly::ini::trim;

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (auto done = set_shared(key, value)) return *done;
  if (key == "density") return assign(density, parse_number<double>(key, value));
  if (key == "box") return assign(box, parse_number<double>(key, value));
  if (key == "procedure") {
    const auto p = parse_procedure(lower(trim(value)));
    if (!p) return std::unexpected("config key 'procedure': expected 'ends' or 'strands'");
    linking.procedure = *p;
    return {};
  }
  if (key == "conversion") return assign(linking.conversion, parse_number<double>(key, value));
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

double Config::box_edge(std::size_t beads) const { return box > 0 ? box : std::cbrt(static_cast<double>(beads) / density); }

void Config::validate() const {
  const auto fail = [](std::string_view msg) { throw std::invalid_argument(std::format("invalid configuration: {}", msg)); };
  validate_shared();
  if (n_crosslinkers < 1) fail("n_crosslinkers must be at least 1");
  if (box < 0) fail("box must be positive, or 0 to derive it from density");
  if (box == 0 && !(density > 0)) fail("density must be positive");
  if (!(linking.conversion > 0 && linking.conversion <= 1)) fail("conversion must be in (0, 1]");
}

std::string Config::to_ini() const {
  return std::format(
      "[system]\n"
      "{}"
      "density = {}\n"
      "box = {}\n\n"
      "[network]\n"
      "procedure = {}\n"
      "conversion = {}\n\n"
      "{}",
      system_ini(), density, box, to_string(linking.procedure), linking.conversion, run_output_ini());
}

}  // namespace gusev
