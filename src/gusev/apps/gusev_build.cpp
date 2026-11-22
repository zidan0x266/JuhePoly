// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// gusev_build: builds an end-linked network by Gusev's Monte Carlo procedure, with the strands as
// Brownian bridges between their cross-linkers, and writes a LAMMPS data file.
// See src/gusev/README.md.
//
//   gusev_build config.ini [key=value ...]
//   gusev_build procedure=strands conversion=0.99
//   gusev_build --print-config         (dump the defaults in INI syntax)
//   gusev_build --help

#include <exception>
#include <iostream>
#include <ostream>  // std::println overload for std::ostream
#include <print>
#include <string>
#include <vector>

#include "gusev/io/config.hpp"
#include "gusev/sim/run.hpp"

namespace {

constexpr std::string_view kUsage =
    "usage: gusev_build [config.ini] [key=value ...]\n"
    "       gusev_build --print-config\n\n"
    "Builds an end-linked network by Gusev's Monte Carlo procedure: cross-linkers placed in\n"
    "a periodic box, strand ends linked to them by Gaussian end-to-end statistics, strands as\n"
    "Brownian bridges. Writes a LAMMPS data file. Keys given on the command line override the\n"
    "file. Run with --print-config to see every key and its default.\n\n"
    "Please cite the papers listed in README.md under 'Acknowledgements'.\n";

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    gusev::Config cfg;
    std::vector<std::string> overrides;
    for (const std::string& a : args) {
      if (a == "-h" || a == "--help") {
        std::print("{}", kUsage);
        return 0;
      }
      if (a == "--print-config") {
        std::print("{}", cfg.to_ini());
        return 0;
      }
      if (a.contains('=')) {
        overrides.push_back(a);
      } else {
        cfg = gusev::Config::from_file(a);
      }
    }
    cfg.apply_overrides(overrides);
    gusev::run(std::move(cfg), std::cout);
    return 0;
  } catch (const std::exception& e) {
    std::println(std::cerr, "gusev_build: {}", e.what());
    return 1;
  }
}
