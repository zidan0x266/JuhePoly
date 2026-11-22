// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The steps that every packing run shares, before and after its method: the log with elapsed
// times, the seed and run_config.ini, the chain lengths drawn from the configured distribution,
// and the output files.
#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <string>
#include <string_view>
#include <vector>

#include "poly/chains/lengths.hpp"
#include "poly/endlink/config.hpp"
#include "poly/endlink/system.hpp"
#include "poly/random.hpp"

namespace poly::endlink {

/// Log lines, each prefixed with the time since the run started.
class Log {
 public:
  explicit Log(std::ostream& out) : out_(out), started_(std::chrono::steady_clock::now()) {}
  void operator()(const std::string& message) {
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started_).count();
    std::println(out_, "[{:9.3f} s] {}", elapsed, message);
    out_.flush();
  }

 private:
  std::ostream& out_;
  std::chrono::steady_clock::time_point started_;
};

/// Check the configuration, draw a fresh seed if it is 0 and write the configuration, seed
/// included, to run_config.ini in the output directory, so that the run can be repeated.
template <class Config>
void start(Config& config) {
  config.validate();
  config.seed = config.seed ? config.seed : Rng::entropy_seed();
  std::filesystem::create_directories(config.output_dir);
  std::ofstream dump(config.output_dir / "run_config.ini");
  std::print(dump, "# effective configuration of this run (seed resolved)\n{}", config.to_ini());
}

/// The chain lengths of a run.
struct Chains {
  chains::ChainLengthDistribution distribution;
  std::vector<std::size_t> lengths;
  chains::LengthStats drawn;
  bool polydisperse;

  /// "200 chains of 10 beads", or with a distribution "200 chains of 9.9965 beads on average
  /// (dispersity 1.4973, schulz_zimm lengths of 2 to 52 beads)"; `noun` names the chains.
  [[nodiscard]] std::string describe(std::string_view noun) const;
};

/// Draw the chain lengths from the distribution of `config`. Monodisperse lengths take no random
/// numbers, so they leave the random sequence of the packing unchanged.
[[nodiscard]] Chains draw_chains(const ConfigBase& config, Rng& rng);

/// Write the LAMMPS data file with `title`, the table of links when `links_file` is given, the
/// table of chain lengths of a polydisperse run and the XYZ file, and log each.
template <class Point>
void write_outputs(const ConfigBase& config, const System<Point>& system, const Chains& chains, std::string_view title, Log& log,
                   std::string_view links_file = {});

}  // namespace poly::endlink
