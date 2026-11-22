// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include <format>
#include <utility>

#include "poly/endlink/driver.hpp"
#include "poly/endlink/writers.hpp"

namespace poly::endlink {

std::string Chains::describe(std::string_view noun) const {
  if (!polydisperse) return std::format("{} {} of {} beads", drawn.chains, noun, distribution.spec().mean);
  return std::format("{} {} of {:.4f} beads on average (dispersity {:.4f}, {} lengths of {} to {} beads)", drawn.chains, noun,
                     drawn.number_average, drawn.dispersity(), to_string(distribution.spec().distribution), drawn.shortest, drawn.longest);
}

Chains draw_chains(const ConfigBase& config, Rng& rng) {
  chains::ChainLengthDistribution distribution(config.length_spec());
  std::vector<std::size_t> lengths = distribution.sample(config.n_chains, rng);
  const chains::LengthStats drawn = chains::length_stats(lengths);
  return {std::move(distribution), std::move(lengths), drawn, !config.lengths.monodisperse()};
}

template <class Point>
void write_outputs(const ConfigBase& config, const System<Point>& system, const Chains& chains, std::string_view title, Log& log,
                   std::string_view links_file) {
  const std::filesystem::path data = config.output_dir / config.data_file;
  write_lammps_data(data, system, config.mass, title);
  log(std::format("wrote {}", data.string()));
  if (!links_file.empty()) {
    const std::filesystem::path links = config.output_dir / links_file;
    write_links(links, system);
    log(std::format("wrote {}", links.string()));
  }
  if (chains.polydisperse) {
    const std::filesystem::path table = config.output_dir / "chain_lengths.dat";
    chains::write_length_table(table, chains.distribution, chains.lengths);
    log(std::format("wrote {}", table.string()));
  }
  if (config.write_xyz) {
    const std::filesystem::path xyz = std::filesystem::path(data).replace_extension(".xyz");
    write_xyz(xyz, system);
    log(std::format("wrote {}", xyz.string()));
  }
}

template void write_outputs(const ConfigBase&, const System<Vec3i>&, const Chains&, std::string_view, Log&, std::string_view);
template void write_outputs(const ConfigBase&, const System<Vec3d>&, const Chains&, std::string_view, Log&, std::string_view);

}  // namespace poly::endlink
