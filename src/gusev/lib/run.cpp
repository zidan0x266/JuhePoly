// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "gusev/sim/run.hpp"

#include <format>
#include <stdexcept>
#include <string>

#include "gusev/sim/linking.hpp"
#include "poly/endlink/driver.hpp"

namespace gusev {

System run(Config config, std::ostream& out) {
  poly::endlink::start(config);
  poly::endlink::Log say(out);
  poly::Rng rng(config.seed);
  const poly::endlink::Chains chains = poly::endlink::draw_chains(config, rng);

  const std::size_t n_beads = chains.drawn.beads + config.n_crosslinkers;
  const double box = config.box_edge(n_beads);
  const double imbalance = static_cast<double>(config.functionality * config.n_crosslinkers) / static_cast<double>(2 * config.n_chains);
  const std::string strands = chains.describe("strands");
  const std::string reaction = config.linking.procedure == Procedure::Ends
                                   ? std::format("target conversion {}", config.linking.conversion)
                                   : std::format("reaction probability {} per end", config.linking.conversion);
  say(std::format("Gusev network ({} procedure, {}): {} and {} cross-linkers of functionality {}, stoichiometric imbalance {:.4f}, "
                  "seed {}",
                  to_string(config.linking.procedure), reaction, strands, config.n_crosslinkers, config.functionality, imbalance,
                  config.seed));
  if (chains.polydisperse) say(std::format("strand lengths: {}", chains.distribution.describe()));
  say(std::format("box edge {:.4f}, {} beads, number density {:.4f}", box, n_beads, static_cast<double>(n_beads) / (box * box * box)));

  System system = place_crosslinkers(chains.lengths, config.n_crosslinkers, config.functionality, box, rng);
  const LinkReport report = link_strands(system, config.linking, config.bond_length, rng);
  const poly::endlink::Statistics stats = poly::endlink::statistics(system, config.bond_length);
  say(std::format("linked {} ends: conversion {:.4f} of the {} groups of the minority component; {} ends found no active site "
                  "within reach",
                  report.reactions, stats.conversion, system.minority_groups(), report.unreachable));
  say(std::format("strands: {} connected, {} dangling, {} free; {} primary loops, {} loops through a periodic image; {} of {} "
                  "cross-linkers fully linked; mean R^2/(N b^2) of the connected strands {:.4f}",
                  stats.connected, stats.dangling, stats.free, stats.primary_loops, stats.image_loops, stats.full_crosslinkers,
                  system.n_crosslinkers(), stats.span_ratio));

  place_beads(system, config.bond_length, rng);
  if (auto ok = poly::endlink::validate(system); !ok) throw std::logic_error("network check failed: " + ok.error());
  say("beads placed and checked");

  poly::endlink::write_outputs(config, system, chains,
                               std::format("LAMMPS data file: Gusev network of {} and {} cross-linkers of functionality {}, "
                                           "conversion {:.4f} (Polymerization gusev, seed {})",
                                           strands, config.n_crosslinkers, config.functionality, stats.conversion, config.seed),
                               say, "strands.dat");
  return system;
}

}  // namespace gusev
