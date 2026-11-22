// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/sim/run.hpp"

#include <cmath>
#include <format>
#include <stdexcept>
#include <string>

#include "poly/endlink/driver.hpp"

namespace rw {

System run(Config config, std::ostream& out) {
  poly::endlink::start(config);
  poly::endlink::Log say(out);
  poly::Rng rng(config.seed);
  const poly::endlink::Chains chains = poly::endlink::draw_chains(config, rng);

  const int cells = config.box_cells(chains.drawn.beads);
  const std::size_t n_beads = chains.drawn.beads + config.n_crosslinkers;
  const double sites = 4.0 * cells * cells * cells;
  const double box = std::sqrt(2.0) * cells * config.bond_length;
  if (chains.polydisperse) {
    say(std::format("random-walk packing: {} chains with {} ({}), and {} cross-linkers, seed {}", config.n_chains,
                    poly::chains::describe_target(config.lengths), chains.distribution.describe(), config.n_crosslinkers,
                    config.seed));
    say(std::format("chain lengths drawn: Ln {:.4f}, Lw {:.4f}, dispersity {:.4f}, {} to {} beads", chains.drawn.number_average,
                    chains.drawn.weight_average, chains.drawn.dispersity(), chains.drawn.shortest, chains.drawn.longest));
  } else {
    say(std::format("random-walk packing: {} and {} cross-linkers, seed {}", chains.describe("chains"), config.n_crosslinkers,
                    config.seed));
  }
  say(std::format("FCC lattice of {}^3 cells, box edge {:.4f}, {:.0f} sites; occupancy {:.4f}, number density {:.4f}", cells,
                  box, sites, static_cast<double>(n_beads) / sites, static_cast<double>(n_beads) / (box * box * box)));

  Packer packer(cells, config.walk, rng);
  System system = packer.pack(chains.lengths, config.n_crosslinkers, config.bond_length, config.functionality);
  if (auto ok = validate(system, config.walk.angles); !ok) throw std::logic_error("packing check failed: " + ok.error());
  say(std::format("packed and checked; {} trapped chain attempts were grown again", packer.dead_ends()));

  poly::endlink::write_outputs(config, system, chains,
                               std::format("LAMMPS data file: {} and {} cross-linkers, FCC random walks (JuhePoly "
                                           "randomwalk, seed {})",
                                           chains.describe("chains"), config.n_crosslinkers, config.seed),
                               say);
  return system;
}

}  // namespace rw
