// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Driver: the run steps of every packing module (poly/endlink/driver.hpp) around this module's
// method: draw the chain lengths, place the cross-linkers, link the chains by Gusev's procedure,
// place the beads, check the result and write the output files, with the table of links as
// strands.dat.
#pragma once

#include <ostream>

#include "gusev/io/config.hpp"
#include "gusev/sim/packing.hpp"

namespace gusev {

System run(Config config, std::ostream& log);

}  // namespace gusev
