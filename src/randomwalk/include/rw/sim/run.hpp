// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Driver: the run steps of every packing module (poly/endlink/driver.hpp) around this module's
// packing: draw the chain lengths, pack the chains and cross-linkers on the lattice, check the
// result and write the output files. The chains stay unlinked; LAMMPS links them.
#pragma once

#include <ostream>

#include "rw/io/config.hpp"
#include "rw/sim/packing.hpp"

namespace rw {

System run(Config config, std::ostream& log);

}  // namespace rw
