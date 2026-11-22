// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output of an end-linking system, the same for every packing module: a LAMMPS data file (atom
// style molecular, with image flags and atom type labels), an extended XYZ file for
// visualization and a table of the chains and their links. Defined for lattice (Vec3i) and
// continuous (Vec3d) positions.
#pragma once

#include <filesystem>
#include <ostream>
#include <string_view>

#include "poly/endlink/system.hpp"

namespace poly::endlink {

/// LAMMPS data file of a system with placed beads: the chain beads, chain after chain (chain c is
/// molecule c + 1), then the cross-linkers, each a molecule of its own; bonds along the chains
/// and from every linked end to its cross-linker. `title` becomes the first line; every type
/// gets mass `mass`.
template <class Point>
void write_lammps_data(std::ostream& out, const System<Point>& system, double mass, std::string_view title);
template <class Point>
void write_lammps_data(const std::filesystem::path& path, const System<Point>& system, double mass, std::string_view title);

/// Extended XYZ: atom type as species, wrapped positions, box in the comment line.
template <class Point>
void write_xyz(std::ostream& out, const System<Point>& system);
template <class Point>
void write_xyz(const std::filesystem::path& path, const System<Point>& system);

/// One line per chain: its length, the atom IDs of the cross-linkers of its ends (0 for an
/// unlinked end) and, with both ends linked, its span.
template <class Point>
void write_links(std::ostream& out, const System<Point>& system);
template <class Point>
void write_links(const std::filesystem::path& path, const System<Point>& system);

}  // namespace poly::endlink
