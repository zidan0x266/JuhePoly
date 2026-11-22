// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency) of the code shared by the modules: lattice
// geometry, molecule and bond bookkeeping, the bond tally and the output writers.

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "poly/chem/progress.hpp"
#include "poly/endlink/system.hpp"
#include "poly/endlink/writers.hpp"
#include "poly/geometry/fcc.hpp"
#include "poly/geometry/vec3.hpp"
#include "poly/io/output_file.hpp"
#include "poly/io/pdb.hpp"
#include "poly/io/psf.hpp"
#include "poly/random.hpp"
#include "poly/topology/bond_table.hpp"
#include "poly/topology/molecules.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond)                                                                                  \
  do {                                                                                               \
    ++g_checks;                                                                                      \
    if (!(cond)) {                                                                                   \
      ++g_failures;                                                                                  \
      std::fprintf(stderr, "CHECK failed: %s  (%s:%d)\n", #cond, __FILE__, __LINE__);                \
    }                                                                                                \
  } while (0)

using poly::geometry::FccGeometry;
using poly::geometry::Site;
using poly::geometry::Vec3i;

std::string read_text(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

void test_geometry() {
  constexpr Vec3i a{1, -2, 3};
  static_assert(a + Vec3i{1, 1, 1} == Vec3i{2, -1, 4});
  static_assert(a.dot(Vec3i{2, 1, 0}) == 0 && a.norm2() == 14 && a[0] == 1 && a[1] == -2 && a[2] == 3);

  // Every site maps to a lattice point and back; the twelve neighbors are distinct sites.
  const FccGeometry lattice(3);
  CHECK(lattice.num_sites() == 108 && lattice.extent() == 6);
  bool round_trip = true, neighbors_distinct = true;
  for (Site s = 0; s < lattice.num_sites(); ++s) {
    const Vec3i p = lattice.position(s);
    round_trip = round_trip && ((p.x + p.y + p.z) & 1) == 0 && lattice.site(p) == s;
    std::set<Site> around;
    for (int k = 0; k < FccGeometry::kCoordination; ++k) around.insert(lattice.neighbor(p, k));
    neighbors_distinct = neighbors_distinct && around.size() == 12 && !around.contains(s);
  }
  CHECK(round_trip);
  CHECK(neighbors_distinct);
  CHECK((lattice.minimum_image(Vec3i{5, -5, 3}) == Vec3i{-1, 1, -3}));
}

void test_endlink() {
  // Two chains of 3 and 2 beads and two cross-linkers on a lattice of period 4 (scale 0.5); the
  // second chain lies across the box boundary.
  const std::vector<std::size_t> lengths{3, 2};
  auto s = poly::endlink::make_system<Vec3i>(lengths, 2, 4, 4, 0.5);
  s.beads = {{0, 0, 0}, {1, 1, 0}, {2, 2, 0}, {-1, -1, 0}, {-2, -2, 0}};
  s.crosslinkers = {{3, 3, 0}, {2, 0, 2}};
  CHECK(s.n_chains() == 2 && s.chain_beads() == 5 && s.atoms() == 7 && s.bonds() == 3 && s.links.size() == 4);
  CHECK(s.type(0) == poly::endlink::AtomType::End && s.type(1) == poly::endlink::AtomType::Bead && s.molecule(4) == 2 && s.molecule(6) == 4);
  s.attach(0, 1, 0);
  s.attach(1, 0, 0, {-1, -1, 0});
  CHECK(s.bonds() == 5 && s.used_sites[0] == 2 && s.type(2) == poly::endlink::AtomType::Bead && validate(s).has_value());
  auto stats = poly::endlink::statistics(s, 1.0);
  CHECK(stats.dangling == 2 && stats.connected == 0 && stats.linked_ends == 2 && stats.conversion == 0.5);

  std::ostringstream data;
  poly::endlink::write_lammps_data(data, s, 1.0, "test");
  CHECK(data.str().starts_with("test\n\n7 atoms\n4 atom types\n5 bonds\n1 bond types\n\n0 2.000000 xlo xhi\n"));
  CHECK(data.str().contains("Atoms # molecular\n\n"
                            "1 1 2 0.000000 0.000000 0.000000 0 0 0\n"
                            "2 1 1 0.500000 0.500000 0.000000 0 0 0\n"
                            "3 1 1 1.000000 1.000000 0.000000 0 0 0\n"
                            "4 2 1 1.500000 1.500000 0.000000 -1 -1 0\n"
                            "5 2 2 1.000000 1.000000 0.000000 -1 -1 0\n"
                            "6 3 3 1.500000 1.500000 0.000000 0 0 0\n"
                            "7 4 3 1.000000 0.000000 1.000000 0 0 0\n"));
  CHECK(data.str().ends_with("Bonds\n\n1 1 1 2\n2 1 2 3\n3 1 3 6\n4 1 6 4\n5 1 4 5\n"));

  // both ends of chain 0 on cross-linker 0 through the same image: a primary loop
  s.attach(0, 0, 0);
  stats = poly::endlink::statistics(s, 1.0);
  CHECK(stats.connected == 1 && stats.primary_loops == 1 && s.span(0) == Vec3i{} && validate(s).has_value());
  s.used_sites[1] = 1;  // a site used by no chain
  CHECK(!validate(s).has_value());

  // continuous positions: wrapped into the box with their images
  const std::vector<std::size_t> one{2};
  auto c = poly::endlink::make_system<poly::geometry::Vec3d>(one, 1, 4, 5.0, 1.0);
  c.beads = {{-0.5, 2.0, 5.25}, {0.25, 2.0, 4.75}};
  c.crosslinkers = {{1.0, 1.0, 1.0}};
  std::ostringstream xyz;
  poly::endlink::write_xyz(xyz, c);
  CHECK(xyz.str().ends_with("2 4.500000 2.000000 0.250000\n2 0.250000 2.000000 4.750000\n3 1.000000 1.000000 1.000000\n"));
  std::ostringstream cdata;
  poly::endlink::write_lammps_data(cdata, c, 1.0, "c");
  CHECK(cdata.str().contains("1 1 2 4.500000 2.000000 0.250000 -1 0 1\n"));
}

void test_normal() {
  // standard normal deviates: mean 0, variance 1, fourth moment 3; the same seed repeats them
  poly::Rng rng(31), again(31);
  const int n = 200000;
  double m1 = 0, m2 = 0, m4 = 0;
  bool repeats = true;
  for (int i = 0; i < n; ++i) {
    const double x = rng.normal();
    repeats = repeats && x == again.normal();
    m1 += x;
    m2 += x * x;
    m4 += x * x * x * x;
  }
  CHECK(repeats);
  CHECK(std::fabs(m1 / n) < 0.01 && std::fabs(m2 / n - 1) < 0.015 && std::fabs(m4 / n - 3) < 0.1);
}

void test_molecules() {
  poly::topology::Molecules<int> m(6);
  CHECK(m.size() == 6 && m.find(4) == 4 && m.units(4) == 1 && m.is_representative(4));
  // On a tie the molecule of the first unit keeps its representative.
  const auto first = m.unite(1, 2);
  CHECK(first && first->kept == 1 && first->absorbed == 2);
  // The larger molecule keeps its representative, whatever the order.
  const auto second = m.unite(3, 2);
  CHECK(second && second->kept == 1 && second->absorbed == 3);
  CHECK(m.find(3) == 1 && m.units(1) == 3 && !m.is_representative(3));
  CHECK(!m.unite(2, 3).has_value());
  CHECK(m.same(2, 3) && !m.same(0, 3));
  m.reset(2);
  CHECK(m.size() == 2 && !m.same(0, 1));
}

void test_bond_table() {
  poly::topology::BondTable<unsigned, 2> bonds(4);
  CHECK(bonds.size() == 4 && bonds.count(0) == 0 && bonds.partners(0).empty());
  bonds.link(0, 1);
  bonds.link(0, 3);
  CHECK(bonds.count(0) == 2 && bonds.full(0) && !bonds.full(1));
  CHECK(bonds.bonded(0, 3) && bonds.bonded(3, 0) && !bonds.bonded(1, 3));
  CHECK(bonds.partners(0)[0] == 1 && bonds.partners(0)[1] == 3 && bonds.partners(1)[0] == 0);
}

void test_tally() {
  poly::chem::ConnectivityTally t;
  int asked = 0;
  const auto same = [&asked](bool answer) {
    return [&asked, answer] {
      ++asked;
      return answer;
    };
  };
  t.count(0, 0, same(true));
  t.count(1, 0, same(true));
  t.count(0, 2, same(true));
  CHECK(asked == 0);  // only bonds between two bonded units need the molecules
  t.count(1, 1, same(true));
  t.count(2, 1, same(false));
  CHECK(asked == 2);
  CHECK(t.monomer_monomer == 1 && t.polymer_monomer == 2 && t.intramolecular == 1 && t.polymer_polymer == 1);
}

void test_writers() {
  const auto dir = std::filesystem::temp_directory_path() / "poly_common_test";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  const std::array<std::string_view, 3> residues{"MON", "ONE", "CHN"};
  const std::vector<std::vector<int>> partners{{1}, {0, 2}, {1}};

  poly::io::write_psf(dir / "a.psf", 3, [&](std::size_t i) { return residues[i]; }, [&](std::size_t i) { return partners[i]; });
  CHECK(read_text(dir / "a.psf") ==
        "       3 !NATOM\n"
        "       1 MAIN 1    MON  P    U            0.010000      0.010000           0\n"
        "       2 MAIN 1    ONE  P    U            0.010000      0.010000           0\n"
        "       3 MAIN 1    CHN  P    U            0.010000      0.010000           0\n"
        "\n"
        "       4 !NBOND: bonds\n"
        "       1       2       2       1       2       3       3       2\n"
        "END\n");

  poly::io::write_pdb(dir / "a.pdb", "TEST SYSTEM", {4.0, 5.0, 6.5}, 2, [&](std::size_t i) { return residues[i]; },
                      [](std::size_t i) { return std::array<double, 3>{1.0 * static_cast<double>(i), 2.5, -0.125}; });
  CHECK(read_text(dir / "a.pdb") ==
        "TITLE     TEST SYSTEM\n"
        "REMARK    SHOW IN SIMULATION BOX\n"
        "CRYST1    4.000    5.000    6.500  90.00  90.00  90.00 P 1\n"
        "ATOM      1   MON mon    1       0.000   2.500  -0.125  0.00  0.00\n"
        "ATOM      2   ONE one    1       1.000   2.500  -0.125  0.00  0.00\n"
        "END\n");

  bool threw = false;
  try {
    [[maybe_unused]] auto out = poly::io::open_output(dir / "missing" / "a.txt");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  std::filesystem::remove_all(dir);
}

}  // namespace

int main() {
  test_geometry();
  test_endlink();
  test_normal();
  test_molecules();
  test_bond_table();
  test_tally();
  test_writers();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
