// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency) of the gusev module: the search for active
// sites through periodic images, the probabilities with which ends choose their sites, the two
// linking procedures, the Brownian bridges and the output files.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "gusev/io/config.hpp"
#include "gusev/sim/linking.hpp"
#include "gusev/sim/packing.hpp"
#include "poly/endlink/writers.hpp"

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

using namespace gusev;
using poly::endlink::Statistics;

System random_network(std::size_t n_strands, std::size_t length, std::size_t n_crosslinkers, int f, double box, std::uint64_t seed) {
  poly::Rng rng(seed);
  const std::vector<std::size_t> lengths(n_strands, length);
  return place_crosslinkers(lengths, n_crosslinkers, f, box, rng);
}

// ------------------------------------------------------------------------------------------------
void test_site_finder() {
  // Every image within the radius, also for radii beyond half the box, against brute force.
  System net = random_network(1, 4, 60, 4, 7.0, 3);
  net.used_sites[5] = 4;  // a cross-linker without active sites is never found
  const SiteFinder finder(net);
  poly::Rng rng(9);
  bool same = true;
  for (int trial = 0; trial < 40; ++trial) {
    const Vec3d p{7.0 * rng.uniform01(), 7.0 * rng.uniform01(), 7.0 * rng.uniform01()};
    const double radius = 0.5 + 6.0 * rng.uniform01();
    std::vector<SiteFinder::Candidate> found;
    finder.within(p, radius, found);
    std::set<std::tuple<std::uint32_t, long, long, long>> fast, brute;
    bool images_ok = true;
    for (const auto& c : found) {
      const Vec3d shift = c.displacement + p - net.crosslinkers[c.crosslinker];
      fast.insert({c.crosslinker, c.image.x, c.image.y, c.image.z});
      images_ok = images_ok && std::lround(shift.x / 7.0) == c.image.x && std::lround(shift.y / 7.0) == c.image.y &&
                  std::lround(shift.z / 7.0) == c.image.z;
    }
    same = same && images_ok;
    for (std::uint32_t k = 0; k < net.n_crosslinkers(); ++k) {
      if (net.free_sites(k) == 0) continue;
      for (long i = -2; i <= 2; ++i) {
        for (long j = -2; j <= 2; ++j) {
          for (long l = -2; l <= 2; ++l) {
            const Vec3d d = net.crosslinkers[k] + Vec3d{7.0 * static_cast<double>(i), 7.0 * static_cast<double>(j), 7.0 * static_cast<double>(l)} - p;
            if (d.norm2() <= radius * radius) brute.insert({k, i, j, l});
          }
        }
      }
    }
    same = same && fast == brute && fast.size() == found.size();
  }
  CHECK(same);
}

// ------------------------------------------------------------------------------------------------
void test_selection_probabilities() {
  // One strand of 2 beads (3 segments: <R^2> = 3 b^2, cutoff 3 b) with end 0 on cross-linker 0;
  // the other end chooses among cross-linker 0 itself (3 active sites, d = 0) and three others
  // at d = 1, 2 and 3 with 4 sites each; one more at d = 3.5 is beyond the cutoff.
  const std::vector<Vec3d> positions{{10, 10, 10}, {11, 10, 10}, {10, 12, 10}, {10, 10, 7}, {13.5, 10, 10}};
  std::vector<double> expected{3.0, 4 * std::exp(-0.5), 4 * std::exp(-2.0), 4 * std::exp(-4.5), 0.0};  // exp(-d^2 / 2)
  double total = 0.0;
  for (const double w : expected) total += w;
  std::vector<double> count(positions.size(), 0.0);
  const int n = 40000;
  bool spans_ok = true;
  for (int trial = 0; trial < n; ++trial) {
    const std::vector<std::size_t> lengths{2};
    System net = poly::endlink::make_system<Vec3d>(lengths, positions.size(), 4, 20.0, 1.0);
    net.crosslinkers = positions;
    net.attach(0, 0, 0);
    poly::Rng rng(static_cast<std::uint64_t>(trial) + 1);
    const LinkReport report = link_strands(net, LinkRules{Procedure::Ends, 1.0}, 1.0, rng);
    const std::uint32_t k = net.link(0, 1).crosslinker;
    count[k] += 1;
    spans_ok = spans_ok && report.reactions == 1 && net.span(0) == positions[k] - positions[0];
  }
  CHECK(spans_ok);
  bool probabilities_ok = true;
  for (std::size_t k = 0; k < positions.size(); ++k) {
    const double p = expected[k] / total;
    const double sd = std::sqrt(n * p * (1 - p)) + 1e-9;
    probabilities_ok = probabilities_ok && std::fabs(count[k] - n * p) < 4 * sd + 0.5;
  }
  CHECK(probabilities_ok && count[4] == 0);
}

// ------------------------------------------------------------------------------------------------
void test_procedures() {
  // ends: the conversion of the minority groups is exactly the target
  for (const double p : {0.5, 0.8, 0.95}) {
    System net = random_network(200, 10, 100, 4, 13.5, 5);
    poly::Rng rng(11);
    const LinkReport report = link_strands(net, LinkRules{Procedure::Ends, p}, 1.0, rng);
    const Statistics s = statistics(net, 1.0);
    CHECK(s.linked_ends == report.reactions && report.reactions == static_cast<std::size_t>(std::llround(p * 400)));
    CHECK(validate(net).has_value());
  }
  // ends with too few cross-linker sites: all of them are used
  System few = random_network(200, 10, 40, 4, 13.5, 5);
  poly::Rng rng_few(12);
  (void)link_strands(few, LinkRules{Procedure::Ends, 1.0}, 1.0, rng_few);
  CHECK(statistics(few, 1.0).full_crosslinkers == 40 && statistics(few, 1.0).conversion == 1.0 && validate(few).has_value());

  // strands: with active sites always in reach, the ends react independently with probability xi
  const double xi = 0.7;
  System net = random_network(4000, 10, 8000, 4, 30.0, 6);
  poly::Rng rng(13);
  const LinkReport report = link_strands(net, LinkRules{Procedure::Strands, xi}, 1.0, rng);
  const Statistics s = statistics(net, 1.0);
  const auto close = [](double count, double p) { return std::fabs(count - 4000 * p) < 4 * std::sqrt(4000 * p * (1 - p)); };
  CHECK(report.unreachable == 0 && validate(net).has_value());
  CHECK(close(static_cast<double>(s.free), (1 - xi) * (1 - xi)) && close(static_cast<double>(s.dangling), 2 * xi * (1 - xi)) &&
        close(static_cast<double>(s.connected), xi * xi));
  // Gaussian linking: R^2 / (N b^2) near 1 when sites are plentiful
  CHECK(s.span_ratio > 0.9 && s.span_ratio < 1.1);

  // strands of different lengths: each reaches no farther than its own N b, and short strands
  // close more primary loops than long ones
  std::vector<std::size_t> lengths(3000);
  for (std::size_t i = 0; i < lengths.size(); ++i) lengths[i] = i % 2 ? 3 : 30;
  poly::Rng mixed_rng(14);
  System mixed = place_crosslinkers(lengths, 1500, 4, 25.0, mixed_rng);
  (void)link_strands(mixed, LinkRules{Procedure::Ends, 0.95}, 1.0, mixed_rng);
  bool within_reach = true;
  std::array<double, 2> loops{}, strands{};
  for (std::size_t c = 0; c < mixed.n_chains(); ++c) {
    if (mixed.linked_ends(c) < 2) continue;
    const double segments = static_cast<double>(mixed.chain_length(c) + 1);
    within_reach = within_reach && mixed.span(c).norm2() <= segments * segments * (1 + 1e-12);
    const std::size_t short_strand = mixed.chain_length(c) == 3 ? 0 : 1;
    strands[short_strand] += 1;
    loops[short_strand] += mixed.link(c, 0).crosslinker == mixed.link(c, 1).crosslinker && mixed.span(c).norm2() == 0.0 ? 1 : 0;
  }
  CHECK(within_reach && validate(mixed).has_value());
  CHECK(loops[0] / strands[0] > 3 * loops[1] / strands[1]);
}

// ------------------------------------------------------------------------------------------------
void test_placement() {
  // uniform positions in the box: every octant holds about an eighth of them
  poly::Rng rng(4);
  const auto points = place_crosslinkers({}, 4000, 4, 20.0, rng).crosslinkers;
  std::vector<int> octant(8, 0);
  bool inside = true;
  for (const Vec3d& p : points) {
    inside = inside && p.x >= 0 && p.x < 20 && p.y >= 0 && p.y < 20 && p.z >= 0 && p.z < 20;
    ++octant[static_cast<std::size_t>((p.x >= 10 ? 4 : 0) + (p.y >= 10 ? 2 : 0) + (p.z >= 10 ? 1 : 0))];
  }
  CHECK(points.size() == 4000 && inside);
  CHECK(std::ranges::all_of(octant, [](int n) { return std::abs(n - 500) < 4 * 21; }));  // binomial sd = 21
}

// ------------------------------------------------------------------------------------------------
void test_bridges() {
  // Gaussian bridge of n = 8 steps from 0 to (3, 0, 0): point k has mean 3k/n along x and
  // variance (b^2/3) k (n - k)/n per component; the mean square step is b^2 (1 - 1/n) + R^2/n^2.
  const std::size_t n = 8;
  const Vec3d to{3, 0, 0};
  poly::Rng rng(17);
  std::vector<Vec3d> points(n - 1);
  double mean_x = 0, var_y = 0, step2 = 0;
  const int samples = 40000;
  for (int i = 0; i < samples; ++i) {
    brownian_bridge({0, 0, 0}, to, 1.0, rng, points);
    mean_x += points[3].x;  // k = 4
    var_y += points[3].y * points[3].y;
    Vec3d prev{0, 0, 0};
    for (const Vec3d& p : points) {
      step2 += (p - prev).norm2();
      prev = p;
    }
    step2 += (to - prev).norm2();
  }
  mean_x /= samples;
  var_y /= samples;
  step2 /= samples * static_cast<double>(n);
  CHECK(std::fabs(mean_x - 1.5) < 0.02);
  CHECK(std::fabs(var_y - 4.0 * 4.0 / 8.0 / 3.0) < 0.02);
  CHECK(std::fabs(step2 - (1.0 - 1.0 / 8.0 + 9.0 / 64.0)) < 0.01);

  // beads of a linked network: one per strand bead, starting from the attached cross-linkers
  System net = random_network(100, 6, 60, 4, 9.0, 8);
  poly::Rng link_rng(19);
  (void)link_strands(net, LinkRules{Procedure::Strands, 0.8}, 1.0, link_rng);
  place_beads(net, 1.0, link_rng);
  CHECK(net.has_beads() && net.beads.size() == 600 && validate(net).has_value());
}

// ------------------------------------------------------------------------------------------------
void test_data_file() {
  System net = random_network(80, 5, 40, 4, 8.0, 21);
  poly::Rng rng(23);
  (void)link_strands(net, LinkRules{Procedure::Ends, 0.9}, 1.0, rng);
  place_beads(net, 1.0, rng);
  std::ostringstream data;
  write_lammps_data(data, net, 1.0, "test");
  std::istringstream in(data.str());
  std::string line, section;
  std::map<std::string, std::size_t> header;
  std::vector<std::vector<double>> atoms;
  std::vector<std::vector<long>> bonds;
  std::getline(in, line);
  CHECK(line == "test");
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream words(line);
    if (line.ends_with(" atoms") || line.ends_with(" bonds")) {
      std::size_t n = 0;
      words >> n;
      header[line.substr(line.find(' ') + 1)] = n;
      continue;
    }
    if (std::isalpha(static_cast<unsigned char>(line[0]))) {
      section = line;
      continue;
    }
    if (section == "Atoms # molecular") {
      std::vector<double> v(9);
      for (double& x : v) words >> x;
      atoms.push_back(v);
    } else if (section == "Bonds") {
      std::vector<long> v(4);
      for (long& x : v) words >> x;
      bonds.push_back(v);
    }
  }
  const Statistics s = statistics(net, 1.0);
  CHECK(atoms.size() == 440 && header["atoms"] == 440);
  CHECK(bonds.size() == 80 * 4 + s.linked_ends && header["bonds"] == bonds.size());
  // types: open strand ends are `end`, full cross-linkers `xlink_full`
  std::size_t ends = 0, full = 0;
  for (const auto& a : atoms) {
    ends += a[2] == 2 ? 1 : 0;
    full += a[2] == 4 ? 1 : 0;
  }
  CHECK(ends == 160 - s.linked_ends && full == s.full_crosslinkers);
  // every bond joins atoms a few bond lengths apart at most (minimum image: an end attached to a
  // periodic image of its cross-linker is bonded across the box)
  bool bonds_ok = true;
  std::vector<int> crosslinker_bonds(40, 0);
  for (const auto& b : bonds) {
    const auto& p = atoms[static_cast<std::size_t>(b[2] - 1)];
    const auto& q = atoms[static_cast<std::size_t>(b[3] - 1)];
    double d2 = 0;
    for (int c = 0; c < 3; ++c) {
      double d = q[3 + c] - p[3 + c];
      d -= 8.0 * std::round(d / 8.0);
      d2 += d * d;
    }
    bonds_ok = bonds_ok && d2 < 36.0;
    for (const long id : {b[2], b[3]}) {
      if (id > 400) ++crosslinker_bonds[static_cast<std::size_t>(id - 401)];
    }
  }
  CHECK(bonds_ok);
  bool sites_ok = true;
  for (std::size_t k = 0; k < 40; ++k) sites_ok = sites_ok && crosslinker_bonds[k] == net.used_sites[k];
  CHECK(sites_ok);

  std::ostringstream strands;
  write_links(strands, net);
  CHECK(static_cast<std::size_t>(std::ranges::count(strands.str(), '\n')) == 80 + 3);
}

// ------------------------------------------------------------------------------------------------
void test_config() {
  Config cfg;
  CHECK(!cfg.set("procedure", "bonds").has_value());
  CHECK(!cfg.set("chain_lenght", "5").has_value());
  cfg.apply_overrides({"chain_length=8", "n_chains=10", "n_crosslinkers=5", "procedure=Strands", "conversion=0.99", "box=12.5",
                       "seed=99", "length_distribution=gaussian", "dispersity=1.1"});
  cfg.validate();
  CHECK(cfg.linking.procedure == Procedure::Strands && cfg.box_edge(100) == 12.5);
  const auto dir = std::filesystem::temp_directory_path() / "gusev_test_config";
  std::filesystem::create_directories(dir);
  {
    std::ofstream out(dir / "c.ini");
    out << cfg.to_ini();
  }
  CHECK(Config::from_file(dir / "c.ini").to_ini() == cfg.to_ini());
  std::filesystem::remove_all(dir);

  for (const auto& overrides : {std::vector<std::string>{"conversion=0"}, std::vector<std::string>{"functionality=0"},
                                std::vector<std::string>{"density=-1"}, std::vector<std::string>{"chain_length=1"},
                                std::vector<std::string>{"length_distribution=gaussian", "dispersity=3"}}) {
    bool threw = false;
    try {
      Config bad;
      bad.apply_overrides(overrides);
      bad.validate();
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    CHECK(threw);
  }
}

}  // namespace

int main() {
  test_site_finder();
  test_selection_probabilities();
  test_procedures();
  test_placement();
  test_bridges();
  test_data_file();
  test_config();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
