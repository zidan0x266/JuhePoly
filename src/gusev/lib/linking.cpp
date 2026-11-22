// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "gusev/sim/linking.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "poly/index_set.hpp"

namespace gusev {

namespace {

/// Sites whose Gaussian weight is below e^-kCut of the best candidate's are left out.
constexpr double kCut = 37.0;

[[nodiscard]] long floor_div(long a, long b) noexcept { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }

}  // namespace

std::string_view to_string(Procedure p) noexcept { return p == Procedure::Ends ? "ends" : "strands"; }

std::optional<Procedure> parse_procedure(std::string_view s) noexcept {
  if (s == "ends") return Procedure::Ends;
  if (s == "strands") return Procedure::Strands;
  return std::nullopt;
}

SiteFinder::SiteFinder(const System& system) : system_(&system) {
  const double box = system.period;
  const std::size_t n = system.n_crosslinkers();
  // About two cross-linkers per cell.
  const double spacing = n ? std::cbrt(box * box * box / static_cast<double>(n)) : box;
  cells_ = std::clamp(static_cast<int>(box / (1.26 * spacing)), 1, 128);
  cell_size_ = box / cells_;
  const auto cell_of = [this](const Vec3d& p) {
    const auto index = [this](double x) { return std::min(cells_ - 1, static_cast<int>(x / cell_size_)); };
    return static_cast<std::size_t>((index(p.x) * cells_ + index(p.y)) * cells_ + index(p.z));
  };
  const auto n_cells = static_cast<std::size_t>(cells_) * static_cast<std::size_t>(cells_) * static_cast<std::size_t>(cells_);
  cell_start_.assign(n_cells + 1, 0);
  for (const Vec3d& p : system.crosslinkers) ++cell_start_[cell_of(p) + 1];
  for (std::size_t c = 0; c < n_cells; ++c) cell_start_[c + 1] += cell_start_[c];
  members_.resize(n);
  std::vector<std::uint32_t> next(cell_start_.begin(), cell_start_.end() - 1);
  for (std::uint32_t k = 0; k < n; ++k) members_[next[cell_of(system.crosslinkers[k])]++] = k;
}

void SiteFinder::within(const Vec3d& point, double radius, std::vector<Candidate>& out) const {
  const System& system = *system_;
  const double r2 = radius * radius;
  const auto first = [&](double x) { return static_cast<long>(std::floor((x - radius) / cell_size_)); };
  const auto last = [&](double x) { return static_cast<long>(std::floor((x + radius) / cell_size_)); };
  const long n = cells_;
  // Cell indices outside [0, n) are periodic images: wrap them and shift by whole boxes.
  for (long i = first(point.x); i <= last(point.x); ++i) {
    const long ix = floor_div(i, n);
    const double sx = static_cast<double>(ix) * system.period;
    for (long j = first(point.y); j <= last(point.y); ++j) {
      const long iy = floor_div(j, n);
      const double sy = static_cast<double>(iy) * system.period;
      for (long l = first(point.z); l <= last(point.z); ++l) {
        const long iz = floor_div(l, n);
        const Vec3d shift{sx, sy, static_cast<double>(iz) * system.period};
        const Vec3i image{static_cast<int>(ix), static_cast<int>(iy), static_cast<int>(iz)};
        const auto c = static_cast<std::size_t>(((i - ix * n) * n + (j - iy * n)) * n + (l - iz * n));
        for (std::uint32_t m = cell_start_[c]; m < cell_start_[c + 1]; ++m) {
          const std::uint32_t k = members_[m];
          if (system.free_sites(k) == 0) continue;
          const Vec3d d = system.crosslinkers[k] + shift - point;
          if (d.norm2() <= r2) out.push_back({k, d, image});
        }
      }
    }
  }
}

LinkReport link_strands(System& system, const LinkRules& rules, double bond_length, poly::Rng& rng) {
  const auto f = static_cast<std::uint32_t>(system.functionality);
  const auto n_crosslinkers = static_cast<std::uint32_t>(system.n_crosslinkers());
  // Site j of cross-linker k is k f + j; its first used_sites[k] sites are taken.
  poly::IndexSet<std::uint32_t> sites(static_cast<std::size_t>(n_crosslinkers) * f);
  for (std::uint32_t k = 0; k < n_crosslinkers; ++k) {
    for (std::uint32_t j = system.used_sites[k]; j < f; ++j) sites.insert(k * f + j);
  }
  const SiteFinder finder(system);
  const double box = system.period;
  const double spacing = n_crosslinkers ? std::cbrt(box * box * box / n_crosslinkers) : box;
  std::vector<SiteFinder::Candidate> found;
  LinkReport report;

  // Link end `end` of chain `c` to cross-linker `k` through `image`, and keep the chain seen from
  // end 0 once both ends are linked, or from its only linked end.
  const auto link = [&](std::size_t c, int end, std::uint32_t k, Vec3i image) {
    sites.erase(k * f + system.used_sites[k]);
    system.attach(c, end, k, image);
    ++report.reactions;
    Link& end0 = system.link(c, 0);
    Link& end1 = system.link(c, 1);
    if (end0.linked() && end1.linked()) {
      end1.image = end1.image - end0.image;
      end0.image = {};
    } else {
      system.link(c, end).image = {};
    }
  };
  const auto random_site = [&] { return sites[rng.below(static_cast<std::uint32_t>(sites.size()))] / f; };
  // An active site for the free end of a chain of `segments` segments whose other end is at x:
  // the images within the contour length, each weighted by its active sites times the Gaussian
  // factor. The search radius grows until every image left out weighs less than e^-kCut of the
  // closest one, which keeps the choice exact also for chains forced to reach far.
  const auto site_near = [&](const Vec3d& x, std::size_t segments) -> std::optional<SiteFinder::Candidate> {
    const double s2 = static_cast<double>(segments) * bond_length * bond_length / 3.0;  // per component
    const double cutoff = static_cast<double>(segments) * bond_length;
    double radius = std::min(cutoff, std::sqrt(2 * kCut * s2) + spacing);
    double closest = 0.0;
    for (;;) {
      found.clear();
      finder.within(x, radius, found);
      if (found.empty()) {
        if (radius >= cutoff) return std::nullopt;
        radius = cutoff;
        continue;
      }
      closest = std::numeric_limits<double>::infinity();
      for (const auto& c : found) closest = std::min(closest, c.displacement.norm2());
      const double needed = closest + 2 * kCut * s2;
      if (radius >= cutoff || radius * radius >= needed) break;
      radius = std::min(cutoff, std::sqrt(needed));
    }
    double total = 0.0;
    std::vector<double> weight(found.size());
    for (std::size_t i = 0; i < found.size(); ++i) {
      weight[i] = system.free_sites(found[i].crosslinker) * std::exp(-(found[i].displacement.norm2() - closest) / (2 * s2));
      total += weight[i];
    }
    double u = total * rng.uniform01();
    std::size_t i = 0;
    while (i + 1 < found.size() && u >= weight[i]) u -= weight[i++];
    return found[i];
  };
  const auto segments = [&](std::size_t c) { return system.chain_length(c) + 1; };

  if (rules.procedure == Procedure::Ends) {
    poly::IndexSet<std::uint32_t> ends(2 * system.n_chains());
    std::size_t linked = 0;
    for (std::size_t c = 0; c < system.n_chains(); ++c) {
      for (int end = 0; end < 2; ++end) {
        if (system.linked(c, end)) {
          ++linked;
        } else {
          ends.insert(static_cast<std::uint32_t>(2 * c) + static_cast<std::uint32_t>(end));
        }
      }
    }
    const auto target = static_cast<std::size_t>(std::llround(rules.conversion * static_cast<double>(system.minority_groups())));
    while (linked < target && !ends.empty() && !sites.empty()) {
      const std::uint32_t e = ends[rng.below(static_cast<std::uint32_t>(ends.size()))];
      ends.erase(e);
      const std::size_t c = e / 2;
      const int end = static_cast<int>(e % 2);
      const Link other = system.link(c, 1 - end);
      if (!other.linked()) {
        link(c, end, random_site(), {});
        ++linked;
        continue;
      }
      const auto site = site_near(system.crosslinkers[other.crosslinker], segments(c));
      if (!site) {
        ++report.unreachable;
        continue;
      }
      link(c, end, site->crosslinker, other.image + site->image);
      ++linked;
    }
  } else {
    for (std::size_t c = 0; c < system.n_chains() && !sites.empty(); ++c) {
      const std::uint32_t first = random_site();
      if (rng.uniform01() < rules.conversion) link(c, 0, first, {});
      const auto site = site_near(system.crosslinkers[first], segments(c));
      if (!site) {
        ++report.unreachable;
        continue;
      }
      if (rng.uniform01() < rules.conversion) link(c, 1, site->crosslinker, site->image);
    }
  }
  return report;
}

}  // namespace gusev
