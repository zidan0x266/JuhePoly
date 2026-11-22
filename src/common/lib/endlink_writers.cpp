// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include <cmath>
#include <cstddef>
#include <format>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <string>
#include <vector>

#include "poly/endlink/writers.hpp"
#include "poly/io/lammps_data.hpp"
#include "poly/io/output_file.hpp"

namespace poly::endlink {

namespace {

using io::TextBuffer;

/// Writes positions wrapped into the box, with or without their image flags.
template <class Point>
class Coordinates;

/// Lattice points: every wrapped coordinate is one of `period` multiples of the scale, so each
/// is formatted once.
template <>
class Coordinates<Vec3i> {
 public:
  explicit Coordinates(const System<Vec3i>& system) : period_(system.period), table_(static_cast<std::size_t>(system.period)) {
    for (std::size_t i = 0; i < table_.size(); ++i) table_[i] = std::format("{:.6f}", system.scale * static_cast<double>(i));
  }
  void wrapped(TextBuffer& out, const Vec3i& p) const {
    out << table_[wrap(p.x)] << ' ' << table_[wrap(p.y)] << ' ' << table_[wrap(p.z)];
  }
  void with_images(TextBuffer& out, const Vec3i& p) const {
    wrapped(out, p);
    out << ' ' << io::floor_div(p.x, period_) << ' ' << io::floor_div(p.y, period_) << ' ' << io::floor_div(p.z, period_);
  }

 private:
  [[nodiscard]] std::size_t wrap(int c) const noexcept { return static_cast<std::size_t>(c - io::floor_div(c, period_) * period_); }

  int period_;
  std::vector<std::string> table_;
};

template <>
class Coordinates<Vec3d> {
 public:
  explicit Coordinates(const System<Vec3d>& system) : period_(system.period), scale_(system.scale) {}
  void wrapped(TextBuffer& out, const Vec3d& p) const {
    out << std::format("{:.6f} {:.6f} {:.6f}", scale_ * wrap(p.x).x, scale_ * wrap(p.y).x, scale_ * wrap(p.z).x);
  }
  void with_images(TextBuffer& out, const Vec3d& p) const {
    const Wrapped x = wrap(p.x), y = wrap(p.y), z = wrap(p.z);
    out << std::format("{:.6f} {:.6f} {:.6f}", scale_ * x.x, scale_ * y.x, scale_ * z.x) << ' ' << x.image << ' ' << y.image << ' '
        << z.image;
  }

 private:
  /// A coordinate mapped into [0, period), with the number of periods it was shifted by.
  struct Wrapped {
    double x;
    int image;
  };
  [[nodiscard]] Wrapped wrap(double x) const noexcept {
    double image = std::floor(x / period_);
    double w = x - image * period_;
    if (w >= period_) {
      w -= period_;
      image += 1;
    } else if (w < 0) {
      w += period_;
      image -= 1;
    }
    return {w, static_cast<int>(image)};
  }

  double period_;
  double scale_;
};

}  // namespace

template <class Point>
void write_lammps_data(std::ostream& stream, const System<Point>& system, double mass, std::string_view title) {
  const Coordinates<Point> coordinates(system);
  TextBuffer out(stream);
  io::write_lammps_header(out, {title, system.atoms(), system.bonds(), kTypeLabels, 1, system.box_length(), mass});

  out << "\nAtoms # molecular\n\n";
  std::size_t id = 0;
  const auto atom = [&](const Point& p, std::size_t molecule, AtomType type) {
    out << ++id << ' ' << molecule << ' ' << static_cast<int>(type) << ' ';
    coordinates.with_images(out, p);
    out << '\n';
  };
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    for (std::size_t i = system.chain_start[c]; i < system.chain_start[c + 1]; ++i) atom(system.beads[i], c + 1, system.bead_type(c, i));
  }
  for (std::size_t k = 0; k < system.n_crosslinkers(); ++k) {
    atom(system.crosslinkers[k], system.n_chains() + k + 1, system.crosslinker_type(k));
  }

  out << "\nBonds\n\n";
  std::size_t bond = 0;
  const auto crosslinker_atom = [&](std::uint32_t k) { return system.chain_beads() + k + 1; };
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    // one-based atom IDs: the chain holds atoms chain_start[c] + 1 to chain_start[c + 1]
    const std::size_t first = system.chain_start[c] + 1, last = system.chain_start[c + 1];
    if (system.linked(c, 0)) out << ++bond << " 1 " << crosslinker_atom(system.link(c, 0).crosslinker) << ' ' << first << '\n';
    for (std::size_t i = first; i < last; ++i) out << ++bond << " 1 " << i << ' ' << i + 1 << '\n';
    if (system.linked(c, 1)) out << ++bond << " 1 " << last << ' ' << crosslinker_atom(system.link(c, 1).crosslinker) << '\n';
  }
}

template <class Point>
void write_lammps_data(const std::filesystem::path& path, const System<Point>& system, double mass, std::string_view title) {
  auto out = io::open_output(path, std::ios::binary);
  write_lammps_data(out, system, mass, title);
  io::check_written(out, path);
}

template <class Point>
void write_xyz(std::ostream& stream, const System<Point>& system) {
  const Coordinates<Point> coordinates(system);
  TextBuffer out(stream);
  const std::string box = std::format("{:.6f}", system.box_length());
  out << system.atoms() << '\n';
  out << "Lattice=\"" << box << " 0 0 0 " << box << " 0 0 0 " << box << "\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n";
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    for (std::size_t i = system.chain_start[c]; i < system.chain_start[c + 1]; ++i) {
      out << static_cast<int>(system.bead_type(c, i)) << ' ';
      coordinates.wrapped(out, system.beads[i]);
      out << '\n';
    }
  }
  for (std::size_t k = 0; k < system.n_crosslinkers(); ++k) {
    out << static_cast<int>(system.crosslinker_type(k)) << ' ';
    coordinates.wrapped(out, system.crosslinkers[k]);
    out << '\n';
  }
}

template <class Point>
void write_xyz(const std::filesystem::path& path, const System<Point>& system) {
  auto out = io::open_output(path, std::ios::binary);
  write_xyz(out, system);
  io::check_written(out, path);
}

template <class Point>
void write_links(std::ostream& out, const System<Point>& system) {
  std::println(out, "# chains of the system: number of beads, atom IDs of the cross-linkers of end 0 and end 1");
  std::println(out, "# (0: unlinked end) and, with both ends linked, the span from the first cross-linker to the second");
  std::println(out, "# chain length xlink_0 xlink_1 span_x span_y span_z");
  for (std::size_t c = 0; c < system.n_chains(); ++c) {
    const auto crosslinker_atom = [&](int end) {
      return system.linked(c, end) ? system.chain_beads() + system.link(c, end).crosslinker + 1 : std::size_t{0};
    };
    const Point span = system.linked_ends(c) == 2 ? system.span(c) : Point{};
    std::println(out, "{} {} {} {} {:.6f} {:.6f} {:.6f}", c + 1, system.chain_length(c), crosslinker_atom(0), crosslinker_atom(1),
                 system.scale * static_cast<double>(span.x), system.scale * static_cast<double>(span.y),
                 system.scale * static_cast<double>(span.z));
  }
}

template <class Point>
void write_links(const std::filesystem::path& path, const System<Point>& system) {
  auto out = io::open_output(path);
  write_links(out, system);
  io::check_written(out, path);
}

template void write_lammps_data(std::ostream&, const System<Vec3i>&, double, std::string_view);
template void write_lammps_data(std::ostream&, const System<Vec3d>&, double, std::string_view);
template void write_lammps_data(const std::filesystem::path&, const System<Vec3i>&, double, std::string_view);
template void write_lammps_data(const std::filesystem::path&, const System<Vec3d>&, double, std::string_view);
template void write_xyz(std::ostream&, const System<Vec3i>&);
template void write_xyz(std::ostream&, const System<Vec3d>&);
template void write_xyz(const std::filesystem::path&, const System<Vec3i>&);
template void write_xyz(const std::filesystem::path&, const System<Vec3d>&);
template void write_links(std::ostream&, const System<Vec3i>&);
template void write_links(std::ostream&, const System<Vec3d>&);
template void write_links(const std::filesystem::path&, const System<Vec3i>&);
template void write_links(const std::filesystem::path&, const System<Vec3d>&);

}  // namespace poly::endlink
