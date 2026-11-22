// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// LAMMPS data files of the modules that write starting structures (randomwalk, builder): text
// output through a buffer and the header of the file, from the title to the Masses section.
#pragma once

#include <charconv>
#include <concepts>
#include <cstddef>
#include <format>
#include <ostream>
#include <span>
#include <string>
#include <string_view>

namespace poly::io {

/// Text output through a fixed-size buffer, with integers formatted by std::to_chars.
class TextBuffer {
 public:
  explicit TextBuffer(std::ostream& out) : out_(out) { buf_.reserve(kChunk + 256); }
  ~TextBuffer() { flush(); }
  TextBuffer(const TextBuffer&) = delete;
  TextBuffer& operator=(const TextBuffer&) = delete;

  TextBuffer& operator<<(std::string_view s) {
    buf_.append(s);
    if (buf_.size() >= kChunk) flush();
    return *this;
  }
  TextBuffer& operator<<(char c) {
    buf_.push_back(c);
    return *this;
  }
  TextBuffer& operator<<(std::integral auto v) {
    char tmp[24];
    const auto end = std::to_chars(tmp, tmp + sizeof tmp, v).ptr;
    return *this << std::string_view(tmp, static_cast<std::size_t>(end - tmp));
  }
  void flush() {
    out_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
    buf_.clear();
  }

 private:
  static constexpr std::size_t kChunk = std::size_t{1} << 20;
  std::ostream& out_;
  std::string buf_;
};

/// Floor division, for the image flags of unwrapped coordinates.
[[nodiscard]] constexpr int floor_div(int a, int b) noexcept { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }

struct LammpsHeader {
  std::string_view title;                        ///< first (comment) line
  std::size_t atoms{0};
  std::size_t bonds{0};
  std::span<const std::string_view> atom_types;  ///< labels, by type - 1
  int bond_types{1};
  double box{0.0};   ///< edge of the cubic box [0, box)
  double mass{1.0};  ///< mass of every atom type
};

/// Everything before the Atoms section: counts, box, atom type labels and masses.
inline void write_lammps_header(TextBuffer& out, const LammpsHeader& h) {
  out << h.title << "\n\n";
  out << h.atoms << " atoms\n" << h.atom_types.size() << " atom types\n";
  out << h.bonds << " bonds\n" << h.bond_types << " bond types\n\n";
  const std::string box = std::format("{:.6f}", h.box);
  for (const char* axis : {"x", "y", "z"}) out << "0 " << box << ' ' << axis << "lo " << axis << "hi\n";

  out << "\nAtom Type Labels\n\n";
  for (std::size_t t = 0; t < h.atom_types.size(); ++t) out << t + 1 << ' ' << h.atom_types[t] << '\n';
  out << "\nMasses\n\n";
  for (std::size_t t = 0; t < h.atom_types.size(); ++t) out << t + 1 << ' ' << std::format("{}", h.mass) << '\n';
}

}  // namespace poly::io
