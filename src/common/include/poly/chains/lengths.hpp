// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Distributions of the precursor chain lengths. A distribution is discrete: a probability P(N)
// for every integer length N between a lower and an optional upper bound. The two-parameter
// shapes (Gaussian, Schulz-Zimm, log-normal, Weibull-Nycander-Gold) are set by the number-average
// chain length Ln and the dispersity Đ = Lw/Ln, the one-parameter shapes (Flory, Poisson) by Ln
// alone.
// The parameters are solved so that this truncated, discrete distribution has exactly the
// requested Ln (and Đ); the textbook relations (such as Đ = 1 + 1/k for Schulz-Zimm) hold only
// for the unbounded shapes and drift for short chains. Lengths are drawn by stratified sampling,
// so that even a few thousand chains follow the distribution closely. The packing modules
// (randomwalk, gusev) share these distributions and their configuration keys.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <filesystem>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "poly/random.hpp"

namespace poly::chains {

enum class LengthDistribution {
  Monodisperse,        ///< every chain has `mean` beads
  Gaussian,            ///< P(N) ∝ exp[-(N - mu)^2 / (2 sigma^2)]
  SchulzZimm,          ///< P(N) ∝ N^(k-1) exp(-N / theta)
  LogNormal,           ///< P(N) ∝ N^-1 exp[-(ln N - mu)^2 / (2 sigma^2)]
  Flory,               ///< P(N) ∝ a^N, the most probable distribution
  Poisson,             ///< P(N) ∝ a^(N-1) / (N-1)!
  WeibullNycanderGold  ///< living polymerization with slow initiation; see weibull_nycander_gold()
};

[[nodiscard]] std::string_view to_string(LengthDistribution d) noexcept;
[[nodiscard]] std::optional<LengthDistribution> parse_length_distribution(std::string_view s) noexcept;
/// Whether the shape has one parameter, so that the number average alone sets it.
[[nodiscard]] constexpr bool one_parameter(LengthDistribution d) noexcept {
  return d == LengthDistribution::Flory || d == LengthDistribution::Poisson;
}

struct LengthSpec {
  LengthDistribution distribution{LengthDistribution::Monodisperse};
  double mean{18.0};          ///< number-average chain length Ln, in beads
  /// Đ = Lw/Ln of a two-parameter shape; 1 gives monodisperse chains whatever the shape. Not
  /// used by the one-parameter shapes, whose dispersity follows from Ln, and left at 1 for them.
  double dispersity{1.0};
  std::size_t min_length{2};  ///< shortest chain
  std::size_t max_length{0};  ///< longest chain; 0: no upper bound

  [[nodiscard]] bool monodisperse() const noexcept {
    return distribution == LengthDistribution::Monodisperse || (dispersity == 1.0 && !one_parameter(distribution));
  }
};

/// Averages of a set of chain lengths.
struct LengthStats {
  std::size_t chains{0};
  std::size_t beads{0};
  std::size_t shortest{0};
  std::size_t longest{0};
  double number_average{0.0};  ///< Ln
  double weight_average{0.0};  ///< Lw
  [[nodiscard]] double dispersity() const noexcept { return number_average > 0 ? weight_average / number_average : 0.0; }
};

[[nodiscard]] LengthStats length_stats(std::span<const std::size_t> lengths) noexcept;

class ChainLengthDistribution {
 public:
  /// Solve the parameters of the shape for the mean and dispersity of `spec`. Throws
  /// std::invalid_argument if the shape cannot reach them within the bounds, naming the range
  /// it can reach.
  explicit ChainLengthDistribution(const LengthSpec& spec);

  [[nodiscard]] const LengthSpec& spec() const noexcept { return spec_; }
  /// Shortest and longest length with a probability above e^-40 of the most probable one.
  [[nodiscard]] std::size_t shortest() const noexcept { return first_; }
  [[nodiscard]] std::size_t longest() const noexcept { return first_ + p_.size() - 1; }
  /// Probability of a chain of `n` beads.
  [[nodiscard]] double probability(std::size_t n) const noexcept {
    return n >= first_ && n - first_ < p_.size() ? p_[n - first_] : 0.0;
  }
  [[nodiscard]] double number_average() const noexcept { return number_average_; }
  [[nodiscard]] double weight_average() const noexcept { return weight_average_; }
  [[nodiscard]] double dispersity() const noexcept { return weight_average_ / number_average_; }
  /// The solved parameters: mu and sigma of a Gaussian, k and theta of a Schulz-Zimm, mu and
  /// sigma of ln N of a log-normal, a of a Flory or Poisson (second zero), and v and r of a
  /// Weibull-Nycander-Gold distribution (as in weibull_nycander_gold()); zero when monodisperse.
  [[nodiscard]] const std::array<double, 2>& parameters() const noexcept { return parameters_; }
  /// The shape and its parameters in words, e.g. "Schulz-Zimm, k = 2, theta = 9".
  [[nodiscard]] std::string describe() const;

  /// `n` chain lengths in random order: chain i takes the length at the quantile (i + u)/n,
  /// u uniform in [0, 1), and the lengths are then shuffled. Monodisperse lengths use no
  /// random numbers.
  [[nodiscard]] std::vector<std::size_t> sample(std::size_t n, Rng& rng) const;

 private:
  /// Tabulate the solved shape (family, width, location) on [lo, hi] and keep its parameters.
  void tabulate(LengthDistribution family, double width, double location, double lo, double hi);

  LengthSpec spec_;
  std::array<double, 2> parameters_{};
  std::size_t first_{0};
  std::vector<double> p_;    ///< P(first_ + i)
  std::vector<double> cdf_;  ///< cumulative sums of p_
  double number_average_{0.0};
  double weight_average_{0.0};
};

/// Weibull, Nycander and Gold's chain-length distribution for a living polymerization with
/// initiation and propagation only (B. Weibull, B. Nycander, Acta Chem. Scand. 1954, 8, 847;
/// L. Gold, J. Chem. Phys. 1958, 28, 91): the number fractions p(k), k = 0..kmax, of initiator
/// molecules that carry k monomer units (k = 0: not initiated), for the number average v over all
/// of them and r = kp/ki. With tau = kp ∫[M]dt, a molecule
/// initiated a time w before the end has grown by a Poisson number of units of mean w, so
/// p(k) = (1/r) ∫_0^tau e^-(tau-w)/r Poisson(k-1; w) dw; p(k) is computed from exact recurrences
/// between neighboring k, each run in its stable direction, and is accurate for every r > 0.
[[nodiscard]] std::vector<double> weibull_nycander_gold(double v, double r, std::size_t kmax);

/// Set one of the configuration keys of the chain lengths: chain_length, length_distribution,
/// dispersity, min_chain_length and max_chain_length. Returns std::nullopt if `key` (already
/// normalized) is none of them.
[[nodiscard]] std::optional<std::expected<void, std::string>> set_length_key(LengthSpec& spec, std::string_view key,
                                                                            std::string_view value);
/// The same keys in INI syntax, one line each.
[[nodiscard]] std::string length_keys_ini(const LengthSpec& spec);
/// `spec` for a system of `n_chains` chains: with max_chain_length 0, no chain is longer than all
/// the chain beads together, n_chains times the number average.
[[nodiscard]] LengthSpec bounded_by_system(LengthSpec spec, std::size_t n_chains);

/// "schulz_zimm lengths of number-average 18 and dispersity 1.5", or for a one-parameter shape
/// "flory lengths of number-average 18".
[[nodiscard]] std::string describe_target(const LengthSpec& spec);

/// For every chain length: the number of chains drawn with it and the number expected from
/// `distribution`, after a header with the averages of both.
void write_length_table(std::ostream& out, const ChainLengthDistribution& distribution, std::span<const std::size_t> lengths);
void write_length_table(const std::filesystem::path& path, const ChainLengthDistribution& distribution,
                        std::span<const std::size_t> lengths);

}  // namespace poly::chains
