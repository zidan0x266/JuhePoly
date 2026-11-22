// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "poly/chains/lengths.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>
#include <utility>

#include "poly/ini.hpp"
#include "poly/io/output_file.hpp"

namespace poly::chains {

namespace {

/// Lengths whose probability is below e^-kCut of the most probable one are left out.
constexpr double kCut = 40.0;
/// Most lengths one distribution may span.
constexpr double kMaxSupport = 1 << 24;
/// Longest chain considered when there is no upper bound.
constexpr double kNoBound = 1e12;

[[noreturn]] void fail(const std::string& message) { throw std::invalid_argument("invalid configuration: " + message); }

/// Poisson probability e^-tau tau^j / j!.
[[nodiscard]] double poisson(double tau, std::size_t j) noexcept {
  const auto x = static_cast<double>(j);
  return std::exp(x * std::log(tau) - tau - std::lgamma(x + 1));
}

/// Number fractions p(k), k = 0..last, of the Weibull-Nycander-Gold distribution for
/// tau = kp ∫[M]dt and r = kp/ki: p(0) = e^-tau/r and, for k >= 1,
/// p(k) = (1/r) ∫_0^tau e^-(tau-w)/r Poisson(k-1; w) dw. For r >= 1, with c = 1 - 1/r,
/// r p(k) = G(k) = Σ_{j>=k} c^(j-k) π(j) = π(k) + c G(k+1), π(j) = Poisson(j; tau), all terms
/// positive. For r < 1, with b = 1/r - 1 and g(k) = r p(k), integration by parts gives
/// g(k) = (π(k-1) - g(k-1))/b, stable upward for k below b tau and, rearranged,
/// g(k-1) = π(k-1) - b g(k) downward above it; each part runs in its stable direction.
std::vector<double> wng_fractions(double tau, double r, std::size_t last) {
  std::vector<double> p(last + 1, 0.0);
  p[0] = std::exp(-tau / r);
  // Beyond `top`, Poisson(j; tau) and with it every term is negligible.
  const std::size_t top = std::max(last, static_cast<std::size_t>(std::ceil(tau + 12 * std::sqrt(tau) + 45)));
  const auto set = [&](std::size_t k, double g) {
    if (k <= last) p[k] = g / r;
  };
  if (r >= 1) {
    const double c = 1 - 1 / r;
    double g = 0.0;
    for (std::size_t k = top + 1; k-- > 1;) set(k, g = poisson(tau, k) + c * g);
  } else {
    const double b = 1 / r - 1;
    const std::size_t split = std::clamp<std::size_t>(static_cast<std::size_t>(std::ceil(std::min(b * tau, 1e15))), 1, top);
    double g = std::exp(-tau) * -std::expm1(-b * tau) / b;  // g(1)
    set(1, g);
    for (std::size_t k = 2; k <= split; ++k) set(k, g = (poisson(tau, k - 1) - g) / b);
    g = 0.0;
    for (std::size_t k = top; k > split + 1; --k) set(k - 1, g = poisson(tau, k - 1) - b * g);
  }
  return p;
}

/// One member of a family. `width` is sigma (Gaussian, log-normal), k (Schulz-Zimm) or r
/// (Weibull-Nycander-Gold), unused for Flory and Poisson; `location` is mu (Gaussian, log-normal),
/// theta (Schulz-Zimm, Flory: a = e^-1/theta), a (Poisson) or tau (Weibull-Nycander-Gold).
struct Shape {
  LengthDistribution family;
  double width;
  double location;
};

/// ln P(n) - ln P(n0) of a shape, written so that it keeps its precision far from n0.
class LogRatio {
 public:
  LogRatio(const Shape& f, double n0)
      : f_(f), n0_(n0), ln_n0_(std::log(n0)), c_(1 / (2 * f.width * f.width)), lgamma_n0_(std::lgamma(n0)) {}

  double operator()(double n) const {
    if (n == n0_) return 0.0;
    switch (f_.family) {
      case LengthDistribution::Gaussian: return -(n - n0_) * (n + n0_ - 2 * f_.location) * c_;
      case LengthDistribution::SchulzZimm: return (f_.width - 1) * (std::log(n) - ln_n0_) - (n - n0_) / f_.location;
      case LengthDistribution::LogNormal: {
        const double ln_n = std::log(n);
        const double d = ln_n - ln_n0_;
        return -d - d * (ln_n + ln_n0_ - 2 * f_.location) * c_;
      }
      case LengthDistribution::Flory: return -(n - n0_) / f_.location;
      case LengthDistribution::Poisson: return (n - n0_) * std::log(f_.location) - (std::lgamma(n) - lgamma_n0_);
      case LengthDistribution::WeibullNycanderGold:  // tabulated, see scan()
      case LengthDistribution::Monodisperse: break;
    }
    return -std::numeric_limits<double>::infinity();
  }

 private:
  Shape f_;
  double n0_;
  double ln_n0_;
  double c_;          ///< 1 / (2 width^2)
  double lgamma_n0_;  ///< ln Γ(n0)
};

/// Most probable length of the continuous shape.
double mode(const Shape& f) {
  switch (f.family) {
    case LengthDistribution::Gaussian: return f.location;
    case LengthDistribution::SchulzZimm: return f.width > 1 ? (f.width - 1) * f.location : 0.0;
    case LengthDistribution::LogNormal: return std::exp(f.location - f.width * f.width);
    case LengthDistribution::Flory: return 0.0;
    case LengthDistribution::Poisson: return f.location + 1;
    case LengthDistribution::WeibullNycanderGold:
    case LengthDistribution::Monodisperse: break;
  }
  return f.location;
}

/// scan() for the Weibull-Nycander-Gold shape, whose weights are tabulated by recurrence.
template <class Visit>
bool scan_tabulated(const Shape& f, double lo, double hi, Visit&& visit) {
  const double tau = f.location;
  const double last = std::min(hi, std::max(std::ceil(tau + 12 * std::sqrt(tau) + 45), lo + 50));
  if (last - lo > kMaxSupport) return false;
  const std::vector<double> p = wng_fractions(tau, f.width, static_cast<std::size_t>(last));
  const auto first = static_cast<std::size_t>(lo);
  std::size_t peak = first;
  for (std::size_t n = first; n < p.size(); ++n) peak = p[n] > p[peak] ? n : peak;
  if (!(p[peak] > 0)) {  // every length from lo on is negligible: all weight at lo
    visit(lo, 1.0);
    return true;
  }
  const double cut = p[peak] * std::exp(-kCut);
  std::size_t begin = peak, end = peak + 1;
  while (begin > first && p[begin - 1] >= cut) --begin;
  while (end < p.size() && p[end] >= cut) ++end;
  for (std::size_t n = begin; n < end; ++n) visit(static_cast<double>(n), p[n] / p[peak]);
  return true;
}

/// Call visit(n, w) in increasing n for the lengths in [lo, hi] whose weight w, relative to the
/// most probable length, is above e^-kCut. All shapes are unimodal, so the scan walks out from
/// the most probable length and stops on either side once the weight has fallen below the cut.
/// Returns false, having stopped early, if the lengths span more than kMaxSupport values.
template <class Visit>
bool scan(const Shape& f, double lo, double hi, Visit&& visit) {
  if (f.family == LengthDistribution::WeibullNycanderGold) return scan_tabulated(f, lo, hi, visit);
  const double m = mode(f);
  double peak = std::floor(std::isnan(m) ? lo : std::clamp(m, lo, hi));
  if (peak + 1 <= hi && LogRatio(f, peak)(peak + 1) > 0) peak += 1;
  const LogRatio ratio(f, peak);
  // A few probes at doubling distances tell whether the tails fit, before they are walked.
  for (double span = 1; span <= kMaxSupport; span *= 2) {
    const bool up = peak + span <= hi && ratio(peak + span) >= -kCut;
    const bool down = peak - span >= lo && ratio(peak - span) >= -kCut;
    if (!up && !down) break;
    if (2 * span > kMaxSupport) return false;
  }
  double first = peak;
  while (first - 1 >= lo && ratio(first - 1) >= -kCut) {
    first -= 1;
    if (peak - first > kMaxSupport) return false;
  }
  for (double n = first; n <= hi; n += 1) {
    const double r = ratio(n);
    if (n > peak && r < -kCut) break;
    if (n - first > kMaxSupport) return false;
    visit(n, std::exp(r));
  }
  return true;
}

struct Moments {
  double m0{0.0};
  double m1{0.0};
  double m2{0.0};
  bool complete{true};  ///< false if the shape was too broad to scan
  [[nodiscard]] double mean() const noexcept { return m1 / m0; }
  [[nodiscard]] double dispersity() const noexcept { return m2 * m0 / (m1 * m1); }
};

Moments moments(const Shape& f, double lo, double hi) {
  Moments m;
  m.complete = scan(f, lo, hi, [&m](double n, double w) {
    m.m0 += w;
    m.m1 += n * w;
    m.m2 += n * n * w;
  });
  return m;
}

constexpr double kInfinity = std::numeric_limits<double>::infinity();

/// The root of the increasing function g between a and b, where g(a) < 0 < g(b), by regula
/// falsi with the Illinois modification; stops once |g| <= tol or the interval is exhausted.
template <class G>
double solve_increasing(G&& g, double a, double b, double ga, double gb, double tol) {
  int moved = 0;  // the end moved by the previous step: -1 lower, +1 upper
  double c = 0.5 * (a + b);
  for (int i = 0; i < 400; ++i) {
    c = (a * gb - b * ga) / (gb - ga);
    if (!(c > a && c < b)) c = 0.5 * (a + b);  // also when g is infinite at an end
    if (!(c > a && c < b)) break;
    const double gc = g(c);
    if (std::abs(gc) <= tol) break;
    if (gc < 0) {
      a = c;
      ga = gc;
      if (moved == -1) gb *= 0.5;
      moved = -1;
    } else {
      b = c;
      gb = gc;
      if (moved == +1) ga *= 0.5;
      moved = +1;
    }
  }
  return c;
}

/// The solver works on two unbounded variables: the mean grows with x and the dispersity with
/// y. Gaussian and log-normal: mu = x, sigma = e^y. Schulz-Zimm: theta = e^x, k = e^-y.
/// Weibull-Nycander-Gold: tau = e^x, r = e^y. Flory and Poisson: theta or a = e^x, y unused.
Shape shape(LengthDistribution family, double y, double x) {
  switch (family) {
    case LengthDistribution::SchulzZimm: return {family, std::exp(-y), std::exp(x)};
    case LengthDistribution::Flory:
    case LengthDistribution::Poisson: return {family, 1.0, std::exp(x)};
    case LengthDistribution::WeibullNycanderGold: return {family, std::exp(y), std::exp(x)};
    default: return {family, std::exp(y), x};
  }
}

/// The x at which the shape of width y has the number average `mean`, or nothing if no x has:
/// a broad shape bounded above may not reach the mean, as its weights fall off from the start.
std::optional<double> solve_location(LengthDistribution family, double y, double mean, double lo, double hi) {
  const auto g = [&](double x) {
    const Moments m = moments(shape(family, y, x), lo, hi);
    return m.complete ? m.mean() - mean : kInfinity;
  };
  const double width = shape(family, y, 0.0).width;
  double x0 = mean;  // the unbounded shape has this mean
  double step = width + 1;
  if (family == LengthDistribution::SchulzZimm) {
    x0 = std::log(mean / width);
    step = 1;
  } else if (family == LengthDistribution::LogNormal) {
    x0 = std::log(mean) - width * width / 2;
    step = std::max(width, 0.5);
  } else if (family == LengthDistribution::Flory) {
    x0 = std::log(std::max(mean - lo, 0.5));  // the unbounded shape from lo on
    step = 1;
  } else if (family == LengthDistribution::Poisson || family == LengthDistribution::WeibullNycanderGold) {
    x0 = std::log(std::max(mean - 1, 0.5));
    step = 1;
  }
  double below = step, a = x0 - below, ga = g(a);
  for (int i = 0; ga >= 0; ++i) {
    if (i == 100) throw std::logic_error("chain lengths: no lower bracket for the mean");
    below *= 2;
    a = x0 - below;
    ga = g(a);
  }
  double above = step, b = x0 + above, gb = g(b);
  for (int i = 0; gb <= 0; ++i) {
    if (i == 100) return std::nullopt;
    above *= 2;
    b = x0 + above;
    gb = g(b);
  }
  return solve_increasing(g, a, b, ga, gb, 1e-13 * mean);
}

/// Range of the width variable y searched for each shape: sigma of a Gaussian from 1e-3 to
/// 1000 times the mean, k of a Schulz-Zimm from 1e9 down to 1e-3, sigma of ln N from 1e-4 to 4,
/// r of a Weibull-Nycander-Gold distribution from 1e-4 to 1e6.
std::pair<double, double> width_range(LengthDistribution family, double mean) {
  switch (family) {
    case LengthDistribution::Gaussian: return {std::log(1e-3), std::log(1e3 * mean)};
    case LengthDistribution::SchulzZimm: return {-std::log(1e9), -std::log(1e-3)};
    case LengthDistribution::WeibullNycanderGold: return {std::log(1e-4), std::log(1e6)};
    default: return {std::log(1e-4), std::log(4.0)};
  }
}

/// Width variable y of the unbounded shape with dispersity d.
double unbounded_width(LengthDistribution family, double mean, double d) {
  switch (family) {
    case LengthDistribution::Gaussian: return std::log(mean * std::sqrt(d - 1));  // Đ = 1 + sigma^2 / mean^2
    case LengthDistribution::SchulzZimm: return std::log(d - 1);                  // Đ = 1 + 1/k
    case LengthDistribution::WeibullNycanderGold: return 0.0;                     // r = 1
    default: return std::log(std::sqrt(std::log(d)));                             // Đ = exp(sigma^2)
  }
}

}  // namespace

std::string_view to_string(LengthDistribution d) noexcept {
  switch (d) {
    case LengthDistribution::Monodisperse: return "monodisperse";
    case LengthDistribution::Gaussian: return "gaussian";
    case LengthDistribution::SchulzZimm: return "schulz_zimm";
    case LengthDistribution::LogNormal: return "lognormal";
    case LengthDistribution::Flory: return "flory";
    case LengthDistribution::Poisson: return "poisson";
    case LengthDistribution::WeibullNycanderGold: return "weibull_nycander_gold";
  }
  return "?";
}

std::optional<LengthDistribution> parse_length_distribution(std::string_view s) noexcept {
  for (const auto d : {LengthDistribution::Monodisperse, LengthDistribution::Gaussian, LengthDistribution::SchulzZimm,
                       LengthDistribution::LogNormal, LengthDistribution::Flory, LengthDistribution::Poisson,
                       LengthDistribution::WeibullNycanderGold}) {
    if (s == to_string(d)) return d;
  }
  return std::nullopt;
}

LengthStats length_stats(std::span<const std::size_t> lengths) noexcept {
  LengthStats s;
  s.chains = lengths.size();
  if (lengths.empty()) return s;
  double squares = 0.0;
  s.shortest = lengths.front();
  for (const std::size_t n : lengths) {
    s.beads += n;
    squares += static_cast<double>(n) * static_cast<double>(n);
    s.shortest = std::min(s.shortest, n);
    s.longest = std::max(s.longest, n);
  }
  s.number_average = static_cast<double>(s.beads) / static_cast<double>(s.chains);
  s.weight_average = squares / static_cast<double>(s.beads);
  return s;
}

ChainLengthDistribution::ChainLengthDistribution(const LengthSpec& spec) : spec_(spec) {
  const double lo = static_cast<double>(spec.min_length);
  const double hi = spec.max_length ? static_cast<double>(spec.max_length) : kNoBound;
  if (spec.min_length < 2) fail("min_chain_length must be at least 2");
  if (spec.max_length && spec.max_length < spec.min_length) fail("max_chain_length must be 0 or at least min_chain_length");
  if (!(spec.dispersity >= 1)) fail("dispersity must be at least 1");
  if (spec.distribution == LengthDistribution::Monodisperse && spec.dispersity != 1) {
    fail("a dispersity above 1 needs a length_distribution other than monodisperse");
  }
  const LengthDistribution family = spec.distribution;
  if (one_parameter(family) && spec.dispersity != 1) {
    fail(std::format("a {} distribution has one parameter: chain_length sets its dispersity, so leave dispersity out",
                     to_string(family)));
  }

  if (spec.monodisperse()) {
    if (spec.mean != std::round(spec.mean)) fail("chain_length must be a whole number for monodisperse chains");
    if (spec.mean < lo || spec.mean > hi) fail("chain_length must lie between min_chain_length and max_chain_length");
    first_ = static_cast<std::size_t>(spec.mean);
    p_ = {1.0};
    cdf_ = {1.0};
    number_average_ = weight_average_ = spec.mean;
    return;
  }
  if (!(spec.mean > lo && spec.mean < hi)) {
    fail(std::format("for a {} distribution, chain_length must lie strictly between min_chain_length and max_chain_length",
                     to_string(family)));
  }
  if (one_parameter(family)) {
    const auto x = solve_location(family, 0.0, spec.mean, lo, hi);
    if (!x) fail("min_chain_length to max_chain_length is too narrow a range for this distribution");
    const Shape f = shape(family, 0.0, *x);
    tabulate(f.family, f.width, f.location, lo, hi);
    return;
  }

  // Outer solve on the width y; for each y the location x is solved for the mean.
  const double target = spec.dispersity;
  const auto g = [&](double y) {
    const auto x = solve_location(family, y, spec.mean, lo, hi);
    if (!x) return kInfinity;  // too broad to have this mean
    const Moments m = moments(shape(family, y, *x), lo, hi);
    return m.complete ? m.dispersity() - target : kInfinity;
  };
  const auto [y_min, y_max] = width_range(family, spec.mean);
  const auto unreachable = [&](bool too_broad, double reached) {
    const auto bounds = spec.max_length ? std::format("between {} and {}", spec.min_length, spec.max_length)
                                        : std::format("of at least {}", spec.min_length);
    std::string hint;
    if (too_broad && (family == LengthDistribution::Gaussian || family == LengthDistribution::WeibullNycanderGold)) {
      hint = "; schulz_zimm and lognormal reach wider distributions";
    } else if (!too_broad && spec.mean != std::round(spec.mean)) {
      hint = "; a chain_length that is not a whole number cannot be narrower";
    }
    fail(std::format("a {} distribution of number-average {} with lengths {} reaches dispersities {} {:.6g}, not {}{}",
                     to_string(family), spec.mean, bounds, too_broad ? "up to" : "down to", reached, target, hint));
  };

  double a = std::clamp(unbounded_width(family, spec.mean, target), y_min, y_max);
  double ga = g(a);
  double b = a, gb = ga;
  if (ga < 0) {
    for (double step = 1; gb < 0; step *= 2) {
      if (b == y_max) unreachable(true, gb + target);
      b = std::min(b + step, y_max);
      gb = g(b);
    }
  } else {
    for (double step = 1; ga > 0; step *= 2) {
      if (a == y_min) unreachable(false, ga + target);
      a = std::max(a - step, y_min);
      ga = g(a);
    }
  }
  const double y = ga == 0 ? a : gb == 0 ? b : solve_increasing(g, a, b, ga, gb, 1e-12 * target);

  const auto x = solve_location(family, y, spec.mean, lo, hi);
  if (!x) fail("min_chain_length to max_chain_length is too narrow a range for this distribution");
  const Shape f = shape(family, y, *x);
  tabulate(f.family, f.width, f.location, lo, hi);
  // The dispersity jumps where broader shapes become infeasible; a target beyond the jump
  // leaves the solver at the edge.
  if (std::abs(dispersity() - target) > 1e-9 * target) unreachable(true, dispersity());
}

void ChainLengthDistribution::tabulate(LengthDistribution family, double width, double location, double lo, double hi) {
  const Shape f{family, width, location};
  switch (f.family) {
    case LengthDistribution::SchulzZimm: parameters_ = {f.width, f.location}; break;
    case LengthDistribution::Flory: parameters_ = {std::exp(-1 / f.location), 0.0}; break;
    case LengthDistribution::Poisson: parameters_ = {f.location, 0.0}; break;
    case LengthDistribution::WeibullNycanderGold: {
      const double tau = f.location, r = f.width;  // v counts the molecules not initiated as zero
      parameters_ = {tau - (r - 1) * -std::expm1(-tau / r), r};
      break;
    }
    default: parameters_ = {f.location, f.width};
  }
  double total = 0.0;
  const bool complete = scan(f, lo, hi, [&](double n, double w) {
    if (p_.empty()) first_ = static_cast<std::size_t>(n);
    p_.push_back(w);
    total += w;
  });
  if (!complete) fail("the chain length distribution is too broad to tabulate; set max_chain_length");
  double m1 = 0.0, m2 = 0.0, sum = 0.0;
  cdf_.reserve(p_.size());
  for (std::size_t i = 0; i < p_.size(); ++i) {
    p_[i] /= total;
    const auto n = static_cast<double>(first_ + i);
    m1 += n * p_[i];
    m2 += n * n * p_[i];
    sum += p_[i];
    cdf_.push_back(sum);
  }
  cdf_.back() = 1.0;
  number_average_ = m1;
  weight_average_ = m2 / m1;
  if (std::abs(number_average_ - spec_.mean) > 1e-9 * spec_.mean) {
    throw std::logic_error(std::format("chain lengths: the number average {} misses {}", number_average_, spec_.mean));
  }
}

std::string ChainLengthDistribution::describe() const {
  if (spec_.monodisperse()) return std::format("monodisperse, {} beads", first_);
  switch (spec_.distribution) {
    case LengthDistribution::Gaussian: return std::format("Gaussian, mu = {:.6g}, sigma = {:.6g}", parameters_[0], parameters_[1]);
    case LengthDistribution::SchulzZimm:
      return std::format("Schulz-Zimm, k = {:.6g}, theta = {:.6g}", parameters_[0], parameters_[1]);
    case LengthDistribution::LogNormal:
      return std::format("log-normal, mu = {:.6g}, sigma = {:.6g} of ln N", parameters_[0], parameters_[1]);
    case LengthDistribution::Flory: return std::format("Flory, a = {:.6g}", parameters_[0]);
    case LengthDistribution::Poisson: return std::format("Poisson, a = {:.6g}", parameters_[0]);
    case LengthDistribution::WeibullNycanderGold:
      return std::format("Weibull-Nycander-Gold, v = {:.6g}, r = kp/ki = {:.6g}", parameters_[0], parameters_[1]);
    case LengthDistribution::Monodisperse: break;
  }
  return "?";
}

std::vector<std::size_t> ChainLengthDistribution::sample(std::size_t n, Rng& rng) const {
  std::vector<std::size_t> lengths(n, first_);
  if (p_.size() == 1) return lengths;
  for (std::size_t i = 0; i < n; ++i) {
    const double u = (static_cast<double>(i) + rng.uniform01()) / static_cast<double>(n);
    const auto k = static_cast<std::size_t>(std::ranges::upper_bound(cdf_, u) - cdf_.begin());
    lengths[i] = first_ + std::min(k, p_.size() - 1);
  }
  for (std::size_t i = n; i > 1; --i) std::swap(lengths[i - 1], lengths[rng.below(static_cast<std::uint32_t>(i))]);
  return lengths;
}

std::vector<double> weibull_nycander_gold(double v, double r, std::size_t kmax) {
  if (!(v >= 0) || !(r > 0)) throw std::invalid_argument("Weibull-Nycander-Gold distribution: v must be >= 0 and r > 0");
  // tau from v = tau - (r - 1)(1 - e^-tau/r), which grows with tau
  const auto excess = [&](double tau) { return tau - (r - 1) * -std::expm1(-tau / r) - v; };
  double a = 0.0, b = v + std::max(r, 1.0) + 1;
  for (int i = 0; i < 200 && b - a > 1e-15 * b; ++i) (excess(0.5 * (a + b)) < 0 ? a : b) = 0.5 * (a + b);
  return wng_fractions(0.5 * (a + b), r, kmax);
}

std::optional<std::expected<void, std::string>> set_length_key(LengthSpec& spec, std::string_view key, std::string_view value) {
  using ini::assign;
  using ini::parse_number;
  if (key == "chain_length") return assign(spec.mean, parse_number<double>(key, value));
  if (key == "length_distribution") {
    const auto d = parse_length_distribution(ini::lower(ini::trim(value)));
    if (!d) {
      return std::unexpected(
          "config key 'length_distribution': expected 'monodisperse', 'gaussian', 'schulz_zimm', 'lognormal', 'flory', "
          "'poisson' or 'weibull_nycander_gold'");
    }
    spec.distribution = *d;
    return std::expected<void, std::string>{};
  }
  if (key == "dispersity") return assign(spec.dispersity, parse_number<double>(key, value));
  if (key == "min_chain_length") return assign(spec.min_length, parse_number<std::size_t>(key, value));
  if (key == "max_chain_length") return assign(spec.max_length, parse_number<std::size_t>(key, value));
  return std::nullopt;
}

std::string length_keys_ini(const LengthSpec& spec) {
  return std::format("chain_length = {}\nlength_distribution = {}\ndispersity = {}\nmin_chain_length = {}\nmax_chain_length = {}\n",
                     spec.mean, to_string(spec.distribution), spec.dispersity, spec.min_length, spec.max_length);
}

LengthSpec bounded_by_system(LengthSpec spec, std::size_t n_chains) {
  if (spec.max_length == 0) spec.max_length = static_cast<std::size_t>(std::llround(static_cast<double>(n_chains) * spec.mean));
  return spec;
}

std::string describe_target(const LengthSpec& spec) {
  if (one_parameter(spec.distribution)) return std::format("{} lengths of number-average {}", to_string(spec.distribution), spec.mean);
  return std::format("{} lengths of number-average {} and dispersity {}", to_string(spec.distribution), spec.mean, spec.dispersity);
}

void write_length_table(std::ostream& out, const ChainLengthDistribution& distribution, std::span<const std::size_t> lengths) {
  const LengthSpec& spec = distribution.spec();
  const LengthStats drawn = length_stats(lengths);
  if (one_parameter(spec.distribution)) {
    std::println(out, "# chain lengths: {} distribution, number-average {} ({})", to_string(spec.distribution), spec.mean,
                 distribution.describe());
  } else {
    std::println(out, "# chain lengths: {} distribution, number-average {}, dispersity {} ({})", to_string(spec.distribution),
                 spec.mean, spec.dispersity, distribution.describe());
  }
  std::println(out, "# distribution: Ln {:.6f}, Lw {:.6f}, dispersity {:.6f}, lengths {} to {}", distribution.number_average(),
               distribution.weight_average(), distribution.dispersity(), distribution.shortest(), distribution.longest());
  std::println(out, "# drawn: {} chains, Ln {:.6f}, Lw {:.6f}, dispersity {:.6f}, lengths {} to {}", drawn.chains,
               drawn.number_average, drawn.weight_average, drawn.dispersity(), drawn.shortest, drawn.longest);
  std::println(out, "# length chains expected");
  if (lengths.empty()) return;
  const std::size_t first = std::min(distribution.shortest(), drawn.shortest);
  const std::size_t last = std::max(distribution.longest(), drawn.longest);
  std::vector<std::size_t> count(last - first + 1, 0);
  for (const std::size_t n : lengths) ++count[n - first];
  for (std::size_t n = first; n <= last; ++n) {
    const double expected = static_cast<double>(drawn.chains) * distribution.probability(n);
    if (count[n - first] > 0 || expected >= 1e-3) std::println(out, "{} {} {:.4f}", n, count[n - first], expected);
  }
}

void write_length_table(const std::filesystem::path& path, const ChainLengthDistribution& distribution,
                        std::span<const std::size_t> lengths) {
  auto out = io::open_output(path);
  write_length_table(out, distribution, lengths);
  io::check_written(out, path);
}

}  // namespace poly::chains
