/*!
 * STL C++ v0.1.6
 * https://github.com/ankane/stl-cpp
 * Unlicense OR MIT License
 *
 * Ported from https://www.netlib.org/a/stl
 *
 * Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990).
 * STL: A Seasonal-Trend Decomposition Procedure Based on Loess.
 * Journal of Official Statistics, 6(1), 3-33.
 */

#pragma once
#include "stl.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>

#define FLOAT_EPSILON 0.0001f

namespace mstl {
class MstlResult {
public:
  std::vector<std::vector<float>> seasonal;
  std::vector<float> trend;
  std::vector<float> remainder;
  [[nodiscard]] inline std::vector<float> seasonal_strength() const {
    std::vector<float> strengths(seasonal.size());
    std::transform(seasonal.begin(), seasonal.end(), strengths.begin(),
                   [this](const auto &seasonal) {
                     return stl::_::strength(seasonal, remainder);
                   });
    return strengths;
  }
  [[nodiscard]] inline float trend_strength() const {
    return stl::_::strength(trend, remainder);
  }
};

class MstlParams {
  size_t iterate_;
  std::optional<float> lambda_;
  std::optional<std::vector<size_t>> swin;
  stl::StlParams stl_params_;

public:
  MstlParams()
      : iterate_(2), lambda_(std::nullopt), swin(std::nullopt), stl_params_() {}
  inline MstlParams &iterations(size_t iter) {
    iterate_ = iter;
    return *this;
  }
  inline MstlParams &lambda(float lambda) {
    lambda_ = lambda;
    return *this;
  }
  inline MstlParams &seasonal_lengths(const std::vector<size_t> &lengths) {
    swin = lengths;
    return *this;
  }
  inline MstlParams &stl_params(const stl::StlParams &stl_params) {
    stl_params_ = stl_params;
    return *this;
  }

  MstlResult fit(const std::vector<float> &series,
                 const std::vector<size_t> &periods) const;
};

MstlParams params() { return MstlParams(); }

namespace _ {
std::vector<float> box_cox(const std::vector<float> &y, float lambda) {
  std::vector<float> out(y.size());
  // floating points equivalent of being equal to 0.f
  if (lambda < FLOAT_EPSILON) {
    std::transform(y.begin(), y.end(), out.begin(),
                   [lambda](auto value) { return std::log(value); });
  } else {
    std::transform(y.begin(), y.end(), out.begin(), [lambda](auto value) {
      return (std::pow(value, lambda) - 1.f) / lambda;
    });
  }
  return out;
}
MstlResult mstl(const std::vector<float> &x,
                const std::vector<size_t> &seas_ids, size_t iterate,
                const std::optional<float> lambda,
                const std::optional<std::vector<size_t>> &swin,
                const stl::StlParams &stl_params) {
  const auto k = x.size();
  const auto indices = ([&seas_ids]() {
    std::vector<size_t> indexes(seas_ids.size());
    std::iota(indexes.begin(), indexes.end(), 0);
    std::sort(indexes.begin(), indexes.end(), [&seas_ids](auto a, auto b) {
      return seas_ids[a] < seas_ids[b];
    });
    return indexes;
  })();

  size_t iter = iterate;

  if (seas_ids.size() == 1) {
    iter = 1;
  }

  std::vector<std::vector<float>> seasonality(seas_ids.size());
  std::vector<float> trend;

  auto deseas = lambda.has_value() ? box_cox(x, lambda.value())
                                   : std::vector(x.begin(), x.end());

  if (seas_ids.empty()) {
    // TODO use Friedman's Super Smoother for trend
    throw std::invalid_argument("periods must not be empty");
  }

  for (size_t j = 0; j < iter; ++j) {
    for (size_t i = 0; i < indices.size(); ++i) {
      const auto idx = indices[i];
      if (j > 0) {
        std::transform(deseas.begin(), deseas.end(), seasonality[idx].begin(),
                       deseas.begin(), [](auto d, auto s) { return s + d; });
      }
      const auto fit = ([&stl_params, &swin, idx, &deseas, &seas_ids, i]() {
        if (swin.has_value()) {
          const auto sw = swin.value()[idx];
          stl::StlParams params = stl_params;
          return params.seasonal_length(sw).fit(deseas, seas_ids[idx]);
        }
        if (stl_params.ns.has_value()) {
          return stl_params.fit(deseas, seas_ids[idx]);
        }

        stl::StlParams params = stl_params;
        return params.seasonal_length(7 + 4 * (i + 1))
            .fit(deseas, seas_ids[idx]);
      })();

      seasonality[idx] = fit.seasonal;
      trend = fit.trend;

      std::transform(deseas.begin(), deseas.end(), seasonality[idx].begin(),
                     deseas.begin(), [](auto d, auto s) { return d - s; });
    }
  }

  std::vector<float> remainder(k);
  std::transform(deseas.begin(), deseas.end(), trend.begin(), remainder.begin(),
                 [](auto deseas_part, auto trend_part) {
                   return deseas_part - trend_part;
                 });
  return {seasonality, trend, remainder};
}
} // namespace _
} // namespace mstl

mstl::MstlResult
mstl::MstlParams::fit(const std::vector<float> &series,
                      const std::vector<size_t> &periods) const {
  for (auto period : periods) {
    if (period < 2) {
      throw std::invalid_argument("each period must be at least 2");
    }
    if (series.size() < period * 2) {
      throw std::invalid_argument("series is shorter than twice the period");
    }
  }

  if (lambda_.has_value()) {
    const auto value = lambda_.value();
    if (value > 1.f || value < 0.f) {
      throw std::invalid_argument("lambda must be between 0 and 1");
    }
  }

  if (swin.has_value() && swin.value().size() != periods.size()) {
    throw std::invalid_argument(
        "seasonal_lengths must have the same length as periods");
  }

  return _::mstl(series, periods, iterate_, lambda_, swin, stl_params_);
}
