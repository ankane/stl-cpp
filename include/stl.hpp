/*
 * STL C++ v0.2.0
 * https://github.com/ankane/stl-cpp
 * Unlicense OR MIT License
 *
 * Ported from https://www.netlib.org/a/stl
 *
 * Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990).
 * STL: A Seasonal-Trend Decomposition Procedure Based on Loess.
 * Journal of Official Statistics, 6(1), 3-33.
 *
 * Bandara, K., Hyndman, R. J., & Bergmeir, C. (2021).
 * MSTL: A Seasonal-Trend Decomposition Algorithm for Time Series with Multiple Seasonal Patterns.
 * arXiv:2107.13462 [stat.AP]. https://doi.org/10.48550/arXiv.2107.13462
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace stl {

namespace detail {

template<typename T>
bool est(
    const std::vector<T>& y,
    size_t n,
    size_t len,
    int ideg,
    T xs,
    T* ys,
    size_t nleft,
    size_t nright,
    std::vector<T>& w,
    bool userw,
    const std::vector<T>& rw
) {
    T range = static_cast<T>(n) - static_cast<T>(1.0);
    T h = std::max(xs - static_cast<T>(nleft), static_cast<T>(nright) - xs);

    if (len > n) {
        h += static_cast<T>((len - n) / 2);
    }

    T h9 = static_cast<T>(0.999) * h;
    T h1 = static_cast<T>(0.001) * h;

    // compute weights
    T a = 0.0;
    for (size_t j = nleft; j <= nright; j++) {
        w.at(j - 1) = 0.0;
        T r = std::abs(static_cast<T>(j) - xs);
        if (r <= h9) {
            if (r <= h1) {
                w.at(j - 1) = 1.0;
            } else {
                w.at(j - 1) = static_cast<T>(std::pow(1.0 - std::pow(r / h, 3.0), 3.0));
            }
            if (userw) {
                w.at(j - 1) *= rw.at(j - 1);
            }
            a += w.at(j - 1);
        }
    }

    if (a <= 0.0) {
        return false;
    } else { // weighted least squares
        for (size_t j = nleft; j <= nright; j++) { // make sum of w(j) == 1
            w.at(j - 1) /= a;
        }

        if (h > 0.0 && ideg > 0) { // use linear fit
            T a = 0.0;
            for (size_t j = nleft; j <= nright; j++) { // weighted center of x values
                a += w.at(j - 1) * static_cast<T>(j);
            }
            T b = xs - a;
            T c = 0.0;
            for (size_t j = nleft; j <= nright; j++) {
                c += w.at(j - 1) * std::pow(static_cast<T>(j) - a, static_cast<T>(2.0));
            }
            if (std::sqrt(c) > 0.001 * range) {
                b /= c;

                // points are spread out enough to compute slope
                for (size_t j = nleft; j <= nright; j++) {
                    w.at(j - 1) *= b * (static_cast<T>(j) - a) + static_cast<T>(1.0);
                }
            }
        }

        *ys = 0.0;
        for (size_t j = nleft; j <= nright; j++) {
            *ys += w.at(j - 1) * y.at(j - 1);
        }

        return true;
    }
}

template<typename T>
void ess(
    const std::vector<T>& y,
    size_t n,
    size_t len,
    int ideg,
    size_t njump,
    bool userw,
    const std::vector<T>& rw,
    std::span<T> ys,
    std::vector<T>& res
) {
    if (n < 2) {
        ys[0] = y.at(0);
        return;
    }

    size_t nleft = 0;
    size_t nright = 0;

    size_t newnj = std::min(njump, n - 1);
    if (len >= n) {
        nleft = 1;
        nright = n;
        for (size_t i = 1; i <= n; i += newnj) {
            bool ok = est(y, n, len, ideg, static_cast<T>(i), &ys[i - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[i - 1] = y.at(i - 1);
            }
        }
    } else if (newnj == 1) { // newnj equal to one, len less than n
        size_t nsh = (len + 1) / 2;
        nleft = 1;
        nright = len;
        for (size_t i = 1; i <= n; i++) { // fitted value at i
            if (i > nsh && nright != n) {
                nleft += 1;
                nright += 1;
            }
            bool ok = est(y, n, len, ideg, static_cast<T>(i), &ys[i - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[i - 1] = y.at(i - 1);
            }
        }
    } else { // newnj greater than one, len less than n
        size_t nsh = (len + 1) / 2;
        for (size_t i = 1; i <= n; i += newnj) { // fitted value at i
            if (i < nsh) {
                nleft = 1;
                nright = len;
            } else if (i >= n - nsh + 1) {
                nleft = n - len + 1;
                nright = n;
            } else {
                nleft = i - nsh + 1;
                nright = len + i - nsh;
            }
            bool ok = est(y, n, len, ideg, static_cast<T>(i), &ys[i - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[i - 1] = y.at(i - 1);
            }
        }
    }

    if (newnj != 1) {
        for (size_t i = 1; i <= n - newnj; i += newnj) {
            T delta = (ys[i + newnj - 1] - ys[i - 1]) / static_cast<T>(newnj);
            for (size_t j = i + 1; j <= i + newnj - 1; j++) {
                ys[j - 1] = ys[i - 1] + delta * static_cast<T>(j - i);
            }
        }
        size_t k = ((n - 1) / newnj) * newnj + 1;
        if (k != n) {
            bool ok = est(y, n, len, ideg, static_cast<T>(n), &ys[n - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[n - 1] = y.at(n - 1);
            }
            if (k != n - 1) {
                T delta = (ys[n - 1] - ys[k - 1]) / static_cast<T>(n - k);
                for (size_t j = k + 1; j <= n - 1; j++) {
                    ys[j - 1] = ys[k - 1] + delta * static_cast<T>(j - k);
                }
            }
        }
    }
}

template<typename T>
void ma(const std::vector<T>& x, size_t n, size_t len, std::vector<T>& ave) {
    size_t newn = n - len + 1;
    auto flen = static_cast<double>(len);
    double v = 0.0;

    // get the first average
    for (size_t i = 0; i < len; i++) {
        v += x.at(i);
    }

    ave.at(0) = static_cast<T>(v / flen);
    if (newn > 1) {
        size_t k = len;
        size_t m = 0;
        for (size_t j = 1; j < newn; j++) {
            // window down the array
            v = v - x.at(m) + x.at(k);
            ave.at(j) = static_cast<T>(v / flen);
            k += 1;
            m += 1;
        }
    }
}

template<typename T>
void fts(
    const std::vector<T>& x,
    size_t n,
    size_t np,
    std::vector<T>& trend,
    std::vector<T>& work
) {
    ma(x, n, np, trend);
    ma(trend, n - np + 1, np, work);
    ma(work, n - 2 * np + 2, 3, trend);
}

template<typename T>
void rwts(std::span<const T> y, const std::vector<T>& fit, std::vector<T>& rw) {
    size_t n = y.size();

    for (size_t i = 0; i < n; i++) {
        rw.at(i) = std::abs(y[i] - fit.at(i));
    }

    size_t mid1 = (n - 1) / 2;
    size_t mid2 = n / 2;

    // sort
    // TODO checked cast
    std::sort(rw.begin(), rw.begin() + static_cast<std::ptrdiff_t>(n));

    T cmad = static_cast<T>(3.0) * (rw.at(mid1) + rw.at(mid2)); // 6 * median abs resid
    T c9 = static_cast<T>(0.999) * cmad;
    T c1 = static_cast<T>(0.001) * cmad;

    for (size_t i = 0; i < n; i++) {
        T r = std::abs(y[i] - fit.at(i));
        if (r <= c1) {
            rw.at(i) = 1.0;
        } else if (r <= c9) {
            rw.at(i) = static_cast<T>(std::pow(1.0 - std::pow(r / cmad, 2.0), 2.0));
        } else {
            rw.at(i) = 0.0;
        }
    }
}

template<typename T>
void ss(
    const std::vector<T>& y,
    size_t n,
    size_t np,
    size_t ns,
    int isdeg,
    size_t nsjump,
    bool userw,
    std::vector<T>& rw,
    std::vector<T>& season,
    std::vector<T>& work1,
    std::vector<T>& work2,
    std::vector<T>& work3,
    std::vector<T>& work4
) {
    for (size_t j = 1; j <= np; j++) {
        size_t k = (n - j) / np + 1;

        for (size_t i = 1; i <= k; i++) {
            work1.at(i - 1) = y.at((i - 1) * np + j - 1);
        }
        if (userw) {
            for (size_t i = 1; i <= k; i++) {
                work3.at(i - 1) = rw.at((i - 1) * np + j - 1);
            }
        }
        ess(work1, k, ns, isdeg, nsjump, userw, work3, std::span{work2}.subspan(1), work4);
        T xs = 0.0;
        size_t nright = std::min(ns, k);
        bool ok = est(work1, k, ns, isdeg, xs, &work2.at(0), 1, nright, work4, userw, work3);
        if (!ok) {
            work2.at(0) = work2.at(1);
        }
        xs = static_cast<T>(k + 1);
        size_t nleft = static_cast<size_t>(std::max(1, static_cast<int>(k) - static_cast<int>(ns) + 1));
        ok = est(work1, k, ns, isdeg, xs, &work2.at(k + 1), nleft, k, work4, userw, work3);
        if (!ok) {
            work2.at(k + 1) = work2.at(k);
        }
        for (size_t m = 1; m <= k + 2; m++) {
            season.at((m - 1) * np + j - 1) = work2.at(m - 1);
        }
    }
}

template<typename T>
void onestp(
    std::span<const T> y,
    size_t np,
    size_t ns,
    size_t nt,
    size_t nl,
    int isdeg,
    int itdeg,
    int ildeg,
    size_t nsjump,
    size_t ntjump,
    size_t nljump,
    size_t ni,
    bool userw,
    std::vector<T>& rw,
    std::vector<T>& season,
    std::vector<T>& trend,
    std::vector<T>& work1,
    std::vector<T>& work2,
    std::vector<T>& work3,
    std::vector<T>& work4,
    std::vector<T>& work5
) {
    size_t n = y.size();

    for (size_t j = 0; j < ni; j++) {
        for (size_t i = 0; i < n; i++) {
            work1.at(i) = y[i] - trend.at(i);
        }

        ss(work1, n, np, ns, isdeg, nsjump, userw, rw, work2, work3, work4, work5, season);
        fts(work2, n + 2 * np, np, work3, work1);
        ess(work3, n, nl, ildeg, nljump, false, work4, std::span{work1}, work5);
        for (size_t i = 0; i < n; i++) {
            season.at(i) = work2.at(np + i) - work1.at(i);
        }
        for (size_t i = 0; i < n; i++) {
            work1.at(i) = y[i] - season.at(i);
        }
        ess(work1, n, nt, itdeg, ntjump, userw, rw, std::span{trend}, work3);
    }
}

template<typename T>
void stl(
    std::span<const T> y,
    size_t np,
    size_t ns,
    size_t nt,
    size_t nl,
    int isdeg,
    int itdeg,
    int ildeg,
    size_t nsjump,
    size_t ntjump,
    size_t nljump,
    size_t ni,
    size_t no,
    std::vector<T>& rw,
    std::vector<T>& season,
    std::vector<T>& trend
) {
    size_t n = y.size();

    if (ns < 3) {
        throw std::invalid_argument{"seasonal_length must be at least 3"};
    }
    if (nt < 3) {
        throw std::invalid_argument{"trend_length must be at least 3"};
    }
    if (nl < 3) {
        throw std::invalid_argument{"low_pass_length must be at least 3"};
    }
    if (np < 2) {
        throw std::invalid_argument{"period must be at least 2"};
    }

    if (isdeg != 0 && isdeg != 1) {
        throw std::invalid_argument{"seasonal_degree must be 0 or 1"};
    }
    if (itdeg != 0 && itdeg != 1) {
        throw std::invalid_argument{"trend_degree must be 0 or 1"};
    }
    if (ildeg != 0 && ildeg != 1) {
        throw std::invalid_argument{"low_pass_degree must be 0 or 1"};
    }

    if (ns % 2 != 1) {
        throw std::invalid_argument{"seasonal_length must be odd"};
    }
    if (nt % 2 != 1) {
        throw std::invalid_argument{"trend_length must be odd"};
    }
    if (nl % 2 != 1) {
        throw std::invalid_argument{"low_pass_length must be odd"};
    }

    std::vector<T> work1(n + 2 * np);
    std::vector<T> work2(n + 2 * np);
    std::vector<T> work3(n + 2 * np);
    std::vector<T> work4(n + 2 * np);
    std::vector<T> work5(n + 2 * np);

    bool userw = false;
    size_t k = 0;

    while (true) {
        onestp(y, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, work1, work2, work3, work4, work5);
        k += 1;
        if (k > no) {
            break;
        }
        for (size_t i = 0; i < n; i++) {
            work1.at(i) = trend.at(i) + season.at(i);
        }
        rwts(y, work1, rw);
        userw = true;
    }

    if (no <= 0) {
        for (size_t i = 0; i < n; i++) {
            rw.at(i) = 1.0;
        }
    }
}

template<typename T>
double var(const std::vector<T>& series) {
    double mean = std::accumulate(series.begin(), series.end(), 0.0) / static_cast<double>(series.size());
    double sum = 0.0;
    for (auto v : series) {
        double diff = v - mean;
        sum += diff * diff;
    }
    return sum / static_cast<double>(series.size() - 1);
}

template<typename T>
double strength(const std::vector<T>& component, const std::vector<T>& remainder) {
    std::vector<T> sr;
    sr.reserve(remainder.size());
    for (size_t i = 0; i < remainder.size(); i++) {
        sr.push_back(component.at(i) + remainder.at(i));
    }
    return std::max(0.0, 1.0 - var(remainder) / var(sr));
}

} // namespace detail

/// A STL result.
template<typename T = float>
class StlResult {
  public:
    /// Returns the seasonal component.
    std::vector<T> seasonal;

    /// Returns the trend component.
    std::vector<T> trend;

    /// Returns the remainder.
    std::vector<T> remainder;

    /// Returns the weights.
    std::vector<T> weights;

    /// Returns the seasonal strength.
    double seasonal_strength() const {
        return detail::strength(seasonal, remainder);
    }

    /// Returns the trend strength.
    double trend_strength() const {
        return detail::strength(trend, remainder);
    }
};

/// A set of STL parameters.
class StlParams {
  public:
    /// @private
    std::optional<size_t> ns_ = std::nullopt;

  private:
    std::optional<size_t> nt_ = std::nullopt;
    std::optional<size_t> nl_ = std::nullopt;
    int isdeg_ = 0;
    int itdeg_ = 1;
    std::optional<int> ildeg_ = std::nullopt;
    std::optional<size_t> nsjump_ = std::nullopt;
    std::optional<size_t> ntjump_ = std::nullopt;
    std::optional<size_t> nljump_ = std::nullopt;
    std::optional<size_t> ni_ = std::nullopt;
    std::optional<size_t> no_ = std::nullopt;
    bool robust_ = false;

  public:
    /// Sets the length of the seasonal smoother.
    StlParams seasonal_length(size_t length) {
        this->ns_ = length;
        return *this;
    }

    /// Sets the length of the trend smoother.
    StlParams trend_length(size_t length) {
        this->nt_ = length;
        return *this;
    }

    /// Sets the length of the low-pass filter.
    StlParams low_pass_length(size_t length) {
        this->nl_ = length;
        return *this;
    }

    /// Sets the degree of locally-fitted polynomial in seasonal smoothing.
    StlParams seasonal_degree(int degree) {
        this->isdeg_ = degree;
        return *this;
    }

    /// Sets the degree of locally-fitted polynomial in trend smoothing.
    StlParams trend_degree(int degree) {
        this->itdeg_ = degree;
        return *this;
    }

    /// Sets the degree of locally-fitted polynomial in low-pass smoothing.
    StlParams low_pass_degree(int degree) {
        this->ildeg_ = degree;
        return *this;
    }

    /// Sets the skipping value for seasonal smoothing.
    StlParams seasonal_jump(size_t jump) {
        this->nsjump_ = jump;
        return *this;
    }

    /// Sets the skipping value for trend smoothing.
    StlParams trend_jump(size_t jump) {
        this->ntjump_ = jump;
        return *this;
    }

    /// Sets the skipping value for low-pass smoothing.
    StlParams low_pass_jump(size_t jump) {
        this->nljump_ = jump;
        return *this;
    }

    /// Sets the number of loops for updating the seasonal and trend components.
    StlParams inner_loops(size_t loops) {
        this->ni_ = loops;
        return *this;
    }

    /// Sets the number of iterations of robust fitting.
    StlParams outer_loops(size_t loops) {
        this->no_ = loops;
        return *this;
    }

    /// Sets whether robustness iterations are to be used.
    StlParams robust(bool robust) {
        this->robust_ = robust;
        return *this;
    }

    /// Decomposes a time series from a vector.
    template<typename T>
    StlResult<T> fit(const std::vector<T>& series, size_t period) const;

    /// Decomposes a time series from a span.
    template<typename T>
    StlResult<T> fit(std::span<const T> series, size_t period) const;
};

/// Creates a new set of STL parameters.
inline StlParams params() {
    return StlParams{};
}

template<typename T>
StlResult<T> StlParams::fit(std::span<const T> series, size_t period) const {
    std::span<const T> y = series;
    size_t np = period;
    size_t n = series.size();

    if (n < 2 * np) {
        throw std::invalid_argument{"series has less than two periods"};
    }

    size_t ns = this->ns_.value_or(np);

    int isdeg = this->isdeg_;
    int itdeg = this->itdeg_;

    StlResult<T> res{
        std::vector<T>(n),
        std::vector<T>(n),
        std::vector<T>(),
        std::vector<T>(n)
    };

    int ildeg = this->ildeg_.value_or(itdeg);
    size_t newns = std::max(ns, static_cast<size_t>(3));
    if (newns % 2 == 0) {
        newns += 1;
    }

    size_t newnp = std::max(np, static_cast<size_t>(2));
    auto nt = static_cast<size_t>(std::ceil((1.5 * static_cast<float>(newnp)) / (1.0 - 1.5 / static_cast<float>(newns))));
    nt = this->nt_.value_or(nt);
    nt = std::max(nt, static_cast<size_t>(3));
    if (nt % 2 == 0) {
        nt += 1;
    }

    size_t nl = this->nl_.value_or(newnp);
    if (nl % 2 == 0 && !this->nl_.has_value()) {
        nl += 1;
    }

    size_t ni = this->ni_.value_or(this->robust_ ? 1 : 2);
    size_t no = this->no_.value_or(this->robust_ ? 15 : 0);

    size_t nsjump = this->nsjump_.value_or(static_cast<size_t>(std::ceil(static_cast<float>(newns) / 10.0)));
    size_t ntjump = this->ntjump_.value_or(static_cast<size_t>(std::ceil(static_cast<float>(nt) / 10.0)));
    size_t nljump = this->nljump_.value_or(static_cast<size_t>(std::ceil(static_cast<float>(nl) / 10.0)));

    detail::stl(y, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, no, res.weights, res.seasonal, res.trend);

    res.remainder.reserve(n);
    for (size_t i = 0; i < n; i++) {
        res.remainder.push_back(y[i] - res.seasonal.at(i) - res.trend.at(i));
    }

    return res;
}

template<typename T>
StlResult<T> StlParams::fit(const std::vector<T>& series, size_t period) const {
    return StlParams::fit(std::span{series}, period);
}

/// A MSTL result.
template<typename T = float>
class MstlResult {
  public:
    /// Returns the seasonal component.
    std::vector<std::vector<T>> seasonal;

    /// Returns the trend component.
    std::vector<T> trend;

    /// Returns the remainder.
    std::vector<T> remainder;

    /// Returns the seasonal strength.
    std::vector<double> seasonal_strength() const {
        std::vector<double> res;
        for (const auto& s : seasonal) {
            res.push_back(detail::strength(s, remainder));
        }
        return res;
    }

    /// Returns the trend strength.
    double trend_strength() const {
        return detail::strength(trend, remainder);
    }

  private:
    MstlResult(std::vector<std::vector<T>>&& seasonal, std::vector<T>&& trend, std::vector<T>&& remainder) : seasonal{std::move(seasonal)}, trend{std::move(trend)}, remainder{std::move(remainder)} { }

    friend class MstlParams;
};

/// A set of MSTL parameters.
class MstlParams {
    size_t iterate_ = 2;
    std::optional<float> lambda_ = std::nullopt;
    std::optional<std::vector<size_t>> swin_ = std::nullopt;
    StlParams stl_params_;

  public:
    /// Sets the number of iterations.
    MstlParams iterations(size_t iterations) {
        this->iterate_ = iterations;
        return *this;
    }

    /// Sets lambda for Box-Cox transformation.
    MstlParams lambda(float lambda) {
        this->lambda_ = lambda;
        return *this;
    }

    /// Sets the lengths of the seasonal smoothers.
    MstlParams seasonal_lengths(const std::vector<size_t>& lengths) {
        this->swin_ = lengths;
        return *this;
    }

    /// Sets the STL parameters.
    MstlParams stl_params(const StlParams& stl_params) {
        this->stl_params_ = stl_params;
        return *this;
    }

    /// Decomposes a time series from a vector.
    template<typename T>
    MstlResult<T> fit(const std::vector<T>& series, const std::vector<size_t>& periods) const;

    /// Decomposes a time series from a span.
    template<typename T>
    MstlResult<T> fit(std::span<const T> series, std::span<const size_t> periods) const;
};

/// Creates a new set of MSTL parameters.
inline MstlParams mstl_params() {
    return MstlParams{};
}

namespace detail {

template<typename T>
std::vector<T> box_cox(std::span<const T> y, float lambda) {
    std::vector<T> res;
    res.reserve(y.size());
    if (lambda != 0.0) {
        for (size_t i = 0; i < y.size(); i++) {
            res.push_back(static_cast<T>(std::pow(y[i], lambda) - 1.0) / lambda);
        }
    } else {
        for (size_t i = 0; i < y.size(); i++) {
            res.push_back(std::log(y[i]));
        }
    }
    return res;
}

template<typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<std::vector<T>>> mstl(
    std::span<const T> x,
    const std::vector<size_t>& seas_ids,
    size_t iterate,
    std::optional<float> lambda,
    const std::optional<std::vector<size_t>>& swin,
    const StlParams& stl_params
) {
    // keep track of indices instead of sorting seas_ids
    // so order is preserved with seasonality
    std::vector<size_t> indices;
    indices.reserve(seas_ids.size());
    for (size_t i = 0; i < seas_ids.size(); i++) {
        indices.push_back(i);
    }
    std::sort(indices.begin(), indices.end(), [&seas_ids](size_t a, size_t b) {
        return seas_ids.at(a) < seas_ids.at(b);
    });

    if (seas_ids.size() == 1) {
        iterate = 1;
    }

    std::vector<std::vector<T>> seasonality;
    seasonality.reserve(seas_ids.size());
    std::vector<T> trend;

    std::vector<T> deseas = lambda.has_value() ? box_cox(x, lambda.value()) : std::vector<T>(x.begin(), x.end());

    if (!seas_ids.empty()) {
        for (size_t i = 0; i < seas_ids.size(); i++) {
            seasonality.push_back(std::vector<T>());
        }

        for (size_t j = 0; j < iterate; j++) {
            for (size_t i = 0; i < indices.size(); i++) {
                size_t idx = indices.at(i);

                if (j > 0) {
                    for (size_t ii = 0; ii < deseas.size(); ii++) {
                        deseas.at(ii) += seasonality.at(idx).at(ii);
                    }
                }

                StlResult<T> fit;
                if (swin) {
                    StlParams clone = stl_params;
                    fit = clone.seasonal_length(swin.value().at(idx)).fit(deseas, seas_ids.at(idx));
                } else if (stl_params.ns_.has_value()) {
                    fit = stl_params.fit(deseas, seas_ids.at(idx));
                } else {
                    StlParams clone = stl_params;
                    fit = clone.seasonal_length(7 + 4 * (i + 1)).fit(deseas, seas_ids.at(idx));
                }

                seasonality.at(idx) = fit.seasonal;
                trend = fit.trend;

                for (size_t ii = 0; ii < deseas.size(); ii++) {
                    deseas.at(ii) -= seasonality.at(idx).at(ii);
                }
            }
        }
    } else {
        // TODO use Friedman's Super Smoother for trend
        throw std::invalid_argument{"periods must not be empty"};
    }

    std::vector<T> remainder;
    remainder.reserve(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        remainder.push_back(deseas.at(i) - trend.at(i));
    }

    return std::make_tuple(trend, remainder, seasonality);
}

} // namespace detail

template<typename T>
MstlResult<T> MstlParams::fit(std::span<const T> series, std::span<const size_t> periods) const {
    // return error to be consistent with stl
    // and ensure seasonal is always same length as periods
    for (auto v : periods) {
        if (v < 2) {
            throw std::invalid_argument{"periods must be at least 2"};
        }
    }

    // return error to be consistent with stl
    // and ensure seasonal is always same length as periods
    for (auto v : periods) {
        if (series.size() < v * 2) {
            throw std::invalid_argument{"series has less than two periods"};
        }
    }

    if (lambda_.has_value()) {
        float lambda = lambda_.value();
        if (lambda < 0 || lambda > 1) {
            throw std::invalid_argument{"lambda must be between 0 and 1"};
        }
    }

    if (swin_.has_value()) {
        const std::vector<size_t>& swin = swin_.value();
        if (swin.size() != periods.size()) {
            throw std::invalid_argument{"seasonal_lengths must have the same length as periods"};
        }
    }

    auto [trend, remainder, seasonal] = detail::mstl(
        series,
        // copy to support bounds checking before C++26
        std::vector(periods.begin(), periods.end()),
        iterate_,
        lambda_,
        swin_,
        stl_params_
    );

    return MstlResult<T> {
        std::move(seasonal),
        std::move(trend),
        std::move(remainder)
    };
}

template<typename T>
MstlResult<T> MstlParams::fit(
    const std::vector<T>& series,
    const std::vector<size_t>& periods
) const {
    return MstlParams::fit(std::span{series}, std::span{periods});
}

} // namespace stl
