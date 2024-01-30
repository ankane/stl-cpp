/*!
 * STL C++ v0.1.4
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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace stl {

template <typename T>
bool est(const T* y, size_t n, size_t len, int ideg, T xs, T* ys, size_t nleft, size_t nright, T* w, bool userw, const T* rw) {
    auto range = ((T) n) - 1.0;
    auto h = std::max(xs - ((T) nleft), ((T) nright) - xs);

    if (len > n) {
        h += (T) ((len - n) / 2);
    }

    auto h9 = 0.999 * h;
    auto h1 = 0.001 * h;

    // compute weights
    auto a = 0.0;
    for (auto j = nleft; j <= nright; j++) {
        w[j - 1] = 0.0;
        auto r = fabs(((T) j) - xs);
        if (r <= h9) {
            if (r <= h1) {
                w[j - 1] = 1.0;
            } else {
                w[j - 1] = pow(1.0 - pow(r / h, 3), 3);
            }
            if (userw) {
                w[j - 1] *= rw[j - 1];
            }
            a += w[j - 1];
        }
    }

    if (a <= 0.0) {
        return false;
    } else { // weighted least squares
        for (auto j = nleft; j <= nright; j++) { // make sum of w(j) == 1
            w[j - 1] /= a;
        }

        if (h > 0.0 && ideg > 0) { // use linear fit
            auto a = 0.0;
            for (auto j = nleft; j <= nright; j++) { // weighted center of x values
                a += w[j - 1] * ((T) j);
            }
            auto b = xs - a;
            auto c = 0.0;
            for (auto j = nleft; j <= nright; j++) {
                c += w[j - 1] * pow(((T) j) - a, 2);
            }
            if (sqrt(c) > 0.001 * range) {
                b /= c;

                // points are spread out enough to compute slope
                for (auto j = nleft; j <= nright; j++) {
                    w[j - 1] *= b * (((T) j) - a) + 1.0;
                }
            }
        }

        *ys = 0.0;
        for (auto j = nleft; j <= nright; j++) {
            *ys += w[j - 1] * y[j - 1];
        }

        return true;
    }
}

template <typename T>
void ess(const T* y, size_t n, size_t len, int ideg, size_t njump, bool userw, const T* rw, T* ys, T* res) {
    if (n < 2) {
        ys[0] = y[0];
        return;
    }

    size_t nleft = 0;
    size_t nright = 0;

    auto newnj = std::min(njump, n - 1);
    if (len >= n) {
        nleft = 1;
        nright = n;
        for (size_t i = 1; i <= n; i += newnj) {
            auto ok = est(y, n, len, ideg, (T) i, &ys[i - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[i - 1] = y[i - 1];
            }
        }
    } else if (newnj == 1) { // newnj equal to one, len less than n
        auto nsh = (len + 1) / 2;
        nleft = 1;
        nright = len;
        for (size_t i = 1; i <= n; i++) { // fitted value at i
            if (i > nsh && nright != n) {
                nleft += 1;
                nright += 1;
            }
            auto ok = est(y, n, len, ideg, (T) i, &ys[i - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[i - 1] = y[i - 1];
            }
        }
    } else { // newnj greater than one, len less than n
        auto nsh = (len + 1) / 2;
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
            auto ok = est(y, n, len, ideg, (T) i, &ys[i - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[i - 1] = y[i - 1];
            }
        }
    }

    if (newnj != 1) {
        for (size_t i = 1; i <= n - newnj; i += newnj) {
            auto delta = (ys[i + newnj - 1] - ys[i - 1]) / ((T) newnj);
            for (auto j = i + 1; j <= i + newnj - 1; j++) {
                ys[j - 1] = ys[i - 1] + delta * ((T) (j - i));
            }
        }
        auto k = ((n - 1) / newnj) * newnj + 1;
        if (k != n) {
            auto ok = est(y, n, len, ideg, (T) n, &ys[n - 1], nleft, nright, res, userw, rw);
            if (!ok) {
                ys[n - 1] = y[n - 1];
            }
            if (k != n - 1) {
                auto delta = (ys[n - 1] - ys[k - 1]) / ((T) (n - k));
                for (auto j = k + 1; j <= n - 1; j++) {
                    ys[j - 1] = ys[k - 1] + delta * ((T) (j - k));
                }
            }
        }
    }
}

template <typename T>
void ma(const T* x, size_t n, size_t len, T* ave) {
    auto newn = n - len + 1;
    auto flen = (T) len;
    auto v = 0.0;

    // get the first average
    for (size_t i = 0; i < len; i++) {
        v += x[i];
    }

    ave[0] = v / flen;
    if (newn > 1) {
        auto k = len;
        auto m = 0;
        for (size_t j = 1; j < newn; j++) {
            // window down the array
            v = v - x[m] + x[k];
            ave[j] = v / flen;
            k += 1;
            m += 1;
        }
    }
}

template <typename T>
void fts(const T* x, size_t n, size_t np, T* trend, T* work) {
    ma(x, n, np, trend);
    ma(trend, n - np + 1, np, work);
    ma(work, n - 2 * np + 2, 3, trend);
}

template <typename T>
void rwts(const T* y, size_t n, const T* fit, T* rw) {
    for (size_t i = 0; i < n; i++) {
        rw[i] = fabs(y[i] - fit[i]);
    }

    auto mid1 = (n - 1) / 2;
    auto mid2 = n / 2;

    // sort
    std::sort(rw, rw + n);

    auto cmad = 3.0 * (rw[mid1] + rw[mid2]); // 6 * median abs resid
    auto c9 = 0.999 * cmad;
    auto c1 = 0.001 * cmad;

    for (size_t i = 0; i < n; i++) {
        auto r = fabs(y[i] - fit[i]);
        if (r <= c1) {
            rw[i] = 1.0;
        } else if (r <= c9) {
            rw[i] = pow(1.0 - pow(r / cmad, 2), 2);
        } else {
            rw[i] = 0.0;
        }
    }
}

template <typename T>
void ss(const T* y, size_t n, size_t np, size_t ns, int isdeg, size_t nsjump, bool userw, T* rw, T* season, T* work1, T* work2, T* work3, T* work4) {
    for (size_t j = 1; j <= np; j++) {
        size_t k = (n - j) / np + 1;

        for (size_t i = 1; i <= k; i++) {
            work1[i - 1] = y[(i - 1) * np + j - 1];
        }
        if (userw) {
            for (size_t i = 1; i <= k; i++) {
                work3[i - 1] = rw[(i - 1) * np + j - 1];
            }
        }
        ess(work1, k, ns, isdeg, nsjump, userw, work3, work2 + 1, work4);
        T xs = 0.0;
        auto nright = std::min(ns, k);
        auto ok = est(work1, k, ns, isdeg, xs, &work2[0], 1, nright, work4, userw, work3);
        if (!ok) {
            work2[0] = work2[1];
        }
        xs = k + 1;
        size_t nleft = std::max(1, (int) k - (int) ns + 1);
        ok = est(work1, k, ns, isdeg, xs, &work2[k + 1], nleft, k, work4, userw, work3);
        if (!ok) {
            work2[k + 1] = work2[k];
        }
        for (size_t m = 1; m <= k + 2; m++) {
            season[(m - 1) * np + j - 1] = work2[m - 1];
        }
    }
}

template <typename T>
void onestp(const T* y, size_t n, size_t np, size_t ns, size_t nt, size_t nl, int isdeg, int itdeg, int ildeg, size_t nsjump, size_t ntjump, size_t nljump, size_t ni, bool userw, T* rw, T* season, T* trend, T* work1, T* work2, T* work3, T* work4, T* work5) {
    for (size_t j = 0; j < ni; j++) {
        for (size_t i = 0; i < n; i++) {
            work1[i] = y[i] - trend[i];
        }

        ss(work1, n, np, ns, isdeg, nsjump, userw, rw, work2, work3, work4, work5, season);
        fts(work2, n + 2 * np, np, work3, work1);
        ess(work3, n, nl, ildeg, nljump, false, work4, work1, work5);
        for (size_t i = 0; i < n; i++) {
            season[i] = work2[np + i] - work1[i];
        }
        for (size_t i = 0; i < n; i++) {
            work1[i] = y[i] - season[i];
        }
        ess(work1, n, nt, itdeg, ntjump, userw, rw, trend, work3);
    }
}

template <typename T>
void stl(const T* y, size_t n, size_t np, size_t ns, size_t nt, size_t nl, int isdeg, int itdeg, int ildeg, size_t nsjump, size_t ntjump, size_t nljump, size_t ni, size_t no, T* rw, T* season, T* trend) {
    if (ns < 3) {
        throw std::invalid_argument("seasonal_length must be at least 3");
    }
    if (nt < 3) {
        throw std::invalid_argument("trend_length must be at least 3");
    }
    if (nl < 3) {
        throw std::invalid_argument("low_pass_length must be at least 3");
    }
    if (np < 2) {
        throw std::invalid_argument("period must be at least 2");
    }

    if (isdeg != 0 && isdeg != 1) {
        throw std::invalid_argument("seasonal_degree must be 0 or 1");
    }
    if (itdeg != 0 && itdeg != 1) {
        throw std::invalid_argument("trend_degree must be 0 or 1");
    }
    if (ildeg != 0 && ildeg != 1) {
        throw std::invalid_argument("low_pass_degree must be 0 or 1");
    }

    if (ns % 2 != 1) {
        throw std::invalid_argument("seasonal_length must be odd");
    }
    if (nt % 2 != 1) {
        throw std::invalid_argument("trend_length must be odd");
    }
    if (nl % 2 != 1) {
        throw std::invalid_argument("low_pass_length must be odd");
    }

    auto work1 = std::vector<T>(n + 2 * np);
    auto work2 = std::vector<T>(n + 2 * np);
    auto work3 = std::vector<T>(n + 2 * np);
    auto work4 = std::vector<T>(n + 2 * np);
    auto work5 = std::vector<T>(n + 2 * np);

    auto userw = false;
    size_t k = 0;

    while (true) {
        onestp(y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, work1.data(), work2.data(), work3.data(), work4.data(), work5.data());
        k += 1;
        if (k > no) {
            break;
        }
        for (size_t i = 0; i < n; i++) {
            work1[i] = trend[i] + season[i];
        }
        rwts(y, n, work1.data(), rw);
        userw = true;
    }

    if (no <= 0) {
        for (size_t i = 0; i < n; i++) {
            rw[i] = 1.0;
        }
    }
}

template <typename T>
T var(const std::vector<T>& series) {
    auto mean = std::accumulate(series.begin(), series.end(), 0.0) / series.size();
    std::vector<T> tmp;
    tmp.reserve(series.size());
    for (auto v : series) {
        tmp.push_back(pow(v - mean, 2));
    }
    return std::accumulate(tmp.begin(), tmp.end(), 0.0) / (series.size() - 1);
}

template <typename T>
T strength(const std::vector<T>& component, const std::vector<T>& remainder) {
    std::vector<T> sr;
    sr.reserve(remainder.size());
    for (size_t i = 0; i < remainder.size(); i++) {
        sr.push_back(component[i] + remainder[i]);
    }
    return std::max(0.0, 1.0 - var(remainder) / var(sr));
}

template <class T>
class StlResult {
public:
    std::vector<T> seasonal;
    std::vector<T> trend;
    std::vector<T> remainder;
    std::vector<T> weights;

    inline T seasonal_strength() {
        return strength(seasonal, remainder);
    }

    inline T trend_strength() {
        return strength(trend, remainder);
    }
};

class StlParams {
    std::optional<size_t> ns_ = std::nullopt;
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
    inline StlParams seasonal_length(size_t ns) {
        this->ns_ = ns;
        return *this;
    }

    inline StlParams trend_length(size_t nt) {
        this->nt_ = nt;
        return *this;
    }

    inline StlParams low_pass_length(size_t nl) {
        this->nl_ = nl;
        return *this;
    }

    inline StlParams seasonal_degree(int isdeg) {
        this->isdeg_ = isdeg;
        return *this;
    }

    inline StlParams trend_degree(int itdeg) {
        this->itdeg_ = itdeg;
        return *this;
    }

    inline StlParams low_pass_degree(int ildeg) {
        this->ildeg_ = ildeg;
        return *this;
    }

    inline StlParams seasonal_jump(size_t nsjump) {
        this->nsjump_ = nsjump;
        return *this;
    }

    inline StlParams trend_jump(size_t ntjump) {
        this->ntjump_ = ntjump;
        return *this;
    }

    inline StlParams low_pass_jump(size_t nljump) {
        this->nljump_ = nljump;
        return *this;
    }

    inline StlParams inner_loops(size_t ni) {
        this->ni_ = ni;
        return *this;
    }

    inline StlParams outer_loops(size_t no) {
        this->no_ = no;
        return *this;
    }

    inline StlParams robust(bool robust) {
        this->robust_ = robust;
        return *this;
    }

    StlResult<float> fit(const float* y, size_t n, size_t np);
    StlResult<double> fit(const double* y, size_t n, size_t np);
    StlResult<float> fit(const std::vector<float>& y, size_t np);
    StlResult<double> fit(const std::vector<double>& y, size_t np);
};

StlParams params() {
    return StlParams();
}

template <typename T>
StlResult<T> stl_fit(const T* y, size_t n, size_t np, std::optional<size_t> ns_, std::optional<size_t> nt_, std::optional<size_t> nl_, int isdeg_, int itdeg_, std::optional<int> ildeg_, std::optional<size_t> nsjump_, std::optional<size_t> ntjump_, std::optional<size_t> nljump_, std::optional<size_t> ni_, std::optional<size_t> no_, bool robust_) {
    if (n < 2 * np) {
        throw std::invalid_argument("series has less than two periods");
    }

    auto ns = ns_.value_or(np);

    auto isdeg = isdeg_;
    auto itdeg = itdeg_;

    auto res = StlResult<T> {
        std::vector<T>(n),
        std::vector<T>(n),
        std::vector<T>(),
        std::vector<T>(n)
    };

    auto ildeg = ildeg_.value_or(itdeg);
    auto newns = std::max(ns, (size_t) 3);
    if (newns % 2 == 0) {
        newns += 1;
    }

    auto newnp = std::max(np, (size_t) 2);
    auto nt = (size_t) ceil((1.5 * newnp) / (1.0 - 1.5 / (T) newns));
    nt = nt_.value_or(nt);
    nt = std::max(nt, (size_t) 3);
    if (nt % 2 == 0) {
        nt += 1;
    }

    auto nl = nl_.value_or(newnp);
    if (nl % 2 == 0 && !nl_.has_value()) {
        nl += 1;
    }

    auto ni = ni_.value_or(robust_ ? 1 : 2);
    auto no = no_.value_or(robust_ ? 15 : 0);

    auto nsjump = nsjump_.value_or((size_t) ceil(((T) newns) / 10.0));
    auto ntjump = ntjump_.value_or((size_t) ceil(((T) nt) / 10.0));
    auto nljump = nljump_.value_or((size_t) ceil(((T) nl) / 10.0));

    stl(y, n, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, no, res.weights.data(), res.seasonal.data(), res.trend.data());

    res.remainder.reserve(n);
    for (size_t i = 0; i < n; i++) {
        res.remainder.push_back(y[i] - res.seasonal[i] - res.trend[i]);
    }

    return res;
}

StlResult<float> StlParams::fit(const float* y, size_t n, size_t np) {
    return stl_fit(y, n, np, ns_, nt_, nl_, isdeg_, itdeg_, ildeg_, nsjump_, ntjump_, nljump_, ni_, no_, robust_);
}
StlResult<double> StlParams::fit(const double* y, size_t n, size_t np) {
    return stl_fit(y, n, np, ns_, nt_, nl_, isdeg_, itdeg_, ildeg_, nsjump_, ntjump_, nljump_, ni_, no_, robust_);
}
StlResult<float> StlParams::fit(const std::vector<float>& y, size_t np) {
    return stl_fit(y.data(), y.size(), np, ns_, nt_, nl_, isdeg_, itdeg_, ildeg_, nsjump_, ntjump_, nljump_, ni_, no_, robust_);
}
StlResult<double> StlParams::fit(const std::vector<double>& y, size_t np) {
    return stl_fit(y.data(), y.size(), np, ns_, nt_, nl_, isdeg_, itdeg_, ildeg_, nsjump_, ntjump_, nljump_, ni_, no_, robust_);
}

}
