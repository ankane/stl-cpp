#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#include "../include/stl.hpp"

#define ASSERT_EXCEPTION(code, type, message) { \
    try {                                       \
        code;                                   \
        assert(false);                          \
    } catch (const type &e) {                   \
        assert(strcmp(e.what(), message) == 0); \
    }                                           \
}

void print_vector(const std::vector<float>& x) {
    for (auto v : x) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
}

std::vector<float> first(const std::vector<float>& x, size_t n) {
    n = std::min(n, x.size());
    return std::vector<float>(x.begin(), x.begin() + n);
}

void assert_in_delta(float exp, float act) {
    assert(fabs(exp - act) < 0.001);
}

void assert_elements_in_delta(const std::vector<float>& exp, const std::vector<float>& act) {
    assert(exp.size() == act.size());
    for (auto i = 0; i < exp.size(); i++) {
        assert_in_delta(exp[i], act[i]);
    }
}

std::vector<float> generate_series() {
    std::vector<float> series = {
        5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
        7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
        3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
    };
    return series;
}

void test_works() {
    auto series = generate_series();
    auto result = stl::params().fit(series, 7);
    assert_elements_in_delta({0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802}, first(result.seasonal, 5));
    assert_elements_in_delta({4.804099, 4.9097075, 5.015316, 5.16045, 5.305584}, first(result.trend, 5));
    assert_elements_in_delta({-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037}, first(result.remainder, 5));
    assert_elements_in_delta({1.0, 1.0, 1.0, 1.0, 1.0}, first(result.weights, 5));
}

void test_robust() {
    auto series = generate_series();
    auto result = stl::params().robust(true).fit(series, 7);
    assert_elements_in_delta({0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711}, first(result.seasonal, 5));
    assert_elements_in_delta({5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114}, first(result.trend, 5));
    assert_elements_in_delta({-0.5465884, 3.0460663, -1.7184906, 1.6089439, -6.5681853}, first(result.remainder, 5));
    assert_elements_in_delta({0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217}, first(result.weights, 5));
}

void test_too_few_periods() {
    ASSERT_EXCEPTION(
        stl::params().fit(generate_series(), 16),
        std::invalid_argument,
        "series has less than two periods"
    );
}

void test_bad_seasonal_degree() {
    ASSERT_EXCEPTION(
        stl::params().seasonal_degree(2).fit(generate_series(), 7),
        std::invalid_argument,
        "seasonal_degree must be 0 or 1"
    );
}

void test_seasonal_strength() {
    auto result = stl::params().fit(generate_series(), 7);
    assert_in_delta(0.284111676315015, result.seasonal_strength());
}

void test_trend_strength() {
    auto result = stl::params().fit(generate_series(), 7);
    assert_in_delta(0.16384245231864702, result.trend_strength());
}

int main() {
    test_works();
    test_robust();
    test_too_few_periods();
    test_bad_seasonal_degree();
    test_seasonal_strength();
    test_trend_strength();
    return 0;
}
