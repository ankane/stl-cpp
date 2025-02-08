#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#if __cplusplus >= 202002L
#include <span>
#endif

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
    for (size_t i = 0; i < exp.size(); i++) {
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

std::vector<float> max_seasonal_series() {
    std::vector<float> series;
    for (size_t i = 0; i < 30; i++) {
        series.push_back(i % 7);
    }
    return series;
}

std::vector<float> max_trend_series() {
    std::vector<float> series;
    for (size_t i = 0; i < 30; i++) {
        series.push_back(i);
    }
    return series;
}

void test_stl_works() {
    auto series = generate_series();
    auto result = stl::params().fit(series, 7);
    assert_elements_in_delta({0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802}, first(result.seasonal, 5));
    assert_elements_in_delta({4.804099, 4.9097075, 5.015316, 5.16045, 5.305584}, first(result.trend, 5));
    assert_elements_in_delta({-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037}, first(result.remainder, 5));
    assert_elements_in_delta({1.0, 1.0, 1.0, 1.0, 1.0}, first(result.weights, 5));
}

#if __cplusplus >= 202002L
void test_stl_span() {
    auto series = generate_series();
    auto result = stl::params().fit(std::span{series}, 7);
    assert_elements_in_delta({0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802}, first(result.seasonal, 5));
    assert_elements_in_delta({4.804099, 4.9097075, 5.015316, 5.16045, 5.305584}, first(result.trend, 5));
    assert_elements_in_delta({-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037}, first(result.remainder, 5));
    assert_elements_in_delta({1.0, 1.0, 1.0, 1.0, 1.0}, first(result.weights, 5));
}
#endif

void test_stl_robust() {
    auto series = generate_series();
    auto result = stl::params().robust(true).fit(series, 7);
    assert_elements_in_delta({0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711}, first(result.seasonal, 5));
    assert_elements_in_delta({5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114}, first(result.trend, 5));
    assert_elements_in_delta({-0.5465884, 3.0460663, -1.7184906, 1.6089439, -6.5681853}, first(result.remainder, 5));
    assert_elements_in_delta({0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217}, first(result.weights, 5));
}

void test_stl_too_few_periods() {
    ASSERT_EXCEPTION(
        stl::params().fit(generate_series(), 16),
        std::invalid_argument,
        "series has less than two periods"
    );
}

void test_stl_bad_seasonal_degree() {
    ASSERT_EXCEPTION(
        stl::params().seasonal_degree(2).fit(generate_series(), 7),
        std::invalid_argument,
        "seasonal_degree must be 0 or 1"
    );
}

void test_stl_seasonal_strength() {
    auto result = stl::params().fit(generate_series(), 7);
    assert_in_delta(0.284111676315015, result.seasonal_strength());
}

void test_stl_seasonal_strength_max() {
    auto series = max_seasonal_series();
    auto result = stl::params().fit(series, 7);
    assert_in_delta(1.0, result.seasonal_strength());
}

void test_stl_trend_strength() {
    auto result = stl::params().fit(generate_series(), 7);
    assert_in_delta(0.16384245231864702, result.trend_strength());
}

void test_stl_trend_strength_max() {
    auto series = max_trend_series();
    auto result = stl::params().fit(series, 7);
    assert_in_delta(1.0, result.trend_strength());
}

void test_mstl_works() {
    auto result = stl::mstl_params().fit(generate_series(), {6, 10});
    assert_elements_in_delta({0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874}, first(result.seasonal[0], 5));
    assert_elements_in_delta({1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514}, first(result.seasonal[1], 5));
    assert_elements_in_delta({5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862}, first(result.trend, 5));
    assert_elements_in_delta({-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475}, first(result.remainder, 5));
}

#if __cplusplus >= 202002L
void test_mstl_span() {
    auto series = generate_series();
    auto result = stl::mstl_params().fit(std::span{series}, {6, 10});
    assert_elements_in_delta({0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874}, first(result.seasonal[0], 5));
    assert_elements_in_delta({1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514}, first(result.seasonal[1], 5));
    assert_elements_in_delta({5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862}, first(result.trend, 5));
    assert_elements_in_delta({-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475}, first(result.remainder, 5));
}
#endif

void test_mstl_unsorted_periods() {
    auto result = stl::mstl_params().fit(generate_series(), {10, 6});
    assert_elements_in_delta({1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514}, first(result.seasonal[0], 5));
    assert_elements_in_delta({0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874}, first(result.seasonal[1], 5));
    assert_elements_in_delta({5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862}, first(result.trend, 5));
    assert_elements_in_delta({-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475}, first(result.remainder, 5));
}

void test_mstl_lambda() {
    auto result = stl::mstl_params().lambda(0.5).fit(generate_series(), {6, 10});
    assert_elements_in_delta({0.43371448, 0.10503793, -0.7178911, 1.2356076, -1.8253292}, first(result.seasonal[0], 5));
    assert_elements_in_delta({1.0437742, 0.8650516, 0.07303603, -1.428663, -1.1990008}, first(result.seasonal[1], 5));
    assert_elements_in_delta({2.0748303, 2.1291165, 2.1834028, 2.2330272, 2.2826517}, first(result.trend, 5));
    assert_elements_in_delta({-1.0801829, 0.900794, -0.7101207, 1.9600279, -1.2583216}, first(result.remainder, 5));
}

void test_mstl_lambda_zero() {
    std::vector<float> series;
    for (auto& v : generate_series()) {
        series.push_back(v + 1);
    }
    auto result = stl::mstl_params().lambda(0.0).fit(series, {6, 10});
    assert_elements_in_delta({0.18727916, 0.029921893, -0.2716494, 0.47748315, -0.7320051}, first(result.seasonal[0], 5));
    assert_elements_in_delta({0.42725056, 0.32145387, -0.019030934, -0.56607914, -0.46765903}, first(result.seasonal[1], 5));
    assert_elements_in_delta({1.592807, 1.6144379, 1.6360688, 1.6559447, 1.6758206}, first(result.trend, 5));
    assert_elements_in_delta({-0.41557717, 0.33677137, -0.24677622, 0.7352363, -0.47615635}, first(result.remainder, 5));
}

void test_mstl_lambda_out_of_range() {
    ASSERT_EXCEPTION(
        stl::mstl_params().lambda(2.0).fit(generate_series(), {6, 10}),
        std::invalid_argument,
        "lambda must be between 0 and 1"
    );
}

void test_mstl_empty_periods() {
    ASSERT_EXCEPTION(
        stl::mstl_params().fit(generate_series(), {}),
        std::invalid_argument,
        "periods must not be empty"
    );
}

void test_mstl_period_one() {
    ASSERT_EXCEPTION(
        stl::mstl_params().fit(generate_series(), {1}),
        std::invalid_argument,
        "periods must be at least 2"
    );
}

void test_mstl_too_few_periods() {
    ASSERT_EXCEPTION(
        stl::mstl_params().fit(generate_series(), {16}),
        std::invalid_argument,
        "series has less than two periods"
    );
}

void test_mstl_seasonal_strength() {
    auto result = stl::mstl_params().stl_params(stl::params().seasonal_length(7)).fit(generate_series(), {7});
    assert_in_delta(0.284111676315015, result.seasonal_strength()[0]);
}

void test_mstl_seasonal_strength_max() {
    auto series = max_seasonal_series();
    auto result = stl::mstl_params().stl_params(stl::params().seasonal_length(7)).fit(series, {7});
    assert_in_delta(1.0, result.seasonal_strength()[0]);
}

void test_mstl_trend_strength() {
    auto result = stl::mstl_params().stl_params(stl::params().seasonal_length(7)).fit(generate_series(), {7});
    assert_in_delta(0.16384245231864702, result.trend_strength());
}

void test_mstl_trend_strength_max() {
    auto series = max_trend_series();
    auto result = stl::mstl_params().stl_params(stl::params().seasonal_length(7)).fit(series, {7});
    assert_in_delta(1.0, result.trend_strength());
}

int main() {
    test_stl_works();
#if __cplusplus >= 202002L
    test_stl_span();
#endif
    test_stl_robust();
    test_stl_too_few_periods();
    test_stl_bad_seasonal_degree();
    test_stl_seasonal_strength();
    test_stl_seasonal_strength_max();
    test_stl_trend_strength();
    test_stl_trend_strength_max();

    test_mstl_works();
#if __cplusplus >= 202002L
    test_mstl_span();
#endif
    test_mstl_unsorted_periods();
    test_mstl_lambda();
    test_mstl_lambda_zero();
    test_mstl_lambda_out_of_range();
    test_mstl_empty_periods();
    test_mstl_period_one();
    test_mstl_too_few_periods();
    test_mstl_seasonal_strength();
    test_mstl_seasonal_strength_max();
    test_mstl_trend_strength();
    test_mstl_trend_strength_max();

    return 0;
}
