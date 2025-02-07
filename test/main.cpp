#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#include "../include/stl.hpp"
#include "../include/mstl.hpp"

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

void stl_test_works() {
    auto series = generate_series();
    auto result = stl::params().fit(series, 7);
    assert_elements_in_delta({0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802}, first(result.seasonal, 5));
    assert_elements_in_delta({4.804099, 4.9097075, 5.015316, 5.16045, 5.305584}, first(result.trend, 5));
    assert_elements_in_delta({-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037}, first(result.remainder, 5));
    assert_elements_in_delta({1.0, 1.0, 1.0, 1.0, 1.0}, first(result.weights, 5));
}

void stl_test_robust() {
    auto series = generate_series();
    auto result = stl::params().robust(true).fit(series, 7);
    assert_elements_in_delta({0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711}, first(result.seasonal, 5));
    assert_elements_in_delta({5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114}, first(result.trend, 5));
    assert_elements_in_delta({-0.5465884, 3.0460663, -1.7184906, 1.6089439, -6.5681853}, first(result.remainder, 5));
    assert_elements_in_delta({0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217}, first(result.weights, 5));
}

void stl_test_too_few_periods() {
    ASSERT_EXCEPTION(
        stl::params().fit(generate_series(), 16),
        std::invalid_argument,
        "series is shorter than twice the period"
    );
}

void stl_test_bad_seasonal_degree() {
    ASSERT_EXCEPTION(
        stl::params().seasonal_degree(2).fit(generate_series(), 7),
        std::invalid_argument,
        "seasonal_degree must be 0 or 1"
    );
}

void stl_test_seasonal_strength() {
    auto result = stl::params().fit(generate_series(), 7);
    assert_in_delta(0.284111676315015, result.seasonal_strength());
}

void stl_test_trend_strength() {
    auto result = stl::params().fit(generate_series(), 7);
    assert_in_delta(0.16384245231864702, result.trend_strength());
}

void test_stl() {
    stl_test_works();
    stl_test_robust();
    stl_test_too_few_periods();
    stl_test_bad_seasonal_degree();
    stl_test_seasonal_strength();
    stl_test_trend_strength();
}

void mstl_handles_single_period() {
    auto series = generate_series();
    auto stl_params = stl::params()
        .trend_length(13)
        .low_pass_length(9)
        .trend_jump(1)
        .low_pass_jump(1)
        .seasonal_jump(1)
        .robust(false)
        .low_pass_degree(1)
        .trend_degree(1)
        .seasonal_degree(1)
        .seasonal_length(11)
        .inner_loops(2)
        .outer_loops(0);
    auto stl_result = stl_params.fit(series, 7);
    auto mstl_result = mstl::params().stl_params(stl_params).fit(series, {7});
    assert_elements_in_delta(stl_result.seasonal, mstl_result.seasonal[0]);
    assert_elements_in_delta(stl_result.remainder, mstl_result.remainder);
    assert_elements_in_delta(stl_result.trend, mstl_result.trend);
}

void mstl_handles_multiple_periods() {
    auto series = generate_series();
    auto stl_params = stl::params()
        .low_pass_length(13)
        .seasonal_length(13)
        .trend_length(13)
        .trend_jump(1)
        .low_pass_jump(1)
        .seasonal_jump(1)
        .robust(false)
        .low_pass_degree(1)
        .trend_degree(1)
        .seasonal_degree(1)
        .inner_loops(2)
        .outer_loops(0);
    auto result = mstl::params()
        .iterations(2)
        .seasonal_lengths({11, 13})
        .stl_params(stl_params)
        .fit(series, std::vector<size_t>{7, 10});
    // Values taken from parallel implementation of MSTL (python statsmodels)
    assert_elements_in_delta({1.02957645,  1.58052462, -2.58504053,3.82336372,  -1.37414519}, first(result.seasonal[0], 5));
    assert_elements_in_delta({-1.130680493627964, 2.4459641040455704, 0.3115169691001893, -0.9364803464881937, -4.19763814690413}, first(result.seasonal[1], 5));
    assert_elements_in_delta({4.899, 5.027 , 5.151, 5.270, 5.387}, first(result.trend, 5));
    assert_elements_in_delta({0.20186475, -0.05349705, -0.8779612 ,  0.84224536,  0.18390715}, first(result.remainder, 5));
}

void mstl_handles_cox() {
    auto series = generate_series();
    auto stl_params = stl::params()
        .low_pass_length(13)
        .seasonal_length(13)
        .trend_length(13)
        .trend_jump(1)
        .low_pass_jump(1)
        .seasonal_jump(1)
        .robust(false)
        .low_pass_degree(1)
        .trend_degree(1)
        .seasonal_degree(1)
        .inner_loops(2)
        .outer_loops(0);
    auto result = mstl::params()
        .iterations(2)
        .seasonal_lengths({11, 13})
        .lambda(0.5f)
        .stl_params(stl_params)
        .fit(series, std::vector<size_t>{7, 10});
    // Values taken from parallel implementation of MSTL (python statsmodels)
    assert_elements_in_delta({1.0345437330165619, 1.002305016841231, -1.2867553566909664, 2.365208882252409, -1.5555646550017448}, first(result.seasonal[0], 5));
    assert_elements_in_delta({-0.7310726692318952, 1.2115820999320608, 0.453518109999968, -1.3655355589288307, -2.6226520547233756}, first(result.seasonal[1], 5));
    assert_elements_in_delta({1.97986303, 2.05898726, 2.13443353, 2.20569258, 2.27523968}, first(result.trend, 5));
    assert_elements_in_delta({0.18880186, -0.27287437, -0.47276916,  0.79463409, -0.09702297}, first(result.remainder, 5));
}

void test_mstl() {
    mstl_handles_single_period();
    mstl_handles_multiple_periods();
    mstl_handles_cox();
}

int main() {
    test_stl();
    test_mstl();
    return 0;
}
