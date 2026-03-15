#include <cassert>
#include <cmath>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

#include <stl.hpp>

#include "helper.hpp"

using stl::Mstl;

template<typename T>
void test_mstl_works() {
    Mstl<T> fit{generate_series<T>(), {6, 10}};
    assert_elements_in_delta(
        {0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874},
        first(fit.seasonal().at(0), 5)
    );
    assert_elements_in_delta(
        {1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514},
        first(fit.seasonal().at(1), 5)
    );
    assert_elements_in_delta(
        {5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475},
        first(fit.remainder(), 5)
    );
}

template<typename T>
void test_mstl_span() {
    std::vector<T> series = generate_series<T>();
    Mstl<T> fit{std::span<const T>(series), {{6, 10}}};
    assert_elements_in_delta(
        {0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874},
        first(fit.seasonal().at(0), 5)
    );
    assert_elements_in_delta(
        {1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514},
        first(fit.seasonal().at(1), 5)
    );
    assert_elements_in_delta(
        {5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475},
        first(fit.remainder(), 5)
    );
}

template<typename T>
void test_mstl_unsorted_periods() {
    Mstl<T> fit{generate_series<T>(), {10, 6}};
    assert_elements_in_delta(
        {1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514},
        first(fit.seasonal().at(0), 5)
    );
    assert_elements_in_delta(
        {0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874},
        first(fit.seasonal().at(1), 5)
    );
    assert_elements_in_delta(
        {5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475},
        first(fit.remainder(), 5)
    );
}

template<typename T>
void test_mstl_lambda() {
    Mstl<T> fit{generate_series<T>(), {6, 10}, { .lambda = 0.5 }};
    assert_elements_in_delta(
        {0.43371448, 0.10503793, -0.7178911, 1.2356076, -1.8253292},
        first(fit.seasonal().at(0), 5)
    );
    assert_elements_in_delta(
        {1.0437742, 0.8650516, 0.07303603, -1.428663, -1.1990008},
        first(fit.seasonal().at(1), 5)
    );
    assert_elements_in_delta(
        {2.0748303, 2.1291165, 2.1834028, 2.2330272, 2.2826517},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-1.0801829, 0.900794, -0.7101207, 1.9600279, -1.2583216},
        first(fit.remainder(), 5)
    );
}

template<typename T>
void test_mstl_lambda_zero() {
    std::vector<T> series;
    for (auto& v : generate_series<T>()) {
        series.push_back(v + 1);
    }
    Mstl<T> fit{series, {6, 10}, { .lambda = 0.0 }};
    assert_elements_in_delta(
        {0.18727916, 0.029921893, -0.2716494, 0.47748315, -0.7320051},
        first(fit.seasonal().at(0), 5)
    );
    assert_elements_in_delta(
        {0.42725056, 0.32145387, -0.019030934, -0.56607914, -0.46765903},
        first(fit.seasonal().at(1), 5)
    );
    assert_elements_in_delta(
        {1.592807, 1.6144379, 1.6360688, 1.6559447, 1.6758206},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-0.41557717, 0.33677137, -0.24677622, 0.7352363, -0.47615635},
        first(fit.remainder(), 5)
    );
}

template<typename T>
void test_mstl_lambda_out_of_range() {
    assert_exception<std::invalid_argument>([]() {
        Mstl<T>{generate_series<T>(), {6, 10}, { .lambda = 2.0 }};
    }, "lambda must be between 0 and 1");
}

template<typename T>
void test_mstl_empty_periods() {
    assert_exception<std::invalid_argument>([]() {
        Mstl<T>{generate_series<T>(), {}};
    }, "periods must not be empty");
}

template<typename T>
void test_mstl_period_one() {
    assert_exception<std::invalid_argument>([]() {
        Mstl<T>{generate_series<T>(), {1}};
    }, "periods must be at least 2");
}

template<typename T>
void test_mstl_too_few_periods() {
    assert_exception<std::invalid_argument>([]() {
        Mstl<T>{generate_series<T>(), {16}};
    }, "series has less than two periods");
}

template<typename T>
void test_mstl_seasonal_strength() {
    Mstl<T> fit{generate_series<T>(), {7}, { .stl_params = { .seasonal_length = 7 } }};
    assert_in_delta(0.284111676315015, fit.seasonal_strength().at(0));
}

template<typename T>
void test_mstl_seasonal_strength_max() {
    std::vector<T> series = max_seasonal_series<T>();
    Mstl<T> fit{series, {7}, { .stl_params = { .seasonal_length = 7 } }};
    assert_in_delta(1.0, fit.seasonal_strength().at(0));
}

template<typename T>
void test_mstl_trend_strength() {
    Mstl<T> fit{generate_series<T>(), {7}, { .stl_params = { .seasonal_length = 7 } }};
    assert_in_delta(0.16384245231864702, fit.trend_strength());
}

template<typename T>
void test_mstl_trend_strength_max() {
    std::vector<T> series = max_trend_series<T>();
    Mstl<T> fit{series, {7}, { .stl_params = { .seasonal_length = 7 } }};
    assert_in_delta(1.0, fit.trend_strength());
}

template<typename T>
void test_type() {
    test_mstl_works<T>();
    test_mstl_span<T>();
    test_mstl_unsorted_periods<T>();
    test_mstl_lambda<T>();
    test_mstl_lambda_zero<T>();
    test_mstl_lambda_out_of_range<T>();
    test_mstl_empty_periods<T>();
    test_mstl_period_one<T>();
    test_mstl_too_few_periods<T>();
    test_mstl_seasonal_strength<T>();
    test_mstl_seasonal_strength_max<T>();
    test_mstl_trend_strength<T>();
    test_mstl_trend_strength_max<T>();
}

void test_mstl() {
    test_type<float>();
    test_type<double>();
}
