#include <cassert>
#include <cmath>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

#include <stl.hpp>

#include "helper.hpp"

using stl::Stl;

template<typename T>
void test_stl_works() {
    std::vector<T> series = generate_series<T>();
    Stl<T> fit{series, 7};
    assert_elements_in_delta(
        {0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802},
        first(fit.seasonal(), 5)
    );
    assert_elements_in_delta(
        {4.804099, 4.9097075, 5.015316, 5.16045, 5.305584},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037},
        first(fit.remainder(), 5)
    );
    assert_elements_in_delta({1.0, 1.0, 1.0, 1.0, 1.0}, first(fit.weights(), 5));
}

template<typename T>
void test_stl_span() {
    std::vector<T> series = generate_series<T>();
    Stl<T> fit{std::span<const T>(series), 7};
    assert_elements_in_delta(
        {0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802},
        first(fit.seasonal(), 5)
    );
    assert_elements_in_delta(
        {4.804099, 4.9097075, 5.015316, 5.16045, 5.305584},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037},
        first(fit.remainder(), 5)
    );
    assert_elements_in_delta({1.0, 1.0, 1.0, 1.0, 1.0}, first(fit.weights(), 5));
}

template<typename T>
void test_stl_robust() {
    std::vector<T> series = generate_series<T>();
    Stl<T> fit{series, 7, { .robust = true }};
    assert_elements_in_delta(
        {0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711},
        first(fit.seasonal(), 5)
    );
    assert_elements_in_delta(
        {5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114},
        first(fit.trend(), 5)
    );
    assert_elements_in_delta(
        {-0.5465884, 3.0460663, -1.7184906, 1.6089439, -6.5681853},
        first(fit.remainder(), 5)
    );
    assert_elements_in_delta(
        {0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217},
        first(fit.weights(), 5)
    );
}

template<typename T>
void test_stl_too_few_periods() {
    assert_exception<std::invalid_argument>([]() {
        Stl<T>{generate_series<T>(), 16};
    }, "series has less than two periods");
}

template<typename T>
void test_stl_bad_seasonal_degree() {
    assert_exception<std::invalid_argument>([]() {
        Stl<T>{generate_series<T>(), 7, { .seasonal_degree = 2 }};
    }, "seasonal_degree must be 0 or 1");
}

template<typename T>
void test_stl_seasonal_strength() {
    Stl<T> fit{generate_series<T>(), 7};
    assert_in_delta(0.284111676315015, fit.seasonal_strength());
}

template<typename T>
void test_stl_seasonal_strength_max() {
    std::vector<T> series = max_seasonal_series<T>();
    Stl<T> fit{series, 7};
    assert_in_delta(1.0, fit.seasonal_strength());
}

template<typename T>
void test_stl_trend_strength() {
    Stl<T> fit{generate_series<T>(), 7};
    assert_in_delta(0.16384245231864702, fit.trend_strength());
}

template<typename T>
void test_stl_trend_strength_max() {
    std::vector<T> series = max_trend_series<T>();
    Stl<T> fit{series, 7};
    assert_in_delta(1.0, fit.trend_strength());
}

template<typename T>
void test_type() {
    test_stl_works<T>();
    test_stl_span<T>();
    test_stl_robust<T>();
    test_stl_too_few_periods<T>();
    test_stl_bad_seasonal_degree<T>();
    test_stl_seasonal_strength<T>();
    test_stl_seasonal_strength_max<T>();
    test_stl_trend_strength<T>();
    test_stl_trend_strength_max<T>();
}

void test_stl() {
    test_type<float>();
    test_type<double>();
}
