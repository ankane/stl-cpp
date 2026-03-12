#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iostream>
#include <optional>
#include <ranges>
#include <string_view>
#include <vector>

template<typename T>
void print_vector(const std::vector<T>& x) {
    for (auto v : x) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
}

template<typename T>
std::vector<T> first(const std::vector<T>& x, ptrdiff_t n) {
    auto view = std::views::take(x, n);
    return std::vector<T>(view.begin(), view.end());
}

template<typename T, typename U>
void assert_in_delta(T exp, U act) {
    assert(std::abs(exp - act) < 0.001);
}

template<typename T>
void assert_elements_in_delta(const std::vector<double>& exp, const std::vector<T>& act) {
    assert(exp.size() == act.size());
    for (size_t i = 0; i < exp.size(); i++) {
        assert_in_delta(exp.at(i), act.at(i));
    }
}

template<typename T>
void assert_exception(const std::function<void(void)>& code, std::optional<std::string_view> message = std::nullopt) {
    std::optional<T> exception;
    try {
        code();
    } catch (const T& e) {
        exception = e;
    }
    assert(exception.has_value());
    if (message) {
        assert(std::string_view{exception.value().what()} == message.value());
    }
}

template<typename T>
std::vector<T> generate_series() {
    std::vector<T> series{
        5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
        7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
        3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
    };
    return series;
}

template<typename T>
std::vector<T> max_seasonal_series() {
    std::vector<T> series;
    series.reserve(30);
    for (size_t i = 0; i < 30; i++) {
        series.push_back(static_cast<T>(i % 7));
    }
    return series;
}

template<typename T>
std::vector<T> max_trend_series() {
    std::vector<T> series;
    series.reserve(30);
    for (size_t i = 0; i < 30; i++) {
        series.push_back(static_cast<T>(i));
    }
    return series;
}
