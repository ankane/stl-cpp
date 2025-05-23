#include <iostream>
#include <vector>

#include "stl.hpp"

int main() {
    std::vector<float> series = {
        5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
        7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
        3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
    };
    size_t period = 7; // period of seasonal component

    auto res = stl::params().fit(series, period);
    for (auto v : res.trend) {
        std::cout << v << std::endl;
    }

    return 0;
}
