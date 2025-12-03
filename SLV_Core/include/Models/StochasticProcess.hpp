#pragma once
#include <vector>
#include <utility>

namespace MyQuantLib {

    // Abstract Base Class for lsv process
    class StochasticProcess {
    public:
        virtual ~StochasticProcess() = default;

        // return (S_next, v_next)
        virtual std::pair<double, double> evolve(
            double s_t, double v_t, double leverage,
            double dt, double dZ_s, double dZ_v,
            double r, double q) const = 0;
    };

}