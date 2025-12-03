#pragma once
#include "StochasticProcess.hpp"
#include <algorithm>
#include <cmath>

namespace MyQuantLib {

    struct HestonParams {
        double v0;    // Initial variance
        double kappa; // Mean reversion
        double theta; // Long-term variance
        double xi;    // Vol of vol
        double rho;   // Correlation
    };

    class HestonProcess : public StochasticProcess {
        HestonParams params_;
    public:
        explicit HestonProcess(const HestonParams& p) : params_(p) {}

		// Getter for Heston parameters
        const HestonParams& getParams() const { return params_; }

        // Implement log-Euleur scheme
        std::pair<double, double> evolve(
            double s_t, double v_t, double leverage,
            double dt, double dZ_s, double dZ_v,
            double r, double q) const override;
    };

}