#pragma once
#include <vector>
#include <cmath>
#include <numeric>

namespace MyQuantLib {

    class KernelDensity {
    public:
        // Nadaraya-Watson estimator for E[L(t,X)^2 * v| X=x]
		// Calculate the conditional expectation of Y given X using kernel density estimation
        static std::vector<double> estimateConditionalExpectation(
            const std::vector<double>& s_particles, // X
            const std::vector<double>& v_particles, // Y
            const std::vector<double>& s_grid,      // Target points
            double bandwidth,                       // h
            double fallback_value                   // Theta (if the sample batch size is not enough)
        );

		// Bandwith selection using Silverman's rule of thumb
        static double silvermanBandwidth(const std::vector<double>& data);
    };

}