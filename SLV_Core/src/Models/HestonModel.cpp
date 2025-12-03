#include "Models/HestonModel.hpp"

namespace MyQuantLib {

    std::pair<double, double> HestonProcess::evolve(
        double s_t, double v_t, double leverage,
        double dt, double dZ_s, double dZ_v,
        double r, double q) const
    {
		// Full Truncation / Reflection for variance
        double v_t_pos = std::max(v_t, 0.0);
        double sqrt_v_dt = std::sqrt(v_t_pos * dt);

        // Evolve Stock (Geometric Brownian Motion scale by leverage)
        // drift = (r - q - 0.5 * (L * sqrt(v))^2) * dt
        // diffusion = L * sqrt(v) * dW
        double vol = leverage * std::sqrt(v_t_pos); // Sigma_SLV
        double drift = (r - q - 0.5 * vol * vol) * dt;
        double diffusion = leverage * sqrt_v_dt * dZ_s;

        double s_next = s_t * std::exp(drift + diffusion);

        // Evolve Variance (Heston)
        double v_next = v_t + params_.kappa * (params_.theta - v_t_pos) * dt
            + params_.xi * sqrt_v_dt * dZ_v;

		// Make sure variance is non-negative
        v_next = std::max(v_next, 0.0);

        return { s_next, v_next };
    }

}