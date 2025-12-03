#include "Models/LocalVolMC.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <iostream>

namespace MyQuantLib {

    std::vector<double> LocalVolMC::priceOptions(
        double T,
        const std::vector<double>& strikes,
        int N_paths,
        int M_steps)
    {
        double dt = T / M_steps;
        double sqrt_dt = std::sqrt(dt);
        double s0 = (*market_).s0;
        double r = market_->r;
        double q = market_->q;

        std::vector<double> final_prices(N_paths);

        // RNG for each thread
        int max_threads = omp_get_max_threads();
        std::vector<std::mt19937> rngs(max_threads);
        for (int i = 0; i < max_threads; ++i) rngs[i].seed(1234 + i);

        // --- 1. MONTE CARLO SIMULATION (Parallel) ---
        #pragma omp parallel for
        for (int i = 0; i < N_paths; ++i) {
            int tid = omp_get_thread_num();
            std::normal_distribution<> d(0.0, 1.0);

            double S = s0;
            double t = 0.0;

            for (int step = 0; step < M_steps; ++step) {
                double k = std::log(S / s0) - r * t;
                double vol = market_->getLocalVol(t, k);

                // Euler Discretization
                // dS = S * ( (r-q)dt + vol*dW )
                // Log-Euler: S_new = S * exp( (r-q - 0.5*v^2)dt + vol*sqrt(dt)*Z )
                double drift = (r - q - 0.5 * vol * vol) * dt;
                double diffusion = vol * sqrt_dt * d(rngs[tid]);

                S *= std::exp(drift + diffusion);
                t += dt;
            }
            final_prices[i] = S;
        }

        // --- 2. CALCULATE OPTION PRICES ---
        std::vector<double> results;
        double df = std::exp(-r * T); // Discount factor

        for (double K : strikes) {
            double sum_payoff = 0.0;

            for (double S_T : final_prices) {
                sum_payoff += std::max(S_T - K, 0.0);
            }

            results.push_back((sum_payoff / N_paths) * df);
        }

        return results;
    }

}