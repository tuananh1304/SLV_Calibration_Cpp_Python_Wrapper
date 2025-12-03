#include "Models/SLVCalibrator.hpp"
#include "Math/Statistics.hpp"
#include <random>
#include <iostream>
#include <omp.h>

namespace MyQuantLib {

    SLVCalibrator::SLVCalibrator(
        std::shared_ptr<MarketData> market,
        std::shared_ptr<HestonProcess> heston,
        double T, int M, int N,
        const std::vector<double>& s_interp_grid)
        : market_(market), heston_(heston), T_(T), M_steps_(M), N_particles_(N), s_interp_grid_(s_interp_grid)
    {
        // Tao luoi thoi gian deu
        time_grid_.resize(M + 1);
        double dt = T / M;
        for (int i = 0; i <= M; ++i) time_grid_[i] = i * dt;
    }

    std::vector<std::vector<double>> SLVCalibrator::calibrate() {
        std::cout << "[C++] Starting SLV Calibration using " << N_particles_ << "particles..." << std::endl;

        double dt = T_ / M_steps_;
        size_t s_grid_size = s_interp_grid_.size();

		// 1. Initialize Leverage Surface
        std::vector<std::vector<double>> leverage_surface(M_steps_ + 1, std::vector<double>(s_grid_size));

		// 2. initialize particles at t=0
        std::vector<double> s_particles(N_particles_, market_->s0);
        std::vector<double> v_particles(N_particles_, heston_->getParams().v0);

		// 3. calculate initial leverage L(0, S)
        // L(0, S) = sigma_LV(0, S) / sqrt(v0)
        for (size_t j = 0; j < s_grid_size; ++j) {
			double k = std::log(s_interp_grid_[j] / market_->s0);
            double lv = market_->getLocalVol(0.0, s_interp_grid_[j]);
            leverage_surface[0][j] = lv / std::sqrt(heston_->getParams().v0);
        }

        // Random Number Generator
        int max_threads = omp_get_max_threads();
        std::vector<std::mt19937> rngs(max_threads);
        for (int i = 0; i < max_threads; ++i) rngs[i].seed(42 + i); // fix the seed to reproduce

        // 4. Time Stepping Loop
        for (int i = 0; i < M_steps_; ++i) {
            double t_curr = time_grid_[i];
            double t_next = time_grid_[i + 1];

			// Create interpolator for current leverage surface
            LinearInterpolator leverage_interp(s_interp_grid_, leverage_surface[i]);

            // Evolve particles (Parallel)
#pragma omp parallel for
            for (int p = 0; p < N_particles_; ++p) {
                int tid = omp_get_thread_num();
                std::normal_distribution<> d(0, 1);

				// Generate correlated Brownian increments
                double Z1 = d(rngs[tid]);
                double Z2 = d(rngs[tid]);
                double rho = heston_->getParams().rho;

                double dZ_s = Z1;
                double dZ_v = rho * Z1 + std::sqrt(1.0 - rho * rho) * Z2;

				// Search for leverage at current S
                double L_val = leverage_interp(s_particles[p]);

				// Heston evolve
                auto res = heston_->evolve(s_particles[p], v_particles[p], L_val,
                    dt, dZ_s, dZ_v, market_->r, market_->q);

                s_particles[p] = res.first;
                v_particles[p] = res.second;
            }

			// Calculate new leverage surface at t_next
            // Bandwidth selection
            double h = KernelDensity::silvermanBandwidth(s_particles);

			// Calculate E[v | S] on interpolation grid
            std::vector<double> cond_exp_v = KernelDensity::estimateConditionalExpectation(
                s_particles, v_particles, s_interp_grid_, h, heston_->getParams().theta
            );

			// Update leverage surface for t_next
            // L(t, S) = sigma_LV(t, S) / sqrt(E[v|S])
            for (size_t j = 0; j < s_grid_size; ++j) {
                double s_val = s_interp_grid_[j];
				double k_val = std::log(s_val / market_->s0) - market_->r * t_next;
                double local_vol = market_->getLocalVol(t_next, k_val);

				double exp_v = std::max(cond_exp_v[j], 1e-8); // Avoid division by zero
                leverage_surface[i + 1][j] = local_vol / std::sqrt(exp_v);
            }
        }

        std::cout << "[C++] Calibration Finished." << std::endl;
        return leverage_surface;
    }

}