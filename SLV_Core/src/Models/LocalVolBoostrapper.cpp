#include "Models/LocalVolBootstrapper.hpp"
#include "Math/Analytics.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

namespace MyQuantLib {

    SVIParams LocalVolBootstrapper::interpolateParams(double t) {
        if (svi_surface_.empty()) return { 0,0,0,0,0 };

        auto it = svi_surface_.lower_bound(t);
        if (it == svi_surface_.begin()) return it->second;
        if (it == svi_surface_.end()) return svi_surface_.rbegin()->second;

        double t2 = it->first;
        double t1 = std::prev(it)->first;
        const auto& p1 = std::prev(it)->second;
        const auto& p2 = it->second;

        double w = (t - t1) / (t2 - t1);

        // Linear interpolation of parameters (Simplified)
        return {
            p1.a * (1 - w) + p2.a * w,
            p1.b * (1 - w) + p2.b * w,
            p1.rho * (1 - w) + p2.rho * w,
            p1.m * (1 - w) + p2.m * w,
            p1.sigma * (1 - w) + p2.sigma * w
        };
    }

    void LocalVolBootstrapper::applySmoothing(std::vector<std::vector<double>>& surface) {
        //Binominal smoothing, can add other type of smoothing here, for future version maybe
        auto temp = surface;
        auto rows = surface.size();
        auto cols = surface[0].size();

        for (int i = 1; i < rows - 1; ++i) {
            for (int j = 1; j < cols - 1; ++j) {
                surface[i][j] = 0.25 * temp[i][j] +
                    0.125 * (temp[i - 1][j] + temp[i + 1][j] + temp[i][j - 1] + temp[i][j + 1]) +
                    0.0625 * (temp[i - 1][j - 1] + temp[i - 1][j + 1] + temp[i + 1][j - 1] + temp[i + 1][j + 1]);
            }
        }
    }

    std::vector<std::vector<double>> LocalVolBootstrapper::bootstrapBySVI() {
        std::cout << "[C++] Bootstrapping Local Volatility using Gatheral's SVI method..." << std::endl;
        const auto& t_grid = market_->tGrid;
        const auto& k_grid = market_->kGrid;
        double s0 = market_->s0;
        double r = market_->r;
        double q = market_->q;
        auto n_t = t_grid.size();
        auto n_k = k_grid.size();

        std::vector<std::vector<double>> loc_var(n_t, std::vector<double>(n_k));
        for (int i = 0; i < n_t; ++i) {
            double t = t_grid[i];
            for (int j = 0; j < n_k; ++j) {
                double k = k_grid[j]; // Log-moneyness
                //double iv = iv_surface[i][j];

                auto [w, dw_dk, d2w_dk2] = SVIModel::get_total_variance_and_greeks_K(k, interpolateParams(t));
                w = std::max(w, 1e-8);

                // dw/dt

                double dt_bump = 0.0001; // Bước nhảy nhỏ 1e-3

                SVIParams p_up = interpolateParams(t + dt_bump);
                SVIParams p_down = interpolateParams(std::max(0.0001, t - dt_bump));

                double w_up = SVIModel::total_variance(k, p_up);
                double w_down = SVIModel::total_variance(k, p_down);

                double dw_dt = (w_up - w_down) / (2 * dt_bump);

                if (dw_dt <= 0.0) {
                    throw std::runtime_error("SVI no-arbitrage violation: dw/dt must be > 0");
                }

                // Nếu ở sát biên thời gian 0, dùng sai phân tiến (Forward diff)
                //if (t < dt_bump) {
                //    dw_dt = (w_up - SVIModel::total_variance(k, interpolateParams(0.0001))) / dt_bump;
                //}

                // Gatheral's Local Vol Formula

                double numerator = dw_dt;
                // Denominator
                double term1 = 1 - (k / w) * dw_dk;
                double term2 = 0.25 * (-0.25 - (1.0 / w) + (k * k) / (w * w)) * (dw_dk * dw_dk);
                double term3 = 0.5 * d2w_dk2;

                double denominator = term1 + term2 + term3;

                if (denominator > 1e-8) {
                    double local_var = numerator / denominator;
                    loc_var[i][j] = std::max(local_var, 1e-6);
                }
                else {
                    loc_var[i][j] = w / t; // Fallback to implied vol if denominator is too small*
                }
            }
        }
        applySmoothing(loc_var);

		// Convert local variance to local volatility

		std::vector<std::vector<double>> loc_vol(n_t, std::vector<double>(n_k));    
		for (int i = 0; i < n_t; ++i) {
			for (int j = 0; j < n_k; ++j) {
				loc_vol[i][j] = std::sqrt(loc_var[i][j]);
			}
		}
        return loc_vol;

    }
}