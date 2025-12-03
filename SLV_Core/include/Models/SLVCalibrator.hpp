#pragma once
#include "Instruments/MarketData.hpp"
#include "Models/HestonModel.hpp"
#include <vector>
#include <memory>

namespace MyQuantLib {

    class SLVCalibrator {
        std::shared_ptr<MarketData> market_;
        std::shared_ptr<HestonProcess> heston_;

        // Simulation parameters
        double T_;
        int M_steps_;
        int N_particles_;

        // Grids
        std::vector<double> time_grid_;
        std::vector<double> s_interp_grid_;

    public:
        SLVCalibrator(
            std::shared_ptr<MarketData> market,
            std::shared_ptr<HestonProcess> heston,
            double T, int M, int N,
            const std::vector<double>& s_interp_grid
        );

		// Return the calibrated leverage surface
        std::vector<std::vector<double>> calibrate();
    };

}