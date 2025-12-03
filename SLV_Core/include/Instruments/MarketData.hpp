#pragma once
#include <vector>
#include "Math/Interpolation.hpp"

namespace MyQuantLib {

    struct MarketData {
        double s0;      // Spot price
        double r;       // Risk-free rate
        double q;       // Dividend yield

        std::vector<double> tGrid; // Time grid
        std::vector<double> kGrid; // log Moneyness grid

		// Local Volatility Surface: rows=time, cols=log moneyness
        std::vector<std::vector<double>> localVolSurface;

        BilinearInterpolator lvInterpolator;

        MarketData(double spot, double rate, double div,
            const std::vector<double>& t,
            const std::vector<double>& k,
            const std::vector<std::vector<double>>& lv)
            : s0(spot), r(rate), q(div), tGrid(t), kGrid(k), localVolSurface(lv) {
            lvInterpolator = BilinearInterpolator(tGrid, kGrid, localVolSurface);
        }

        double getLocalVol(double t, double k) const {
            return lvInterpolator.getValue(t, k);
        }
    };

}