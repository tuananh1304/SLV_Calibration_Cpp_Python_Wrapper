#pragma once
#include <vector>
#include <cmath>
#include <map>

namespace MyQuantLib {

    struct SVIParams {
        double a, b, rho, m, sigma;
    };

    struct SVITermStructure {
		std::vector<double> maturities;
        std::vector<SVIParams> params;
    };

    //SVIParams interpolateParams(double t, const SVITermStructure& svi_ts) {
    //    if (t <= svi_ts.maturities.front()) return svi_ts.params.front(); // Extrapolate first but keeping the initial params
    //    if (t >= svi_ts.maturities.back()) return svi_ts.params.back();   // Extrapolate last but keeping the last params
    //    // Linear interpolation

    //    auto it = std::upper_bound(svi_ts.maturities.begin(), svi_ts.maturities.end(), t);
    //    size_t idx = it - svi_ts.maturities.begin();
    //    double t2 = svi_ts.maturities[idx];
    //    double t1 = svi_ts.maturities[idx - 1];
    //    double w = (t - t1) / (t2 - t1);
    //    const SVIParams& p1 = svi_ts.params[idx - 1];
    //    const SVIParams& p2 = svi_ts.params[idx];
    //    return {
    //            p1.a * (1 - w) + p2.a * w,
    //            p1.b * (1 - w) + p2.b * w,
    //            p1.rho * (1 - w) + p2.rho * w,
    //            p1.m * (1 - w) + p2.m * w,
    //            p1.sigma * (1 - w) + p2.sigma * w
    //    };
    //}

    class SVIModel {
    public:
        // Raw SVI formula: w(k) = a + b * (rho*(k-m) + sqrt((k-m)^2 + sigma^2))
        static double total_variance(double k, const SVIParams& p) {
            return p.a + p.b * (p.rho * (k - p.m) + std::sqrt((k - p.m) * (k - p.m) + p.sigma * p.sigma));
        }

        // Lấy Implied Vol từ Log-Moneyness k
        static double get_iv(double k, double T, const SVIParams& p) {
            double w = total_variance(k, p);
            if (w < 0) w = 0;
            if (T <= 1e-5) return 0.0;
            return std::sqrt(w / T); // sigma = sqrt(w / T)
        }

		static std::tuple<double, double, double> get_total_variance_and_greeks_K(double k, const SVIParams& p) {
			double w = total_variance(k, p);
            double discr = std::sqrt((k - p.m) * (k - p.m) + p.sigma * p.sigma);
			if (w < 0) w = 0;
			// first order greeks by k: dw/dk
			double dw_dk = p.b * (p.rho + (k - p.m) / discr);
            double d2w_dk2 = p.b * (p.sigma * p.sigma) / (discr * discr * discr);
            return std::make_tuple(w, dw_dk, d2w_dk2);
			
		}
    };

}