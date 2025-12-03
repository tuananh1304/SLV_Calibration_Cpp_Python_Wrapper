#pragma once
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace MyQuantLib {
	// Linear interpolation class
	class LinearInterpolator {
		std::vector<double> x_, y_;
	public:
		LinearInterpolator() = default;
		LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y) :x_(x), y_(y) {
			if (x.size() != y.size() || x.size() < 2) {
				throw std::runtime_error("Input vectors must have the same size and contain at least two points.");
			}

			if (!std::is_sorted(x.begin(), x.end())) {
				std::vector<std::pair<double, double>> xyPairs;
				auto size = x.size();
				for (size_t i = 0; i < size; ++i) {
					xyPairs.emplace_back(x[i], y[i]);
				}

				std::sort(xyPairs.begin(), xyPairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
				x_.resize(size);
				y_.resize(size);
				for (size_t i = 0; i < size; ++i) {
					x_[i] = xyPairs[i].first;
					y_[i] = xyPairs[i].second;
				}
			}

		}

		double operator()(double xVal) const {
			if (xVal <= x_.front()) {
				return y_.front();
			}

			if (xVal >= x_.back()) {
				return y_.back();
			}

			auto it = std::lower_bound(x_.begin(), x_.end(), xVal);
			size_t idx = std::distance(x_.begin(), it) - 1;
			double t = (xVal - x_[idx]) / (x_[idx + 1] - x_[idx]);
			return y_[idx] + t * (y_[idx + 1] - y_[idx]);
		}
	};

	class BilinearInterpolator {
		std::vector<double> x_, y_;
		std::vector<std::vector<double>> z_; // Grid z[x][y]
	public:
		BilinearInterpolator() = default;
		BilinearInterpolator(const std::vector<double>& x,
			const std::vector<double>& y,
			const std::vector<std::vector<double>>& z)
			: x_(x), y_(y), z_(z) {
		}

		double getValue(double xVal, double yVal) const {
			// clamp to grid bounds
			double x = std::max(x_.front(), std::min(xVal, x_.back()));
			double y = std::max(y_.front(), std::min(yVal, y_.back()));

			// find grid cell
			auto itx = std::lower_bound(x_.begin(), x_.end(), x);
			size_t i = (itx == x_.begin()) ? 0 : std::distance(x_.begin(), itx) - 1;
			if (i >= x_.size() - 1) i = x_.size() - 2;

			auto ity = std::lower_bound(y_.begin(), y_.end(), y);
			size_t j = (ity == y_.begin()) ? 0 : std::distance(y_.begin(), ity) - 1;
			if (j >= y_.size() - 1) j = y_.size() - 2;

			// interpolation
			double x1 = x_[i], x2 = x_[i + 1];
			double y1 = y_[j], y2 = y_[j + 1];

			double q11 = z_[i][j];
			double q21 = z_[i + 1][j];
			double q12 = z_[i][j + 1];
			double q22 = z_[i + 1][j + 1];

			double denom = (x2 - x1) * (y2 - y1);

			return (q11 * (x2 - x) * (y2 - y) +
				q21 * (x - x1) * (y2 - y) +
				q12 * (x2 - x) * (y - y1) +
				q22 * (x - x1) * (y - y1)) / denom;
		}
	};
}