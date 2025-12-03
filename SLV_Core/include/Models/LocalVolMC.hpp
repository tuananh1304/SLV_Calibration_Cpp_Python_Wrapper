#pragma once
#include <Instruments/MarketData.hpp>
#include <vector>
#include <memory>
#include <random>

namespace MyQuantLib {

	class LocalVolMC {
		std::shared_ptr<MarketData> market_;
	public:
		explicit LocalVolMC(std::shared_ptr<MarketData> market) : market_(market) {}
		// Simulate one step of local volatility model
		// Hàm chính: Trả về danh sách giá Option tương ứng với danh sách Strikes đầu vào
		std::vector<double> priceOptions(
			double T,
			const std::vector<double>& strikes,
			int N_paths,
			int M_steps = 100
		);
	};
}