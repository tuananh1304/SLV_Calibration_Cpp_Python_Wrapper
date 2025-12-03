#pragma once
#include "Instruments/MarketData.hpp"
#include "SVI.hpp"
#include <vector>
#include <memory>
#include <map>
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

namespace MyQuantLib {

    class LocalVolBootstrapper {
        std::shared_ptr<MarketData> market_;
        std::map<double, SVIParams> svi_surface_; // Key = Maturity Time

    public:
        LocalVolBootstrapper(std::shared_ptr<MarketData> market) : market_(market) {}

        // Thêm tham số SVI đã fit cho một kỳ hạn
        void addSlice(double t, const SVIParams& params) {
            svi_surface_[t] = params;
        }

        // Main function: Local Vol Surface (Dupire)
        /*std::vector<std::vector<double>> bootstrap();*/

        // Main function: Local Vol Surface (Gatheral by SVI parameterization)
        std::vector<std::vector<double>> bootstrapBySVI();

        void exportToCSV(const std::vector<std::vector<double>>& surface,
            const std::vector<double>& t_grid,
            const std::vector<double>& k_grid,
            const std::string& filename) {

            std::ofstream file(filename);

            if (!file.is_open()) {
                std::cerr << "Error: Could not open file " << filename << std::endl;
                return;
            }

            // 1. Viết Header (Hàng đầu tiên là danh sách k_grid/Strikes)
            file << "Time/Strike"; // Ô góc trên cùng bên trái
            for (double k : k_grid) {
                file << "," << k;
            }
            file << "\n";

            // 2. Viết Dữ liệu (Từng dòng một: Time -> Vol1, Vol2...)
            for (size_t i = 0; i < t_grid.size(); ++i) {
                // Cột đầu tiên là Time
                file << t_grid[i];

                // Các cột tiếp theo là Volatility tại t đó
                for (size_t j = 0; j < k_grid.size(); ++j) {
                    // Quan trọng: Set precision cao để debug lỗi số học nhỏ
                    // std::fixed giúp tránh format dạng 1e-05 khó đọc trong Excel
                    file << "," << std::fixed << std::setprecision(8) << surface[i][j];
                }
                file << "\n";
            }

            file.close();
            std::cout << ">>> Exported debug data to: " << filename << std::endl;
        }

    private:
        // Helper: Interpolate SVI parameters theo thời gian
        SVIParams interpolateParams(double t);

        // Helper: Gaussian Smoothing đơn giản (thay cho scipy.ndimage.gaussian_filter)
        void applySmoothing(std::vector<std::vector<double>>& surface);
    };

}