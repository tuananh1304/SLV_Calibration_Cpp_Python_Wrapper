//#include <iostream>
//#include <vector>
//#include <memory>
//#include <cmath>
//#include <map>
//#include <iomanip>
//#include <chrono>
//
//// INCLUDE CÁC HEADER CỦA BẠN (Đảm bảo đường dẫn đúng)
//#include "Instruments/MarketData.hpp"
//#include "Models/LocalVolBootstrapper.hpp"
//#include "Models/HestonModel.hpp"
//#include "Models/SLVCalibrator.hpp"
//#include "Models/SVI.hpp"
//
//// Lưu ý: Đảm bảo namespace trùng với code trong thư viện của bạn (QuantLib hay MyQuantLib)
//using namespace MyQuantLib;
//
//// --- 1. Helper Function mô phỏng numpy.linspace ---
//std::vector<double> linspace(double start, double end, int num) {
//    std::vector<double> result;
//    if (num <= 0) return result;
//    if (num == 1) {
//        result.push_back(start);
//        return result;
//    }
//    double step = (end - start) / (num - 1);
//    for (int i = 0; i < num; ++i) {
//        result.push_back(start + i * step);
//    }
//    return result;
//}
//
//int main() {
//    try {
//        std::cout << "=== SLV C++ Reproduction of Python Workflow ===" << std::endl;
//
//        // ==========================================================
//        // STEP 1: PARAMETERS SETUP (Tương ứng Python dictionary)
//        // ==========================================================
//
//        // Market Params
//        double s0 = 100.0;
//        double r = 0.02;
//        double q = 0.01;
//
//        // Grids: np.linspace(0, 1.0, 21) & np.linspace(60, 180, 51)
//        auto t_grid = linspace(0.0, 1.0, 21);
//        auto s_grid = linspace(60.0, 180.0, 51);
//
//        // Heston Params
//        HestonParams heston_p;
//        heston_p.v0 = 0.04;
//        heston_p.kappa = 2.0;
//        heston_p.theta = 0.04;
//        heston_p.xi = 0.5;
//        heston_p.rho = -0.7;
//
//        // Numerical Params
//        double T_calib = 1.0;
//        int M_steps = 20;
//        int N_particles = 50000;
//
//        std::cout << "[1] Parameters initialized." << std::endl;
//
//        // ==========================================================
//        // STEP 2 & 3: BOOTSTRAP LOCAL VOL FROM SVI
//        // (Bỏ qua bước tạo giá giả và fit lại, dùng luôn True Params)
//        // ==========================================================
//
//        // Tạo MarketData tạm thời (với dummy local vol) để Bootstrapper dùng các grid/rate
//        // Dummy LV = 0.2
//        std::vector<std::vector<double>> dummy_lv(t_grid.size(), std::vector<double>(s_grid.size(), 0.2));
//        auto temp_market = std::make_shared<MarketData>(s0, r, q, t_grid, s_grid, dummy_lv);
//
//        // Khởi tạo Bootstrapper
//        LocalVolBootstrapper bootstrapper(temp_market);
//
//        // Nạp True SVI Params (Python: true_svi_params)
//        // Format SVIParams struct: {a, b, rho, m, sigma}
//        // Python: 0.25: (0.04, 0.2, -0.8, 0.1, 0.3)
//        bootstrapper.addSlice(0.25, { 0.04, 0.2, -0.8, 0.1, 0.3 });
//
//        // Python: 0.5: (0.05, 0.18, -0.75, 0.05, 0.25)
//        bootstrapper.addSlice(0.50, { 0.05, 0.18, -0.75, 0.05, 0.25 });
//
//        // Python: 1.0: (0.06, 0.15, -0.7, 0.0, 0.2)
//        bootstrapper.addSlice(1.00, { 0.06, 0.15, -0.7, 0.0, 0.2 });
//
//        std::cout << "[2] Bootstrapping Local Volatility Surface..." << std::endl;
//        auto start_boot = std::chrono::high_resolution_clock::now();
//
//        // Chạy Dupire Formula
//        auto target_lv_surface = bootstrapper.bootstrap();
//
//        auto end_boot = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double> diff_boot = end_boot - start_boot;
//        std::cout << "    -> Finished in " << diff_boot.count() << " s" << std::endl;
//
//        // ==========================================================
//        // STEP 4: PREPARE FINAL ENVIRONMENT
//        // ==========================================================
//
//        // Tạo MarketData CHÍNH THỨC với Local Vol vừa tính được
//        auto market_env = std::make_shared<MarketData>(
//            s0, r, q, t_grid, s_grid, target_lv_surface
//        );
//
//        auto heston_process = std::make_shared<HestonProcess>(heston_p);
//
//        // ==========================================================
//        // STEP 5: RUN PARTICLE METHOD CALIBRATION
//        // ==========================================================
//
//        std::cout << "[3] Running SLV Particle Calibration (N=" << N_particles << ")..." << std::endl;
//
//        SLVCalibrator calibrator(
//            market_env,
//            heston_process,
//            T_calib,
//            M_steps,
//            N_particles,
//            s_grid
//        );
//
//        auto start_calib = std::chrono::high_resolution_clock::now();
//
//        auto leverage_surface = calibrator.calibrate();
//
//        auto end_calib = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double> diff_calib = end_calib - start_calib;
//
//        std::cout << "    -> Calibration DONE in " << diff_calib.count() << " s" << std::endl;
//
//        // ==========================================================
//        // STEP 6: VERIFICATION OUTPUT
//        // ==========================================================
//        std::cout << "\n=== RESULTS PREVIEW ===" << std::endl;
//        std::cout << "Format: [Time Index] -> Spot=100 (At-The-Money)" << std::endl;
//
//        // Tìm index của Spot ~ 100
//        size_t s_idx = 0;
//        for (size_t i = 0; i < s_grid.size(); ++i) {
//            if (s_grid[i] >= 100.0) { s_idx = i; break; }
//        }
//
//        std::cout << std::left << std::setw(10) << "Time(yr)"
//            << std::setw(15) << "Target LV(%)"
//            << std::setw(15) << "Leverage L(t,S)" << std::endl;
//        std::cout << "----------------------------------------" << std::endl;
//
//        for (int i = 0; i <= M_steps; i += 5) { // In mỗi 5 bước
//            double t = t_grid[i]; // Giả sử lưới thời gian khớp nhau
//            double lv = target_lv_surface[i][s_idx];
//            double lev = leverage_surface[i][s_idx];
//
//            std::cout << std::left << std::setw(10) << t
//                << std::setw(15) << lv
//                << std::setw(15) << lev << std::endl;
//        }
//
//    }
//    catch (const std::exception& e) {
//        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
//        return 1;
//    }
//
//    return 0;
//}