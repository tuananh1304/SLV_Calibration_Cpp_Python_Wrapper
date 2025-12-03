#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Quan trọng: Để chuyển đổi tự động vector <-> list

// Include đầy đủ các file header
#include "Instruments/MarketData.hpp"
#include "Models/HestonModel.hpp"
#include "Models/SLVCalibrator.hpp"
#include "Models/SVI.hpp"
#include "Models/LocalVolBootstrapper.hpp"
#include "Models/LocalVolMC.hpp"
#include "Math/Analytics.hpp"

namespace py = pybind11;
using namespace MyQuantLib;

PYBIND11_MODULE(slv_engine, m) {
    m.doc() = "Stochastic Local Volatility Engine (Python Bindings)";

    // --- 1. DATA STRUCTURES ---

    // MarketData (Giữ nguyên)
    py::class_<MarketData, std::shared_ptr<MarketData>>(m, "MarketData")
        .def(py::init<double, double, double,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<std::vector<double>>&>())
        .def_readonly("s0", &MarketData::s0) // Cho phép Python đọc s0 nếu cần
        .def_readonly("t_grid", &MarketData::tGrid)
        .def_readonly("k_grid", &MarketData::kGrid);

    // HestonParams (Giữ nguyên)
    py::class_<HestonParams>(m, "HestonParams")
        .def(py::init<>()) // Constructor mặc định
        .def(py::init<double, double, double, double, double>())
        .def_readwrite("v0", &HestonParams::v0)
        .def_readwrite("kappa", &HestonParams::kappa)
        .def_readwrite("theta", &HestonParams::theta)
        .def_readwrite("xi", &HestonParams::xi)
        .def_readwrite("rho", &HestonParams::rho);

    // SVIParams (MỚI THÊM)
    py::class_<SVIParams>(m, "SVIParams")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double>())
        .def_readwrite("a", &SVIParams::a)
        .def_readwrite("b", &SVIParams::b)
        .def_readwrite("rho", &SVIParams::rho)
        .def_readwrite("m", &SVIParams::m)
        .def_readwrite("sigma", &SVIParams::sigma);

	// SVIModel
	m.def("get_iv", SVIModel::get_iv, "Get Implied Volatility from Log-Moneyness",
		    py::arg("k"),
		    py::arg("T"),
		    py::arg("params"));

    // --- 2. MODELS & PROCESSES ---

    // HestonProcess
    py::class_<HestonProcess, std::shared_ptr<HestonProcess>>(m, "HestonProcess")
        .def(py::init<const HestonParams&>());

	py::class_<LocalVolMC, std::shared_ptr<LocalVolMC>>(m, "LocalVolMC")
		.def(py::init<std::shared_ptr<MarketData>>())
		.def("price_options", &LocalVolMC::priceOptions, "Price options under Local Vol Model",
			py::arg("T"),
			py::arg("strikes"),
			py::arg("N_paths"),
			py::arg("M_steps") = 100);

    // --- 3. ENGINES & CALCULATORS ---

    // LocalVolBootstrapper
    // Calculate Surface from SVI Params
    py::class_<LocalVolBootstrapper>(m, "LocalVolBootstrapper")
        .def(py::init<std::shared_ptr<MarketData>>())
        .def("add_slice", &LocalVolBootstrapper::addSlice, "Add SVI params for a specific maturity")
        .def("bootstrapBySVI", &LocalVolBootstrapper::bootstrapBySVI, "Calculate and return Local Vol Surface (2D Array)");

    // SLVCalibrator 
    py::class_<SLVCalibrator>(m, "SLVCalibrator")
        .def(py::init<std::shared_ptr<MarketData>,
            std::shared_ptr<HestonProcess>,
            double, int, int,
            const std::vector<double>&>())
        .def("calibrate", &SLVCalibrator::calibrate, "Run Monte Carlo Particle Method");

    // Analytical Black Schole Prices
	py::class_<Analytics>(m, "Analytics")
		.def_static("black_scholes_call", &Analytics::black_scholes_call, "Black-Scholes Call Price",
			py::arg("s0"),
			py::arg("K"),
			py::arg("T"),
			py::arg("r"),
			py::arg("q"),
			py::arg("sigma"))
		.def_static("implied_volatility", &Analytics::implied_volatility, "Calculate Implied Volatility",
			py::arg("price"),
			py::arg("s0"),
			py::arg("K"),
			py::arg("T"),
			py::arg("r"),
			py::arg("q"));

}