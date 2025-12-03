#include "Math/Statistics.hpp"
#include <algorithm>
#include <omp.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MyQuantLib {

   double KernelDensity::silvermanBandwidth(const std::vector<double>& data) {
       if (data.empty()) return 1.0;

       double sum = 0.0, sumSq = 0.0;
       for (double x : data) {
           sum += x;
           sumSq += x * x;
       }
       double mean = sum / data.size();
       double var = (sumSq / data.size()) - (mean * mean);
       double stdDev = std::sqrt(var);

       // h = 1.06 * sigma * n^(-1/5)
       return 1.06 * stdDev * std::pow((double)data.size(), -0.2);
   }

   std::vector<double> KernelDensity::estimateConditionalExpectation(
       const std::vector<double>& s_particles,
       const std::vector<double>& v_particles,
       const std::vector<double>& s_grid,
       double bandwidth,
       double fallback_value)
   {
       size_t N = s_particles.size();
       size_t M = s_grid.size();
       std::vector<double> result(M);

       // Avoid dividing by 0 by using a floor here
       double h = std::max(bandwidth, 1e-6);
       double inv_h_sq_2 = 1.0 / (2.0 * h * h);
       double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

       // use OpenMP to accelerate looping though grid_s
       #pragma omp parallel for schedule(static)
       for (int i = 0; i < (int)M; ++i) {
           double s_target = s_grid[i];
           double num = 0.0;
           double den = 0.0;

           for (size_t j = 0; j < N; ++j) {
               double diff = s_target - s_particles[j];
               // Gaussian Kernel: K(u) = exp(-0.5 * u^2)
               // Bỏ qua hằng số chuẩn hóa 1/(h*sqrt(2pi)) ở cả tử và mẫu để tối ưu
               double kernel = std::exp(-(diff * diff) * inv_h_sq_2);

               num += kernel * v_particles[j];
               den += kernel;
           }

           if (den > 1e-10) {
               result[i] = num / den;
           }
           else {
               result[i] = fallback_value; // fall back for boundary values
           }
       }
       return result;
   }

}