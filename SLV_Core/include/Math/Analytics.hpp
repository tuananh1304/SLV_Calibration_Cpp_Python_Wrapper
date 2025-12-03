#pragma once  
#include <cmath>  
#include <algorithm>  
#include <stdexcept>  

#ifndef M_SQRT1_2  
#define M_SQRT1_2 0.7071067811865476 // 1/sqrt(2)  
#endif  

namespace MyQuantLib {  

   class Analytics {  
   public:  
       // Cumulative function (Standard Normal CDF)  
       static double norm_cdf(double x) {  
           return 0.5 * std::erfc(-x * M_SQRT1_2);  
       }  

       // density function (PDF)  
       static double norm_pdf(double x) {  
           return 0.3989422804014327 * std::exp(-0.5 * x * x);  
       }  

       // Black Scholes Call Price  
       static double black_scholes_call(double S, double K, double T, double r, double q, double sigma) {  
           if (T <= 0) return std::max(S - K, 0.0);  
           double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));  
           double d2 = d1 - sigma * std::sqrt(T);  
           return S * std::exp(-q * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);  
       }  

       // Vega for newton-raphson  
       static double black_scholes_vega(double S, double K, double T, double r, double q, double sigma) {  
           if (T <= 0) return 0.0;  
           double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));  
           return S * std::exp(-q * T) * norm_pdf(d1) * std::sqrt(T);  
       }  

       // Find implied volatility using Newton-Raphson  
       static double implied_volatility(double price, double S, double K, double T, double r, double q) {  
           double sigma = 0.3; // Initial guess  
           for (int i = 0; i < 10; ++i) {  
               double p = black_scholes_call(S, K, T, r, q, sigma);  
               double diff = price - p;  
               if (std::abs(diff) < 1e-8) return sigma;  

               double vega = black_scholes_vega(S, K, T, r, q, sigma);  
               if (std::abs(vega) < 1e-12) break; // Avoid division by zero  

               sigma += diff / vega;  
           }  
           return sigma; // Fallback  
       }  
   };  

}