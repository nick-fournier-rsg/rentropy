// shared_utils.h
#ifndef UTILS_H
#define UTILS_H
#include <Rcpp.h>
using namespace Rcpp;

const double ALT_MAX_DELTA = 1.0e-14;
const double MAX_GAMMA_DIFF = 1.0e-5;
const double MAX_GAMMA = 1.0e-8;
const double MIN_GAMMA = 1.0e-10;
const double IMPORTANCE_ADJUST = 2.0;
const int IMPORTANCE_ADJUST_COUNT = 100;
const double MIN_IMPORTANCE = 1.0;
const double MAX_RELAXATION_FACTOR = 1.0e6;
const double MIN_CONTROL_VALUE = 0.1;
const double MAX_GAMMA_CLAMP = 10.0;


// For single-zone version (vectors)
template <typename T>
inline List make_result(const T& weights,
                        const T& relaxation,
                        bool converged,
                        int iterations,
                        double delta,
                        double max_gamma_diff,
                        double elapsed_time
                    ) {
  return List::create(
    _["weights"] = weights,
    _["relaxation_factors"] = relaxation,
    _["converged"] = converged,
    _["iterations"] = iterations,
    _["delta"] = delta,
    _["max_gamma_diff"] = max_gamma_diff,
    _["elapsed_time"] = elapsed_time
  );
}

#endif