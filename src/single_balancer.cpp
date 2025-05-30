#include <Rcpp.h>
#include <cmath>
#include <chrono>
#include <vector>
#include "utils.h"
using namespace Rcpp;

//' @title Entropy Balancer
//' @description Perform entropy balancing on a given incidence matrix with specified control totals and weights.
//' @param incidence Numeric matrix of incidence data (samples x controls).
//' @param initial_weights Numeric vector of initial weights for each sample.
//' @param control_totals Numeric vector of target totals for each control.
//' @param weights_lb Numeric vector of lower bounds for weights (optional).
//' @param weights_ub Numeric vector of upper bounds for weights (optional).
//' @param controls_importance Numeric vector of importance weights for each control (optional).
//' @param master_control_index Integer index of the master control (optional, default -1).
//' @param max_iterations Maximum number of iterations to run (default 10000).
//' @param max_delta Maximum allowed change in weights for convergence (default 1e-9).
//' @param print_every_n Print progress every n iterations (default 10).
//' @return A list containing the final weights, relaxation factors, convergence status, number of iterations, delta, and max gamma difference.
//' @export
// [[Rcpp::export]]
List entropy_balancer(
  const NumericMatrix& incidence,
  const NumericVector& initial_weights,
  const NumericVector& control_totals,
  const NumericVector& weights_lb = NumericVector::create(),
  const NumericVector& weights_ub = NumericVector::create(),
  const NumericVector& controls_importance = NumericVector::create(),
  int master_control_index = -1,
  int max_iterations = 10000,
  double max_delta = 1e-8,
  int print_every_n = 10
) {
  const int control_count = incidence.ncol();
  const int sample_count = incidence.nrow();
  const double* inc_ptr = incidence.begin(); // column-major layout

  if (!Rf_isMatrix(incidence)) stop("incidence must be a numeric matrix");
  if (initial_weights.size() != sample_count)
    stop("initial_weights must match number of rows in incidence matrix");
  if (control_totals.size() != control_count)
    stop("control_totals must match number of columns in incidence matrix");
  if (weights_lb.size() && weights_lb.size() != sample_count)
    stop("weights_lb must match initial_weights");
  if (weights_ub.size() && weights_ub.size() != sample_count)
    stop("weights_ub must match initial_weights");
  if (controls_importance.size() && controls_importance.size() != control_count)
    stop("controls_importance must match control_count");

  NumericVector lb = weights_lb.size() ? clone(weights_lb) : rep(0.0, sample_count);
  NumericVector ub = weights_ub.size() ? clone(weights_ub) : rep(R_PosInf, sample_count);
  NumericVector importance = controls_importance.size() ? clone(controls_importance) : rep(1.0, control_count);
  NumericVector weights = clone(initial_weights);
  NumericVector relaxation(control_count, 1.0);

  std::vector<double> incidence2(control_count * sample_count);
  for (int j = 0; j < control_count; ++j)
    for (int i = 0; i < sample_count; ++i) {
      double val = inc_ptr[j * sample_count + i];
      incidence2[j * sample_count + i] = val * val;
    }

  IntegerVector control_indexes(control_count);
  int k = 0;
  for (int i = 0; i < control_count; ++i)
    if (i != master_control_index) control_indexes[k++] = i;
  if (master_control_index >= 0) control_indexes[k++] = master_control_index;

  double* weights_ptr = weights.begin();
  const double* lb_ptr = lb.begin();
  const double* ub_ptr = ub.begin();

  double delta = 0.0, max_gamma_diff = 0.0, importance_adjustment = 1.0;
  bool converged = false;
  int iter = 0;

  auto start = std::chrono::high_resolution_clock::now();

  for (; iter < max_iterations; ++iter) {
    delta = 0.0;
    NumericVector gamma(control_count, 1.0);  // Fresh allocation each iteration

    if (iter > 0 && iter % IMPORTANCE_ADJUST_COUNT == 0)
      importance_adjustment /= IMPORTANCE_ADJUST;

    for (int i = 0; i < control_count; ++i) {
      int c = control_indexes[i];
      int offset_c = c * sample_count;
      double xx = 0.0, yy = 0.0;

      for (int j = 0; j < sample_count; ++j) {
        double w = weights_ptr[j];
        double inc = inc_ptr[offset_c + j];
        xx += w * inc;
        yy += w * incidence2[offset_c + j];
      }

      double imp = (c == master_control_index)
        ? importance[c]
        : std::max(importance[c] * importance_adjustment, MIN_IMPORTANCE);

      if (xx <= 0.0 || !R_finite(xx) || !R_finite(yy)) continue;

      double relaxed = std::max(control_totals[c] * relaxation[c], MIN_CONTROL_VALUE);
      double gamma_val = 1.0 - (xx - relaxed) / (yy + relaxed / imp);
      gamma_val = std::max(gamma_val, MIN_GAMMA);
      if (!R_finite(gamma_val)) gamma_val = MIN_GAMMA;

      gamma[c] = gamma_val;
      double log_gamma = (gamma_val == 1.0) ? 0.0 : std::log(gamma_val);

      for (int j = 0; j < sample_count; ++j) {
        double old_w = weights_ptr[j];
        double inc = inc_ptr[offset_c + j];
        double new_w = old_w * std::exp(log_gamma * inc);

        // Clamp
        new_w = std::min(std::max(new_w, lb_ptr[j]), ub_ptr[j]);

        double diff = new_w - old_w;
        delta += (diff >= 0.0) ? diff : -diff;

        weights_ptr[j] = new_w;
      }

      // Always apply relaxation update (matches stable baseline)
      double relax_factor = relaxation[c] * std::pow(1.0 / gamma_val, 1.0 / imp);
      relaxation[c] = std::min(relax_factor, MAX_RELAXATION_FACTOR);
    }

    delta /= sample_count;

    max_gamma_diff = 0.0;
    for (int c = 0; c < control_count; ++c)
      max_gamma_diff = std::max(max_gamma_diff, fabs(gamma[c] - 1.0));

    if (print_every_n > 0 && (iter + 1) % print_every_n == 0)
      Rcpp::Rcout << "[Iter " << (iter + 1) << "] delta = " << delta << ", max_gamma_diff = " << max_gamma_diff << "\n";

    converged = (delta < max_delta && max_gamma_diff < MAX_GAMMA_DIFF);
    bool stalled = (delta < ALT_MAX_DELTA);
    if (converged || stalled) break;
  }

  auto end = std::chrono::high_resolution_clock::now();
  double elapsed_time = std::chrono::duration<double>(end - start).count();

  Rcpp::Rcout << (converged ? "Converged" : "Stalled") << " after " << (iter + 1) << " iterations.\n";
  Rcpp::Rcout << "Elapsed time: " << elapsed_time << " seconds.\n";

  return make_result(weights, relaxation, converged, iter, delta, max_gamma_diff, elapsed_time);
}
