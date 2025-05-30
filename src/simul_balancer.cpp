#include <Rcpp.h>
#include <cmath>
#include <chrono>
#include "utils.h"
using namespace Rcpp;

//' @title Simultaneos Entropy Balancer
//' @description This function performs simultaneous entropy balancing on a set of samples across multiple zones, i.e., 2D balancing.
//' @param incidence A matrix of incidence values with dimensions (control_count × sample_count).
//' @param sub_weights A matrix of weights with dimensions (zone_count × sample_count) (mutable).
//' @param sub_controls A matrix of control values with dimensions (zone_count × control_count).
//' @param parent_controls A vector of parent control values with length (control_count).
//' @param parent_weights A vector of parent weights with length (sample_count).
//' @param weights_lb A vector of lower bounds for weights with length (sample_count).
//' @param weights_ub A vector of upper bounds for weights with length (sample_count).
//' @param controls_importance A vector of importance values for each control with length (control_count).
//' @param master_control_index An integer index for the master control, -1 means no master control.
//' @param max_iterations An integer specifying the maximum number of iterations to perform.
//' @param max_delta A double specifying the maximum allowed change in weights for convergence.
//' @param print_every_n An integer specifying how often to print progress information.
//' @return A list containing the updated weights, relaxation matrix, convergence status, number of iterations, delta, max gamma difference, and elapsed time.
//' @export
// [[Rcpp::export]]
List simul_entropy_balancer(
  NumericMatrix incidence,                  // [control_count × sample_count]  
  NumericMatrix sub_weights,                // [zone_count × sample_count] (mutable)  
  NumericMatrix sub_controls,               // [zone_count × control_count]  
  NumericVector parent_controls,            // [control_count]
  NumericVector parent_weights,             // [sample_count]
  NumericVector weights_lb,                 // [sample_count]
  NumericVector weights_ub,                 // [sample_count]
  NumericVector controls_importance,        // [control_count]  
  int master_control_index = -1,            // -1 means no master control
  int max_iterations = 10000,
  double max_delta = 1e-8,
  int print_every_n = 10
) {
  auto start = std::chrono::high_resolution_clock::now();

  // infer counts from dimensions
  int sample_count = incidence.ncol();
  int control_count = incidence.nrow();
  int zone_count = sub_weights.nrow();
  double delta;
  double max_gamma_diff;
  bool converged = false;
  int iter = 0;

  if (print_every_n == 0) print_every_n = -1;

  // Precompute incidence^2
  NumericMatrix incidence2(control_count, sample_count);
  for (int c = 0; c < control_count; ++c)
    for (int j = 0; j < sample_count; ++j)
      incidence2(c, j) = incidence(c, j) * incidence(c, j);

  // Prepare relaxation_matrix and gamma
  NumericMatrix relaxation_matrix(zone_count, control_count);
  std::fill(relaxation_matrix.begin(), relaxation_matrix.end(), 1.0);

  NumericMatrix gamma(zone_count, control_count);
  std::fill(gamma.begin(), gamma.end(), 1.0);

  // Prepare control reordering
  IntegerVector control_indexes(control_count);
  int k = 0;
  for (int i = 0; i < control_count; ++i)
    if (i != master_control_index) control_indexes[k++] = i;
  if (master_control_index >= 0) control_indexes[k] = master_control_index;

  NumericMatrix weights_previous(zone_count, sample_count);
  double importance_adjustment = 1.0;

  for (; iter < max_iterations; ++iter) {
    std::copy(sub_weights.begin(), sub_weights.end(), weights_previous.begin());

    if (iter > 0 && iter % IMPORTANCE_ADJUST_COUNT == 0)
      importance_adjustment /= IMPORTANCE_ADJUST;

    for (int idx = 0; idx < control_count; ++idx) {
      int c = control_indexes[idx];

      double imp = (c == master_control_index)
        ? controls_importance[c]
        : std::max(controls_importance[c] * importance_adjustment, MIN_IMPORTANCE);

      for (int z = 0; z < zone_count; ++z) {
        double xx = 0.0, yy = 0.0;
        for (int j = 0; j < sample_count; ++j) {
          double w = sub_weights(z, j);
          double inc = incidence(c, j);
          xx += w * inc;
          yy += w * incidence2(c, j);
        }

        if (xx <= 0.0) continue;

        double relaxed = sub_controls(z, c) * relaxation_matrix(z, c);
        relaxed = std::max(relaxed, MIN_CONTROL_VALUE);

        double gamma_val = 1.0 - (xx - relaxed) / (yy + relaxed / imp);
        gamma_val = std::max(gamma_val, MIN_GAMMA);
        gamma(z, c) = gamma_val;

        double log_gamma = std::log(gamma_val);
        for (int j = 0; j < sample_count; ++j) {
          double w_old = sub_weights(z, j);
          double inc = incidence(c, j);
          double new_w = w_old * std::exp(log_gamma * inc);
          new_w = std::max(weights_lb[j], std::min(new_w, weights_ub[j]));
          sub_weights(z, j) = new_w;
        }

        double rf = relaxation_matrix(z, c) * std::pow(1.0 / gamma_val, 1.0 / imp);
        relaxation_matrix(z, c) = std::min(rf, MAX_RELAXATION_FACTOR);
      }
    }

    // Rescale each column to match parent_weights
    for (int j = 0; j < sample_count; ++j) {
      double total = 0.0;
      for (int z = 0; z < zone_count; ++z)
        total += sub_weights(z, j);
      if (total > 0.0) {
        double scale = parent_weights[j] / total;
        for (int z = 0; z < zone_count; ++z)
          sub_weights(z, j) *= scale;
      }
    }

    // Convergence check
    delta = 0.0;
    for (int z = 0; z < zone_count; ++z)
      for (int j = 0; j < sample_count; ++j)
        delta += std::abs(sub_weights(z, j) - weights_previous(z, j));
    delta /= sample_count;

    max_gamma_diff = 0.0;
    for (int z = 0; z < zone_count; ++z)
      for (int c = 0; c < control_count; ++c)
        max_gamma_diff = std::max(max_gamma_diff, std::abs(gamma(z, c) - 1.0));

    if (print_every_n > 0 && (iter + 1) % print_every_n == 0) {
      Rcpp::Rcout << "[Iter " << (iter + 1) << "] delta = " << delta
                  << ", max_gamma_diff = " << max_gamma_diff << "\n";
    }

    converged = (delta < max_delta) && (max_gamma_diff < MAX_GAMMA);
    bool no_progress = delta < ALT_MAX_DELTA;

    if (converged || no_progress) break;
  }

  auto end = std::chrono::high_resolution_clock::now();
  double elapsed_time = std::chrono::duration<double>(end - start).count();
  Rcpp::Rcout << (converged ? "Converged" : "Stalled") << " after " << (iter + 1) << " iterations.\n";
  Rcpp::Rcout << "Elapsed time: " << elapsed_time << " seconds.\n";

  return make_result(sub_weights, relaxation_matrix, converged, iter, delta, max_gamma_diff, elapsed_time);
}
