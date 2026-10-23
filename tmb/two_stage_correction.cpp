// Second-stage geostatistical correction to the dynamical model (issue #21).
//
// Beta-binomial model for the bioassay counts, died_i of n_i, with mean
// ilogit(lambda_i) and intra-class correlation rho_i:
//
//   lambda_i = m_i + omega(s_i) + xi(s_i, t_i) + sum_k b_k[j_k(i)]
//
// where m_i is the dynamical model's logit prediction (an offset, never
// updated), omega is a Matern field correcting the initial conditions, xi is
// the accumulated sum since t0 of an AR(1)-in-time, Matern-in-space field of
// annual selection anomalies eta, and the b_k are iid noise effects: in the
// final model u, per pixel-year, p, per pixel, and sy, per study-year
// (R/two_stage_correction.R). The prediction target is m + omega + xi; the
// b_k enter predictions of new assays only as fresh draws. TMB's Laplace
// approximation integrates the random effects out (#34).
//
// Random effects:
//   w_omega  node values of omega                        (n_nodes)
//   x        node values of xi in years t0+1, ..., T     (n_nodes_xi x n_years)
//   iid      the iid effects, every term's levels stacked (n_levels)
//
// xi(., t0) = 0 is not a parameter: it enters only as the zero that the first
// time difference is taken from. x rather than eta is the random effect so
// that each observation touches the nodes of a single year; the map
// x -> eta (eta_t = x_t - x_{t-1}) is unit lower triangular, so its Jacobian
// is 1 and the prior on eta is the prior on x.

#include <TMB.hpp>

using namespace density;
using namespace R_inla;
using namespace Eigen;

// log density of log sigma when sigma is half-normal with the given scale,
// including the Jacobian |d sigma / d log sigma| = sigma
template<class Type>
Type log_half_normal_log_sd(Type log_sigma, Type scale) {
  Type sigma = exp(log_sigma);
  return log(Type(2.0)) + dnorm(sigma, Type(0.0), scale, true) + log_sigma;
}

// SPDE precision (alpha = 2, d = 2) scaled to have marginal SD sigma
template<class Type>
SparseMatrix<Type> matern_precision(spde_t<Type> spde, Type kappa,
                                    Type sigma) {
  Type tau_spde = Type(1.0) / (sigma * kappa * sqrt(Type(4.0 * M_PI)));
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  return Q * (tau_spde * tau_spde);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(m);                 // dynamical model logit prediction (offset)
  DATA_SPARSE_MATRIX(A_omega);    // observations -> omega mesh nodes
  DATA_SPARSE_MATRIX(A_xi);       // observations -> vec(x) (node x year)
  DATA_IMATRIX(iid_index);        // 0-based level in iid, observation x term
  DATA_IVECTOR(iid_term);         // 0-based term of each level in iid
  DATA_STRUCT(spde, spde_t);      // FEM matrices of the omega mesh
  DATA_STRUCT(spde_xi, spde_t);   // and of the xi mesh
  DATA_VECTOR(range_prior);       // (mean, sd) of the normal prior on log range (km)
  DATA_VECTOR(phi_prior);         // (mean, sd) of the normal prior on logit phi
  DATA_VECTOR(sd_scale);          // half-normal scales of sigma_omega, sigma_eta and every iid SD
  DATA_VECTOR(died);              // the counts
  DATA_VECTOR(n_mosquito);
  DATA_VECTOR(rho);               // beta-binomial intra-class correlation

  PARAMETER_VECTOR(w_omega);
  PARAMETER_ARRAY(x);
  PARAMETER_VECTOR(iid);
  PARAMETER(log_sigma_omega);
  PARAMETER(log_range_omega);
  PARAMETER(log_sigma_eta);
  PARAMETER(log_range_eta);
  PARAMETER(logit_phi);
  PARAMETER_VECTOR(log_sigma_iid); // one per iid term

  // range r = sqrt(8) / kappa
  Type log_sqrt8 = log(sqrt(Type(8.0)));
  Type kappa_omega = exp(log_sqrt8 - log_range_omega);
  Type kappa_eta = exp(log_sqrt8 - log_range_eta);
  Type phi = invlogit(logit_phi);
  vector<Type> sigma_iid = exp(log_sigma_iid);
  Type nll = 0.0;

  // omega: Matern field at the nodes
  nll += GMRF(matern_precision(spde, kappa_omega,
                               exp(log_sigma_omega)))(w_omega);
  vector<Type> lambda = m + A_omega * w_omega;

  // xi: accumulated AR(1) of Matern fields. eta holds the time differences of
  // x (with x_t0 = 0). SEPARABLE(f, g) applies f to the last (year) dimension
  // and g to the first (node); TMB's AR1 is the unit-variance stationary
  // process, so eta has stationary marginal SD sigma_eta, from the GMRF
  int n_nodes = x.dim(0);
  int n_years = x.dim(1);
  array<Type> eta(n_nodes, n_years);
  for (int i = 0; i < n_nodes; i++) {
    eta(i, 0) = x(i, 0);
    for (int t = 1; t < n_years; t++) eta(i, t) = x(i, t) - x(i, t - 1);
  }
  nll += SEPARABLE(AR1(phi), GMRF(matern_precision(spde_xi,
                                                   kappa_eta,
                                                   exp(log_sigma_eta))))(eta);
  vector<Type> x_vec = x.vec();
  lambda += A_xi * x_vec;

  // iid effects
  for (int j = 0; j < iid.size(); j++) {
    nll -= dnorm(iid(j), Type(0.0), sigma_iid(iid_term(j)), true);
  }
  for (int i = 0; i < died.size(); i++) {
    for (int k = 0; k < iid_index.cols(); k++) lambda(i) += iid(iid_index(i, k));
  }

  // log p and log (1 - p) without overflow; the binomial coefficient is a
  // constant and left out
  for (int i = 0; i < died.size(); i++) {
    Type log_p = -logspace_add(Type(0.0), -lambda(i));
    Type log_q = -logspace_add(Type(0.0), lambda(i));
    Type c = (Type(1.0) - rho(i)) / rho(i);
    Type a = c * exp(log_p);
    Type b = c * exp(log_q);
    nll -= lgamma(died(i) + a) + lgamma(n_mosquito(i) - died(i) + b) -
      lgamma(n_mosquito(i) + c) - lgamma(a) - lgamma(b) + lgamma(c);
  }

  // priors on the hyperparameters, on the optimiser's scale: normal on log
  // range and on logit phi (no Jacobian: the priors are defined on those
  // scales), and half-normal on every SD
  nll -= dnorm(log_range_omega, range_prior(0), range_prior(1), true);
  nll -= dnorm(log_range_eta, range_prior(0), range_prior(1), true);
  nll -= dnorm(logit_phi, phi_prior(0), phi_prior(1), true);
  nll -= log_half_normal_log_sd(log_sigma_omega, sd_scale(0));
  nll -= log_half_normal_log_sd(log_sigma_eta, sd_scale(1));
  for (int k = 0; k < log_sigma_iid.size(); k++) {
    nll -= log_half_normal_log_sd(log_sigma_iid(k), sd_scale(2));
  }

  return nll;
}
