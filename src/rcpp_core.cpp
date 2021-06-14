// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <mvt.h>
#include <4beta.h>

  //////////////////////
 // SGT-distribution //
//////////////////////

// R::beta and R::lbeta are risky...
// [[Rcpp::export]]
double log_sgt0(double x, double sigma, double skew, double p, double q, const bool mean_cent, const bool var_adj) {
  double m = 0;
  double v = 1;
  double exp_p = exp(p);
  double exp_q = exp(q);
  if(var_adj == true) {
    if((exp_p * exp_q) < 2) return -arma::datum::inf;
    v = pow(exp_q, -(1/exp_p)) * pow(
      (3 * pow(skew, 2) + 1)
      * (R::beta((3/exp_p), (exp_q - (2/exp_p))) / R::beta((1/exp_p), exp_q))
      - 4 * pow(skew, 2)
      * pow((R::beta((2/exp_p), (exp_q - (1/exp_p))) / R::beta((1/exp_p), exp_q)), 2),
      -(0.5));
  }
  if(mean_cent == true) {
    if((exp_p * exp_q) < 1) return -arma::datum::inf;
    m = (2 * v * sigma * skew * pow(exp_q, (1/exp_p)) * R::beta((2/exp_p), (exp_q-(1/exp_p)))) / R::beta((1/exp_p), exp_q);
  }
  double sgn;
  if((x+m) < 0) sgn = -1.0;
  if((x+m) >= 0) sgn = 1.0;
  double ret = log(exp_p)
    - log(2)
    - log(v)
    - log(sigma)
    - q/exp_p
    - R::lbeta(1/exp_p, exp_q)
    - (1/exp_p + exp_q) * log(
        1 + pow(sgn*(x+m), exp_p) / (exp_q * pow((v * sigma), exp_p) * pow((skew*sgn + 1), exp_p))
    );
  return ret;
}
  //////////////////////////////////////
 // Multivariate normal distribution //
//////////////////////////////////////

// The slower but more general multinormal implementation below is a modification of: https://gallery.rcpp.org/articles/dmvnorm_arma/
static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
double dmvnrm_arma(arma::vec const &x,
                   arma::vec const &mean,
                   arma::mat const &sigma,
                   bool const logd = true) {

  arma::uword const xdim = x.size();
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  z = (x - mean).t();
  inplace_tri_mat_mult(z, rooti);
  double out = other_terms - 0.5 * arma::dot(z, z);

  if(logd) {
    return out;
  }
  return exp(out);
}

// [[Rcpp::export]]
double dmvnrm_arma_diagonal(arma::vec const &x,
                            arma::vec const &mean,
                            arma::mat const &sds,
                            bool const logd = true) {

  arma::vec log_dens(x.size());
  for(int i = 0; i != x.size(); i++) {
    log_dens(i) = arma::log_normpdf(x(i), mean(i), sds(i));
  }

  double out = accu(log_dens);
  if(logd) {
    return out;
  }
  return exp(out);
}

  /////////////////////////////////////
 // Likelihood function and friends //
/////////////////////////////////////

struct LikelihoodParallel : public RcppParallel::Worker {

  //Output
  arma::mat& log_likes;

  //Inputs
  arma::mat& SGT;
  const arma::mat& E;
  const int t;
  const bool mean_cent;
  const bool var_adj;

  LikelihoodParallel(arma::mat& log_likes, arma::mat& SGT, const arma::mat& E, const int t, const bool mean_cent, const bool var_adj)
    : log_likes(log_likes), SGT(SGT), E(E), t(t), mean_cent(mean_cent), var_adj(var_adj) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(int ij = begin; ij != end; ++ij) {
      arma::vec sgt_param = vectorise(SGT.row(ij / t));
      log_likes(ij % t, ij / t) = log_sgt0(E(ij), 1.0, sgt_param(0), sgt_param(1), sgt_param(2), mean_cent, var_adj);
    }
  }

};

// [[Rcpp::export]]
arma::mat fill_xx(arma::mat xx, arma::mat yy, const int m, const int t) {

  bool constant = true;
  if(xx.n_cols % m == 0) constant = false;
  int lags = xx.n_cols / m;

  bool carry_on = true;
  int count = 0;
  while(carry_on) {
    count = count + 1;
    int row = t - count;
    for(int i = 0; i != lags; ++i) {
      if(row + i + 1 == t) break;
      for(int j = 0; j != m; ++j) {
        xx.row(row + i + 1)(constant + m*i + j) = yy.row(row)(j);
      }
    }
    carry_on = xx.row(row).has_nan();
  }
  return(xx);
}

// NOT USED!
// [[Rcpp::export]]
Rcpp::List garch_out(arma::mat yy, arma::mat fit, arma::mat B, arma::mat GARCH,
                     int t, int m) {
  arma::mat V_diags(t, m);
  arma::mat E(t, m, arma::fill::zeros);
  arma::mat VE = (yy - fit) * B.t();
  E.row(0) = VE.row(0);
  arma::vec v_last(m, arma::fill::zeros);
  V_diags.row(0) = arma::exp(v_last).t();
  arma::vec r_exp = arma::exp(GARCH.col(2));
  for(int i = 1; i != VE.n_rows; i++) {
    arma::vec last_abse = arma::abs(E.row(i-1)).t();
    for(int j = 0; j != m; j++) {
      v_last(j) = GARCH(j,0) * v_last(j) + GARCH(j,1) * pow(last_abse(j), r_exp(j));
    }
    V_diags.row(i) = arma::exp(v_last).t();
    E.row(i) = VE.row(i) / arma::exp(v_last.t());
  }
  Rcpp::List out;
  out["E"] = E;
  out["V_diags"] = V_diags;
  return out;
}

// [[Rcpp::export]]
double log_like(arma::vec state, arma::mat yy, arma::mat xx,
                const int first_b, const int first_sgt, const int first_garch, const int first_regime, const int first_yna,
                const int m, const int A_rows, const int t, const arma::ivec regimes, const arma::ivec yna_indices,
                const bool B_inverse, const bool mean_cent, const bool var_adj, const bool parallel_likelihood) {

  //Construct A matrix
  arma::mat A(A_rows, m);
  if(A_rows > 0) {
    for(int i = 0; i != first_b; ++i) A(i) = state(i);
  }

  //Construct B matrix
  arma::mat B(m, m);

  //Modify this if first_garch != -1, i.e. unit diagonal in B (2/4)
  //SVAR or VAR (recursive)
  if((first_sgt - first_b) == (m * m)) {
    for(int i = first_b; i != first_sgt; ++i) {
      B(i - first_b) = state(i);
    }
  } else {
    int row = 0;
    int col = 1;
    int index = 0;
    for(int i = first_b; i != first_sgt; ++i) {
      row = row + 1;
      if(row > m) {
        row = 1;
        col = col + 1;
      }
      while(col > row) {
        B(index) = 0;
        index = index + 1;
        row = row + 1;
        if(row > m) {
          row = 1;
          col = col + 1;
        }
      }
      B(index) = state(i);
      index = index + 1;
    }
  }
  if(!B_inverse) B = arma::inv(B);

  //Construct SGT matrix
  arma::mat SGT(m, 3);
  for(int i = first_sgt; i != first_sgt + m * 3; ++i) SGT(i - first_sgt) = state(i);

  //Construct GARCH matrix
  arma::mat GARCH(m, 2);
  if(first_garch != first_regime) {
    for(int i = first_garch; i != first_garch + m * 2; ++i) GARCH(i - first_garch) = state(i);
  }

  //Construct REGIME matrix ( < 0 not allowed)
  arma::mat REGIME(m, regimes.size() + 1, arma::fill::ones);
  if(first_regime != first_yna) {
    for(int i = first_regime; i != first_yna; ++i) {
      REGIME(i - first_regime) = state(i);
      if(REGIME(i - first_regime) <= 0) return -arma::datum::inf;
    }
    for(int i = 0; i != m; i++) {
      double last_regime = REGIME.n_cols - arma::sum(REGIME.row(i).subvec(0, REGIME.n_cols - 2));
      if(last_regime <= 0) return -arma::datum::inf;
      REGIME.col(REGIME.n_cols - 1)(i) = last_regime;
    }
  }

  //Rebuild yy and xx if there are missing values in the data
  if(yna_indices(0) != -1) {
    for(int i = 0; i != yna_indices.size(); ++i) {
      yy(yna_indices(i) - 1) = state(first_yna + i);
    }
    xx = fill_xx(xx, yy, m, t);
  }

  //Recover structural shocks
  arma::mat fit(t, m, arma::fill::zeros);
  if(A_rows > 0) {
    fit = xx * A;
  }
  arma::mat E(t, m, arma::fill::zeros);

  // (garch == true OR regimes == true)
  arma::mat RVE = (yy - fit) * B.t();
  arma::mat RV_diags(t, m, arma::fill::ones);
  if(first_garch != first_regime || first_regime != first_yna) {

    E = RVE;
    arma::vec v_last(m, arma::fill::ones);
    arma::vec r_last(m, arma::fill::ones);
    RV_diags.row(0) = v_last.t() % r_last.t(); //not necessary

    for(int i = 1; i != RVE.n_rows; i++) {

      // (garch == true)
      if(first_garch != first_regime) {

        // Fix the unconditional expected value
        if(var_adj == true) {
          arma::vec last_dev = E.row(i-1).t() % E.row(i-1).t(); //^2
          for(int j = 0; j != m; j++) {
            v_last(j) = (1 - GARCH(j,0) - GARCH(j,1)) + GARCH(j,0) * v_last(j) + GARCH(j,1) * last_dev(j);
          }
          RV_diags.row(i) = arma::sqrt(v_last).t(); //sqrt
          E.row(i) = RVE.row(i) / arma::sqrt(v_last).t(); //sqrt
        }

        // Do NOT fix the unconditional expected value
        if(var_adj == false) {
          arma::vec last_dev = E.row(i-1).t() % E.row(i-1).t(); //^2
          for(int j = 0; j != m; j++) {
            v_last(j) = 1 + GARCH(j,0) * v_last(j) + GARCH(j,1) * last_dev(j);
          }
          RV_diags.row(i) = arma::sqrt(v_last).t(); //sqrt
          E.row(i) = RVE.row(i) / arma::sqrt(v_last).t(); //sqrt
        }
      }

      // (regimes == true)
      if(first_regime != first_yna) {
        int REGIME_col = 0;
        for(int j = 1; j != REGIME.n_cols; j++) {
          if(i < regimes(j-1)) break;
          REGIME_col = j;
        }
        arma::vec R_diag = REGIME.col(REGIME_col);
        RV_diags.row(i) = RV_diags.row(i) % R_diag.t();
        E.row(i) = E.row(i) / R_diag.t(); //needs to be E not RVE if garch == true
      }
    }
  } else {
    E = RVE;
  }

  //Compute the log-likelihood
  arma::mat log_likes(t, m);

  //Sequential
  if(parallel_likelihood == false) {
    for(int ij = 0; ij != log_likes.size(); ++ij) {
      arma::vec sgt_param = vectorise(SGT.row(ij / t));
      log_likes(ij % t, ij / t) = log_sgt0(E(ij), 1, sgt_param(0), sgt_param(1), sgt_param(2), mean_cent, var_adj);
    }
  }

  //Parallel
  if(parallel_likelihood == true) {
    LikelihoodParallel Wrkr(log_likes, SGT, E, t, mean_cent, var_adj);
    parallelFor(0, log_likes.size(), Wrkr);
  }

  //The log-likelihood to be returned
  double ret;
  double log_abs_det;

  // (garch == true OR regimes == true)
  if(first_garch != first_regime || first_regime != first_yna) {
    arma::vec BRV_inv_dets(t);
    for(int i = 0; i != t; i++) {
      arma::mat RV = arma::diagmat(RV_diags.row(i).t());

      // With extremely bad parameter values the diagonal of RV may underflow
      // (e.g. some numerical optimization algorithms might cause such behavior)
      for(int j = 0; j != m; j++) {
        if(RV(j,j) == 0) return -arma::datum::inf;
      }
      arma::mat BRV = arma::inv( arma::diagmat(RV) ) * B;

      // Log-determinant
      double det_val;
      double det_sign;
      arma::log_det(det_val, det_sign, BRV);
      BRV_inv_dets(i) = det_val;
    }
    log_abs_det = accu(BRV_inv_dets);

  } else {
    double det_val;
    double det_sign;
    arma::log_det(det_val, det_sign, B);
    log_abs_det = t * det_val;
  }

  ret = log_abs_det + accu(log_likes);
  if(std::isnan(ret)) return -arma::datum::inf;
  return ret;
}

  ///////////
 // Prior //
///////////

// [[Rcpp::export]]
bool check_permutation(arma::mat B) {
  arma::vec euclidean_norms(B.n_cols);
  for(int j = 0; j != B.n_cols; j++) euclidean_norms(j) = arma::sum(B.col(j) % B.col(j));
  euclidean_norms = arma::sqrt(euclidean_norms);
  for(int j = 0; j != B.n_cols; j++) B.col(j) = B.col(j) / euclidean_norms(j);
  for(int i = 0; i != B.n_rows; i++) {
    for(int j = 0; j != B.n_cols; j++) {
      if(i == j) {
        if(B(i,j) <= 0) return false;
      }
      if(i < j) {
        if(B(i,i) <= B(i,j)) return false;
      }
    }
  }
  return true;
}

// [[Rcpp::export]]
double log_dirichlet(arma::vec x, const double alpha = 2) {
  if(arma::sum(x) != 1) return -arma::datum::inf;
  for(int i = 0; i < x.size(); i++) if(x(i) < 0) return -arma::datum::inf;
  return arma::sum( (alpha - 1) * arma::log(x) );
}

// [[Rcpp::export]]
double log_prior(arma::vec state, const arma::mat yy, const arma::mat xx,
                 const int first_b, const int first_sgt, const int first_garch,  const int first_regime, const int first_yna,
                 const int m, const int A_rows, const int t, const arma::ivec regimes,
                 const arma::vec a_mean, const arma::mat a_cov, const bool prior_A_diagonal,
                 const arma::vec b_mean, const arma::mat b_cov,
                 const double p_prior_mode, const double p_prior_scale,
                 const double q_prior_mode, const double q_prior_scale,
                 const double dirichlet_alpha,
                 const bool B_inverse) {

  //Prior on A (e.g. Minnesota prior)
  double log_a_prior = 0;
  if(a_cov(0) != -1) {

    //Construct A matrix
    arma::mat A(A_rows, m);
    if(A_rows > 0) {
      for(int i = 0; i != first_b; ++i) A(i) = state(i);
    }

    //Evaluate log prior density at A
    if(prior_A_diagonal == true) {
      log_a_prior = dmvnrm_arma_diagonal(arma::vectorise(A), a_mean, a_cov, prior_A_diagonal);
    } else {
      log_a_prior = dmvnrm_arma(arma::vectorise(A), a_mean, a_cov, prior_A_diagonal);
    }
  }

  // Prior on B
  double log_b_prior = 0;

  // Construct B matrix
  arma::mat B(m, m);

  //SVAR or VAR (recursive)
  if((first_sgt - first_b) == (m * m)) {
    for(int i = first_b; i != first_sgt; ++i) {
      B(i - first_b) = state(i);
    }
  } else {
    int row = 0;
    int col = 1;
    int index = 0;
    for(int i = first_b; i != first_sgt; ++i) {
      row = row + 1;
      if(row > m) {
        row = 1;
        col = col + 1;
      }
      while(col > row) {
        B(index) = 0;
        index = index + 1;
        row = row + 1;
        if(row > m) {
          row = 1;
          col = col + 1;
        }
      }
      B(index) = state(i);
      index = index + 1;
    }
  }

  // Permutation check if 'type == "svar"'
  if((first_sgt - first_b) == (m * m)) {
    if(B_inverse) {
      if(!check_permutation(arma::inv(B))) {
        return -arma::datum::inf;
      }
    } else {
      if(!check_permutation(B)) {
        return -arma::datum::inf;
      }
    }
  }

  // Evaluate log prior density at B
  if(b_cov(0) != -1) {
    log_b_prior = dmvnrm_arma(arma::vectorise(B), b_mean, b_cov);
  }

  //Construct SGT matrix
  arma:: mat SGT(m, 3);
  for(int i = first_sgt; i != first_sgt + m * 3; ++i) SGT(i - first_sgt) = state(i);

  //Construct GARCH matrix and check feasibility of alpha and beta
  arma:: mat GARCH(m, 2);
  if(first_garch != first_regime) {
    for(int i = first_garch; i != first_garch + m * 2; ++i) GARCH(i - first_garch) = state(i);
    for(int i = 0; i != m; ++i) {
      if(GARCH.row(i)(0) < 0) return -arma::datum::inf;
      if(GARCH.row(i)(0) > 1) return -arma::datum::inf;
      if(GARCH.row(i)(1) < 0) return -arma::datum::inf;
      if(GARCH.row(i)(1) > 1) return -arma::datum::inf;
    }
  }

  //Bounds on skewness parameter
  for(int i = 0; i != m; ++i) {
    if(SGT.row(i)(0) > 0.99) return -arma::datum::inf;
    if(SGT.row(i)(0) < -0.99) return -arma::datum::inf;
  }

  //Compute the log-prior on p and q parameters
  arma::vec log_pq_prior(m);
  for(int i = 0; i != m; ++i) {
    arma::vec sgt_param = vectorise(SGT.row(i));
    arma::vec sgt_param_exp = arma::exp(sgt_param);
    if(sgt_param_exp(1) < 0.01) return -arma::datum::inf;
    if(sgt_param_exp(1) > 3.99) return -arma::datum::inf;
    log_pq_prior(i) =
      arma::log_normpdf(sgt_param_exp(1), p_prior_mode, p_prior_scale) +
      arma::log_normpdf(sgt_param(2), q_prior_mode, q_prior_scale);
  }

  //Construct REGIME matrix ( < 0 not allowed)
  arma::mat REGIME(m, regimes.size() + 1, arma::fill::ones);
  if(first_regime != first_yna) {
    for(int i = first_regime; i != first_yna; ++i) {
      REGIME(i - first_regime) = state(i);
      if(REGIME(i - first_regime) <= 0) return -arma::datum::inf;
    }
    for(int i = 0; i != m; i++) {
      double last_regime = REGIME.n_cols - arma::sum(REGIME.row(i).subvec(0, REGIME.n_cols - 2));
      if(last_regime <= 0) return -arma::datum::inf;
      REGIME.col(REGIME.n_cols - 1)(i) = last_regime;
    }
  }

  //Dirichlet prior on regime parameters
  double log_diri_prior = 0;
  if(first_regime != first_yna) {
    for(int i = 0; i != REGIME.n_rows; i++) {
      arma::vec dirichlet_vec = (REGIME.row(0).t() / REGIME.n_cols);
      log_diri_prior = log_diri_prior + log_dirichlet(dirichlet_vec, dirichlet_alpha);
    }
  }

  //Return the log-prior
  double ret = accu(log_pq_prior) + log_a_prior + log_b_prior + log_diri_prior;
  if(std::isnan(ret)) ret = -arma::datum::inf;
  return ret;
}

  ///////////////
 // Posterior //
///////////////

struct Model {

  // Likelihood arguments
  arma::mat yy;
  arma::mat xx;
  int first_b;
  int first_sgt;
  int first_garch;
  int first_regime;
  int first_yna;
  int m;
  int A_rows;
  int t;
  arma::ivec regimes;
  arma::ivec yna_indices;
  bool mean_cent;
  bool var_adj;
  bool B_inverse;

  // Prior arguments
  arma::vec a_mean;
  arma::mat a_cov;
  bool prior_A_diagonal;
  arma::vec b_mean;
  arma::mat b_cov;
  arma::vec p_prior;
  arma::vec q_prior;
  double dirichlet_alpha;

  void init(Rcpp::List model_R) {
    arma::mat yy_ = model_R["yy"]; yy = yy_;
    arma::mat xx_ = model_R["xx"]; xx = xx_;
    int first_b_ = model_R["first_b"]; first_b = first_b_;
    int first_sgt_ = model_R["first_sgt"]; first_sgt = first_sgt_;
    int first_garch_ = model_R["first_garch"]; first_garch = first_garch_;
    int first_regime_ = model_R["first_regime"]; first_regime =first_regime_;
    int first_yna_ = model_R["first_yna"]; first_yna = first_yna_;
    int m_ = model_R["m"]; m = m_;
    int A_rows_ = model_R["A_rows"]; A_rows = A_rows_;
    int t_ = model_R["t"]; t = t_;
    arma::ivec regimes_ = model_R["regimes"]; regimes = regimes_;
    arma::ivec yna_indices_ = model_R["yna_indices"]; yna_indices = yna_indices_;
    bool mean_cent_ = model_R["mean_cent"]; mean_cent = mean_cent_;
    bool var_adj_ = model_R["var_adj"]; var_adj = var_adj_;
    bool B_inverse_ = model_R["B_inverse"]; B_inverse = B_inverse_;

    arma::vec a_mean_ = model_R["a_mean"]; a_mean = a_mean_;
    arma::mat a_cov_ = model_R["a_cov"]; a_cov = a_cov_;
    bool prior_A_diagonal_ = model_R["prior_A_diagonal"]; prior_A_diagonal = prior_A_diagonal_;
    arma::vec b_mean_ = model_R["b_mean"]; b_mean = b_mean_;
    arma::mat b_cov_ = model_R["b_cov"]; b_cov = b_cov_;
    arma::vec p_prior_ = model_R["p_prior"]; p_prior = p_prior_;
    arma::vec q_prior_ = model_R["q_prior"]; q_prior = q_prior_;
    double dirichlet_alpha_ = model_R["dirichlet_alpha"]; dirichlet_alpha = dirichlet_alpha_;
  }
};

double log_posterior(const arma::vec par, const Model M, const bool parallel_likelihood) {

  double ll = log_like(par, M.yy, M.xx,
                       M.first_b, M.first_sgt, M.first_garch, M.first_regime, M.first_yna,
                       M.m, M.A_rows, M.t,
                       M.regimes, M.yna_indices, M.B_inverse,
                       M.mean_cent, M.var_adj, parallel_likelihood);
  double lp = log_prior(par, M.yy, M.xx,
                        M.first_b, M.first_sgt, M.first_garch, M.first_regime, M.first_yna,
                        M.m, M.A_rows, M.t, M.regimes,
                        M.a_mean, M.a_cov, M.prior_A_diagonal, M.b_mean, M.b_cov,
                        M.p_prior(0), M.p_prior(1),
                        M.q_prior(0), M.q_prior(1),
                        M.dirichlet_alpha,
                        M.B_inverse);
  return ll + lp;
}

  ///////////////////////////////
 // DE-MC sampler and friends //
///////////////////////////////

void draw(arma::mat& draws, arma::vec& densities, arma::vec& asums, int state_row, int last_row,
          const double gamma, const int K, const int n, const bool parallel_likelihood,
          Model M) {

  //Initialize the state with current state of the chain
  arma::vec state = vectorise(draws.row(state_row));

  //Initialize last density
  double last_density = densities(state_row);

  //Acceptance vector to compute acceptance rates
  arma::vec accepted(K);
  accepted.fill(0);

  //Proposal density
  double proposal_density;

  //Inner loop of iterations
  for(int i = 0; i != K; ++i) {

    //Draw proposal
    double rnd1 = arma::randu() * last_row ;
    double rnd2 = arma::randu() * (last_row - 1);
    int rnd1_int = rnd1;
    int rnd2_int = rnd2;
    if(rnd1_int >= rnd2_int) ++rnd2_int; //Makes sure the draw happens without replacement
    arma::vec r1 = vectorise(draws.row(rnd1_int));
    arma::vec r2 = vectorise(draws.row(rnd2_int));
    arma::vec proposal = state + gamma * (r1 - r2);

    //Compute proposal density
    proposal_density = log_posterior(proposal, M, parallel_likelihood);

    //Accept or reject proposal
    double log_ratio = proposal_density - last_density;
    bool accept;
    if(log_ratio > 0) {
      accept = true;
    } else if(log_ratio > std::log(arma::randu())) {
      accept = true;
    } else {
      accept = false;
    }

    //Update state and density
    if(accept == true) {
      state = proposal;
      last_density = proposal_density;
      accepted(i) = 1;
    }
  }

  //Save the state, density and acceptance rate
  draws.row(state_row + n) = state.t();
  densities(state_row + n) = last_density;
  asums(state_row + n) = arma::sum(accepted);

  return void();
}

struct DrawParallel : public RcppParallel::Worker {

  arma::mat& draws;
  arma::vec& densities;
  arma::vec& asums;
  const double gamma;
  const int K;
  const int n;
  const bool parallel_likelihood;
  const Model M;

  DrawParallel(arma::mat& draws, arma::vec& densities, arma::vec& asums,
               const double gamma, const int K, const int n, const bool parallel_likelihood, const Model M)
    : draws(draws), densities(densities), asums(asums),
      gamma(gamma), K(K), n(n), parallel_likelihood(parallel_likelihood), M(M) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(int j = begin; j != end; ++j) {
      draw(draws, densities, asums, j, j - n,
           gamma, K, n, parallel_likelihood,
           M);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List sampler(const int N, const int n, const int m0, const int K, const double gamma, const arma::mat init_draws,
                   const bool output_as_input, const arma::mat old_chain, const bool new_chain,
                   const bool parallel, const bool parallel_likelihood,
                   const Rcpp::List model_R, const bool progress_bar) {

  // Initialize the model object
  Model M;
  M.init(model_R);

  // Initialize the chains
  arma::mat draws(N * n + m0, init_draws.n_cols);

  // Initial chains/draws from multivariate normal...
  if(!output_as_input) {
    for(int i = 0; i != m0; ++i) {
      draws.row(i) = init_draws.row(i);
    }

  //...or from an existing sample
  } else {

    // Draw random rows from an existing sample to initialize new chains...
    if(new_chain) {
      for(int i = 0; i != m0; ++i) {
        double rnd_double = arma::randu() * old_chain.n_rows;
        int rnd_int = rnd_double;
        draws.row(i) = old_chain.row(rnd_int);
      }

    //...or copy the entire sample to continue the existing chains
    } else {
      for(int i = 0; i != old_chain.n_rows; ++i) {
        draws.row(i) = old_chain.row(i);
      }
    }
  }

  // Initialize the loop variables
  arma::vec densities(N * n + m0);
  arma::vec asums(N * n + m0);
  int last_row = m0 - 1;
  if(new_chain == false && output_as_input == true) last_row = old_chain.n_rows - 1;

  // Initial densities
  for(int j = 0; j != n; ++j) {
    arma::vec state = vectorise(draws.row(last_row - j));
    double density = log_posterior(state, M, parallel_likelihood);
    densities(last_row - j) = density;
  }

  // Progress bar variables (Should be corrected)
  double progress = 0.0;
  int bar_width = 70;
  int update_every = N/100;

  // Sequential implementation:
  if(parallel == false) {
    for(int i = 0; i != N; ++i) {

      // If an existing sample is continued we skip iterations already carried
      if(i == 0 && output_as_input == true && new_chain == false) {
        i = (old_chain.n_rows - m0)/n;
        progress += double(i)/double(N);
      }

      // Update progress bar
      if(i % update_every == 0) {
        if(progress_bar) {
          progress += 0.01;
          Rcpp::Rcout << "[";
          int pos = bar_width * progress;
          for (int i = 0; i < bar_width; ++i) {
            if (i < pos) Rcpp::Rcout << "=";
            else if (i == pos) Rcpp::Rcout << ">";
            else Rcpp::Rcout << " ";
          }
          Rcpp::Rcout << "] " << int(progress * 100.0) << " %\r";
          Rcpp::Rcout.flush();
        }
      }

      // Update states of the chains
      for(int j = 0; j != n; ++j) {
        int state_row = last_row - j;
        draw(draws, densities, asums, state_row, last_row,
             gamma, K, n, parallel_likelihood, M);
      }
      last_row = last_row + n;
    }
  }

  // Parallel implementation:
  if(parallel == true) {
    DrawParallel Wrkr(draws, densities, asums,
                      gamma, K, n, parallel_likelihood, M);
    for(int i = 0; i != N; ++i) {

      //If an existing sample is continued we skip iterations already carried
      if(i == 0 && output_as_input == true && new_chain == false) {
        i = (old_chain.n_rows - m0)/n;
        progress += double(i)/double(N);
      }

      // Update progress bar
      if(i % update_every == 0) {
        if(progress_bar) {
          progress += 0.01;
          Rcpp::Rcout << "[";
          int pos = bar_width * progress;
          for (int i = 0; i < bar_width; ++i) {
            if (i < pos) Rcpp::Rcout << "=";
            else if (i == pos) Rcpp::Rcout << ">";
            else Rcpp::Rcout << " ";
          }
          Rcpp::Rcout << "] " << int(progress * 100.0) << " %\r";
          Rcpp::Rcout.flush();
        }
      }

      // Update states of the chains
      parallelFor(last_row - n + 1, last_row + 1, Wrkr);
      last_row = last_row + n;
    }
  }
  if(progress_bar) Rcpp::Rcout << std::endl;

  // Return what needs to be returned
  Rcpp::List ret_list(3);
  ret_list[0] = draws;
  ret_list[1] = densities;
  ret_list[2] = asums;
  return ret_list;
}

  /////////////////////
 // IRF and friends // (Not ready...)
/////////////////////

// [[Rcpp::export]]
arma::mat stackA_cpp(arma::mat A, const bool constant = true) {
  A = A.t();
  if(constant) A = A.cols(1, A.n_cols - 1);
  const int m = A.n_rows;
  const int lags = A.n_cols / m;
  const int eye_dim = m * lags - m;
  arma::mat AA(A.n_rows + eye_dim, A.n_cols, arma::fill::zeros);
  for(int i = 0; i < A.n_rows; i++) {
    AA.row(i) = A.row(i);
  }
  for(int i = A.n_rows; i < AA.n_rows; i++) {
    for(int j = 0; j < AA.n_cols; j++) {
      if((i - A.n_rows) == j) AA(i,j) = 1;
    }
  }
  return(AA);
}

// [[Rcpp::export]]
Rcpp::List irf_cpp(const arma::mat s, const int horizon, const arma::vec cumulate,
                   const arma::vec shock_sizes, const arma::vec shocks,
                   const int A_rows, const int first_b, const int first_sgt, const int m,
                   const bool B_inverse, const bool parallel) {

  Rcpp::List ret(shocks.size());

  // Sequential
  if(parallel == false) {

    // Loop over shocks
    for(int j = 0; j < shocks.size(); j++) {
      arma::vec e(m, arma::fill::zeros);
      e(shocks(j)) = shock_sizes(shocks(j));
      arma::cube irfs(shocks.size(), horizon + 1, s.n_rows);

      // Loop over posterior sample
      for(int i = 0; i < s.n_rows; i++) {

        // Construct A matrix
        arma::mat A(A_rows, m);
        if(A_rows > 0) {
          for(int j = 0; j != first_b; ++j) A(j) = s.row(i)(j);
        }
        arma::mat AA = stackA_cpp(A);

        // Construct B matrix
        arma::mat B(m, m);

        //SVAR or VAR (recursive)
        if((first_sgt - first_b) == (m * m)) {
          for(int j = first_b; j != first_sgt; ++j) {
            B(j - first_b) = s.row(i)(j);
          }
        } else {
          int row = 0;
          int col = 1;
          int index = 0;
          for(int j = first_b; j != first_sgt; ++j) {
            row = row + 1;
            if(row > m) {
              row = 1;
              col = col + 1;
            }
            while(col > row) {
              B(index) = 0;
              index = index + 1;
              row = row + 1;
              if(row > m) {
                row = 1;
                col = col + 1;
              }
            }
            B(index) = s.row(i)(j);
            index = index + 1;
          }
        }
        if(B_inverse) B = arma::inv(B);

        // Loop over horizons
        for(int h = 0; h < (horizon + 1); h++) {

          arma::vec zero = B * e;
          arma::vec zero_long(AA.n_cols, arma::fill::zeros);
          if(h == 0) {
            for(int hh = 0; hh < zero.size(); hh++) zero_long(hh) = zero(hh);
            irfs.slice(i).col(h) = zero;
          } else {
            arma::vec response = arma::powmat(AA, (h-1)) * zero_long;
            for(int hh = 0; hh < zero.size(); hh++) irfs.slice(i).col(h)(hh) = response(hh);
          }
        }
      }

      if(cumulate(0) != -1) {
        for(int i = 0; i < s.n_rows; i++) {
          for(int hh = 0; hh < cumulate.size(); hh++) {
            for(int h = 1; h < (horizon + 1); h++) {
              irfs.slice(i).row(cumulate(hh))(h) = arma::sum(irfs.slice(i).row(cumulate(hh)).subvec(0, h));
            }
          }
        }
      }
      ret(j) = irfs;
    }
  }

  // Parallel
  if(parallel == true) {
    //...
  }

  return(ret);
}

  /////////////////////////
 // Marginal likelihood //
/////////////////////////

double log_alpha(const double numerator_density, const double denumerator_density) {
  double ret = numerator_density - denumerator_density;
  if(ret > 0) {
    return 0;
  } else {
    return ret;
  }
}

double dproposal(const arma::vec theta, const arma::vec theta_star, const arma::mat sigma_star) {
  arma::mat x = theta.t();
  arma::vec retvec = dmvt(x, theta_star, sigma_star, 1, true);
  double ret = retvec(0);
  return ret;
}

struct MlParallel : public RcppParallel::Worker {

  // Output
  arma::vec& denumerator_log_vec;

  // Inputs
  arma::mat& proposal_mat;
  const double logden_star;
  const bool parallel_likelihood;
  const Model M;

  // Constructor
  MlParallel(arma::vec& denumerator_log_vec, arma::mat& proposal_mat, const double logden_star,
             const bool parallel_likelihood, const Model M)
    : denumerator_log_vec(denumerator_log_vec), proposal_mat(proposal_mat), logden_star(logden_star),
      parallel_likelihood(parallel_likelihood), M(M) {}

  //  Instructions
  void operator()(std::size_t begin, std::size_t end) {
    for(int j = begin; j != end; ++j) {
      arma::vec proposal = proposal_mat.row(j).t();
      double proposal_density = log_posterior(proposal, M, parallel_likelihood);
      denumerator_log_vec(j) = log_alpha(proposal_density, logden_star);
    }
  }
};

// For Debugging
// [[Rcpp::export]]
void log_text(std::string text, int iter) {
  std::ofstream myfile;
  myfile.open("log.txt");
  myfile << text << " / iter = " << iter << " \n";
  myfile.close();
  return void();
}

// [[Rcpp::export]]
double logSumExp(arma::vec x) {
  double x_max = max(x);
  arma::vec x_tilde = x - x_max;
  return(x_max + log(sum(exp(x_tilde))));
}

// [[Rcpp::export]]
Rcpp::List log_ml_cpp(const arma::vec proposal_densities, const arma::vec posterior_densities,
                      const arma::vec theta_star, const arma::mat sigma_star, const double logden_star,
                      const arma::uword J, const bool parallel, const bool parallel_likelihood,
                      const Rcpp::List model_R) {

  // Initialize the model object
  Model M;
  M.init(model_R);

  log_text("log_ml_cpp starts", 0);

  arma::vec numerator_log_vec(proposal_densities.size());
  for(int i = 0; i < numerator_log_vec.size(); i++) {
    double theta_density = posterior_densities(i);
    double proposal_density = proposal_densities(i);
    numerator_log_vec(i) = log_alpha(logden_star, theta_density) + proposal_density;
    log_text("First loop, ", i);
  }

  log_text("First loop has ended", 0);

  arma::vec denumerator_log_vec(J);
  arma::mat proposal_mat = rmvt(J, theta_star, sigma_star, 1);

  // Will be used in the future with zero restrictions
  //if(false) {
  //  for(int i = 0; i != no_dev.size(); i++) {
  //    arma::vec theta_star_col(proposal_mat.n_rows);
  //    theta_star_col.fill(theta_star(no_dev(i)));
  //    proposal_mat.col(no_dev(i)) = theta_star_col;
  //  }
  //}

  log_text("Main loop starts", 0);

  // Sequential
  if(parallel == false) {
    for(int j = 0; j < J; j++) {

      log_text("Main loop - Drawing proposal...", j);
      arma::vec proposal = proposal_mat.row(j).t();

      log_text("Main loop - Computing posterior density...", j);
      double proposal_density = log_posterior(proposal, M, parallel_likelihood);

      log_text("Main loop - Finishing...", j);
      denumerator_log_vec(j) = log_alpha(proposal_density, logden_star);
    }
  }

  log_text("Main loop has ended", 0);

  // Parallel
  if(parallel== true) {
    MlParallel Wrkr(denumerator_log_vec, proposal_mat, logden_star,
                    parallel_likelihood, M);
    parallelFor(0, denumerator_log_vec.size(), Wrkr);
  }

  double log_numerator = -log(numerator_log_vec.size()) + logSumExp(numerator_log_vec);
  double log_denumerator = -log(denumerator_log_vec.size()) + logSumExp(denumerator_log_vec);
  double log_posterior_ordinate = log_numerator - log_denumerator;
  double log_marginal_likelihood = logden_star - log_posterior_ordinate;
  Rcpp::List ret_list(4);
  ret_list[0] = log_marginal_likelihood;
  ret_list[1] = log_posterior_ordinate;
  ret_list[2] = numerator_log_vec;
  ret_list[3] = denumerator_log_vec;
  return ret_list;
}

