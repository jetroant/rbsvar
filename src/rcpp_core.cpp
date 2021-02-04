// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>

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
        1 + pow(abs(x+m), exp_p) / (exp_q * pow((v * sigma), exp_p) * pow((skew*sgn + 1), exp_p))
    );
  return ret;
}

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
// Multinormal implementation ends.

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
                const int first_b, const int first_sgt, const int first_garch, const int first_yna,
                const int m, const int A_rows, const int t, const arma::vec yna_indices,
                const bool mean_cent, const bool var_adj, const bool parallel_likelihood) {

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

  //Construct SGT matrix
  arma:: mat SGT(m, 3);
  for(int i = first_sgt; i != first_sgt + m * 3; ++i) SGT(i - first_sgt) = state(i);

  //Construct GARCH matrix
  arma:: mat GARCH(m, 3);
  if(first_garch != first_yna) {
    for(int i = first_garch; i != first_garch + m * 3; ++i) GARCH(i - first_garch) = state(i);
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

  // (garch == true)
  arma::mat V_diags(t, m);
  if(first_garch != first_yna) {
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

  // (garch == false)
  } else {
    E = (yy - fit) * B.t();
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
  // (garch == true)
  if(first_garch != first_yna) {
    arma::vec BV_inv_dets(t);
    for(int i = 0; i != t; i++) {
      arma::mat V = arma::diagmat(V_diags.row(i).t());

      //With extremely bad parameter values the diagonal of V underflows
      //(e.g. some numerical optimization algorithms might cause such behavior)
      for(int j = 0; j != m; j++) {
        if(V(j,j) == 0) return -arma::datum::inf;
      }
      arma::mat BV = arma::inv( arma::diagmat(V) ) * B;
      BV_inv_dets(i) = log(abs(arma::det(BV)));
    }
    ret = accu(BV_inv_dets) + accu(log_likes);

  // (garch == false)
  } else {
    ret = t * log(abs(arma::det(B))) + accu(log_likes);
  }
  if(std::isnan(ret)) return -arma::datum::inf;
  return ret;
}

// [[Rcpp::export]]
double log_prior(arma::vec state, const arma::mat yy, const arma::mat xx,
                 const int first_b, const int first_sgt, const int first_garch, const int first_yna,
                 const int m, const int A_rows, const int t,
                 const arma::vec a_mean, const arma::mat a_cov, const bool prior_A_diagonal,
                 const arma::vec b_mean, const arma::mat b_cov,
                 const double p_prior_mode, const double p_prior_scale,
                 const double q_prior_mode, const double q_prior_scale) {

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

  //Prior on B
  double log_b_prior = 0;
  if(b_cov(0) != -1) {

    //Construct B matrix
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

    //Evaluate log prior density at B
    log_b_prior = dmvnrm_arma(arma::vectorise(B), b_mean, b_cov);
  }

  //Construct SGT matrix
  arma:: mat SGT(m, 3);
  for(int i = first_sgt; i != first_sgt + m * 3; ++i) SGT(i - first_sgt) = state(i);

  //Construct GARCH matrix and check feasibility of gammas and r
  arma:: mat GARCH(m, 3);
  if(first_garch != first_yna) {
    for(int i = first_garch; i != first_garch + m * 3; ++i) GARCH(i - first_garch) = state(i);
    for(int i = 0; i != m; ++i) {
      if(GARCH.row(i)(0) < 0) return -arma::datum::inf;
      if(GARCH.row(i)(1) < 0) return -arma::datum::inf;
      if(GARCH.row(i)(2) > SGT.row(i)(1) + SGT.row(i)(2)) return -arma::datum::inf;
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
    log_pq_prior(i) =
      log_sgt0(sgt_param(1) - p_prior_mode, p_prior_scale, -0.4, 2, 5, false, false) +
      arma::log_normpdf(sgt_param(2), q_prior_mode, q_prior_scale);
  }

  //Return the log-prior
  double ret = accu(log_pq_prior) + log_a_prior + log_b_prior;
  if(std::isnan(ret)) ret = -arma::datum::inf;
  return ret;
}

// [[Rcpp::export]]
void draw(arma::mat& draws, arma::vec& densities, arma::vec& asums, int state_row, int last_row,
          const double gamma, const int K, const int n,
          const arma::mat& yy, const arma::mat& xx,
          const int first_b, const int first_sgt, const int first_garch, const int first_yna,
          const int m, const int A_rows, const int t, const arma::vec yna_indices,
          const arma::vec a_mean, const arma::mat a_cov, const bool prior_A_diagonal,
          const arma::vec b_mean, const arma::mat b_cov,
          const double p_prior_mode, const double p_prior_scale,
          const double q_prior_mode, const double q_prior_scale,
          const bool mean_cent, const bool var_adj, const bool parallel_likelihood) {

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
    proposal_density =
      log_like(proposal, yy, xx,
               first_b, first_sgt, first_garch, first_yna,
               m, A_rows, t, yna_indices,
               mean_cent, var_adj, parallel_likelihood) +
      log_prior(proposal, yy, xx,
                first_b, first_sgt, first_garch, first_yna,
                m, A_rows, t,
                a_mean, a_cov, prior_A_diagonal, b_mean, b_cov,
                p_prior_mode, p_prior_scale, q_prior_mode, q_prior_scale);

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
  const arma::mat& yy;
  const arma::mat& xx;
  const int first_b;
  const int first_sgt;
  const int first_garch;
  const int first_yna;
  const bool mean_cent;
  const bool var_adj;
  const arma::vec a_mean;
  const arma::mat a_cov;
  const bool prior_A_diagonal;
  const arma::vec b_mean;
  const arma::mat b_cov;
  const arma::vec yna_indices;
  const RcppParallel::RVector<double> par_vec;

  DrawParallel(arma::mat& draws, arma::vec& densities, arma::vec& asums, const arma::mat& yy, const arma::mat& xx,
               const int first_b, const int first_sgt, const int first_garch, const int first_yna,
               const arma::vec a_mean, const arma::mat a_cov, const bool prior_A_diagonal, const arma::vec b_mean, const arma::mat b_cov,
               Rcpp::NumericVector par_vec, const arma::vec yna_indices, const bool mean_cent, const bool var_adj)
    : draws(draws), densities(densities), asums(asums), yy(yy), xx(xx),
      first_b(first_b), first_sgt(first_sgt), first_garch(first_garch), first_yna(first_yna),
      a_mean(a_mean), a_cov(a_cov), prior_A_diagonal(prior_A_diagonal), b_mean(b_mean), b_cov(b_cov),
      par_vec(par_vec), yna_indices(yna_indices), mean_cent(mean_cent), var_adj(var_adj) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(int j = begin; j != end; ++j) {
      draw(draws, densities, asums, j, j - par_vec[1],
           par_vec[4], par_vec[3], par_vec[1], yy, xx,
           first_b, first_sgt, first_garch, first_yna,
           par_vec[5], par_vec[6], par_vec[7], yna_indices,
           a_mean, a_cov, prior_A_diagonal, b_mean, b_cov,
           par_vec[8], par_vec[9], par_vec[10], par_vec[11], mean_cent, var_adj, par_vec[12]);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List sampler(const int N, const int n, const int m0, const int K, const double gamma,
                   const arma::vec init_mode, const arma::mat init_scale,
                   const bool output_as_input, const arma::mat old_chain, const bool new_chain,
                   const bool parallel, const bool parallel_likelihood,
                   const arma::mat yy, const arma::mat xx,
                   const int m, const int A_rows, const int t, const arma::vec yna_indices,
                   const bool mean_cent, const bool var_adj,
                   const int first_b, const int first_sgt, const int first_garch, const int first_yna,
                   const arma::vec a_mean, const arma::mat a_cov, const bool prior_A_diagonal,
                   const arma::vec b_mean, const arma::mat b_cov,
                   const double p_prior_mode, const double p_prior_scale,
                   const double q_prior_mode, const double q_prior_scale,
                   const bool progress_bar) {

  //Collect parameters
  Rcpp::NumericVector par_vec(12);
  par_vec[0] = N;
  par_vec[1] = n;
  par_vec[2] = m0;
  par_vec[3] = K;
  par_vec[4] = gamma;
  par_vec[5] = m;
  par_vec[6] = A_rows;
  par_vec[7] = t;
  par_vec[8] = p_prior_mode;
  par_vec[9] = p_prior_scale;
  par_vec[10] = q_prior_mode;
  par_vec[11] = q_prior_scale;
  par_vec[12] = parallel_likelihood;

  //Initialize the chains
  arma::mat draws(N * n + m0, init_mode.n_elem);

  //Initial chains/draws from multivariate normal...
  if(!output_as_input) {
    for(int i = 0; i != m0; ++i) {
      arma::vec xn = arma::randn(init_mode.n_elem);
      if(init_scale.size() > 1) {
        draws.row(i) = init_mode.t() + xn.t() * arma::chol(init_scale);
      } else {
        arma::mat init_scale_root = init_mode.t();
        for(int j = 0; j != init_mode.size(); ++j) {
          init_scale_root(j) = sqrt(init_scale(0));
        }
        draws.row(i) = init_mode.t() + xn.t() % init_scale_root;
      }
    }

  //...or from an existing sample
  } else {

    //Draw random rows from an existing sample to initialize new chains...
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

  //Initialize the loop variables
  arma::vec densities(N * n + m0);
  arma::vec asums(N * n + m0);
  int last_row = m0 - 1;
  if(new_chain == false && output_as_input == true) last_row = old_chain.n_rows - 1;

  //Initial densities
  for(int j = 0; j != n; ++j) {
    arma::vec state = vectorise(draws.row(last_row - j));
    double density =
      log_like(state, yy, xx,
               first_b, first_sgt, first_garch, first_yna,
               m, A_rows, t, yna_indices,
               mean_cent, var_adj, parallel_likelihood) +
      log_prior(state, yy, xx,
                first_b, first_sgt, first_garch, first_yna,
                m, A_rows, t,
                a_mean, a_cov, prior_A_diagonal, b_mean, b_cov,
                p_prior_mode, p_prior_scale, q_prior_mode, q_prior_scale);
    densities(last_row - j) = density;
  }

  //Progress bar variables
  double progress = 0.0;
  int bar_width = 70;
  int update_every = N/100;

  //Sequential
  if(parallel == false) {
    for(int i = 0; i != N; ++i) {

      //If an existing sample is continued we skip iterations already carried
      if(i == 0 && output_as_input == true && new_chain == false) {
        i = (old_chain.n_rows - m0)/n;
        progress += double(i)/double(N);
      }

      //Update progress bar
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

      //Update states of the chains
      for(int j = 0; j != n; ++j) {
        int state_row = last_row - j;
        draw(draws, densities, asums, state_row, last_row,
             gamma, K, n, yy, xx,
             first_b, first_sgt, first_garch, first_yna,
             m, A_rows, t, yna_indices,
             a_mean, a_cov, prior_A_diagonal, b_mean, b_cov,
             p_prior_mode, p_prior_scale, q_prior_mode, q_prior_scale,
             mean_cent, var_adj, parallel_likelihood);
      }
      last_row = last_row + n;
    }
  }

  //Parallel
  if(parallel == true) {
    DrawParallel Wrkr(draws, densities, asums,
                      yy, xx, first_b, first_sgt, first_garch, first_yna,
                      a_mean, a_cov, prior_A_diagonal, b_mean, b_cov,
                      par_vec, yna_indices, mean_cent, var_adj);
    for(int i = 0; i != N; ++i) {

      //If an existing sample is continued we skip iterations already carried
      if(i == 0 && output_as_input == true && new_chain == false) {
        i = (old_chain.n_rows - m0)/n;
        progress += double(i)/double(N);
      }

      //Update progress bar
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

      //Update states of the chains
      parallelFor(last_row - n + 1, last_row + 1, Wrkr);
      last_row = last_row + n;
    }
  }
  if(progress_bar) Rcpp::Rcout << std::endl;

  //Return what is needed to be returned
  Rcpp::List ret_list(3);
  ret_list[0] = draws;
  ret_list[1] = densities;
  ret_list[2] = asums;
  return ret_list;
}





