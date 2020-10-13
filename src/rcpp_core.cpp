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
arma::vec sgt_bounds(arma::vec state, int first_sgt, int m) {

  //Construct SGT matrix
  arma:: mat SGT(m, 3);
  for(int i = first_sgt; i != state.size(); ++i) SGT(i - first_sgt) = state(i);

  //Check bounds
  for(int i = 0; i != m; ++i) {
    if(SGT.row(i)(0) > 0.99) SGT.row(i)(0) = 0.99;
    if(SGT.row(i)(0) < -0.99) SGT.row(i)(0) = -0.99;
    if(SGT.row(i)(1) > 10) SGT.row(i)(1) = 10;
    if(SGT.row(i)(1) < -2.3) SGT.row(i)(1) = -2.3;
    if(SGT.row(i)(2) > 10) SGT.row(i)(2) = 10;
    if(SGT.row(i)(2) < -2.3) SGT.row(i)(2) = -2.3;
  }

  //Return the bounded values
  for(int i = first_sgt; i != state.size(); ++i) state(i) = SGT(i - first_sgt);
  return state;
}

// [[Rcpp::export]]
double log_like(arma::vec state, const arma::mat yy, const arma::mat xx,
                const int first_b, const int first_sgt,
                const int m, const int A_rows, const int t, const bool mean_cent, const bool var_adj,
                const bool parallel_likelihood, bool bounds = true) {

  //Makes sure the likelihood is well defined in some extreme cases
  if(bounds) {
    state = sgt_bounds(state, first_sgt, m);
  }

  //Construct A matrix
  arma::mat A(A_rows, m);
  if(A_rows > 0) {
    for(int i = 0; i != first_b; ++i) A(i) = state(i);
  }

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

  //Construct SGT matrix
  arma:: mat SGT(m, 3);
  for(int i = first_sgt; i != state.size(); ++i) SGT(i - first_sgt) = state(i);

  //Recover structural shocks
  arma::mat fit(t, m, arma::fill::zeros);
  if(A_rows > 0) {
    fit = xx * A;
  }
  const arma::mat E = (yy - fit) * B.t();

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

  //Return the log-likelihood
  double ret = t * log(abs(arma::det(B))) + accu(log_likes);
  if(std::isnan(ret)) ret = -arma::datum::inf;
  return ret;
}

// [[Rcpp::export]]
double log_prior(arma::vec state, const arma::mat yy, const arma::mat xx,
                 const int first_b, const int first_sgt,
                 const int m, const int A_rows, const int t,
                 const double p_prior_mode, const double p_prior_scale,
                 const double q_prior_mode, const double q_prior_scale,
                 const double shrinkage) {

  if(!std::isinf(shrinkage)) {

    //Construct the prior in init, just eval here (t(df = 2)?)

    //Construct A matrix
    //arma::mat A(A_rows, m);
    //if(A_rows > 0) {
    //  for(int i = 0; i != first_b; ++i) A(i) = state(i);
    //}

    //Impose Minnesota prior...
  }

  //Construct SGT matrix
  arma:: mat SGT(m, 3);
  for(int i = first_sgt; i != state.size(); ++i) SGT(i - first_sgt) = state(i);

  //Compute the log-prior on p and q parameters
  arma::vec log_pq_prior(m);
  for(int i = 0; i != m; ++i) {
    arma::vec sgt_param = vectorise(SGT.row(i));
    log_pq_prior(i) =
      log_sgt0(sgt_param(1) - p_prior_mode, p_prior_scale, -0.4, 2, 5, false, false) +
      arma::log_normpdf(sgt_param(2), q_prior_mode, q_prior_scale);
  }

  //Return the log-prior
  double ret = accu(log_pq_prior);
  if(std::isnan(ret)) ret = -arma::datum::inf;
  return ret;
}

// [[Rcpp::export]]
void draw(arma::mat& draws, arma::vec& densities, arma::vec& asums, int state_row, int last_row,
          const double gamma, const int K, const int n,
          const arma::mat& yy, const arma::mat& xx,
          const int first_b, const int first_sgt,
          const int m, const int A_rows, const int t, const bool mean_cent, const bool var_adj,
          const double p_prior_mode, const double p_prior_scale,
          const double q_prior_mode, const double q_prior_scale,
          const double shrinkage, const bool parallel_likelihood) {

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
      log_like(proposal, yy, xx, first_b, first_sgt, m, A_rows, t, mean_cent, var_adj, parallel_likelihood) +
      log_prior(proposal, yy, xx, first_b, first_sgt, m, A_rows, t,
                p_prior_mode, p_prior_scale, q_prior_mode, q_prior_scale,
                shrinkage);

    //Rcpp::Rcout << "last den: "<< last_density << std::endl;
    //Rcpp::Rcout << "den: "<< proposal_density << std::endl;

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
  const bool mean_cent;
  const bool var_adj;
  const RcppParallel::RVector<double> par_vec;

  DrawParallel(arma::mat& draws, arma::vec& densities, arma::vec& asums,
               const arma::mat& yy, const arma::mat& xx, const int first_b, const int first_sgt, Rcpp::NumericVector par_vec, const bool mean_cent, const bool var_adj)
    : draws(draws), densities(densities), asums(asums), yy(yy), xx(xx), first_b(first_b), first_sgt(first_sgt), par_vec(par_vec), mean_cent(mean_cent), var_adj(var_adj) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(int j = begin; j != end; ++j) {
      draw(draws, densities, asums, j, j - par_vec[1],
           par_vec[4], par_vec[3], par_vec[1], yy, xx, first_b, first_sgt, par_vec[5], par_vec[6], par_vec[7], mean_cent, var_adj,
           par_vec[8], par_vec[9], par_vec[10], par_vec[11], par_vec[12], par_vec[13]);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List sampler(const int N, const int n, const int m0, const int K, const double gamma,
                   const arma::vec init_mode, const arma::mat init_scale,
                   const bool output_as_input, const arma::mat old_chain, const bool new_chain,
                   const bool parallel, const bool parallel_likelihood,
                   const arma::mat yy, const arma::mat xx,
                   const int m, const int A_rows, const int t, const bool mean_cent, const bool var_adj,
                   const int first_b, const int first_sgt,
                   const double p_prior_mode, const double p_prior_scale,
                   const double q_prior_mode, const double q_prior_scale,
                   const double shrinkage, const bool progress_bar) {

  //Collect parameters
  Rcpp::NumericVector par_vec(13);
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
  par_vec[12] = shrinkage;
  par_vec[13] = parallel_likelihood;

  //Initialize the chains
  arma::mat draws(N * n + m0, init_mode.n_elem);

  //Initial chains/draws from multivariate normal...
  if(!output_as_input) {
    for(int i = 0; i != m0; ++i) {
      arma::vec xn = arma::randn(init_mode.n_elem);
      draws.row(i) = init_mode.t() + xn.t() * arma::chol(init_scale);
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
      log_like(state, yy, xx, first_b, first_sgt, m, A_rows, t, mean_cent, var_adj, parallel_likelihood) +
      log_prior(state, yy, xx, first_b, first_sgt, m, A_rows, t,
                p_prior_mode, p_prior_scale, q_prior_mode, q_prior_scale,
                shrinkage);
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
             gamma, K, n, yy, xx, first_b, first_sgt, m, A_rows, t, mean_cent, var_adj,
             p_prior_mode, p_prior_scale, q_prior_mode, q_prior_scale,
             shrinkage, parallel_likelihood);
      }
      last_row = last_row + n;
    }
  }

  //Parallel
  if(parallel == true) {
    DrawParallel Wrkr(draws, densities, asums,
                      yy, xx, first_b, first_sgt, par_vec, mean_cent, var_adj);
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

// [[Rcpp::export]]
void test() {
  float progress = 0.0;
  while (progress < 1.0) {
    int barWidth = 70;

    Rcpp::Rcout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) Rcpp::Rcout << "=";
      else if (i == pos) Rcpp::Rcout << ">";
      else Rcpp::Rcout << " ";
    }
    Rcpp::Rcout << "] " << int(progress * 100.0) << " %\r";
    Rcpp::Rcout.flush();

    progress += 0.16; // for demonstration only
  }
  Rcpp::Rcout << std::endl;
}




