
#Finds optimal initial values, helper for 'init_rbsvar()' - NO EXPORT
init_optim <- function(pre_init,
                       type,
                       B_inverse,
                       xy,
                       prior,
                       ols_cov,
                       min_density,
                       cpp_args,
                       method,
                       skip_hessian,
                       parallel_likelihood,
                       max_cores,
                       maxit,
                       trace,
                       REPORT,
                       verbose) {

  obj_full <- function(state) {
    ret <- -(

      log_like(state, xy$yy, xy$xx,
               cpp_args$first_b, cpp_args$first_sgt, cpp_args$first_garch, cpp_args$first_regime, cpp_args$first_yna,
               cpp_args$m, cpp_args$A_rows, cpp_args$t, cpp_args$regimes, cpp_args$yna_indices, cpp_args$B_inverse,
               cpp_args$mean_cent, cpp_args$var_adj, parallel_likelihood) +
        log_prior(state, xy$yy, xy$xx,
                  cpp_args$first_b, cpp_args$first_sgt, cpp_args$first_garch, cpp_args$first_regime, cpp_args$first_yna,
                  cpp_args$m, cpp_args$A_rows, cpp_args$t, cpp_args$regimes,
                  prior$A$mean, prior$A$cov, cpp_args$prior_A_diagonal, prior$B$mean, prior$B$cov,
                  prior$p[1], prior$p[2], prior$q[1], prior$q[2], cpp_args$dirichlet_alpha,
                  cpp_args$B_inverse)
    )
    if(ret == Inf) ret <- -min_density
    ret
  }

  # Number of cores used is restricted if 'max_cores' is not NA
  if(maxit > 0) {
    if(verbose == TRUE) cat("Searching for optimal initial parameters values...\n")
    if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = max_cores)
    opt_full <- optim(pre_init,
                      obj_full,
                      method = method,
                      control = list(maxit = maxit, trace = trace, REPORT = REPORT),
                      hessian = !skip_hessian)
    if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
  } else {
    if(verbose == TRUE) cat("Skipping the search for optimal initial parameters values...\n")
    opt_full <- list("par" = pre_init)
    skip_hessian <- TRUE
  }

  # Approximate the scale matrix based on numerically approximated Hessian
  # and collect the results
  if(!skip_hessian) {
    init_scale <- tryCatch({
      if(min(eigen(opt_full$hessian)$values) < 0) {
        loader <- function(x) {
          hes <- opt_full$hessian
          diag(hes) <- diag(hes) + x
          abs(0.001 - min(eigen(hes)$values))
        }
        load <- optimize(loader, interval = c(0,1))$minimum
        diag(opt_full$hessian) <- diag(opt_full$hessian) + load
      }
      chol2inv(chol(opt_full$hessian))
    }, error = function(e) {
      if(verbose == TRUE) cat("Inversion of numerical Hessian failed, but not to worry, initial covariance matrix will be acquired some other way. \n")
      NULL
    })
  } else {
    init_scale <- NULL
  }

  if(is.null(init_scale)) {
    init_scale <- diag(length(pre_init))*0.001

    # OLS covariance estimates can be used as initial covariance of autoregressive parameters
    if(!is.null(ols_cov)) {
      init_scale[1:nrow(ols_cov), 1:ncol(ols_cov)] <- ols_cov
    }

    # Coarse initial estimate for the covariance of the elements in B may be obtained via bootstrap
    if(verbose == TRUE) cat("Computing coarse bootstrap estimate for the initial covariance of the elements in B. \n")
    uu <- xy$yy - xy$xx %*% matrix(opt_full$par[1:ncol(ols_cov)], ncol = ncol(xy$yy))
    b_boot <- matrix(NA, nrow = 1000, ncol = cpp_args$first_sgt - cpp_args$first_b)
    for(i in 1:1000) {
      uu_rnd <- uu[sample.int(nrow(uu), nrow(uu), replace = TRUE),]
      sigma_rnd <- t(uu_rnd) %*% uu_rnd / nrow(uu_rnd)
      if(type == "svar") {
        B_init <- expm::sqrtm(sigma_rnd)
        Binv_init <- solve(B_init)
        if(B_inverse) b_init <- c(Binv_init) else b_init <- c(B_init)
      }
      if(type == "var") {
        B_init <- t(chol(sigma_rnd))
        Binv_init <- solve(B_init)
        if(B_inverse) b_init <- c(Binv_init)[which(!upper.tri(B_init))] else b_init <- c(B_init)[which(!upper.tri(B_init))]
      }
      b_boot[i,] <- b_init
    }
    init_b_cov <- crossprod(scale(b_boot, scale = FALSE)) / nrow(b_boot)
    if(min(eigen(init_b_cov)$val) < 0) diag(init_b_cov) <- diag(init_b_cov) + abs(min(eigen(init_b_cov)$val))
    if(min(eigen(init_b_cov)$val) < 1e-6) diag(init_b_cov) <- diag(init_b_cov) + 1e-6
    b_indices <- (cpp_args$first_b + 1):(cpp_args$first_sgt)
    init_scale[b_indices, b_indices] <- init_b_cov
  }

  # Collect and return the initial parameter values and the scale matrix
  init <- list("init_mode" = opt_full$par,
               "init_scale" = init_scale,
               "opt" = opt_full)
  init
}

#' Initializes rbsvar model (Documentation incomplete)
#'
#' @param y The data.
#' @param lags ...
#' @param constant ...
#' @param type ...
#' @param garch ...
#' @param mean_cent ...
#' @param var_adj ...
#' @param p_prior ...
#' @param q_prior ...
#' @param r_prior ...
#' @param B_prior ...
#' @param B_inverse ...
#' @param shrinkage ...
#' @param minnesota_means ...
#' @param prior ...
#' @param init ...
#' @param maxit ...
#' @param method ...
#' @param skip_hessian ...
#' @param trace ...
#' @param REPORT ...
#' @param min_density ...
#' @param parallel_likelihood ...
#' @param max_cores ...
#' @param verbose ...
#' @return A list that may be passed for example to \code{est_rbsvar()} or \code{eval_rbsvar()}.
#'
#' @export
init_rbsvar <- function(y,
                        lags = 1,
                        constant = TRUE,
                        type = c("svar", "var")[1],
                        garch = FALSE,
                        regimes = c(),
                        mean_cent = TRUE,
                        var_adj = FALSE,
                        p_prior = c(2, 1),
                        q_prior = c(1, 1),
                        dirichlet_alpha = 2,
                        B_prior = c(NA, NA),
                        B_inverse = TRUE,
                        shrinkage = Inf,
                        minnesota_means = NULL,
                        prior = list(),
                        init = NULL,
                        maxit = 0,
                        method = "L-BFGS-B",
                        skip_hessian = NULL,
                        trace = 1,
                        REPORT = 100,
                        min_density = -1000000,
                        parallel_likelihood = FALSE,
                        max_cores = NA,
                        verbose = TRUE) {

  start_time <- Sys.time()
  xy <- build_xy(y, lags, constant)
  xx <- xy$xx
  yy <- xy$yy
  m <- ncol(xy$yy)
  if(lags == 0 & shrinkage != Inf) {
    shrinkage <- Inf
    if(verbose == TRUE) cat("As 'lags = 0', 'shrinkage' parameter is omitted. \n")
  }
  if(var_adj == TRUE & mean_cent == FALSE) stop("'mean_cent == TRUE' required if 'var_adj == TRUE'")
  if(garch == TRUE & var_adj == FALSE) {
    if(verbose) cat("NOTE: 'var_adj == FALSE' and 'garch == TRUE' --> Unconditional expected value of the GARCH process cannot be fixed. \n")
  }

  # Missing values in data: (Only type = "var" allowed)

  if(sum(is.na(yy)) != 0) {
    if(type == "svar") stop("'type = svar' only supports balanced data for now (no missing values in 'y').")

    # Save the indices of y corresponding to missing values
    yna_indices <- which(is.na(yy))
    yna_init <- rep(NA, length(yna_indices))
    for(i in 1:length(yna_indices)) {
      the_col <- ceiling(yna_indices[i] / nrow(yy))
      yna_init[i] <- mean(yy[,the_col], na.rm = T)
    }

    # Remove NA rows from yy and xx to be used for calculation of initial values
    na_rows <- which(is.na(apply(yy, 1, sum)))
    xx <- xx[-na_rows,]
    yy <- yy[-na_rows,]
  } else {
    yna_indices <- c(-1)
    yna_init <- c()
  }


  # Prior:
  prior$p <- p_prior
  prior$q <- q_prior
  if(sum(is.na(B_prior)) == 0) {
    if(type == "svar") {
      prior$B <- list("mean" = c(diag(ncol(y)) * B_prior[1]),
                      "cov" = diag(ncol(y)^2) * B_prior[2])
    } else {
      prior$B <- list("mean" = diag(ncol(y))[which(!upper.tri(diag(ncol(y))))] * B_prior[1],
                      "cov" = diag((ncol(y) * (ncol(y) + 1)) / 2) * B_prior[2])
    }
  }
  if(is.null(prior$B)) prior$B <- list("mean" = matrix(-1), "cov" = matrix(-1))
  if(is.null(prior$A)) prior$A <- list("mean" = matrix(-1), "cov" = matrix(-1))

  # Initial values and Minnesota-prior:

  # A and Minnesota-prior
  if(!is.null(xx)) {
    if(ncol(xx) < nrow(xx) & shrinkage == Inf) {
      ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
      uu <- yy - xx %*% ols_est
      sigma <- t(uu) %*% uu / nrow(uu)
      ols_cov <- kronecker(sigma, chol2inv(chol(crossprod(xx))))
    } else {
      if(shrinkage == Inf) stop("Non-positive degrees of freedom. Prior required for identification.")
      if(prior$A$cov[1] != -1) {
        if(verbose == TRUE) cat("As 'prior$A' is spesified, 'shrinkage' parameter is omitted. \n")
      } else {
        xy0 <- build_xy0(y, lags, shrinkage, minnesota_means, constant)
        xx0 <- xy0$xx0
        yy0 <- xy0$yy0
        ols_est_prior <- chol2inv(chol(crossprod(xx0))) %*% t(xx0) %*% yy0
        sigma_prior <- diag(xy0$arsigmas)^2
        ols_cov_prior <- kronecker(sigma_prior, chol2inv(chol(crossprod(xx0))))
        prior$A$mean <- c(ols_est_prior)
        prior$A$cov <- ols_cov_prior

        xx_post <- rbind(xx, xx0)
        yy_post <- rbind(yy, yy0)
        ols_est_post <- chol2inv(chol(crossprod(xx_post))) %*% t(xx_post) %*% yy_post
        uu_post <- yy_post - xx_post %*% ols_est_post
        sigma_post <- t(uu_post) %*% uu_post / nrow(uu_post)
        ols_cov_post <- kronecker(sigma_post, chol2inv(chol(crossprod(xx_post))))
        ols_est <- ols_est_post
        ols_cov <- ols_cov_post
        sigma <- sigma_post
      }
    }
  } else {
    ols_est <- c()
    ols_cov <- NULL
    sigma <- t(yy) %*% yy / nrow(yy)
    xy$xx <- matrix(0)
  }

  if(verbose) {
    if(prior$A$mean[1] == -1 | prior$B$mean[1] == -1) cat("Note: Prior is improper. \n")
  }

  # B
  if(type == "svar") {
    B_init <- expm::sqrtm(sigma)
    right_permutation <- check_permutation(B_init)
    if(!right_permutation) B_init <- t(chol(B_init))
    Binv_init <- solve(B_init)
    if(B_inverse) b_init <- c(Binv_init) else b_init <- c(B_init)
  }
  if(type == "var") {
    B_init <- t(chol(sigma))
    Binv_init <- solve(B_init)
    if(B_inverse) b_init <- c(Binv_init)[which(!upper.tri(B_init))] else b_init <- c(B_init)[which(!upper.tri(B_init))]
  }

  # Sgt-parameters
  SGT_param <- matrix(0, nrow = ncol(y), ncol = 3)
  colnames(SGT_param) <- c("lambda", "p", "q")
  SGT_param[,"p"] <- log(p_prior[1])
  SGT_param[,"q"] <- q_prior[1]

  if(is.null(init)) {
    if(verbose) cat("Looking for initial sgt-parameters... \n")
    neg_ll_sgt <- function(par, e, max_ret = 1000000) {
      if(abs(par[1]) > 0.99) return(max_ret)
      ret <- -sum(
        sgt::dsgt(e, mu = 0, sigma = 1, lambda = par[1],
                  p = exp(par[2]), q = exp(par[3]), mean.cent = mean_cent, var.adj = var_adj, log = T)
      )
      if(abs(ret) == Inf | is.nan(ret)) return(max_ret)
      ret
    }
    fit_sgt <- function(e, init_sgt) {
      optim(init_sgt, neg_ll_sgt, method = "L-BFGS-B", e = e, upper = c(0.99, 1, 4))
    }
    E <- (xy$yy - xy$xx %*% matrix(ols_est, ncol = ncol(y))) %*% t(Binv_init)
    for(i in 1:ncol(E)) {
      try(SGT_param[i,] <- fit_sgt(e = E[,i], init_sgt = SGT_param[i,])$par,
          silent = TRUE)
    }
  }

  # Garch-parameters
  if(garch == TRUE) {
    GARCH_param <- matrix(0, nrow = m, ncol = 2)
    colnames(GARCH_param) <- c("alpha", "beta")
    GARCH_param[,"alpha"] <- 0.1
    GARCH_param[,"beta"] <- 0.1
  }

  # Regime-parameters (relative scales between regimes)
  if(length(regimes) != 0) {
    REGIME_param <- matrix(1, nrow = m, ncol = length(regimes))
    regime_init <- c(REGIME_param)
  } else {
    regimes <- c(-1)
    regime_init <- c()
  }

  # User provided initial values
  if(!is.null(init)) {
    SGT_param[] <- init[(length(ols_est) + length(b_init) + 1):(length(ols_est) + length(b_init) + length(SGT_param))]
    if(garch == TRUE) {
      GARCH_param[] <- init[(length(ols_est) + length(b_init) + length(SGT_param)):(length(ols_est) + length(b_init) + length(SGT_param) + length(GARCH_param))]
    }
  }

  # Check the feasibility of the initial values (p, q and alpha and beta)
  count <- 0
  if(var_adj == FALSE) {
    for(i in 1:m) {
      if(exp(SGT_param[i, "p"]) * exp(SGT_param[i, "q"]) <= 1) {
        if(count == 0 & verbose == TRUE) cat("Initial parameter values of 'p' and/or 'q' altered as the ones provided were not feasible. \n")
        count <- count + 1
        SGT_param[i, 2:3] <- c(log(2), log(1))
      }
    }
  } else {
    for(i in 1:m) {
      if(exp(SGT_param[i, "p"]) * exp(SGT_param[i, "q"]) <= 2) {
        if(count == 0 & verbose == TRUE) cat("Initial parameter values of 'p' and/or 'q' altered as the ones provided were not feasible. \n")
        count <- count + 1
        SGT_param[i, 2:3] <- c(log(2), log(2))
      }
    }
  }
  sgt_init <- c(SGT_param)
  if(garch == TRUE) {
    count <- 0
    for(i in 1:m) {
      if(GARCH_param[i, "alpha"] < 0 | GARCH_param[i, "alpha"] >= 1 | GARCH_param[i, "beta"] < 0 | GARCH_param[i, "beta"] >= 1) {
        if(count == 0 & verbose == TRUE) cat("Initial parameter values of 'alpha' and/or 'beta' altered as the ones provided were not feasible. \n")
        count <- count + 1
        GARCH_param[i, "alpha"] <- 0.1
        GARCH_param[i, "beta"] <- 0.1
      }
    }
    garch_init <- c(GARCH_param)
  } else {
    garch_init <- c()
  }

  if(is.null(init)) {
    if(lags > 0) {
      pre_init <- c(ols_est, b_init, sgt_init, garch_init, regime_init, yna_init)
    } else {
      pre_init <- c(b_init, sgt_init, garch_init, regime_init, yna_init)
    }
  } else {
    pre_init <- init
  }

  # Heuristic for skipping evaluation of Hessian in case of large parameter space
  if(is.null(skip_hessian)) {
    skip_hessian <- ifelse(length(pre_init) > 100, TRUE, FALSE)
    if(verbose == TRUE & skip_hessian == TRUE & maxit > 0) cat("Approximation of full Hessian to be skipped... \n")
  }
  if(maxit == 0) skip_hessian <- TRUE

  # Check whether prior covariance matrix of A is diagonal
  if(prior$A$cov[1] != -1) {
    if(sum(prior$A$cov[upper.tri(prior$A$cov)] != 0) + sum(prior$A$cov[lower.tri(prior$A$cov)] != 0) == 0) {
      prior_A_diagonal <- TRUE
      prior$A$cov <- matrix(sqrt(diag(prior$A$cov)), ncol = 1)
    } else {
      prior_A_diagonal <- FALSE
      warning("For now, NON-DIAGONAL prior covariance matrix for A leads to inefficient posterior evaluation, especially when A is large.")
    }
  } else {
    prior_A_diagonal <- FALSE
  }

  # Parameters required by cpp functions
  cpp_args <- list("m" = ncol(xy$yy),
                   "first_b" = length(pre_init) - length(b_init) - length(sgt_init) - length(garch_init) - length(regime_init) - length(yna_init),
                   "first_sgt" = length(pre_init) - length(sgt_init) - length(garch_init) - length(regime_init) - length(yna_init),
                   "first_garch" = length(pre_init) - length(garch_init) - length(regime_init) - length(yna_init),
                   "first_regime" = length(pre_init) - length(regime_init) - length(yna_init),
                   "first_yna" = length(pre_init) - length(yna_init),
                   "regimes" = regimes,
                   "yna_indices" = yna_indices,
                   "A_rows" = ncol(xy$xx),
                   "t" = nrow(xy$yy),
                   "prior_A_diagonal" = prior_A_diagonal,
                   "B_inverse" = B_inverse,
                   "mean_cent" = mean_cent,
                   "var_adj" = var_adj,
                   "dirichlet_alpha" = dirichlet_alpha)
  if(lags == 0) cpp_args$A_rows <- 0

  # Optimal initial values
  init <- init_optim(pre_init = pre_init,
                     type = type,
                     B_inverse = B_inverse,
                     xy = xy,
                     prior = prior,
                     ols_cov = ols_cov,
                     cpp_args = cpp_args,
                     method = method,
                     skip_hessian = skip_hessian,
                     min_density = min_density,
                     parallel_likelihood = parallel_likelihood,
                     max_cores = max_cores,
                     maxit = maxit,
                     trace = trace,
                     REPORT = REPORT,
                     verbose = verbose)

  if(verbose == TRUE) cat(paste0("DONE. Robust bayesian ", type, " model initialized.\n"))
  model <- list("y" = y,
                "xy" = xy,
                "lags" = lags,
                "type" = type,
                "constant" = constant,
                "garch" = garch,
                "prior" = prior,
                "init" = init,
                "cpp_args" = cpp_args,
                "time_init" = Sys.time() - start_time)
  model
}

#' Estimates rbsvar model by sampling from the posterior using DE-MC algorithm (Documentation incomplete)
#'
#' @param model A list outputted by \code{init_rbsvar()}.
#' @param N ...
#' @param n ...
#' @param K ...
#' @param m0 ...
#' @param rel_gamma ...
#' @param output ...
#' @param new_chain ...
#' @param parallel_chains ...
#' @param parallel_likelihood ...
#' @param max_cores ...
#' @param progress_bar ...
#' @param verbose ...
#' @return A list that containing (among other things) the chain(s) sampled from the posterior. May be passed for example to \code{irf()}.
#'
#' @export
est_rbsvar <- function(model,
                       N,
                       n = 2,
                       K = 2,
                       m0 = NA,
                       rel_gamma = 1,
                       output = NULL,
                       new_chain = FALSE,
                       parallel_chains = FALSE,
                       parallel_likelihood = FALSE,
                       max_cores = NA,
                       progress_bar = TRUE,
                       verbose = TRUE) {

  start_time <- Sys.time()
  d <- length(model$init$init_mode)
  gamma <- (2.38/sqrt(2*d)) * rel_gamma
  if(is.na(m0)) m0 <- n*d

  # Informational messages
  if(verbose == TRUE) cat("Running differential evolution Markov chains... \n")
  if(verbose == TRUE) cat(paste0("N = ", N, "\n"))
  if(verbose == TRUE) cat(paste0("pop. size (n) = ", n,"\n"))
  if(verbose == TRUE) cat(paste0("gamma = ", round(gamma, 2),"\n"))
  if(verbose == TRUE) cat(paste0("K = ", K, "\n"))
  if(verbose == TRUE) cat(paste0("m0 = ", m0, "\n"))

  # Check whether old output is provided as an input and act accordingly
  if(!is.null(output)) {
    if(new_chain == TRUE) {
      if(is.null(output$burn)) stop("If 'new_chain = TRUE' then 'output$burn' argument must exist.")
      old_chain <- output$chains[-c(1:(output$args$m0 + output$burn*output$args$n)),]
    } else {
      old_chain <- output$chains
    }
    output_as_input <- TRUE
  } else {
    old_chain <- diag(2)
    output_as_input <- FALSE
  }

  if(!output_as_input) {
    init_draws <- mvtnorm::rmvnorm(n = m0, mean = model$init$init_mode, sigma = model$init$init_scale)
  } else {
    init_draws <- matrix(-1, ncol = length(model$init$init_mode))
  }

  # Max number of threads used
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = max_cores)

  # Sampler runs here
  model_R <- build_model_R(model)
  sampler_out <- sampler(N = N, n = n, m0 = m0, K = K, gamma = gamma,
                         init_draws = init_draws,
                         output_as_input = output_as_input, old_chain = old_chain, new_chain = new_chain,
                         parallel = parallel_chains, parallel_likelihood = parallel_likelihood,
                         model_R = model_R, progress_bar = progress_bar)
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())

  # Clean up the output
  names(sampler_out) <- c("draws", "densities", "asums")
  sampler_out$densities[1:(m0 - n)] <- NA
  sampler_out$asums[1:m0] <- NA
  if(output_as_input == TRUE && new_chain == FALSE) {
    sampler_out$densities[1:nrow(old_chain)] <- output$densities
    sampler_out$asums[1:nrow(old_chain)] <- output$asums
  }

  # Average acceptance rate
  arate <- mean(na.omit(sampler_out$asums)/K)
  if(verbose == TRUE) cat(paste0("Average acceptance rate: ", round(arate, 3), "\n"))

  # Name the variables
  if(model$lags > 0) {
    colnames(sampler_out$draws) <- c(paste0("a", 1:model$cpp_args$first_b),
                                     paste0("b", (model$cpp_args$first_b + 1):model$cpp_args$first_sgt - model$cpp_args$first_b),
                                     paste0("sgt", (model$cpp_args$first_sgt + 1):length(model$init$init_mode) - model$cpp_args$first_sgt))
  } else {
    colnames(sampler_out$draws) <- c(paste0("b", 1:model$cpp_args$first_sgt),
                                     paste0("sgt", (model$cpp_args$first_sgt + 1):length(model$init$init_mode) - model$cpp_args$first_sgt))
  }
  if(model$garch == TRUE) {
    garch_indices <- (model$cpp_args$first_garch + 1):ncol(sampler_out$draws)
    colnames(sampler_out$draws)[garch_indices] <- paste0("garch", 1:length(garch_indices))
  }
  if(model$cpp_args$regimes[1] != -1) {
    regime_indices <- (model$cpp_args$first_regime + 1):ncol(sampler_out$draws)
    colnames(sampler_out$draws)[regime_indices] <- paste0("reg", 1:length(regime_indices))
  }
  if(model$cpp_args$yna_indices[1] != -1) {
    yna_indices <- (model$cpp_args$first_yna + 1):ncol(sampler_out$draws)
    colnames(sampler_out$draws)[yna_indices] <- paste0("y", 1:length(yna_indices))
  }

  time <- Sys.time() - start_time
  if(verbose == TRUE) cat("Time:", round(time, 2), attributes(time)$units, "\n")

  note <- ""
  if(output_as_input == TRUE & new_chain == FALSE) note <- paste0("New N = ", (N - output$args$N))

  # Collect and return everything
  ret <- list("chains" = sampler_out$draws,
              "densities" = sampler_out$densities,
              "asums" = sampler_out$asums,
              "arate" = arate,
              "args" = list("N" = N, "n" = n, "gamma" = gamma, "K" = K, "m0" = m0,
                            "parallel_chains" = parallel_chains, "parallel_likelihood" = parallel_likelihood, "max_cores" = max_cores),
              "model" = model,
              "time" = time,
              "note" = note)
  ret
}

#' C++ wrapper for evaluating the posterior density (Documentation incomplete)
#'
#' @param model A list outputted by \code{init_rbsvar()}.
#' @param par ...
#' @param parallel_likelihood ...
#' @param max_cores ...
#' @return Proportional posterior density at the given point.
#'
#' @export
eval_rbsvar <- function(model,
                        par = NULL,
                        parallel_likelihood = FALSE,
                        max_cores = NA) {
  if(is.null(par)) par <- model$init$init_mode
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = max_cores)
  xy <- model$xy
  cpp_args <- model$cpp_args
  prior <- model$prior
  ret <- (
    log_like(par, xy$yy, xy$xx,
                      cpp_args$first_b, cpp_args$first_sgt, cpp_args$first_garch, cpp_args$first_regime, cpp_args$first_yna,
                      cpp_args$m, cpp_args$A_rows, cpp_args$t, cpp_args$regimes, cpp_args$yna_indices, cpp_args$B_inverse,
                      cpp_args$mean_cent, cpp_args$var_adj, parallel_likelihood) +
      log_prior(par, xy$yy, xy$xx,
                         cpp_args$first_b, cpp_args$first_sgt, cpp_args$first_garch, cpp_args$first_regime, cpp_args$first_yna,
                         cpp_args$m, cpp_args$A_rows, cpp_args$t, cpp_args$regimes,
                         prior$A$mean, prior$A$cov, cpp_args$prior_A_diagonal, prior$B$mean, prior$B$cov,
                         prior$p[1], prior$p[2], prior$q[1], prior$q[2], cpp_args$dirichlet_alpha,
                         cpp_args$B_inverse)
  )
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
  ret
}




