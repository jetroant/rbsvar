
#Finds optimal initial values NO EXPORT
init_optim <- function(pre_init_state,
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
      log_like(state, xy$yy, xy$xx, cpp_args$first_b, cpp_args$first_sgt, cpp_args$m, cpp_args$A_rows, cpp_args$t, cpp_args$mean_cent, cpp_args$var_adj, parallel_likelihood) +
        log_prior(state, xy$yy, xy$xx, cpp_args$first_b, cpp_args$first_sgt, cpp_args$m, cpp_args$A_rows, cpp_args$t,
                  prior$p[1], prior$p[2], prior$q[1], prior$q[2], prior$shrinkage)
    )
    if(ret == Inf) ret <- -min_density
    ret
  }


  #Number of cores used is restricted if 'max_cores' is not NA
  if(maxit > 0) {
    if(verbose == TRUE) cat("Searching for optimal initial parameters values...\n")
    if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = max_cores)
    opt_full <- optim(pre_init_state,
                      obj_full,
                      method = method,
                      control = list(maxit = maxit, trace = trace, REPORT = REPORT),
                      hessian = !skip_hessian)
    if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
  } else {
    if(verbose == TRUE) cat("Skipping the search for optimal initial parameters values...\n")
    opt_full <- list("par" = pre_init_state)
    skip_hessian <- TRUE
  }

  #Approximate the scale matrix based on numerically approximated Hessian
  #and collect the results
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
      if(verbose == TRUE) cat("Inversion of numerical Hessian failed, initial scale and orientation might be inefficient.\n")
      NULL
    })
  } else {
    init_scale <- NULL
  }

  #If normal approximation fails, OLS covariance estimates are used for
  #autoregressive parameters
  if(is.null(init_scale) & maxit > 0) {
    init_scale <- diag(length(pre_init_state))*0.001
    if(!is.null(ols_cov)) {
      init_scale[1:nrow(ols_cov), 1:ncol(ols_cov)] <- ols_cov
    }
  } else {
    init_scale <- matrix(0.001)
  }

  #Collect and return the initial parameter values and the scale matrix
  init <- list("init_mode" = opt_full$par,
               "init_scale" = init_scale,
               "opt" = opt_full)
  init
}

#Initializes rbsvar model
init_rbsvar <- function(y,
                        lags = 1,
                        constant = TRUE,
                        type = c("svar", "var")[1],
                        mean_cent = TRUE,
                        var_adj = FALSE,
                        shrinkage = Inf,
                        p_prior = c(log(2), 1),
                        q_prior = c(log(1), 2),
                        init = NULL,
                        method = "L-BFGS-B",
                        skip_hessian = NULL,
                        min_density = -1000000,
                        parallel_likelihood = FALSE,
                        max_cores = NA,
                        maxit = 500,
                        trace = 1,
                        REPORT = 10,
                        verbose = TRUE) {

  start_time <- Sys.time()
  xy <- build_xy(y, lags, constant)
  xx <- xy$xx
  yy <- xy$yy
  m <- ncol(xy$yy)

  ### Prior ###
  prior <- list("p" = p_prior, "q" = q_prior, "shrinkage" = shrinkage)

  ### Pre-initial values ###

  # 1) Autoregressive parameters
  if(!is.null(xx)) {
    if(ncol(xx) < nrow(xx)) {
      ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
      uu <- yy - xx %*% ols_est
      sigma <- t(uu) %*% uu / nrow(uu)
      ols_cov <- kronecker(sigma, chol2inv(chol(crossprod(xx))))
    } else {
      if(is.null(prior)) stop("Non-positive degrees of freedom. Prior required for identification.")
      if(verbose == TRUE) cat("OLS estimates not available, prior mean used for initialization of autoregressive parameters.")
      A0 <- rbind(matrix(0, ncol = m, nrow = 1),
                  diag(m),
                  matrix(0, nrow = (lags-1)*m, ncol = m))
      if(constant == FALSE) A0 <- A0[-1,]
      ols_est <- A0
      uu <- yy - xx %*% ols_est
      sigma <- t(uu) %*% uu / nrow(uu)
      ols_cov <- NULL
    }
  } else {
    ols_cov <- NULL
    sigma <- t(yy) %*% yy / nrow(yy)
    xy$xx <- matrix(0)
  }

  # 2) B-matrix
  if(type == "svar") {
    B_init <- solve(diag(sqrt(diag(sigma))))
    b_init <- c(B_init)
  }
  if(type == "var") {
    B_init <- solve(chol(sigma))
    b_init <- B_init[which(!upper.tri(B_init))]
  }

  # 3) Sgt-parameters
  SGT_param <- matrix(0, nrow = ncol(y), ncol = 3)
  colnames(SGT_param) <- c("lambda", "p", "q")
  SGT_param[,"p"] <- p_prior[1]
  SGT_param[,"q"] <- q_prior[1]
  sgt_init <- c(SGT_param)

  if(is.null(init)) {
    if(lags > 0) {
      pre_init <- c(ols_est, b_init, sgt_init)
    } else {
      pre_init <- c(b_init, sgt_init)
    }
  } else {
    pre_init <- init
  }

  if(is.null(skip_hessian)) {
    skip_hessian <- ifelse(length(pre_init) > 100, TRUE, FALSE)
    if(verbose == TRUE & skip_hessian == TRUE) cat("Approximation of full Hessian to be skipped... \n")
  }

  #Parameters required by cpp functions
  cpp_args <- list("m" = ncol(xy$yy),
                   "first_b" = length(pre_init) - length(b_init) - 3*ncol(xy$yy),
                   "first_sgt" = length(pre_init) - 3*ncol(xy$yy),
                   "A_rows" = ncol(xy$xx),
                   "t" = nrow(xy$yy),
                   "mean_cent" = mean_cent,
                   "var_adj" = var_adj)
  if(lags == 0) cpp_args$A_rows <- 0

  ### Optimal initial values ###

  init <- init_optim(pre_init_state = pre_init,
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

  if(verbose == TRUE) cat(paste0("Robust bayesian ", type, " model initialized.\n"))
  model <- list("y" = y,
                "xy" = xy,
                "lags" = lags,
                "type" = type,
                "prior" = prior,
                "init" = init,
                "cpp_args" = cpp_args,
                "time_init" = Sys.time() - start_time)
  model
}

#Estimates rbsvar model
est_rbsvar <- function(model,
                       N,
                       n = 4,
                       K = 2,
                       m0 = NA,
                       gamma = NA,
                       output = NULL,
                       new_chain = FALSE,
                       parallel_chains = FALSE,
                       parallel_likelihood = FALSE,
                       max_cores = NA,
                       progress_bar = TRUE,
                       verbose = TRUE) {

  start_time <- Sys.time()
  d <- length(model$init$init_mode)
  if(is.na(gamma)) gamma <- 2.38/sqrt(2*d)
  if(is.na(m0)) m0 <- n*d

  ### Informational messages
  if(verbose == TRUE) cat("Running differential evolution Markov chains... \n")
  if(verbose == TRUE) cat(paste0("N = ", N, "\n"))
  if(verbose == TRUE) cat(paste0("pop. size (n) = ", n,"\n"))
  if(verbose == TRUE) cat(paste0("gamma = ", round(gamma, 2),"\n"))
  if(verbose == TRUE) cat(paste0("K = ", K, "\n"))
  if(verbose == TRUE) cat(paste0("m0 = ", m0, "\n"))

  ### Check whether old output is provided as an input and act accordingly
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

  ### Max number of threads used
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = max_cores)

  ### Sampler runs here
  sampler_out <- sampler(N = N, n = n, m0 = m0, K = K, gamma = gamma,
                         init_mode = model$init$init_mode, init_scale = model$init$init_scale,
                         output_as_input = output_as_input, old_chain = old_chain, new_chain = new_chain,
                         parallel = parallel_chains, parallel_likelihood = parallel_likelihood,
                         yy = model$xy$yy, xx = model$xy$xx,
                         m = model$cpp_args$m, A_rows = model$cpp_args$A_rows, t = model$cpp_args$t, mean_cent = model$cpp_args$mean_cent, var_adj = model$cpp_args$var_adj,
                         first_b = model$cpp_args$first_b, first_sgt = model$cpp_args$first_sgt,
                         p_prior_mode = model$prior$p[1], p_prior_scale = model$prior$p[2],
                         q_prior_mode = model$prior$q[1], q_prior_scale = model$prior$q[2],
                         shrinkage = model$prior$shrinkage, progress_bar = progress_bar)
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())

  ### Clean up the output
  names(sampler_out) <- c("draws", "densities", "asums")
  sampler_out$densities[1:(m0 - n)] <- NA
  sampler_out$asums[1:m0] <- NA
  if(output_as_input == TRUE && new_chain == FALSE) {
    sampler_out$densities[1:nrow(old_chain)] <- output$densities
    sampler_out$asums[1:nrow(old_chain)] <- output$asums
  }

  ### Average acceptance rate
  arate <- mean(na.omit(sampler_out$asums)/K)
  if(verbose == TRUE) cat(paste0("Average acceptance rate: ", round(arate, 3), "\n"))

  ### Name the variables
  if(model$lags > 0) {
    colnames(sampler_out$draws) <- c(paste0("a", 1:model$cpp_args$first_b),
                                     paste0("b", (model$cpp_args$first_b + 1):model$cpp_args$first_sgt - model$cpp_args$first_b),
                                     paste0("sgt", (model$cpp_args$first_sgt + 1):ncol(sampler_out$draws) - model$cpp_args$first_sgt))
  } else {
    colnames(sampler_out$draws) <- c(paste0("b", 1:model$cpp_args$first_sgt),
                                     paste0("sgt", (model$cpp_args$first_sgt + 1):ncol(sampler_out$draws) - model$cpp_args$first_sgt))

  }
    time <- Sys.time() - start_time
  if(verbose == TRUE) cat("Time:", round(time, 2), attributes(time)$units, "\n")

  note <- ""
  if(output_as_input == TRUE & new_chain == FALSE) note <- paste0("New N = ", (N - output$args$N))

  #Collect and return everything
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

#c++ wrapper for evaluating the posterior density
eval_rbsvar <- function(model, par = NULL, parallel_likelihood = FALSE, max_cores = NA) {
  if(is.null(par)) par <- model$init$init_mode
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = max_cores)
  xy <- model$xy
  cpp_args <- model$cpp_args
  prior <- model$prior
  ret <- (
    log_like(par, xy$yy, xy$xx, cpp_args$first_b, cpp_args$first_sgt, cpp_args$m, cpp_args$A_rows, cpp_args$t, cpp_args$mean_cent, cpp_args$var_adj, parallel_likelihood) +
      log_prior(par, xy$yy, xy$xx, cpp_args$first_b, cpp_args$first_sgt, cpp_args$m, cpp_args$A_rows, cpp_args$t,
                prior$p[1], prior$p[2], prior$q[1], prior$q[2], prior$shrinkage)
  )
  if(!is.na(max_cores)) RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
  ret
}


