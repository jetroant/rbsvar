
# Helper function - NO EXPORT
list2conStat <- function(chain_list,
                         variable,
                         index_split = NULL,
                         n_eff = FALSE) {

  if(is.null(nrow(chain_list[[1]]))) {
    for(i in 1:length(chain_list)) chain_list[[i]] <- matrix(chain_list[[i]], ncol = 1)
  }

  if(is.null(index_split)) index_split <- parallel::splitIndices(nrow(chain_list[[1]]), 2)
  sub_chain_list <- list()
  for(i in 1:length(chain_list)) {
    sub_chain_list[[(i*2-1)]] <- chain_list[[i]][index_split[[1]], variable]
    sub_chain_list[[(i*2)]] <- chain_list[[i]][index_split[[2]], variable]
  }
  ret <- rep(NA, 2)
  names(ret) <- c("R_hat", "n_eff")

  #R_hat
  m <- length(sub_chain_list)
  n <- length(sub_chain_list[[1]])
  sub_chain_means <- sapply(sub_chain_list, mean)
  overall_mean <- mean(sub_chain_means)
  B <- (n/(m-1)) * sum((sub_chain_means - overall_mean)^2)
  W <- mean(sapply(sub_chain_list, var))
  var_hat <- ((n-1)/n) * W + (1/n) * B
  ret["R_hat"] <- sqrt(var_hat/W)

  #n_eff
  if(n_eff == TRUE) {
    acf_mod <- function(x) acf(x, lag.max = (n-2), plot = FALSE)$acf[,1,1]
    roo_hat <- apply(sapply(sub_chain_list, acf_mod), 1, mean)
    lag <- 1
    keep_going <- TRUE
    while(keep_going == TRUE) {
      if(sum(c(roo_hat[lag], roo_hat[lag + 1]) < 0) == 2) keep_going <- FALSE
      lag <- lag + 1
    }
    roo_hat <- roo_hat[1:(lag-2)]
    ret["n_eff"] <- (m*n)/(1+2*sum(roo_hat))
  }

  ret
}

###########################
### Misc. res functions ###
###########################

#' R-hat convergence statistic (Documentation incomplete)
#'
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param burn ...
#' @param n_eff ...
#' @param graph ...
#' @param verbose ...
#' @return A matrix of element-wise convergence statistics in the first column.
#'
#' @export
convergence <- function(output,
                        burn,
                        n_eff = FALSE,
                        graph = TRUE,
                        verbose = TRUE) {

  n <- output$args$n
  N <- output$args$N
  m0 <- output$args$m0
  chains <- output$chains
  burn <- burn + 1

  if(verbose == TRUE) cat("Computing convergence statistics...")

  #Collect chains
  chain_list <- list()
  for(i in 1:n) {
    row_indices <- seq(from = m0 + i - n, by = n, length.out = N + 1)
    chain_list[[i]] <- chains[row_indices,][-c(1:burn),]
  }

  #Compute R_hat and n_eff
  ret_mat <- matrix(NA, ncol = 2, nrow = ncol(chains))
  colnames(ret_mat) <- c("R_hat", "n_eff")
  rownames(ret_mat) <- colnames(chains)
  index_split <- parallel::splitIndices(nrow(chain_list[[1]]), 2)
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = ncol(chains), style = 3)
  for(i in 1:ncol(chains)) {
    ret_mat[i,] <- list2conStat(chain_list = chain_list,
                                variable = i,
                                index_split = index_split,
                                n_eff = n_eff)
    if(verbose == TRUE) setTxtProgressBar(pb, i)
  }
  if(verbose == TRUE) close(pb)

  if(graph == TRUE) {
    plot(ts(ret_mat[,"R_hat"]), main = "", ylab = "R_hat", xlab = "Variable")
    abline(h = 1.1, lty = 2)
  }
  ret_mat
}

#' Plots the chains and marginal distirbution w.r.t one parameter (Documentation incomplete)
#'
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param variable ...
#' @param burn ...
#' @param breaks ...
#' @param real_val ...
#' @param n_eff ...
#'
#' @export
plot_chains <- function(output,
                        variable = "b1",
                        burn = 0,
                        breaks = "Sturges",
                        real_val = NA,
                        n_eff = FALSE) {

  n <- output$args$n
  N <- output$args$N
  m0 <- output$args$m0
  chains <- output$chains
  burn <- burn + 1

  #Is B_inv needed?
  if(grepl("_", variable) == TRUE) {
    variable <- strsplit(variable, "_")[[1]][1]
    B_inv <- TRUE
  } else {
    B_inv <- FALSE
  }

  #Collect chains
  chain_list <- list()
  for(i in 1:n) {

    #Pick one chain
    row_indices <- seq(from = m0 + i - n, by = n, length.out = N + 1)
    one_chain <- chains[row_indices,]

    #Pick the variable
    if(B_inv == TRUE) {
      b_chain <- one_chain[,output$indices$b_indices]
      b_inv_chain <- matrix(NA, ncol = output$args$m^2, nrow = nrow(b_chain))
      for(j in 1:nrow(b_chain)) {
        B <- diag(output$args$m)
        B[B != 1] <- b_chain[j,]
        b_inv_chain[j,] <- c(solve(B))
      }
      colnames(b_inv_chain) <- paste0("b", 1:ncol(b_inv_chain))
      chain_list[[i]] <- b_inv_chain[,variable]

    } else {
      chain_list[[i]] <- one_chain[,variable]
    }
  }

  #Convergence diagnostics
  chain_list_burned <- chain_list
  for(i in 1:length(chain_list_burned)) {
    chain_list_burned[[i]] <- chain_list_burned[[i]][-c(1:burn)]
  }
  conv <- list2conStat(chain_list = chain_list_burned,
                       variable = 1,
                       n_eff = n_eff)

  #Plot the chains
  title <- ifelse(B_inv == FALSE, variable, paste0(variable, "_inv"))
  par(mfrow = c(2,1))
  par(mar = c(3,4.3,3,1))
  ylims <- c(Inf, -Inf)
  for(i in 1:n) {
    if(min(chain_list[[i]], na.rm = T) < ylims[1]) ylims[1] <- min(chain_list[[i]], na.rm = T)
    if(max(chain_list[[i]], na.rm = T) > ylims[2]) ylims[2] <- max(chain_list[[i]], na.rm = T)
  }
  for(i in 1:length(chain_list)) {
    if(i == 1) {
      plot(ts(chain_list[[i]]), main = title,
           xlab = "", ylab = "Parameter value", ylim = ylims)
    } else {
      lines(chain_list[[i]])
    }
  }
  if(burn > 1) abline(v = burn, lty = 2, lwd = 2)

  #Plot histogram
  if(n_eff == TRUE) title <- paste0("R_hat = ", round(conv["R_hat"], 2), " / n_eff = ", round(conv["n_eff"], 2))
  if(n_eff == FALSE) title <- paste0("R_hat = ", round(conv["R_hat"], 2))
  full_post <- c()
  for(i in 1:length(chain_list_burned)) {
    full_post <- c(full_post, chain_list_burned[[i]])
  }
  hist(full_post, probability = TRUE, breaks = breaks,
       main = title, xlab = "", ylab = "Density")
  if(!is.na(real_val)) abline(v = real_val, lty = 2, lwd = 3)
  par(mfrow = c(1,1))
}

#' Plots the densities of the chains as a function of iterations (Documentation incomplete)
#'
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param verbose ...
#' @param drop ...
#'
#' @export
plot_densities <- function(output,
                           verbose = TRUE,
                           drop = 1) {
  model <- output$model
  d <- output$densities[-c(1:output$args$m0)]
  n <- output$args$n
  init_den <- eval_rbsvar(output$model)
  par(mfrow = c(n,1))
  for(i in 1:n) {
    the_d <- d[seq(from = i, by = n, length.out = length(d)/n)][-c(1:drop)]
    plot(ts(the_d, start = drop), ylab = paste0("Chain ", i))
    abline(h = init_den, lty = 2)
    if(i == 1) {
      if(min(d[-c(1:(drop*n))]) < init_den) legend("topleft", c("Initial density"), lty = 2, bty = "n")
    }
    if(verbose) cat(paste0("Chain ", i, " - Last density: ", the_d[length(the_d)], " / Max density: ", max(the_d), "\n"))
  }
  if(verbose) cat(paste0("Initial density: ", init_den, "\n"))
  par(mfrow = c(1,1))
}

#' Convenience function returning a sample from the posterior (Documentation incomplete)
#'
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param burn ...
#' @param N ...
#' @return A numerical matrix.
#'
#' @export
post_sample <- function(output, burn = 0, N = NA) {
  s <- output$chains[-c(1:output$args$m0),]
  d <- output$densities[-c(1:output$args$m0)]
  if(burn != 0) {
    s <- s[-c(1:(burn*output$args$n)),]
    d <- d[-c(1:(burn*output$args$n))]
  }
  if(!is.na(N)) {
    rndi <- sample.int(nrow(s), N, replace = T)
    s <- s[rndi,]
    d <- d[rndi]
  }
  ret <- list("s" = s,
              "d" = d)
}

############
### IRFs ###
############

#' Computes impulse response functions (Documentation incomplete)
#'
#' @param model A list outputted by \code{init_rbsvar()}.
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param burn ...
#' @param horizon ...
#' @param N ...
#' @param cumulate ...
#' @param shock_sizes ...
#' @param shocks ...
#' @param verbose ...
#' @param parallel ...
#' @return A list that may be passed to \code{irf_plot()}.
#'
#' @export
irf <- function(model,
                output,
                burn = 0,
                horizon = 40,
                N = 1000,
                cumulate = c(),
                shock_sizes = NULL,
                shocks = NULL,
                verbose = TRUE,
                parallel = FALSE) {

  if(requireNamespace("expm", quietly = TRUE)) {
    #Do nothing
  } else {
    stop("Package 'expm' needed for computation of impulse response functions. \n 'expm' was not found. The package can be installed by 'install.packages('expm')'. \n")
  }

  # Posterior sample
  post <- post_sample(output, burn, N)
  s <- post$s

  # Collect inputs
  if(model$type != "svar") stop("model$type must be 'svar'")
  m <- ncol(model$y)
  p <- model$lags
  ret <- list()
  if(is.null(shock_sizes)) shock_sizes <- rep(1, m)
  if(is.null(shocks)) shocks <- 1:m
  b_indices <- (model$cpp_args$first_b + 1):model$cpp_args$first_sgt
  a_indices <- 1:model$cpp_args$first_b
  constant <- model$constant
  B_inverse <- model$cpp_args$B_inverse

  # For Rcpp implementation to come:
  #if(length(cumulate) == 0) cumulate <- c(-1)
  #ret <- irf_cpp(s, horizon = 3, cumulate,
  #               shock_sizes, shocks,
  #               model$cpp_args$A_rows, model$cpp_args$first_b, model$cpp_args$first_sgt, m,
  #               B_inverse, parallel)

  #Compute and collect B and A sample
  a_sample <- matrix(NA, nrow = nrow(s), ncol = length(a_indices))
  b_sample <- matrix(NA, nrow = nrow(s), ncol = length(b_indices))
  for(i in 1:nrow(s)) {
    a_sample[i,] <- s[i, a_indices]
    b_sample[i,] <- s[i, b_indices]
    if(B_inverse) b_sample[i,] <- c(solve(matrix(b_sample[i,], ncol = m)))
  }

  for(shock_index in 1:m) {

    if(!(shock_index %in% shocks)) {
      ret[[shock_index]] <- NULL
      next
    }

    e <- rep(0, m)
    e[shock_index] <- shock_sizes[shock_index]
    irfs <- array(NA, dim = c(m, horizon+1, nrow(s)))

    if(verbose == TRUE) print(paste0("Computing impulse responses... (", which(shocks == shock_index), "/", length(shocks), ")"))
    if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = nrow(s), style = 3)
    for(row_index in 1:nrow(s)) {

      B <- diag(m)
      B[,] <- b_sample[row_index,]
      A <- matrix(a_sample[row_index,], ncol = m)
      AA <- stackA(A, constant = constant)

      for(h in 1:(horizon+1)) {
        if(h == 1) {
          zero <- B %*% e
          zero_long <- c(zero, rep(0, (ncol(AA)-length(zero))))
          irfs[,h,row_index] <- zero
        } else {
          irfs[,h,row_index] <- (expm::`%^%`(AA, (h-1)) %*% zero_long)[1:length(zero)]
        }
      }

      if(verbose == TRUE) setTxtProgressBar(pb, row_index)
    }
    if(verbose == TRUE) close(pb)

    if(length(cumulate) > 0) {
      for(i in 1:length(cumulate)) {
        for(j in 1:N) {
          irfs[cumulate[i],,j] <- cumsum(irfs[cumulate[i],,j])
        }
      }
    }
    ret[[shock_index]] <- irfs
  }

  ret
}

#' Plots impulse response functions (Documentation incomplete)
#'
#' @param irf_obj A list outputted by \code{irf()}.
#' @param varnames ...
#' @param probs ...
#' @param mar ...
#' @param mfrow ...
#' @param color ...
#' @param leg ...
#' @param plot_median ...
#' @param normalize ...
#'
#' @export
irf_plot <- function(irf_obj,
                     varnames = NULL,
                     probs = c(0.05, 0.16, 0.84, 0.95),
                     mar = c(2,2,2,1),
                     mfrow = NULL,
                     color = "tomato",
                     leg = TRUE,
                     plot_median = FALSE,
                     normalize = NULL) {

  if(max(probs) > 0.99 | min(probs) < 0.01) stop("Values of 'probs' need to be between 0.01 and 0.99.")
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)

  if(requireNamespace("fanplot", quietly = TRUE)) {
    #Do nothing
  } else {
    stop("Package 'fanplot' needed for plotting the impulse response functions. \n 'fanplot' was not found. The package can be installed by 'install.packages('fanplot')'. \n")
  }

  m <- length(irf_obj)
  for(i in 1:length(irf_obj)) if(nrow(irf_obj[[i]]) > m) m <- nrow(irf_obj[[i]])
  if(is.null(varnames)) varnames <- paste0("var", 1:m)
  if(is.null(mfrow)) mfrow <- c(m, m)
  par(mar = mar)
  par(mfrow = mfrow)
  indexmat <- matrix(1:m^2, ncol = m)
  row <- 0
  col <- 0
  for(fig_index in 1:m^2) {

    if((fig_index-1) %% m == 0) {
      row <- row + 1
      col <- 1
    } else {
      col <- col + 1
    }

    if(length(irf_obj) < row) next
    if(is.null(irf_obj[[row]])) next

    sub_irfs <- t(irf_obj[[row]][col,,])
    if(!is.null(normalize)) sub_irfs <- sub_irfs * normalize[col]
    median_sub_irfs <- ts(apply(sub_irfs, 2, median), start = 0)

    quant <- function(column) quantile(column, probs = probs)
    quantiles_sub_irfs <- apply(sub_irfs, 2, quant)

    plot(median_sub_irfs, lwd = 2, lty = 2, col = color, ylab = "", xlab = "",
         main = paste0("Shock ", row, " on ", varnames[col]),
         ylim = c(min(quantiles_sub_irfs), max(quantiles_sub_irfs)))
    grid()
    fanplot::fan(data = quantiles_sub_irfs, data.type = "values", probs = probs,
                 start = 0, fan.col = colorRampPalette(c(color, "white")),
                 rlab = NULL, ln = NULL)
    abline(h = 0, lwd = 2, lty = 2)
    if(plot_median) lines(median_sub_irfs, lwd = 2, lty = 1)

    post_mass <- (max(probs[-length(probs)]) - min(probs[-1]))*100
    if(col == 1 & row == 1 & leg == TRUE) legend("topleft", c(paste0(post_mass,"% of post. prob. mass")), lwd = 0, bty = "n", col = "tomato")
  }
  par(mfrow = c(1,1))
}

#######################################################
### Shock decompositions and narrative restrictions ###
#######################################################

shock_decomp <- function(model,
                         output,
                         burn = 0,
                         N = 1000,
                         verbose = TRUE) {

  if(model$type != "svar") stop("model$type must be 'svar'")
  m <- ncol(model$y)
  p <- model$lags
  yy <- model$xy$yy
  xx <- model$xy$xx
  n <- output$args$n
  m0 <- output$args$m0
  chains <- output$chains
  ret_E <- array(NA, dim = c(N, nrow(yy), m))
  ret_U <- list()
  for(i in 1:m) ret_U[[i]] <- array(NA, dim = c(N, nrow(yy), m))
  b_indices <- (model$cpp_args$first_b + 1):model$cpp_args$first_sgt
  a_indices <- 1:model$cpp_args$first_b

  #Collect A and B sample
  burn <- burn + 1
  b_list <- list()
  a_list <- list()
  for(i in 1:n) {

    #Pick one chain
    row_indices <- seq(from = m0 + i - n, by = n, length.out = output$args$N + 1)
    one_chain <- chains[row_indices,][-c(1:burn),]

    #Pick A
    a_list[[i]] <- one_chain[,a_indices]

    #Pick B
    b_list[[i]] <- one_chain[,b_indices]
  }
  for(i in 1:n) {
    if(i == 1) {
      full_b <- b_list[[i]]
      full_a <- a_list[[i]]
    } else {
      full_b <- rbind(full_b, b_list[[i]])
      full_a <- rbind(full_a, a_list[[i]])
    }
  }

  #Sample N draws and compute B_inv
  rows <- sample.int(nrow(full_a), N, replace = TRUE)
  full_a <- full_a[rows,]
  full_b <- full_b[rows,]
  full_binv <- matrix(NA, ncol = m^2, nrow = nrow(full_b))
  for(j in 1:nrow(full_b)) {
    B <- diag(m)
    B[,] <- full_b[j,]
    full_binv[j,] <- c(solve(B))
  }

  if(verbose == TRUE) print(paste0("Computing shock decompositions..."))
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {

    B <- diag(m)
    B_inv <- diag(m)
    B[,] <- full_b[i,]
    B_inv[,] <- full_binv[i,]
    A <- matrix(full_a[i,], ncol = m)
    U <- yy - xx %*% A
    E <- U %*% t(B)
    ret_E[i,,] <- E
    for(j in 1:m) {
      ret_U[[j]][i,,] <- matrix(c(t(E)) * B_inv[j,], byrow = TRUE, ncol = m)
    }
    if(verbose == TRUE) setTxtProgressBar(pb, i)
  }
  if(verbose == TRUE) close(pb)

  ret <- list("E" = ret_E,
              "U" = ret_U,
              "par" = list("a" = full_a, "b" = full_b, "b_inv" = full_binv))
  ret
}

shock_decomp_pick <- function(decomp_obj,
                              show_date,
                              start_date,
                              freq) {

  row_indices <- ts(1:dim(decomp_obj$E)[2],
                    start = start_date,
                    frequency = freq)
  row_index <- window(row_indices,
                      start = show_date,
                      end = show_date)

  E <- decomp_obj$E[,row_index,]

  m <- dim(decomp_obj$U[[1]])[3]
  U <- list()
  for(i in 1:m) {
    U[[i]] <- decomp_obj$U[[i]][, row_index, ]
  }

  ret <- list("E" = E,
              "U" = U)
  ret
}

narrative_check <- function(decomp_picked,
                            shock = 1,
                            variable = 1,
                            total = TRUE,
                            retvec = FALSE,
                            verbose = TRUE) {

  if(!is.null(names(decomp_picked))) {
    N <- dim(decomp_picked$U[[variable]])[1]
    agree <- rep(NA, N)
    for(i in 1:N) {
      own_contribution <- abs(decomp_picked$U[[variable]][i, shock])
      other_contributions <- abs(decomp_picked$U[[variable]][i, -shock])
      if(total == TRUE) agree[i] <- ifelse(own_contribution > sum(other_contributions), TRUE, FALSE)
      if(total == FALSE) agree[i] <- ifelse(own_contribution > max(other_contributions), TRUE, FALSE)
    }

    if(verbose == TRUE) cat(round(mean(agree)*100, 2),
                            "% of the posterior sample agrees with the restriction. \n")
    if(retvec == TRUE) return(agree)
    if(retvec == FALSE) return(mean(agree))

  } else {
    M <- length(decomp_picked)
    N <- dim(decomp_picked[[1]]$U[[variable]])[1]
    agree_mat <- matrix(NA, nrow = N, ncol = M)
    for(j in 1:M) {
      agree_mat[,j] <- narrative_check(decomp_picked[[j]],
                                       shock = shock,
                                       variable = variable,
                                       total = total,
                                       verbose = FALSE,
                                       retvec = TRUE)
    }
    agreed <- mean(apply(agree_mat, 1, sum) == M)
    if(verbose == TRUE) cat(round(agreed*100, 2),
                            "% of the posterior sample agrees with the restrictions. \n")
    if(retvec == TRUE) return(apply(agree_mat, 1, sum) == M)
    if(retvec == FALSE) return(agreed)
  }
}

narrative_sign_probs <- function(decomp_obj,
                                 start_date,
                                 freq,
                                 dates,
                                 signs) {

  N <- dim(decomp_obj$E)[1]
  m <- dim(decomp_obj$E)[3]
  probs <- matrix(NA, nrow = length(dates), ncol = m)
  total_probs <- rep(NA, m)
  for(i in 1:m) {
    for(j in 1:length(dates)) {
      picked <- shock_decomp_pick(decomp_obj, show_date = dates[[j]], start_date = start_date, freq = freq)
      if(signs[j] == 1) agree_dummy <- picked$E[,i] > 0
      if(signs[j] == -1) agree_dummy <- picked$E[,i] < 0
      if(j == 1) {
        agree_dummies <- agree_dummy
      } else {
        agree_dummies <- cbind(agree_dummies, agree_dummy)
      }
    }
    probs[,i] <- apply(agree_dummies, 2, mean)
    total_probs[i] <- mean(apply(agree_dummies, 1, function(x) ifelse(sum(x) == length(x), TRUE, FALSE)))
  }
  ret <- rbind(probs, total_probs)
  rownames(ret) <- c(1:length(dates), "Total")
  ret
}

###########################
### Marginal likelihood ###
###########################

#' Computes a numerical estimate for the log-marginal-likelihood (Documentation incomplete)
#'
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param model A list outputted by \code{init_rbsvar()}.
#' @param burn ...
#' @param M ...
#' @param J ...
#' @param rel_tune ...
#' @param parallel ...
#' @param parallel_likelihood ...
#' @return A list.
#'
#' @export
marginal_likelihood <- function(output,
                                model,
                                burn,
                                M = NA,
                                J = NULL,
                                rel_tune = 1,
                                parallel = FALSE,
                                parallel_likelihood = FALSE) {

  start_time <- Sys.time()
  post <- post_sample(output, burn, M)
  s <- post$s

  # Too high J seems to induce irregular problems with C stack limit...
  # J = 10 000 seems to provide reasonable accuracy at least in some cases...
  if(is.null(J)) J <- 10000

  tune <- (2.38 / sqrt(2 * ncol(s))) * rel_tune
  s_demeaned <- scale(s, scale = FALSE)
  sigma_star <- (crossprod(s_demeaned) / nrow(s_demeaned)) * tune^2

  # To be used when zero restrictions are allowed
  #if(FALSE) {
  #  nodev <- which(diag(sigma_star) == 0) - 1
  #} else {
  #  nodev <- c(-1)
  #}
  if(min(eigen(sigma_star)$val) < 0) diag(sigma_star) <- diag(sigma_star) + abs(min(eigen(sigma_star)$val)) + 0.0001

  theta_star <- apply(s, 2, mean)
  logden_star <- eval_rbsvar(model, par = theta_star)
  theta_maxden <- s[which(post$d == max(post$d))[1],]

  # If density of theta_star is higher than point in the sample, algorithm fails
  count <- 10
  while(logden_star > eval_rbsvar(model, par = theta_maxden)) {
    theta_star <- theta_star + rnorm(length(theta_star), 0, 1) * sqrt(diag(sigma_star))
    logden_star <- eval_rbsvar(model, par = theta_star)
    count <- count - 1
    if(count < 0) stop("Something went wrong. (while-loop never ended)")
  }

  # If density is not well defined at theta_star, algorithm fails
  count <- 10
  while(logden_star == -Inf) {
    theta_star <- theta_star + (theta_maxden - theta_star) / 2
    logden_star <- eval_rbsvar(model, par = theta_star)
    count <- count - 1
    if(count < 0) stop("Something went wrong. (while-loop never ended)")
  }

  # Compute proposal densities (R-implementation is more robust than RcppDist implementation)
  proposal_densities <- mvtnorm::dmvt(s, delta = theta_star, sigma = sigma_star, df = 1, log = TRUE)

  model_R <- build_model_R(model)
  ret <- log_ml_cpp(proposal_densities = proposal_densities, posterior_densities = post$d,
                    theta_star = theta_star, sigma_star = sigma_star, logden_star = logden_star, J = J,
                    parallel = parallel, parallel_likelihood = parallel_likelihood,
                    model_R = model_R)

  names(ret) <- c("log_marginal_likelihood", "log_mean_posterior_ordinate",
                  "numerator_log_vec", "denumerator_log_vec")
  ret$time <- Sys.time() - start_time
  ret
}

###################
### Forecasting ###
###################

forecast_rbsvar <- function(model,
                            output,
                            burn = 0,
                            horizon = 1,
                            N = 1000,
                            cumulate = c(),
                            plots = TRUE,
                            verbose = TRUE) {

 # Later...

}

##############################################
### Sign restriction probability algorithm ###
##############################################

#' Computes sign restriction probabilities a'la Lanne & Luoto (2020) (Documentation incomplete)
#'
#' @param model A list outputted by \code{init_rbsvar()}.
#' @param output A list outputted by \code{est_rbsvar()}.
#' @param burn ...
#' @param signs ...
#' @param horizons ...
#' @param limit ...
#' @param N ...
#' @param cumulate ...
#' @param verbose ...
#' @return A list.
#'
#' @export
sign_probabilities <- function(model,
                               output,
                               burn,
                               signs,
                               horizons,
                               limit = 3.2,
                               N = 10000,
                               cumulate = c(),
                               verbose = TRUE) {

  if(model$prior$A$mean[1] == -1 | model$prior$B$mean[1] == -1) stop("Prior is improper. Proper prior needed.")
  if(model$type != "svar") stop("model$type == 'svar' required.")

  if(is.null(ncol(signs))) signs <- matrix(signs, ncol = 1)
  if(!is.list(horizons)) horizons <- list(horizons)
  if(length(horizons) != ncol(signs)) stop("'horizons' misspecified.")
  if(nrow(signs) != ncol(model$y)) stop("'signs' misspecified.")

  # Function specific helper functions
  simulate_prior <- function(model, output, N = NA) {
    if(is.null(model$prior$A) | is.null(model$prior$B)) stop("Proper prior needed.")
    if(is.na(N)) N <- output$args$N * output$args$n
    A_prior_cov <- model$prior$A$cov
    if(ncol(A_prior_cov) == 1) A_prior_cov <- diag(c(A_prior_cov))
    A_sample <- mvtnorm::rmvnorm(n = N, mean = model$prior$A$mean, sigma = A_prior_cov)
    B_prior_cov <- model$prior$B$cov
    if(ncol(B_prior_cov) == 1) B_prior_cov <- diag(c(B_prior_cov))
    B_sample <- mvtnorm::rmvnorm(n = N, mean = model$prior$B$mean, sigma = B_prior_cov)
    prior_sample <- cbind(A_sample, B_sample)
    colnames(prior_sample) <- colnames(output$chains)[1:ncol(prior_sample)]
    prior_sample
  }
  check_signs <- function(x, signs) {
    if(sum(sign(x)[signs != 0] != signs[signs != 0]) == 0) {
      1
    } else {
      0
    }
  }

  # Compute posterior impulse responses
  if(verbose) cat("Posterior: \n")
  irf_obj <- irf(model = model,
                 output = output,
                 burn = burn,
                 N = N,
                 horizon = max(unlist(horizons)),
                 cumulate = cumulate,
                 verbose = verbose)
  if(verbose) cat("------ \n")

  # Simulate from the prior and compute prior impulse responses
  if(verbose) cat("Prior: \n")
  output_prior <- list()
  output_prior$chains <- simulate_prior(model, output)
  output_prior$args$m0 <- output$args$m0
  output_prior$args$n <- output$args$n
  irf_obj_prior <- irf(model = model,
                       output = output_prior,
                       burn = 0,
                       N = N,
                       horizon = max(unlist(horizons)),
                       cumulate = cumulate,
                       verbose = verbose)
  if(verbose) cat("------ \n")

  shock_permutations <- gtools::permutations(nrow(signs), ncol(signs))
  sign_permutations <- gtools::permutations(2, ncol(signs), v = c(-1, 1), repeats.allowed = TRUE)

  # Go through different permutations
  if(verbose) cat(paste0("Going through: ", nrow(shock_permutations), " x ", nrow(sign_permutations),
                         " = ", nrow(shock_permutations) * nrow(sign_permutations)," permutations... \n"))
  significant_permutations <- list()
  count <- 0
  for(shock_perm_index in 1:nrow(shock_permutations)) {
    shock_perm_now <- shock_permutations[shock_perm_index,]
    for(sign_perm_index in 1:nrow(sign_permutations)) {
      signs_perm_now <- sign_permutations[sign_perm_index,]

      prob_mat_total <- matrix(NA, ncol = length(shock_perm_now), nrow = N)
      prob_mat_total_prior <- matrix(NA, ncol = length(shock_perm_now), nrow = N)
      for(shock_index in 1:length(shock_perm_now)) {

        shock_now <- shock_perm_now[shock_index]
        horizons_now <- horizons[[shock_index]]
        signs_now <- signs[,shock_index] * signs_perm_now[shock_index]

        # Posterior
        irf_array <- irf_obj[[shock_now]]
        prob_mat <- matrix(NA, ncol = length(horizons_now), nrow = N)
        for(horizon_index in 1:ncol(prob_mat)) {
          this_horizon <- horizons_now[horizon_index] + 1
          irf_mat <- irf_array[,this_horizon,]
          prob_mat[,horizon_index] <- apply(irf_mat, 2, check_signs, signs = signs_now)
        }
        prob_vec <- apply(prob_mat, 1, sum)
        prob_vec <- ifelse(prob_vec == ncol(prob_mat), 1, 0)
        prob_mat_total[,shock_index] <- prob_vec

        # Prior
        irf_array_prior <- irf_obj_prior[[shock_now]]
        prob_mat_prior <- matrix(NA, ncol = length(horizons_now), nrow = N)
        for(horizon_index in 1:ncol(prob_mat_prior)) {
          this_horizon <- horizons_now[horizon_index] + 1
          irf_mat_prior <- irf_array_prior[,this_horizon,]
          prob_mat_prior[,horizon_index] <- apply(irf_mat_prior, 2, check_signs, signs = signs_now)
        }
        prob_vec_prior <- apply(prob_mat_prior, 1, sum)
        prob_vec_prior <- ifelse(prob_vec_prior == ncol(prob_mat_prior), 1, 0)
        prob_mat_total_prior[,shock_index] <- prob_vec_prior
      }

      prob_vec_total <- apply(prob_mat_total, 1, sum)
      prob_vec_total <- ifelse(prob_vec_total == ncol(prob_mat_total), 1, 0)
      prob_vec_total_prior <- apply(prob_mat_total_prior, 1, sum)
      prob_vec_total_prior <- ifelse(prob_vec_total_prior == ncol(prob_mat_total_prior), 1, 0)

      # Bayes factor
      bayes_factor <- mean(prob_vec_total) / mean(prob_vec_total_prior)
      if(bayes_factor > limit) {
        count <- count + 1
        significant_permutations[[count]] <- list("bayes_factor" = bayes_factor,
                                                  "shocks" = shock_perm_now,
                                                  "sign_flips" = signs_perm_now,
                                                  "prob" = mean(prob_vec_total),
                                                  "prob_prior" = mean(prob_vec_total_prior))
      }
    }
  }

  if(verbose) cat("DONE. \n")
  if(verbose) cat("------ \n")

  if(length(significant_permutations) == 0) {
    if(verbose) cat(paste0("Sign restrictions not supported by the data, i.e. Bayes factor < ",
                           limit, " for any permutation of the shocks. \n"))
    return(NULL)
  } else {
    if(length(significant_permutations) == 1) {
      if(verbose) cat(paste0("One permutation of the shocks supported by the the data, i.e. Bayes factor > ",
                             limit, " for one permutation. \n"))
    } else {
      if(verbose) cat(paste0(length(significant_permutations), " permutations of the shocks potentially supported by the data, i.e. Bayes factor < ",
                             limit, " for ", length(significant_permutations),  " permutations. \n"))
    }
  }

  # Matrix of Bayes factors
  bayes_factors <- diag(length(significant_permutations))
  for(i in 1:length(significant_permutations)) {
    for(j in 1:length(significant_permutations)) {
      if(i == j) bayes_factors[i,j] <- significant_permutations[[i]]$bayes_factor
      if(i != j) bayes_factors[i,j] <- significant_permutations[[i]]$prob / significant_permutations[[j]]$prob
    }
  }

  ret <- list("permutations" = significant_permutations,
              "bayes_factors" = bayes_factors,
              "args" = list("signs" = signs, "horizons" = horizons, "limit" = limit))
  ret
}



