
#NO EXPORT
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

  # Collect inputs
  if(model$type != "svar") stop("model$type must be 'svar'")
  m <- ncol(model$y)
  p <- model$lags
  n <- output$args$n
  m0 <- output$args$m0
  ret <- list()
  if(is.null(shock_sizes)) shock_sizes <- rep(1, m)
  if(is.null(shocks)) shocks <- 1:m
  b_indices <- (model$cpp_args$first_b + 1):model$cpp_args$first_sgt
  a_indices <- 1:model$cpp_args$first_b

  # Posterior sample
  #s <- output$chains[-c(1:(m0 + burn)),]
  #s <- s[sample.int(nrow(s), N, replace = T),]

  #B_inverse <- model$cpp_args$B_inverse
  #if(length(cumulate) == 0) cumulate <- c(-1)
  #ret <- irf_cpp(s, horizon = 3, cumulate,
  #               shock_sizes, shocks,
  #               model$cpp_args$A_rows, model$cpp_args$first_b, model$cpp_args$first_sgt, m,
  #               B_inverse, parallel)

  #Compute and collect B_inv and A sample
  burn <- burn + 1
  b_inv_list <- list()
  a_list <- list()
  for(i in 1:n) {

    #Pick one chain
    row_indices <- seq(from = m0 + i - n, by = n, length.out = output$args$N + 1)
    one_chain <- output$chains[row_indices,][-c(1:burn),]

    #Compute B_inv
    b_chain <- one_chain[,b_indices]
    b_inv_chain <- matrix(NA, ncol = m^2, nrow = nrow(b_chain))
    for(j in 1:nrow(b_chain)) {
      B <- diag(m)
      B[,] <- b_chain[j,]
      b_inv_chain[j,] <- c(solve(B))
    }
    b_inv_list[[i]] <- b_inv_chain

    #Pick A
    a_list[[i]] <- one_chain[,a_indices]
  }
  for(i in 1:n) {
    if(i == 1) {
      full_binv <- b_inv_list[[i]]
      full_a <- a_list[[i]]
    } else {
      full_binv <- rbind(full_binv, b_inv_list[[i]])
      full_a <- rbind(full_a, a_list[[i]])
    }
  }

  rows <- sample.int(nrow(full_a), N, replace = TRUE)
  full_a <- full_a[rows,]
  full_binv <- full_binv[rows,]

  for(shock_index in 1:m) {

    if(!(shock_index %in% shocks)) {
      ret[[shock_index]] <- NULL
      next
    }

    e <- rep(0, m)
    e[shock_index] <- shock_sizes[shock_index]
    irfs <- array(NA, dim = c(m, horizon+1, nrow(full_a)))

    if(verbose == TRUE) print(paste0("Computing impulse responses... (", which(shocks == shock_index), "/", length(shocks), ")"))
    if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = nrow(full_a), style = 3)
    for(row_index in 1:nrow(full_a)) {

      B_inv <- diag(m)
      B_inv[,] <- full_binv[row_index,]
      A <- matrix(full_a[row_index,], ncol = m)
      AA <- stackA(A)

      for(h in 1:(horizon+1)) {
        if(h == 1) {
          zero <- B_inv %*% e
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

irf_plot <- function(irf_obj,
                     varnames = NULL,
                     probs = NULL,
                     mar = c(2,2,2,1),
                     mfrow = NULL,
                     color = "tomato",
                     leg = TRUE,
                     normalize = NULL) {

  if(requireNamespace("fanplot", quietly = TRUE)) {
    #Do nothing
  } else {
    stop("Package 'fanplot' needed for plotting the impulse response functions. \n 'fanplot' was not found. The package can be installed by 'install.packages('fanplot')'. \n")
  }

  m <- length(irf_obj)
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

    if(is.null(irf_obj[[row]])) next

    sub_irfs <- t(irf_obj[[row]][col,,])
    if(!is.null(normalize)) sub_irfs <- sub_irfs * normalize[col]
    mean_sub_irfs <- ts(apply(sub_irfs, 2, mean), start = 0)

    if(is.null(probs)) {
      p <- c(0.0249, 0.025, seq(0.1, 0.9, 0.1), 0.975, 0.9751)
    } else {
      p <- probs
    }
    quant <- function(column) quantile(column, probs = p)
    quantiles_sub_irfs <- apply(sub_irfs, 2, quant)

    plot(mean_sub_irfs, lwd = 2, lty = 2, col = color, ylab = "", xlab = "",
         main = paste0("Shock ", row, " on ", varnames[col]),
         ylim = c(min(quantiles_sub_irfs), max(quantiles_sub_irfs)))
    grid()
    fanplot::fan(data = quantiles_sub_irfs, data.type = "values", probs = p,
                 start = 0, fan.col = colorRampPalette(c(color, "white")),
                 rlab = NULL, ln = NULL)
    abline(h = 0, lwd = 2, lty = 2)

    post_mass <- (max(p[-length(p)]) - min(p[-1]))*100
    if(col == 1 & row == 1 & leg == TRUE) legend("topright", c(paste0(post_mass,"% of post. prob. mass")), lwd = 0, bty = "n", col = "tomato")
  }
  par(mfrow = c(1,1))
}

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





