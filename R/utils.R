
#Builds xx and yy from y and p (lags) NO EXPORT
build_xy <- function(y, p, constant = TRUE) {
  if(p == 0) {
    yy <- y
    xx <- NULL
    if(constant == TRUE) xx <- matrix(rep(1, nrow(yy)), ncol = 1)
  } else {
    n <- ncol(y)
    for(i in 1:p) {
      if(i == 1) {
        temp <- rbind(rep(NA, n), y[-(nrow(y)),])
        xx <- temp
      } else {
        temp <- rbind(rep(NA, n), temp[-(nrow(y)),])
        xx <- cbind(xx, temp)
      }
    }
    if(constant == TRUE) xx <- cbind(rep(1, nrow(xx)), xx)
    xx <- xx[-c(1:p),]
    yy <- y[-c(1:p),]
  }
  ret <- list("xx" = xx, "yy" = yy)
  ret
}

#Builds xx and yy dummy prior parts from y and p (lags) NO EXPORT
build_xy0 <- function(y, p, shrinkage, minnesota_means, constant = TRUE) {
  n <- ncol(y)
  if(is.null(minnesota_means)) minnesota_means <- rep(1, n)
  if(length(minnesota_means) == 1) minnesota_means <- rep(minnesota_means, n)
  if(length(minnesota_means) != n) stop("'minnesota_means' not of length 'ncol(y)'.")
  arsigmas <- rep(NA, n)
  for(i in 1:n) arsigmas[i] <- sqrt(ar(y[,i], aic = F, order.max = p)$var.pred)
  if(constant) {
    yy0 <- rbind(diag((minnesota_means * arsigmas) / shrinkage),
                 matrix(0, ncol = n, nrow = n*p - n + 1))
    xx0 <- matrix(0, ncol = nrow(yy0) - 1, nrow = nrow(yy0) - 1)
    for(i in 1:p) {
      indices <- (n*i-n+1):(n*i)
      diag(xx0)[indices] <- i * (arsigmas / shrinkage)
    }
    xx0 <- cbind(0, xx0)
    xx0 <- rbind(xx0,
                 c(0.00001, rep(0, ncol(xx0)-1)))
  } else {

    yy0 <- rbind(diag((minnesota_means * arsigmas) / shrinkage),
                 matrix(0, ncol = n, nrow = n*p - n))
    xx0 <- matrix(0, ncol = nrow(yy0), nrow = nrow(yy0))
    for(i in 1:p) {
      indices <- (n*i-n+1):(n*i)
      diag(xx0)[indices] <- i * (arsigmas / shrinkage)
    }
  }
  ret <- list("xx0" = xx0, "yy0" = yy0, "arsigmas" = arsigmas)
  ret
}

#Stacks VAR(p) coefficient matrix to VAR(1) coefficient matrix
#(companion form, used by irf()) NO EXPORT
stackA <- function(A, constant = TRUE) {
  if(constant == TRUE) A <- t(A)[,-1]
  if(constant == FALSE) A <- t(A)
  m <- nrow(A)
  lags <- ncol(A)/m
  eye <- diag(m*lags-m)
  A <- rbind(A, cbind(eye, matrix(0, ncol = m, nrow= nrow(eye))))
  A
}
