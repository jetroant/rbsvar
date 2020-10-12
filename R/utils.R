
#Builds xx and yy from y and p (VAR) NO EXPORT
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
