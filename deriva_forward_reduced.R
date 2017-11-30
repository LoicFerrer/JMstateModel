deriva_forward_reduced <-
  function (b, funcpa) 
  {
    m <- length(b)
    bh2 <- bh <- rep(0, m)
    v <- rep(0, m)
    fcith <- rep(0, m)
    rl <- funcpa(b)
    for (i in 1:m) {
      bh <- b
      th <- max(1e-06, (1e-04 * abs(b[i])))
      bh[i] <- bh[i] + th
      fcith[i] <- funcpa(bh)
      cat(paste(c("Grad", ":", i, "/", m, sep = "")), "\n")
    }
    th1 <- max(1e-06, (1e-04 * abs(b[1])))
    for (i in 1:m) {
      bh <- b
      thi <- max(1e-06, (1e-04 * abs(b[i])))
      th <- thi * th1
      bh[i] <- bh[i] + thi
      bh[1] <- bh[1] + th1
      temp <- funcpa(bh)
      v[i] <- -(temp - (fcith[1]) - (fcith[i]) + rl)/th
      cat(paste(c("Hess", ":", i, "/", m, sep = "")), "\n")
    }
    result <- list(v = v, rl = rl)
    return(result)
  }