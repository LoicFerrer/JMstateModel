H.longSplinePH <-
  function(betas) # modified
  {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
      Ys <- as.vector(Xs %*% betas) + Zsb
      Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
      eta.s <- {
        if (is.null(transform.value))
          Ws.intF.vl.alph * Ys
        else Ws.intF.vl.alph * transform.value(Ys)
      }
    }
    if (parameterization %in% c("slope", "both")) {
      Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
      Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
      eta.s <- if (parameterization == "both") 
        eta.s + Ws.intF.sl.alph * Ys.deriv
      else Ws.intF.sl.alph * Ys.deriv
    }
    exp.eta.tw.P <- exp(eta.tw1) * P
    H1 <- XtX/sigma^2
    Int <- wk * exp(eta.ws + eta.s)
    if (!is.null(transform.value))
      Y <- as.vector(Xtime %*% betas) + Ztime.b
    H2 <- matrix(0, ncx, ncx)
    for (i in 1:ncx) {
      for (j in i:ncx) {
        XX <- if (parameterization == "value") {
          if (is.null(transform.value))
            Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
          else {
            to.add <- d * c(WintF.vl %*% alpha) * H.transform.value(Xtime[, i], Xtime[, j], Y)
            Ws.intF.vl.alph * H.transform.value(Xs[, i], Xs[, j], Ys) +
              Ws.intF.vl.alph^2 * gr.transform.value(Xs[, i], Ys) * gr.transform.value(Xs[, j], Ys)
          }
        }
        else if (parameterization == "slope") {
          if (i %in% indFixed && j %in% indFixed) {
            ii <- match(i, indFixed)
            jj <- match(j, indFixed)
            Ws.intF.sl.alph^2 * Xs.deriv[, ii] * Xs.deriv[, jj]
          }
          else 0
        }
        else {
          if (i %in% indFixed && j %in% indFixed) {
            ii <- match(i, indFixed)
            jj <- match(j, indFixed)
            if (is.null(transform.value))
              (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * Xs.deriv[, ii]) *
              (Ws.intF.vl.alph * Xs[, j] + Ws.intF.sl.alph * Xs.deriv[, jj])
            else {
              to.add <- d * c(WintF.vl %*% alpha) * H.transform.value(Xtime[, i], Xtime[, j], Y)
              Ws.intF.vl.alph * H.transform.value(Xs[, i], Xs[, j], Ys) +
                (Ws.intF.vl.alph * gr.transform.value(Xs[, i], Ys) + Ws.intF.sl.alph * Xs.deriv[, ii]) *
                (Ws.intF.vl.alph * gr.transform.value(Xs[, j], Ys) + Ws.intF.sl.alph * Xs.deriv[, jj])
            }
          }
          else if (i %in% indFixed && !j %in% indFixed) {
            ii <- match(i, indFixed)
            if (is.null(transform.value))
              (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * Xs.deriv[, ii]) * (Ws.intF.vl.alph * Xs[, j])
            else {
              to.add <- d * c(WintF.vl %*% alpha) * H.transform.value(Xtime[, i], Xtime[, j], Y)
              Ws.intF.vl.alph * H.transform.value(Xs[, i], Xs[, j], Ys) +
                (Ws.intF.vl.alph * gr.transform.value(Xs[, i], Ys) + Ws.intF.sl.alph * Xs.deriv[, ii]) *
                (Ws.intF.vl.alph * gr.transform.value(Xs[, j], Ys))
            } 
          }
          else if (!i %in% indFixed && j %in% indFixed) {
            jj <- match(j, indFixed)
            if (is.null(transform.value))
              (Ws.intF.vl.alph * Xs[, i]) * (Ws.intF.vl.alph * Xs[, j] + Ws.intF.sl.alph * Xs.deriv[, jj])
            else {
              to.add <- d * c(WintF.vl %*% alpha) * H.transform.value(Xtime[, i], Xtime[, j], Y)
              Ws.intF.vl.alph * H.transform.value(Xs[, i], Xs[, j], Ys) +
                (Ws.intF.vl.alph * gr.transform.value(Xs[, i], Ys)) *
                (Ws.intF.vl.alph * gr.transform.value(Xs[, j], Ys) + Ws.intF.sl.alph * Xs.deriv[, jj])
            }
          }
          else {
            if (is.null(transform.value))
              Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
            else {
              to.add <- d * c(WintF.vl %*% alpha) * H.transform.value(Xtime[, i], Xtime[, j], Y)
              Ws.intF.vl.alph * H.transform.value(Xs[, i], Xs[, j], Ys) +
                Ws.intF.vl.alph^2 * gr.transform.value(Xs[, i], Ys) * gr.transform.value(Xs[, j], Ys)
            } 
          }
        }
        ki <- exp.eta.tw.P * rowsum(Int * XX, id.GK, reorder = FALSE)
        ki <- {
          if (is.null(transform.value))
            rowsum(ki, idT, reorder = FALSE)
          else rowsum(ki - to.add, idT, reorder = FALSE)
        } 
        kii <- c((p.byt * ki) %*% wGH)
        H2[i, j] <- sum(kii, na.rm = TRUE)
      }
    }
    H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
    H1 + H2
  }
