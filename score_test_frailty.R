score_test_frailty <- 
  function(object)
  {
    if (object$method != "spline-PH-GH")
      stop("Joint multi-state model is only implemented with 'method = spline-PH-GH'")
    if(options()$digits < 7)
      stop("You have to improve the precision of your machine (with options(digits = k)) or reduce the precision in 'deriva_forward()'")

    transform.value <- object$transform.value
    
    #### Longitudinal sub-part ####
    method <- object$method
    parameterization <- object$parameterization
    logT <- object$y$logT
    id.GK <- rep(seq_along(logT), each = object$control$GKk)
    eta.yx <- as.vector(object$x$X %*% object$coefficients$betas)
    GH <- JM:::gauher(object$control$GHk)
    ncz <- ncol(object$x$Z)
    b <- as.matrix(expand.grid(rep(list(GH$x), ncz)))
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))
    wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    if (object$control$typeGH == "simple") {
      b <- sqrt(2) * t(object$control$inv.chol.VC %*% t(b))
      wGH <- wGH * object$control$det.inv.chol.VC
    }
    else {
      b <- sqrt(2) * b
      VCdets <- object$control$det.inv.chol.VCs
    }
    dimnames(b) <- NULL
    Ztb <- object$x$Z %*% t(b)
    if (parameterization %in% c("value", "both")) {
      Ztime.b <- object$x$Ztime %*% t(b)
      Zsb <- object$x$Zs %*% t(b)
    }
    if (parameterization %in% c("slope", "both")) {
      if (length(object$derivForm$indRandom) > 1 || object$derivForm$indRandom) {
        Ztime.b.deriv <- object$x$Ztime.deriv %*% t(b[, object$derivForm$indRandom, drop = FALSE])
        Zsb.deriv <- object$x$Zs.deriv %*% t(b[, object$derivForm$indRandom, drop = FALSE])
      }
      else {
        Ztime.b.deriv <- matrix(0, nrow(object$x$Ztime.deriv), k)
        Zsb.deriv <- matrix(0, nrow(object$x$Zs.deriv), k)
      }
    }
    if (object$control$typeGH != "simple") {
      lis.b <- vector("list", object$n)
      for (i in 1:object$n) {
        lis.b[[i]] <- t(object$control$inv.chol.VCs[[i]] %*% t(b)) + 
          rep(object$control$ranef[i, ], each = k)
        Ztb[object$id == i, ] <- object$x$Z[object$id == i, , drop = FALSE] %*% 
          t(lis.b[[i]])
      }
      lis.b2 <- lapply(lis.b, function(b) if (ncz == 1) 
        b * b
        else t(apply(b, 1, function(x) x %o% x)))
      for (i in seq_along(logT)) {
        if (parameterization %in% c("value", "both")) {
          bb <- t(lis.b[[object$x$idT[i]]])
          Ztime.b[i, ] <- object$x$Ztime[i, , drop = FALSE] %*% bb
          Zsb[id.GK == i, ] <- object$x$Zs[id.GK == i, ] %*% bb
        }
        if (parameterization %in% c("slope", "both") && 
            (length(object$derivForm$indRandom) > 1 || object$derivForm$indRandom)) {
          bb <- t(lis.b[[object$x$idT[i]]][, object$derivForm$indRandom, drop = FALSE])
          Ztime.b.deriv[i, ] <- object$x$Ztime.deriv[i, , drop = FALSE] %*% bb
          Zsb.deriv[id.GK == i, ] <- object$x$Zs.deriv[id.GK == i, ] %*% bb
        }
      }
    }
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(object$y$y, mu.y, object$coefficients$sigma, TRUE)
    log.p.yb <- rowsum(logNorm, object$id, reorder = FALSE)
    dimnames(log.p.yb) <- NULL
    
    #### Survival sub-part ####
    eta.tw1 <- if (!is.null(object$x$W)) 
      as.vector(object$x$W %*% object$coefficients$gammas)
    else 0
    eta.tw2 <- as.vector(object$x$W2 %*% object$coefficients$gammas.bs)
    if (parameterization %in% c("value", "both")) {
      Y <- as.vector(object$x$Xtime %*% object$coefficients$betas) + Ztime.b
      Ys <- as.vector(object$x$Xs %*% object$coefficients$betas) + Zsb
      eta.t <- {
        if (is.null(transform.value))
          eta.tw2 + eta.tw1 + c(object$x$WintF.vl %*% object$coefficients$alpha) * Y
        else eta.tw2 + eta.tw1 + c(object$x$WintF.vl %*% object$coefficients$alpha) * transform.value(Y)
      }
      eta.s <- {
        if (is.null(transform.value))
          c(object$x$Ws.intF.vl %*% object$coefficients$alpha) * Ys
        else c(object$x$Ws.intF.vl %*% object$coefficients$alpha) * transform.value(Ys)
      }
    }
    if (parameterization %in% c("slope", "both")) {
      Y.deriv <- as.vector(object$x$Xtime.deriv %*% object$coefficients$betas[object$derivForm$indFixed]) + 
        Ztime.b.deriv
      Ys.deriv <- as.vector(object$x$Xs.deriv %*% object$coefficients$betas[object$derivForm$indFixed]) + 
        Zsb.deriv
      eta.t <- if (parameterization == "both") 
        eta.t + c(object$x$WintF.sl %*% object$coefficients$Dalpha) * Y.deriv
      else eta.tw2 + eta.tw1 + c(object$x$WintF.sl %*% object$coefficients$Dalpha) * Y.deriv
      eta.s <- if (parameterization == "both") 
        eta.s + c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv
      else c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv
    }
    eta.ws <- as.vector(object$x$W2s %*% object$coefficients$gammas.bs)
    
    #### Cumulative intensities ####
    log.hazard <- eta.t
    log.survival <- -exp(eta.tw1) * object$x$P * rowsum(object$x$wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- rowsum(object$y$d * log.hazard + log.survival, object$x$idT, reorder = FALSE)
    
    #### Random effects ####
    log.p.b <- if (object$control$typeGH == "simple") {
      rep(JM:::dmvnorm(b, rep(0, ncz), object$coefficients$D, TRUE), each = object$n)
    }
    else {
      matrix(JM:::dmvnorm(do.call(rbind, lis.b), rep(0, ncz), 
                          object$coefficients$D, TRUE), object$n, k, byrow = TRUE)
    }
    p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
    if (object$control$typeGH != "simple") 
      p.ytb <- p.ytb * VCdets
    p.yt <- c(p.ytb %*% wGH)  
    
    # Likelihood function under H1
    func_ll_H1 <- function(par){
      
      sigma2_v <- par[1]
      betas <- par[(1 + 1) : (1 + length(object$coefficients$betas))]
      sigma <- par[(1 + length(c(sigma2_v, betas))) : length(c(sigma2_v, betas, sigma))]
      gammas <- if (!is.null(object$x$W)) 
        par[(1 + length(c(sigma2_v, betas, sigma))) : length(c(sigma2_v, betas, sigma, object$coefficients$gammas))]
      else NULL
      alpha <- if (parameterization %in% c("value", "both"))
        par[(1 + length(c(sigma2_v, betas, sigma, gammas))) :
              length(c(sigma2_v, betas, sigma, gammas, object$coefficients$alpha))]
      else NULL
      Dalpha <- if (parameterization %in% c("slope", "both"))
        par[(1 + length(c(sigma2_v, betas, sigma, gammas, alpha))) :
              length(c(sigma2_v, betas, sigma, gammas, alpha, object$coefficients$Dalpha))]
      else NULL
      gammas.bs <- par[(1 + length(c(sigma2_v, betas, sigma, gammas, alpha, Dalpha))) :
                         length(c(sigma2_v, betas, sigma, gammas, alpha, Dalpha, object$coefficients$gammas.bs))]
      B <- par[(1 + length(c(sigma2_v, betas, sigma, gammas, alpha, Dalpha, gammas.bs))) :
                 npar]
      
      # Extension of corpcor::rebuild.cov to a vecteur of correlation parameters as argument
      rebuild.cov.vect <- function (r, v) {
        if (any(v < 0)) 
          stop("Negative variance encountered!")
        sd <- sqrt(v)
        r.mat <- matrix(1 , ncz, ncz)
        r.mat[upper.tri(r.mat)] <- r
        r.mat <- t(r.mat)
        r.mat[upper.tri(r.mat)] <- r
        m <- sweep(sweep(r.mat, 1, sd, "*"), 2, sd, "*")
        return(m)
      }
      
      B <- rebuild.cov.vect(B[(ncz+1):(ncz*(ncz+1)/2)], B[seq_len(ncz)]^2)
      
      GHk_score <- 50
      GH_score <- JM:::gauher(GHk_score)
      u <- GH_score$x 
      wGH_u <- GH_score$w
      wGH_u <- 1/sqrt(pi) * wGH_u
      u <- sqrt(2) * u
      
      ## Survival sub-part ##
      eta.tw1 <- if (!is.null(object$x$W)) 
        as.vector(object$x$W %*% gammas)
      else 0
      eta.tw2 <- as.vector(object$x$W2 %*% gammas.bs)
      if (parameterization %in% c("value", "both")) {
        Y <- as.vector(object$x$Xtime %*% betas) + 
          Ztime.b
        Ys <- as.vector(object$x$Xs %*% betas) + Zsb
        eta.t <- {
          if (is.null(transform.value))
            eta.tw2 + eta.tw1 + c(object$x$WintF.vl %*% alpha) * Y
          else eta.tw2 + eta.tw1 + c(object$x$WintF.vl %*% alpha) * transform.value(Y)
        }
        eta.s <-  {
          if (is.null(transform.value))
            c(object$x$Ws.intF.vl %*% alpha) * Ys
          else c(object$x$Ws.intF.vl %*% alpha) * transform.value(Ys)
        }
      }
      if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(object$x$Xtime.deriv %*% betas[object$derivForm$indFixed]) + 
          Ztime.b.deriv
        Ys.deriv <- as.vector(object$x$Xs.deriv %*% betas[object$derivForm$indFixed]) + 
          Zsb.deriv
        eta.t <- if (parameterization == "both") 
          eta.t + c(object$x$WintF.sl %*% Dalpha) * Y.deriv
        else eta.tw2 + eta.tw1 + c(object$x$WintF.sl %*% Dalpha) * Y.deriv
        eta.s <- if (parameterization == "both") 
          eta.s + c(object$x$Ws.intF.sl %*% Dalpha) * Ys.deriv
        else c(object$x$Ws.intF.sl %*% Dalpha) * Ys.deriv
      }
      eta.ws <- as.vector(object$x$W2s %*% gammas.bs)
      
      log.hazard.u <- apply(eta.t, 2, function(x) rep(x, each = GHk_score)) + rep(sqrt(sigma2_v) * u, length(logT))
      log.survival.u <- -exp(rep(eta.tw1, each = GHk_score)) *
        rep(exp(sqrt(sigma2_v) * u), length(logT)) *
        rep(object$x$P, each = GHk_score) * 
        apply(rowsum(object$x$wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE), 2, function(x) rep(x, each = GHk_score))
      dimnames(log.survival.u) <- NULL
      
      id.GHu <- c(apply(matrix(c(GHk_score * (object$x$idT - 1) + 1, 
                                 GHk_score * (object$x$idT - 1) + 1 + GHk_score - 1), ncol = 2), 1,
                        function(x) seq(from = x[1], to = x[2])))
      log.p.tbu <- rowsum(rep(object$y$d, each = GHk_score) * log.hazard.u + log.survival.u,
                          id.GHu, reorder = FALSE)
      p.tbu <- exp(log.p.tbu)
      
      p.tb <- matrix( , nrow = object$n, ncol = object$control$GHk^ncz)
      for (i in seq_len(object$n)){
        p.tb[i, ] <- c(t(p.tbu[(GHk_score * (i - 1) + 1) : (GHk_score*i), ]) %*% wGH_u)
      }
      log.p.tb <- log(p.tb)
      
      ## Random effects sub-part ##      
      log.p.b <- if (object$control$typeGH == "simple") {
        rep(JM:::dmvnorm(b, rep(0, ncz), B, TRUE), each = object$n)
      }
      else {
        matrix(JM:::dmvnorm(do.call(rbind, lis.b), rep(0, ncz), 
                            B, TRUE), object$n, k, byrow = TRUE)
      }
      
      ## Longitudinal sub-part ##
      eta.yx <- as.vector(object$x$X %*% betas)
      mu.y <- eta.yx + Ztb
      logNorm <- dnorm(object$y$y, mu.y, sigma, TRUE)
      log.p.yb <- rowsum(logNorm, object$id, reorder = FALSE)
      dimnames(log.p.yb) <- NULL
      
      p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
      if (object$control$typeGH != "simple") 
        p.ytb <- p.ytb * VCdets
      p.yt <- c(p.ytb %*% wGH)
      
      return(sum(log(p.yt)))
    }
    
    sigma2_v <- 0
    betas <- object$coefficients$betas
    sigma <- object$coefficients$sigma
    gammas <- object$coefficients$gammas
    gammas.bs <- object$coefficients$gammas.bs
    alpha <- object$coefficients$alpha
    Dalpha <- object$coefficients$Dalpha
    B <- c(sqrt(diag(object$coefficients$D)),
           cov2cor(object$coefficients$D)[upper.tri(cov2cor(object$coefficients$D))])

    par <- c(sigma2_v, betas, sigma, gammas, alpha, Dalpha, gammas.bs, B)
    npar <- length(par)
    
    deriv_H1 <- deriva_forward_reduced(par, func_ll_H1)
    
    #### Test Statistic ####
    U_i <- 1/(2*p.yt) *
      c( (p.ytb * 
            (rowsum(object$y$d + log.survival, group = object$x$idT, reorder = FALSE)^2 +
               rowsum(log.survival, group = object$x$idT, reorder = FALSE))) %*% wGH)
    U <- sum(U_i)

    var_U <- c(deriv_H1$v[1] - deriv_H1$v[2:length(par)] %*% vcov(object) %*% deriv_H1$v[2:length(par)])
    
    pval <- pchisq(pmax(0,U)^2 / var_U, df = 1, lower.tail = F)/2
    
    list(U_i = U_i,
         U = U,
         var_U = var_U,
         stat = pmax(0,U)^2 / var_U,
         pval = pval,
         conv = object$convergence)
  }
