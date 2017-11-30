JMstateModel <- 
  function (lmeObject, survObject, timeVar, 
            parameterization = c("value", "slope", "both"), 
            method = c("weibull-PH-aGH", "weibull-PH-GH", "weibull-AFT-aGH", "weibull-AFT-GH", "piecewise-PH-aGH", 
                       "piecewise-PH-GH", "Cox-PH-aGH", "Cox-PH-GH", "spline-PH-aGH", "spline-PH-GH", "ch-Laplace"),
            interFact = NULL, derivForm = NULL, lag = 0, scaleWB = NULL, CompRisk = FALSE, init = NULL, control = list(),
            Mstate = FALSE, data.Mstate = NULL, ID.Mstate = NULL, init.type.ranef = "mean", b.postVar.opt = FALSE,
            transform.value = NULL, gr.transform.value = NULL, H.transform.value = NULL,
            # 'Mstate' is TRUE when a joint multi-state model is fitted
            # 'data.Mstate' is the database used in 'coxph'
            # 'ID.Mstate' is the column name of the id in 'data.Mstate' (class(ID.Mstate)=="character")
            # 'init" may be a 'jointModel' object
            # 'init.type.ranef' is the type of posterior random effects used in the computation of the integral over the random effects
            # 'b.postVar.opt' specifies if the hessian of the posterior individual random effects must be computed analytically ("FALSE") or not ("TRUE" --> optim)
            # 'transform.value' is the transformation function of the true current marker value which is considered in the dependence function
            # 'gr.transform.value' is the derivative of the transformation function --> function with, in order: arugments x (fixed effects) and y (true current level)
            # 'H.transform.value' is the second derivative of the transformation function --> function with, in order: arugments x (fixed effects) and y (true current level)
            ...) 
  {
    cl <- match.call()
    #### Errors message ####
    if (!inherits(lmeObject, "lme")) 
      stop("\n'lmeObject' must inherit from class lme.")
    if (length(lmeObject$group) > 1) 
      stop("\nnested random-effects are not allowed in lme().")
    if (!is.null(lmeObject$modelStruct$corStruct)) 
      warning("correlation structure in 'lmeObject' is ignored.\n")
    if (!is.null(lmeObject$modelStruct$varStruct)) 
      warning("variance structure in 'lmeObject' is ignored.\n")
    if (!inherits(survObject, "coxph") && !inherits(survObject, "survreg")) 
      stop("\n'survObject' must inherit from class coxph or class survreg.")
    if (!is.matrix(survObject$x)) 
      stop("\nuse argument 'x = TRUE' in ", 
           if (inherits(survObject, "coxph")) "'coxph()'."
           else "'survreg()'.")
    if (length(timeVar) != 1 || !is.character(timeVar)) 
      stop("\n'timeVar' must be a character string.")
    method. <- match.arg(method)
    method <- switch(method., `weibull-AFT-GH` = , `weibull-AFT-aGH` = "weibull-AFT-GH", 
                     `weibull-PH-GH` = , `weibull-PH-aGH` = "weibull-PH-GH", 
                     `piecewise-PH-GH` = , `piecewise-PH-aGH` = "piecewise-PH-GH", 
                     `Cox-PH-GH` = , `Cox-PH-aGH` = "Cox-PH-GH", `spline-PH-GH` = , 
                     `spline-PH-aGH` = "spline-PH-GH", `ch-Laplace` = "ch-Laplace")
    parameterization <- match.arg(parameterization)
    if (method == "Cox-PH-GH" && !inherits(survObject, "coxph")) 
      stop("\nfor 'method = Cox-PH-GH', 'survObject' must inherit from class coxph.")
    if (parameterization %in% c("slope", "both") && method %in% c("Cox-PH-GH", "ch-Laplace")) 
      stop("\nthe slope parameterization is not currently available for methods 'Cox-PH-GH' & 'ch-Laplace'.")
    if (parameterization %in% c("slope", "both") && is.null(derivForm)) {
      stop("\nwhen parameterization is 'slope' or 'both' you need to specify the 'derivForm' argument.")
    }
    if (parameterization %in% c("slope", "both") && !is.list(derivForm)) {
      stop("\nthe 'derivForm' argument must be a list with components 'fixed' (a formula),\n\t'indFixed'", 
           "(a numeric vector), 'random' (a formula) and 'indRandom' (a numeric vector).")
    }
    if (!is.null(interFact) && !is.list(interFact)) {
      stop("\nthe 'interFact' argument must be a list -- check the help file for more info.")
    }
    if (!is.null(interFact) && method %in% c("Cox-PH-GH", "ch-Laplace")) {
      stop("\nincluding interaction terms is not currently available for methods 'Cox-PH-GH' & 'ch-Laplace'.")
    }
    if ((CompRisk || Mstate) && (method != "spline-PH-GH" || is.null(survObject$strata))) {
      stop("\nto fit a competing risks or joint multi-state model you must choose as method 'spline-PH-GH'", 
           " and include a strata() in the specification of the coxph().")
    }
    if (Mstate && (!is.null(survObject$cluster))) {
      stop("\nto fit a multi-state model with coxph(), you can not use time-dependent covariate.",
           "Thus, do not use 'cluster()' in coxph().")
    }
    if (Mstate && (is.null(ID.Mstate) || !is.vector(data.Mstate[, ID.Mstate]) || is.null(data.Mstate))) {
      stop("\nwhen Mstate is TRUE, you must specify data.Mstate and ID.Mstate",
           " which respectively correspond to the databased used in coxph() and the column name of the id in data.Mstate.")
    }
    if (CompRisk && Mstate) {
      stop("\nyou must choose between a joint model with competing risks or a joint multi-state model.")
    }
    if (!is.null(transform.value) && !is.function(transform.value))
      stop("\nthe 'transform.value' argument must be a function.")
    if (!is.null(transform.value) && (method != "spline-PH-GH"))
      stop("\nto use the 'transform.value' argument, you must choose as method 'spline-PH-GH'.")
    if (!is.null(transform.value) && is.null(gr.transform.value))
      stop("\nwhen you use the 'transform.value' argument, you must define the derivative in 'gr.transform.value'.")
    if (!is.null(gr.transform.value) && !is.function(gr.transform.value))
      stop("\nthe 'gr.transform.value' argument must be a function.")
    if (!is.null(transform.value) && is.null(H.transform.value))
      stop("\nwhen you use the 'transform.value' argument, you must define the second derivative in 'H.transform.value'.")
    if (!is.null(H.transform.value) && !is.function(H.transform.value))
      stop("\nthe 'H.transform.value' argument must be a function.")
    if ((b.postVar.opt && class(init) != "jointModel") | (b.postVar.opt && init.type.ranef != "mode"))
      stop("\n'b.postVar.opt' can be TRUE only when 'class(init) == 'jointModel'' and 'init.type.ranef == 'mode''.")
    #### Event history sub-part ####
    formT <- formula(survObject)
    if (inherits(survObject, "coxph")) {
      W <- survObject$x
      keepW <- suppressWarnings(!is.na(survObject$coefficients))
      W <- W[, keepW, drop = FALSE]
      if (CompRisk || Mstate) {
        nRisks <- length(unique(survObject$strata)) 
      }
      else {
        nRisks <- 1
      }
      surv <- survObject$y 
      if (attr(surv, "type") == "right") {
        LongFormat <- FALSE
        Time <- survObject$y[, 1]
        d <- survObject$y[, 2]
      }      
      else if (attr(surv, "type") == "counting") {
        LongFormat <- TRUE
        if (is.null(survObject$model)) 
          stop("\nplease refit the Cox model including in the ", 
               "call to coxph() the argument 'model = TRUE'.")
        Time <- survObject$y[, 2]
        d <- survObject$y[, 3]
      }
      idT <- if (!is.null(survObject$model$cluster)) {
        as.vector(unclass(survObject$model$cluster))
      }
      else {
        if (!Mstate) {
          if(!CompRisk) {
            seq_along(Time)
          }
          else { rep(seq_len(length(Time)/nRisks), each = nRisks) }
        }
        else { data.Mstate[, ID.Mstate] }
      }
      idT <- match(idT, unique(idT))
    }
    else { 
      W <- survObject$x[, -1, drop = FALSE]
      Time <- exp(survObject$y[, 1])
      d <- survObject$y[, 2]
      idT <- seq_along(Time)
      LongFormat <- FALSE
      nRisks <- 1
    }
    nT <- length(unique(idT))
    if (LongFormat && is.null(survObject$model$cluster) && !Mstate )
      stop("\nuse argument 'model = TRUE' and cluster() in coxph().")
    if (!length(W)) 
      W <- NULL
    if (sum(d) < 5) 
      warning("\nmore than 5 events are required.")
    WintF.vl <- WintF.sl <- as.matrix(rep(1, length(Time)))
    if (!is.null(interFact)) {
      if (!is.null(interFact$value)) 
        WintF.vl <- if (is.null(survObject$model) || !is.null(interFact$data)) {
          model.matrix(interFact$value, data = interFact$data)
        }
      else {
        model.matrix(interFact$value, data = survObject$model)
      }
      if (!is.null(interFact$slope)) 
        WintF.sl <- if (is.null(survObject$model) || !is.null(interFact$data)) {
          model.matrix(interFact$slope, data = interFact$data)
        }
      else {
        model.matrix(interFact$slope, data = survObject$model)
      }
    }
    #### Longitudinal sub-part ####
    id <- as.vector(unclass(lmeObject$groups[[1]]))
    b <- if (class(init) == "jointModel") {
      if (init.type.ranef == "mean")
        data.matrix(ranef(init, type = "mean"))
      else if (init.type.ranef == "mode") {
        if (b.postVar.opt == FALSE)
          data.matrix(modified.ranef.jointModel(init, type = "mode"))
        else {
          b <- modified.ranef.jointModel(init, type = "mode", postVar = TRUE)
          b.postVar <- attr(b, "postVar")
          data.matrix(b)
        }
      }
      else stop("\n'init.type.ranef' must be 'mean' or 'mode'.")
    } 
    else data.matrix(ranef(lmeObject))
    dimnames(b) <- NULL
    nY <- nrow(b)
    if (nY != nT) 
      stop("sample sizes in the longitudinal and event processes differ; ", 
           "maybe you forgot the cluster() argument.\n")
    if (class(init) == "jointModel") {
      init.object <- init
      init <- init$coefficients
    }
    else init.object <- NULL
    TermsX <- lmeObject$terms
    data <- lmeObject$data[all.vars(TermsX)]
    data <- data[complete.cases(data), ]
    formYx <- formula(lmeObject)
    mfX <- model.frame(TermsX, data = data)
    X <- model.matrix(formYx, mfX)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    Z <- model.matrix(formYz, mfZ)
    y.long <- model.response(mfX, "numeric")
    data.id <- data[!duplicated(id), ]
    data.id <- data.id[idT, ]
    if (!timeVar %in% names(data)) 
      stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")
    max.timeY <- tapply(data[[timeVar]], id, max)
    max.timeT <- tapply(Time, idT, max)
    if (!all(max.timeT >= max.timeY)) {
      idnams <- factor(lmeObject$groups[[1]])
      stop("\nit seems that there are longitudinal measurements taken after the event times for some subjects ",
           "(i.e., check subject(s): ", paste(levels(idnams)[(max.timeT < max.timeY)], collapse = ", "), ").")
    }
    data.id[[timeVar]] <- pmax(Time - lag, 0)
    #### Value and slope from the mixed sub-model ####
    if (parameterization %in% c("value", "both")) {
      mfX.id <- model.frame(TermsX, data = data.id)
      mfZ.id <- model.frame(TermsZ, data = data.id)
      Xtime <- model.matrix(formYx, mfX.id)
      Ztime <- model.matrix(formYz, mfZ.id)
      long <- c(X %*% fixef(lmeObject)) + rowSums(Z * b[id, ])
    }
    if (parameterization %in% c("slope", "both")) {
      mfX.deriv <- model.frame(terms(derivForm$fixed), data = data)
      TermsX.deriv <- attr(mfX.deriv, "terms")
      mfZ.deriv <- model.frame(terms(derivForm$random), data = data)
      TermsZ.deriv <- attr(mfZ.deriv, "terms")
      mfX.deriv.id <- model.frame(TermsX.deriv, data = data.id)
      mfZ.deriv.id <- model.frame(TermsZ.deriv, data = data.id)
      Xtime.deriv <- model.matrix(derivForm$fixed, mfX.deriv.id)
      Ztime.deriv <- model.matrix(derivForm$random, mfZ.deriv.id)
      Xderiv <- model.matrix(derivForm$fixed, mfX.deriv)
      Zderiv <- model.matrix(derivForm$random, mfZ.deriv)
      long.deriv <- as.vector(c(Xderiv %*% fixef(lmeObject)[derivForm$indFixed]) + 
                                if (length(derivForm$indRandom) > 1 || derivForm$indRandom) 
                                  rowSums(Zderiv * b[id, derivForm$indRandom, drop = FALSE])
                              else rep(0, nrow(Zderiv)))
    }
    if (parameterization == "value") 
      long.deriv <- NULL
    if (parameterization == "slope") 
      long <- NULL
    y <- list(y = y.long, logT = log(Time), d = d, lag = lag)
    x <- list(X = X, Z = Z, W = W, WintF.vl = WintF.vl, WintF.sl = WintF.sl, 
              idT = idT, nRisks = nRisks)
    x <- switch(parameterization, 
                value = c(x, list(Xtime = Xtime, Ztime = Ztime)),
                slope = c(x, list(Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
                both = c(x, list(Xtime = Xtime, Ztime = Ztime, Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)))
    #### Baseline hazards/intensities ####
    ind.noadapt <- method. %in% c("weibull-AFT-GH", "weibull-PH-GH", "piecewise-PH-GH", "Cox-PH-GH", "spline-PH-GH")
    con <- list(only.EM = FALSE, iter.EM = if (method == "spline-PH-GH") 120 else 50, 
                iter.qN = 350, optimizer = "optim", tol1 = 0.001, tol2 = 1e-04, 
                tol3 = if (!CompRisk | !Mstate) sqrt(.Machine$double.eps) else 1e-09,
                numeriDeriv = "fd", eps.Hes = 1e-06, parscale = NULL, 
                step.max = 0.1, backtrackSteps = 2, knots = NULL, ObsTimes.knots = TRUE, 
                lng.in.kn = if (method == "piecewise-PH-GH") 6 else 5, 
                ord = 4, equal.strata.knots = TRUE, typeGH = if (ind.noadapt) "simple" else "adaptive", 
                GHk = if (ncol(Z) < 3 && nrow(Z) < 2000) 15 else 9, 
                GKk = if (method == "piecewise-PH-GH" || length(Time) > nRisks * nT) 7 else 15, 
                verbose = FALSE)
    if (method == "Cox-PH-GH") {
      con$only.EM <- TRUE
      con$iter.EM <- 200
      con$GHk <- if (ncol(Z) == 1) 
        15
      else if (ncol(Z) == 2) 
        11
      else 9
    }
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (con$typeGH != "simple" && !"GHk" %in% namc) {
      con$GHk <- if (ncol(Z) <= 3 && nrow(Z) < 2000) 
        5
      else 3
    }
    if (length(noNms <- namc[!namc %in% namC]) > 0) 
      warning("unknown names in 'control': ", paste(noNms, collapse = ", "))
    if (method == "Cox-PH-GH" && !con$only.EM) 
      stop("with method 'Cox-PH-GH' only the EM algorithm is used.\n")
    if (method == "Cox-PH-GH" && any(!is.na(match(c("iter.qN", "optimizer"), namc)))) 
      warning("method 'Cox-PH-GH' uses only the EM algorithm.\n")
    if (method %in% c("weibull-AFT-GH", "weibull-PH-GH", "spline-PH-GH", "spline-PH-Laplace")) {
      wk <- JM:::gaussKronrod(con$GKk)$wk
      sk <- JM:::gaussKronrod(con$GKk)$sk
      if (LongFormat) {
        Time0 <- survObject$y[, 1]
        P <- (Time - Time0)/2
        P1 <- (Time + Time0)/2
        st <- outer(P, sk) + P1
      }
      else {
        P <- as.vector(Time)/2
        st <- outer(P, sk + 1)
      }
      dimnames(st) <- names(P) <- NULL
      id.GK <- rep(seq_along(Time), each = con$GKk)
      data.id2 <- data.id[id.GK, , drop = FALSE]
      data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
      if (parameterization %in% c("value", "both")) {
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        Xs <- model.matrix(formYx, mfX)
        Zs <- model.matrix(formYz, mfZ)
      }
      if (parameterization %in% c("slope", "both")) {
        mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
        mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
        Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
        Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
      }
      Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
      Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
      x <- c(x, list(P = P, st = c(t(st)), wk = wk, Ws.intF.vl = Ws.intF.vl, Ws.intF.sl = Ws.intF.sl))
      x <- switch(parameterization, 
                  value = c(x, list(Xs = Xs, Zs = Zs)),
                  slope = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv)),
                  both = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv, Xs = Xs, Zs = Zs)))
      #### method == "spline-PH-GH" || method == "spline-PH-Laplace" ####    
      if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
        strt <- if (is.null(survObject$strata)) 
          gl(1, length(Time))
        else survObject$strata
        nstrt <- length(levels(strt))
        split.Time <- split(Time, strt)
        ind.t <- if (LongFormat) {
          if(!Mstate){
            unlist(tapply(idT, idT, function(x) c(rep(FALSE, length(x) - 1), TRUE)))
          }
          else { unlist(tapply(idT, idT, function(x) c(as.logical(data.Mstate[data.Mstate[, ID.Mstate] %in% x, "status"])))) }
        }
        else {
          rep(TRUE, length(Time))
        }
        kn <- if (con$equal.strata.knots) {
          kk <- if (is.null(con$knots)) {
            pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
            pp <- tail(head(pp, -1), -1)
            quantile(Time[ind.t], pp, names = FALSE)
          }
          else {
            con$knots
          }
          kk <- kk[kk < max(Time)]
          rr <- rep(list(sort(c(rep(range(Time, st), con$ord), kk))), nstrt)
          names(rr) <- names(split.Time)
          rr
        }
        else {
          spt <- if (length(Time) > nT & !CompRisk & !Mstate)
            mapply(function(x, y) {
              x[unlist(tapply(y, y, function(z) c(rep(FALSE, length(z) - 1), TRUE)))]
            }, split.Time, split(idT, strt), SIMPLIFY = FALSE)
          else mapply(function(x, y) {x[y]}, split.Time, split(ind.t, strt))
          lapply(spt, function(t) {
            kk <- if (is.null(con$knots)) {
              pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
              pp <- tail(head(pp, -1), -1)
              quantile(t, pp, names = FALSE)
            }
            else {
              con$knots
            }
            kk <- kk[kk < max(t)]
            sort(c(rep(range(Time, st), con$ord), kk))
          })
        }        
        con$knots <- kn
        W2 <- mapply(function(k, t) splineDesign(k, t, ord = con$ord), 
                     kn, split.Time, SIMPLIFY = FALSE)
        if (any(sapply(W2, colSums) == 0)) 
          stop("\nsome of the knots of the B-splines basis are set outside the range", 
               "\n   of the observed event times for one of the strata; refit the model", 
               "\n   setting the control argument 'equal.strata.knots' to FALSE.")
        W2 <- mapply(function(w2, ind) {
          out <- matrix(0, length(Time), ncol(w2))
          out[strt == ind, ] <- w2
          out
        }, W2, levels(strt), SIMPLIFY = FALSE)
        W2 <- do.call(cbind, W2)
        strt.s <- rep(strt, each = con$GKk)
        split.Time <- split(c(t(st)), strt.s)
        W2s <- mapply(function(k, t) splineDesign(k, t, ord = con$ord), 
                      kn, split.Time, SIMPLIFY = FALSE)
        W2s <- mapply(function(w2s, ind) {
          out <- matrix(0, length(Time) * con$GKk, ncol(w2s))
          out[strt.s == ind, ] <- w2s
          out
        }, W2s, levels(strt), SIMPLIFY = FALSE)
        W2s <- do.call(cbind, W2s)
        y <- c(y, list(strata = strt))
        x <- c(x, list(W2 = W2, W2s = W2s))
      }
    }
    if (method == "piecewise-PH-GH") { 
      wk <- JM:::gaussKronrod(con$GKk)$wk
      sk <- JM:::gaussKronrod(con$GKk)$sk
      nk <- length(sk)
      if (is.null(con$knots) || !is.numeric(con$knots)) {
        Q <- con$lng.in.kn + 1
        qs <- if (con$ObsTimes.knots) {
          unique(quantile(Time, seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)])
        }
        else {
          unique(quantile(Time[d == 1], seq(0, 1, len = Q - 1), names = FALSE))
        }
        qs <- qs + 1e-06
        if (max(qs) > max(Time)) 
          qs[which.max(qs)] <- max(Time) - 1e-06
        con$knots <- qs
        qs <- c(0, qs, max(Time) + 1)
        Q <- length(qs) - 1
      }
      else {
        qs <- c(0, sort(con$knots), max(Time) + 1)
        Q <- length(qs) - 1
      }
      ind <- findInterval(Time, qs, rightmost.closed = TRUE)
      D <- matrix(0, length(ind), Q)
      D[cbind(seq_along(ind), ind)] <- 1
      D <- D * d
      Tiq <- outer(Time, qs, pmin)
      Lo <- Tiq[, 1:Q]
      Up <- Tiq[, 2:(Q + 1)]
      T <- Up - Lo
      P <- T/2
      P[P < con$tol3] <- as.numeric(NA)
      P1 <- (Up + Lo)/2
      st <- matrix(0, nY, nk * Q)
      skQ <- rep(sk, Q)
      for (i in seq_len(nY)) {
        st[i, ] <- rep(P[i, ], each = nk) * skQ + rep(P1[i, ], each = nk)
      }
      y <- c(y, list(ind.D = ind))
      id.GK <- rep(seq_len(nY), rowSums(!is.na(st)))
      P <- c(t(P))
      data.id2 <- data.id[rep(seq_len(nY), each = nk * Q), ]
      data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
      data.id2 <- data.id2[!is.na(data.id2[[timeVar]]), ]
      if (parameterization %in% c("value", "both")) {
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        Xs <- model.matrix(formYx, mfX)
        Zs <- model.matrix(formYz, mfZ)
      }
      if (parameterization %in% c("slope", "both")) {
        mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
        mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
        Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
        Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
      }
      Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
      Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
      x <- c(x, list(P = P[!is.na(P)], st = st[!is.na(st)], 
                     wk = wk, id.GK = id.GK, Q = Q, Ws.intF.vl = Ws.intF.vl, 
                     Ws.intF.sl = Ws.intF.sl))
      x <- switch(parameterization, 
                  value = c(x, list(Xs = Xs, Zs = Zs)),
                  slope = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv)),
                  both = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv, Xs = Xs, Zs = Zs)))
    }
    if (method == "Cox-PH-GH") {
      unqT <- sort(unique(Time[d == 1]))
      times <- lapply(Time, function(t) unqT[t >= unqT])
      ind.len <- sapply(times, length)
      indT <- rep(1:nrow(data.id), ind.len)
      data.id2 <- data.id[indT, ]
      data.id2[timeVar] <- pmax(unlist(times, use.names = FALSE) - lag, 0)
      if (parameterization %in% c("value", "both")) {
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        Xtime2 <- model.matrix(formYx, mfX)
        Ztime2 <- model.matrix(formYz, mfZ)
      }
      if (parameterization %in% c("slope", "both")) {
        mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
        mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
        Xtime2.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
        Ztime2.deriv <- model.matrix(derivForm$random, mfZ.deriv)
      }
      x <- c(x, list(indT = indT))
      x <- switch(parameterization,
                  value = c(x, list(Xtime2 = Xtime2, Ztime2 = Ztime2)),
                  slope = c(x, list(Xtime2.deriv = Xtime2.deriv, Ztime2.deriv = Ztime2.deriv)),
                  both = c(x, list(Xtime2.deriv = Xtime2.deriv, Ztime2.deriv = Ztime2.deriv, Xtime2 = Xtime2, Ztime2 = Ztime2)))
    }
    #### Out ####
    if (class(init.object) == "jointModel" &&
        (!identical(X, init.object$x$X) | !identical(Z, init.object$x$Z) |
         !identical(W, init.object$x$W) | !identical(round(W2, 6), round(init.object$x$W2, 6)) |
         !identical(WintF.vl, init.object$x$WintF.vl) | !identical(WintF.sl, init.object$x$WintF.sl)))
      stop ("The model in 'init' must be the same. Please use the same data, formulas and names of variables.")
    VC <- if (class(init.object) == "jointModel") {
      ncz <- ncol(Z)
      init$D[seq_len(ncz), seq_len(ncz)]
    }
    else lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", lmeObject$sigma^2)[[1]]
    if (con$typeGH != "simple") {
      if (b.postVar.opt == TRUE)
        Vs <- b.postVar
      else {
        Vs <- vector("list", nY)
        inv.VC <- solve(VC)
        if (class(init.object) == "jointModel") {
          idT.GK <- rep(idT, each = init.object$control$GKk)
          eta.yx <- as.vector(X %*% init$betas)
          eta.tw1 <- if (!is.null(W)) 
            as.vector(W %*% init$gammas)
          else rep(0, length(init.object$y$logT))
          eta.tw2 <- as.vector(W2 %*% init$gammas.bs)
          eta.ws <- as.vector(init.object$x$W2s %*% init$gammas.bs)
          if (parameterization %in% c("value", "both")) {
            Y <- as.vector(Xtime %*% init$betas + init.object$EB$Ztimeb)
            Ys <- as.vector(init.object$x$Xs %*% init$betas + rowSums(init.object$x$Zs * b[idT.GK, ]))
            eta.t <- {
              if (is.null(transform.value))
                eta.tw2 + eta.tw1 + c(WintF.vl %*% init$alpha) * Y
              else eta.tw2 + eta.tw1 + c(WintF.vl %*% init$alpha) * transform.value(Y)
            }
            eta.s <- {
              if (is.null(transform.value))
                c(init.object$x$Ws.intF.vl %*% init$alpha) * Ys
              else c(init.object$x$Ws.intF.vl %*% init$alpha) * transform.value(Ys)
            }
          }
          if (parameterization %in% c("slope", "both")) {
            Y.deriv <- as.vector(Xtime.deriv %*% init$betas[init.object$derivForm$indFixed]) +
              if (length(init.object$derivForm$indRandom) > 1) init.object$EB$Ztimeb.deriv
            else 0
            Ys.deriv <- as.vector(init.object$x$Xs.deriv %*% init$betas[init.object$derivForm$indFixed]) +
              if (length(init.object$derivForm$indRandom) > 1) 
                as.vector(rowSums(init.object$x$Zs.deriv * b[idT.GK, init.object$derivForm$indRandom, drop = FALSE]))
            else 0
            eta.t <- if (parameterization %in% "both") 
              eta.t + c(init.object$x$WintF.sl %*% init$Dalpha) * Y.deriv
            else eta.tw2 + eta.tw1 + c(init.object$x$WintF.sl %*% init$Dalpha) * Y.deriv
            eta.s <- if (parameterization %in% "both") 
              eta.s + c(init.object$x$Ws.intF.sl %*% init$Dalpha) * Ys.deriv
            else c(init.object$x$Ws.intF.sl %*% init$Dalpha) * Ys.deriv
          }
          if (parameterization %in% c("value", "both"))
            Ws.intF.vl.alph <- c(init.object$x$Ws.intF.vl %*% init$alpha)
          if (parameterization %in% c("slope", "both"))
            Ws.intF.sl.alph <- c(init.object$x$Ws.intF.sl %*% init$Dalpha)
          Zs.deriv.full <- matrix(0, nrow(init.object$x$Zs), ncol(init.object$x$Zs))
          Zs.deriv.full[ , init.object$derivForm$indRandom] <- init.object$x$Zs.deriv
        }
        if (is.null(transform.value)) {
          for (i in 1:nY) {
            Z.i <- Z[id == i, , drop = FALSE]
            Vs[[i]] <- if (class(init.object) == "jointModel") {
              H.tb <- {
                if (init.object$parameterization == "value") {
                  Zs.i <- Zs[which(i == idT.GK), , drop = FALSE]
                  Ws.intF.vl.alph.i <- Ws.intF.vl.alph[which(i == idT.GK)]
                  Zs.assoc.i <- Zs.i * Ws.intF.vl.alph.i
                  spl.exp.eta.tw1.P.i <- split(exp(eta.tw1[idT==i]) * P[idT==i], seq_len(sum(idT==i)))
                  idT.GK.i <- id.GK[i == idT.GK]
                  spl.crprod.Zs.assoc.i <- split(lapply(split(Zs.assoc.i, seq_len(nrow(Zs.i))),
                                                        function(x) tcrossprod(x)), idT.GK.i)
                  spl.exp.eta.ws.s.i <- split(exp(eta.ws[which(i == idT.GK)] + eta.s[which(i == idT.GK)]), idT.GK.i)
                  spl.res.GK.i <- mapply(function(x, y) mapply(function(x,y) x * y, x, y, SIMPLIFY = F),
                                         spl.exp.eta.ws.s.i, spl.crprod.Zs.assoc.i, SIMPLIFY = F)
                  spl.res.i <- lapply(lapply(spl.res.GK.i, function(x)
                    mapply(function(y, z) y * z, x, split(wk, seq_len(init.object$control$GKk)), SIMPLIFY = F)),
                    function(x) Reduce('+', x))
                  Reduce('+', mapply(function(x, y) x * y, spl.exp.eta.tw1.P.i, spl.res.i, SIMPLIFY = F))
                }
                else if (init.object$parameterization == "slope") {
                  Zs.deriv.full.i <- Zs.deriv.full[which(i == idT.GK), , drop = FALSE]
                  Ws.intF.sl.alph.i <- Ws.intF.sl.alph[which(i == idT.GK)]
                  Zs.assoc.i <- Zs.deriv.full.i * Ws.intF.sl.alph.i
                  spl.exp.eta.tw1.P.i <- split(exp(eta.tw1[idT==i]) * P[idT==i], seq_len(sum(idT==i)))
                  idT.GK.i <- id.GK[i == idT.GK]
                  spl.crprod.Zs.assoc.i <- split(lapply(split(Zs.assoc.i, seq_len(nrow(Zs.i))),
                                                        function(x) tcrossprod(x)), idT.GK.i)
                  spl.exp.eta.ws.s.i <- split(exp(eta.ws[which(i == idT.GK)] + eta.s[which(i == idT.GK)]), idT.GK.i)
                  spl.res.GK.i <- mapply(function(x, y) mapply(function(x,y) x * y, x, y, SIMPLIFY = F),
                                         spl.exp.eta.ws.s.i, spl.crprod.Zs.assoc.i, SIMPLIFY = F)
                  spl.res.i <- lapply(lapply(spl.res.GK.i, function(x)
                    mapply(function(y, z) y * z, x, split(wk, seq_len(init.object$control$GKk)), SIMPLIFY = F)),
                    function(x) Reduce('+', x))
                  Reduce('+', mapply(function(x, y) x * y, spl.exp.eta.tw1.P.i, spl.res.i, SIMPLIFY = F))
                }
                else if (init.object$parameterization == "both") {
                  Zs.i <- Zs[which(i == idT.GK), , drop = FALSE]
                  Zs.deriv.full.i <- Zs.deriv.full[which(i == idT.GK), , drop = FALSE]
                  Ws.intF.vl.alph.i <- Ws.intF.vl.alph[which(i == idT.GK)]
                  Ws.intF.sl.alph.i <- Ws.intF.sl.alph[which(i == idT.GK)]
                  Zs.assoc.i <- Zs.i * Ws.intF.vl.alph.i + Zs.deriv.full.i * Ws.intF.sl.alph.i
                  spl.exp.eta.tw1.P.i <- split(exp(eta.tw1[idT==i]) * P[idT==i], seq_len(sum(idT==i)))
                  idT.GK.i <- id.GK[i == idT.GK]
                  spl.crprod.Zs.assoc.i <- split(lapply(split(Zs.assoc.i, seq_len(nrow(Zs.i))),
                                                        function(x) tcrossprod(x)), idT.GK.i)
                  spl.exp.eta.ws.s.i <- split(exp(eta.ws[which(i == idT.GK)] + eta.s[which(i == idT.GK)]), idT.GK.i)
                  spl.res.GK.i <- mapply(function(x, y) mapply(function(x,y) x * y, x, y, SIMPLIFY = F),
                                         spl.exp.eta.ws.s.i, spl.crprod.Zs.assoc.i, SIMPLIFY = F)
                  spl.res.i <- lapply(lapply(spl.res.GK.i, function(x)
                    mapply(function(y, z) y * z, x, split(wk, seq_len(init.object$control$GKk)), SIMPLIFY = F)),
                    function(x) Reduce('+', x))
                  Reduce('+', mapply(function(x, y) x * y, spl.exp.eta.tw1.P.i, spl.res.i, SIMPLIFY = F))
                }
              }
              solve(crossprod(Z.i)/init$sigma^2 + inv.VC + H.tb)
            }
            else solve(crossprod(Z.i)/lmeObject$sigma^2 + inv.VC)
          }
        }
        else {
          if (class(init.object) != "jointModel"){
            for(i in 1:nY) {
              Z.i <- Z[id == i, , drop = FALSE]
              Vs[[i]] <- solve(crossprod(Z.i)/lmeObject$sigma^2 + inv.VC)
            }
          }
          else {
            exp.eta.tw.P <- exp(eta.tw1) * P
            Int <- wk * exp(eta.ws + eta.s)
            H.tb <- array(0, dim = c(nY, ncz, ncz))
            for (i in 1:ncz) {
              for (j in i:ncz) {
                ZZ <- if (parameterization == "value") {
                  to.add <- d * c(WintF.vl %*% init$alpha) * H.transform.value(Ztime[, i], Ztime[, j], Y)
                  Ws.intF.vl.alph * H.transform.value(Zs[, i], Zs[, j], Ys) +
                    Ws.intF.vl.alph^2 * gr.transform.value(Xs[, i], Ys) * gr.transform.value(Xs[, j], Ys)
                }
                else if (parameterization == "slope") {
                  if (i %in% init.object$derivForm$indRandom && j %in% init.object$derivForm$indRandom) {
                    ii <- match(i, init.object$derivForm$indRandom)
                    jj <- match(j, init.object$derivForm$indRandom)
                    Ws.intF.sl.alph^2 * Zs.deriv[, ii] * Zs.deriv[, jj]
                  }
                  else 0
                }
                else {
                  if (i %in% init.object$derivForm$indRandom && j %in% init.object$derivForm$indRandom) {
                    ii <- match(i, init.object$derivForm$indRandom)
                    jj <- match(j, init.object$derivForm$indRandom)
                    to.add <- d * c(WintF.vl %*% init$alpha) * H.transform.value(Ztime[, i], Ztime[, j], Y)
                    Ws.intF.vl.alph * H.transform.value(Zs[, i], Zs[, j], Ys) +
                      (Ws.intF.vl.alph * gr.transform.value(Zs[, i], Ys) + Ws.intF.sl.alph * Zs.deriv[, ii]) *
                      (Ws.intF.vl.alph * gr.transform.value(Zs[, j], Ys) + Ws.intF.sl.alph * Zs.deriv[, jj])
                  }
                  else if (i %in% init.object$derivForm$indRandom && !j %in% init.object$derivForm$indRandom) {
                    ii <- match(i, init.object$derivForm$indRandom)
                    to.add <- d * c(WintF.vl %*% init$alpha) * H.transform.value(Ztime[, i], Ztime[, j], Y)
                    Ws.intF.vl.alph * H.transform.value(Zs[, i], Zs[, j], Ys) +
                      (Ws.intF.vl.alph * gr.transform.value(Zs[, i], Ys) + Ws.intF.sl.alph * Zs.deriv[, ii]) *
                      (Ws.intF.vl.alph * gr.transform.value(Zs[, j], Ys))
                  }
                  else if (!i %in% init.object$derivForm$indRandom && j %in% init.object$derivForm$indRandom) {
                    jj <- match(j, init.object$derivForm$indRandom)
                    to.add <- d * c(WintF.vl %*% init$alpha) * H.transform.value(Ztime[, i], Ztime[, j], Y)
                    Ws.intF.vl.alph * H.transform.value(Zs[, i], Zs[, j], Ys) +
                      (Ws.intF.vl.alph * gr.transform.value(Zs[, i], Ys)) *
                      (Ws.intF.vl.alph * gr.transform.value(Zs[, j], Ys) + Ws.intF.sl.alph * Zs.deriv[, jj])
                  }
                  else {
                    to.add <- d * c(WintF.vl %*% init$alpha) * H.transform.value(Ztime[, i], Ztime[, j], Y)
                    Ws.intF.vl.alph * H.transform.value(Zs[, i], Zs[, j], Ys) +
                      Ws.intF.vl.alph^2 * gr.transform.value(Zs[, i], Ys) * gr.transform.value(Zs[, j], Ys)
                  }
                }
                ki <- exp.eta.tw.P * rowsum(Int * ZZ, id.GK, reorder = FALSE)
                H.tb[, i, j] <- as.vector(rowsum( ki - to.add, idT, reorder = FALSE))
              }
            }
            for(i in 1:nY) {
              H.tb[i, , ][lower.tri(H.tb[i, , ])] <- t(H.tb[i, , ])[lower.tri(H.tb[i, , ])]
              Z.i <- Z[id == i, , drop = FALSE]
              Vs[[i]] <- solve(crossprod(Z.i)/init$sigma^2 + inv.VC + H.tb[i, , ])
            }
          }
        } 
      }
      con$inv.chol.VCs <- try(lapply(Vs, function(x) solve(chol(solve(x)))), TRUE)
      con$det.inv.chol.VCs <- sapply(con$inv.chol.VCs, det)
    }
    con$inv.chol.VC <- solve(chol(solve(VC)))
    con$det.inv.chol.VC <- det(con$inv.chol.VC)
    con$ranef <- b
    if (all(VC[upper.tri(VC)] == 0)) 
      VC <- diag(VC)
    if (class(init.object) != "jointModel") {
      init.surv <- initial.surv(Time, d, W, WintF.vl, WintF.sl,
                                id, times = data[[timeVar]], method, parameterization, 
                                long = long, long.deriv = long.deriv, 
                                extra = list(W2 = x$W2, control = con, ii = idT, strata = survObject$strata), 
                                LongFormat = CompRisk | Mstate | length(Time) > nT,
                                transform.value = transform.value)
      
      if (method == "Cox-PH-GH" && length(init.surv$lambda0) < 
          length(unqT)) 
        init.surv$lambda0 <- basehaz(survObject)$hazard
      initial.values <- c(list(betas = fixef(lmeObject), sigma = lmeObject$sigma, D = VC), init.surv)
      if (!is.null(init)) {
        nams1 <- names(init)
        nams2 <- names(initial.values)
        if (!is.list(init) || length(noNms <- nams1[!nams1 %in% nams2])) {
          warning("unknown names in 'init': ", paste(noNms, collapse = ", "))
        }
        else {
          initial.values[nams1] <- init
        }
      }
    }
    else initial.values <- init
    rmObjs <- c(names(x), "y.long", "mfX", "mfZ", "data.id2")
    rm(list = rmObjs)
    gc()
    out <- switch(method, 
                  `Cox-PH-GH` = JM:::phGH.fit(x, y, id, initial.values, parameterization, derivForm, con),
                  `weibull-AFT-GH` = JM:::weibullAFTGH.fit(x, y, id, initial.values, scaleWB, parameterization, derivForm, con),
                  `weibull-PH-GH` = JM:::weibullPHGH.fit(x, y, id, initial.values, scaleWB, parameterization, derivForm, con),
                  `piecewise-PH-GH` = JM:::piecewisePHGH.fit(x, y, id, initial.values, parameterization, derivForm, con), 
                  `spline-PH-GH` = splinePHGH.fit(x, y, id, initial.values, parameterization, derivForm, con, transform.value),
                  `ch-Laplace` = JM:::chLaplace.fit(x, y, id, initial.values, b, parameterization, derivForm, con))
    H <- out$Hessian
    if (any(is.na(H) | !is.finite(H))) {
      warning("infinite or missing values in Hessian at convergence.\n")
    }
    else {
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        warning("Hessian matrix at convergence is not positive definite.\n")
    }
    out$coefficients <- out$coefficients[!sapply(out$coefficients, is.null)]
    out$x <- x
    out$y <- y
    out$times <- data[[timeVar]]
    out$data <- data
    out$data.id <- data.id
    out$method <- method
    out$termsYx <- TermsX
    out$termsYz <- TermsZ
    if (parameterization %in% c("slope", "both")) {
      out$termsYx.deriv <- TermsX.deriv
      out$termsYz.deriv <- TermsZ.deriv
    }
    out$termsT <- survObject$terms
    out$formYx <- formYx
    out$formYz <- formYz
    out$formT <- formT
    out$timeVar <- timeVar
    out$control <- con
    out$parameterization <- parameterization
    out$derivForm <- derivForm
    out$interFact <- interFact
    out$CompRisk <- CompRisk
    out$Mstate <- Mstate
    out$data.Mstate <- data.Mstate
    out$ID.Mstate <- data.Mstate[, ID.Mstate]
    out$transform.value <- transform.value
    out$gr.transform.value <- gr.transform.value
    out$H.transform.value <- H.transform.value
    out$LongFormat <- LongFormat
    out$assignY <- attr(lmeObject$fixDF, "assign")[-1]
    out$assignT <- survObject$assign
    out$call <- cl
    class(out) <- "jointModel"
    out
  }

