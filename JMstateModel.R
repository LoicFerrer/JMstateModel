JMstateModel <-
  function (lmeObject, survObject, timeVar, 
            parameterization = c("value", "slope", "both"), 
            method = c("weibull-PH-aGH", "weibull-PH-GH", "weibull-AFT-aGH", "weibull-AFT-GH", "piecewise-PH-aGH", 
                       "piecewise-PH-GH", "Cox-PH-aGH", "Cox-PH-GH", "spline-PH-aGH", "spline-PH-GH", "ch-Laplace"),
            interFact = NULL, derivForm = NULL, lag = 0, scaleWB = NULL, CompRisk = FALSE, init = NULL, control = list(),
            Mstate = FALSE, data.Mstate = NULL, ID.Mstate = NULL, 
            # data.Mstate is the data used in the coxph,
            # ID.Mstate is the name of the column of the id in data.Mstate (class(ID.Mstate)=="character")
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
    if (!inherits(survObject, "coxph") && !inherits(survObject, 
                                                    "survreg")) 
      stop("\n'survObject' must inherit from class coxph or class survreg.")
    if (!is.matrix(survObject$x)) 
      stop("\nuse argument 'x = TRUE' in ", if (inherits(survObject, 
                                                         "coxph")) 
        "'coxph()'."
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
    if (parameterization %in% c("slope", "both") && method %in% 
          c("Cox-PH-GH", "ch-Laplace")) 
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
      stop("\nto fit a competing risks or multi-state joint model you must choose as method 'spline-PH-GH'", 
           " and include a strata() in the specification of the coxph().")
    }
    if (Mstate && (!is.null(survObject$cluster))) {
      stop("\nto fit a multi-state model with coxph(), you can not use time-dependent covariate.",
           "Thus, do not use 'cluster()' in coxph().")
    }
    if (Mstate && (is.null(ID.Mstate) || !is.vector(data.Mstate[ ,ID.Mstate]) || is.null(data.Mstate))) {
      stop("\nwhen Mstate is TRUE, you must specify data.Mstate and ID.Mstate",
           " which respectively correspond to the data used in coxph() and the name of the column of the id in data.Mstate.")
    }
    if (CompRisk && Mstate) {
      stop("\nyou must choose between a competing risks joint model and a multi-state joint model.")
    }
    
    #### Survival sub-part ####
    formT <- formula(survObject)
    if (inherits(survObject, "coxph")) { # si coxph pour le survobject (ce qui est le cas)
      W <- survObject$x # Matrice des facteurs pronostiques
      keepW <- suppressWarnings(!is.na(survObject$coefficients)) # On ne garde que les facteurs pronostiques pr lesquels il n'y a pas de pb ds l'estimation du Cox
      W <- W[, keepW, drop = FALSE]
      
      # Si risques compétitifs, nRisks <- nb de strates (= nb de risques de base à estimer)
      if (CompRisk || Mstate) {
        nRisks <- length(unique(survObject$strata)) 
      }
      else {
        nRisks <- 1
      }
      
      # Vecteur réponse (status) et temps associés (start & stop)
      surv <- survObject$y 
      
      # Processus de comptage
      if (attr(surv, "type") == "right") {
        LongFormat <- FALSE
        Time <- survObject$y[, 1]
        d <- survObject$y[, 2]
      } # Si pas de process de comptage, LongFormat <- FALSE, Time <- temps de survie, d <- indicateur d'événement
      
      else if (attr(surv, "type") == "counting") {
        LongFormat <- TRUE
        if (is.null(survObject$model)) 
          stop("\nplease refit the Cox model including in the ", 
               "call to coxph() the argument 'model = TRUE'.")
        Time <- survObject$y[, 2]
        d <- survObject$y[, 3]
      } # Si on a process de comptage, LongFormat <- TRUE, Time <- temps de survie, d <- indicateur d'événement
      
      
      idT <- if (!is.null(survObject$model$cluster)) {
        as.vector(unclass(survObject$model$cluster)) # Si on utilise des cov. dpdtes du temps (d'où cluster = TRUE),  idT <- ID de l'individu (donné par chaque cluster)
      } # as.vector(unclass) permet de mettre en numérique un facteur
      else {
        if (!Mstate) {
          if(!CompRisk) {
            seq_along(Time) # Si pas de cov dpdtes du temps et pas de risques compét', idT <- 1:nombre_de_temps (c'est logique : chaque valeur est associé à un individu)
          }
          else { rep(seq_len(length(Time)/nRisks), each = nRisks) } # Si on a des cov non dptes du temps mais on a des risques compét, chaque ID est répété X (= nb_de_strates) fois
        }
        else { data.Mstate[ ,ID.Mstate] } # Dans le cas où Mstate = TRUE
      }
      idT <- match(idT, unique(idT))  # match returns a vector of the positions of (first) matches of its first argument in its second.
    }
    
    # si pas coxph pr le survobject (ce qui n'est pas notre cas!)
    else { 
      W <- survObject$x[, -1, drop = FALSE]
      Time <- exp(survObject$y[, 1])
      d <- survObject$y[, 2]
      idT <- seq_along(Time)
      LongFormat <- FALSE
      nRisks <- 1
    }
    
    nT <- length(unique(idT)) # Nombre d'individus
    
    if (LongFormat && is.null(survObject$model$cluster) && !Mstate ) 
        stop("\nuse argument 'model = TRUE' and cluster() in coxph().")
    if (!length(W)) 
      W <- NULL # Si pas de facteurs pronostiques
    if (sum(d) < 5) 
      warning("\nmore than 5 events are required.") # Au moins 5 événements observés sont requis ds le process de survie
    
    # WintF.vl : Matrice de design correspondant à ce qui a été formulé dans interFact$value
    # WintF.sl : Matrice de design correspondant à ce qui a été formulé dans interFact$slope
    WintF.vl <- WintF.sl <- as.matrix(rep(1, length(Time)))
    if (!is.null(interFact)) {
      if (!is.null(interFact$value)) 
        WintF.vl <- if (is.null(survObject$model) || !is.null(interFact$data)) { #
          model.matrix(interFact$value, data = interFact$data) # WintF.vl si model = FALSE ds coxph
        }
      else {
        model.matrix(interFact$value, data = survObject$model)  # WintF.vl si model = TRUE ds coxph 
      }
      if (!is.null(interFact$slope)) 
        WintF.sl <- if (is.null(survObject$model) || !is.null(interFact$data)) {
          model.matrix(interFact$slope, data = interFact$data) # Idem avec WintF.sl
        }
      else {
        model.matrix(interFact$slope, data = survObject$model) # Idem
      }
    }
    
    #### Longitudinal sub-part ####
    id <- as.vector(unclass(lmeObject$groups[[1]])) # ID utilisés ds le process longit (avec des répétitions pr chaque ID puisque plusieurs valeurs par indiv.)
    b <- data.matrix(ranef(lmeObject)) # effets aléatoires pr chaque individu
    dimnames(b) <- NULL
    nY <- nrow(b)
      if (nY != nT) 
        stop("sample sizes in the longitudinal and event processes differ; ", 
             "maybe you forgot the cluster() argument.\n")
    TermsX <- lmeObject$terms # formula et terms, etc... de lme
    data <- lmeObject$data[all.vars(TermsX)] # 
    data <- data[complete.cases(data), ] # Données utilisées ds lme (matrice)
    formYx <- formula(lmeObject)
    mfX <- model.frame(TermsX, data = data) # Données utilisées ds lme + fonction du temps (f2)
    X <- model.matrix(formYx, mfX) # Matrice des effets fixes du modele lineaire mixte
    
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]]) # Idem avec les effets aléatoires
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    Z <- model.matrix(formYz, mfZ) # Matrice des effets aléatoires du modele lineaire mixte
    
    y.long <- model.response(mfX, "numeric") # Vecteur des réponses longitudinales (logPSA)
    
    data.id <- data[!duplicated(id), ]
    data.id <- data.id[idT, ] # a data.frame containing the variables for the linear mixed effects model at the time of the event. (OK)
    
    if (!timeVar %in% names(data)) 
      stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")
    max.timeY <- tapply(data[[timeVar]], id, max) # Dernier temps des données longit de chaque indiv.
    max.timeT <- 
      if(!Mstate) { tapply(Time, idT, max) } # Dernier temps des données de survie de chaque indiv., 
    else { tapply(Time, idT, min) } # ou premier temps de transition en multi-états
    if (!all(max.timeT >= max.timeY)) {
      idnams <- factor(lmeObject$groups[[1]])
      stop("\nit seems that there are longitudinal measurements taken after the event times for some subjects ", 
           "(i.e., check subject(s): ", paste(levels(idnams)[(max.timeT < max.timeY)], collapse = ", "), ").") # On vérifie qu'il n'y a pas de données longit apres l'événement
    }
    data.id[[timeVar]] <- pmax(Time - lag, 0)
        
    #### Current level and slope since the mixed sub-model ####
    if (parameterization %in% c("value", "both")) {
      mfX.id <- model.frame(TermsX, data = data.id)
      mfZ.id <- model.frame(TermsZ, data = data.id)
      Xtime <- model.matrix(formYx, mfX.id)
      Ztime <- model.matrix(formYz, mfZ.id)
      long <- c(X %*% fixef(lmeObject)) + rowSums(Z * b[id, ]) # curent level: Y_i^{*}(t) (from the mixed sub-model)
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
                              else rep(0, nrow(Zderiv))) # slope: dY_i^{*}(t) / dt (since the mixed sub-model)
    }
    if (parameterization == "value") 
      long.deriv <- NULL # si dépendence = niveau courrant, pas besoin de pente
    if (parameterization == "slope") 
      long <- NULL # si dépendence = pente, pas besoin du niveau courrant
    
    y <- list(y = y.long, logT = log(Time), d = d, lag = lag) # Liste avec vecteur des réponses de la partie des longit et surv (réponses longit, log(temps de survie), et indicateur d'événement)
    x <- list(X = X, Z = Z, W = W, WintF.vl = WintF.vl, WintF.sl = WintF.sl, 
              idT = idT, nRisks = nRisks) # Liste avec matrices des effets fixes et des effets aléatoires du modele mixte, Matrice des facteurs pronostiques, 
    # Matrice de design (intercept & cov. demandées) correspondant à ce qui a été formulé dans interFact$value,
    # Matrice de design (intercept & cov. demandées) correspondant à ce qui a été formulé dans interFact$slope,
    # ID associés aux temps de suivi
    # Nombre de risques
    # switch(EXPR, ...) evaluates EXPR and accordingly chooses one of the further arguments (in ...). 
    # switch(1, invisible(pi), pi) # ... (rien)
    # switch(2, invisible(pi), pi) # [1] 3.141593
    x <- switch(parameterization, 
                value = c(x, list(Xtime = Xtime, Ztime = Ztime)),
                slope = c(x, list(Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
                both = c(x, list(Xtime = Xtime, Ztime = Ztime, Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)))
    # x : liste avec les matrices de design pour les process longit et de surv
    
    #### Mise en place des méthodes d'intégration ####
    ind.noadapt <- method. %in% c("weibull-AFT-GH", "weibull-PH-GH", 
                                  "piecewise-PH-GH", "Cox-PH-GH", "spline-PH-GH")
    con <- list(only.EM = FALSE, iter.EM = if (method == "spline-PH-GH") 120 else 50, 
                iter.qN = 350, optimizer = "optim", tol1 = 0.001, tol2 = 1e-04, 
                tol3 = if (!CompRisk | !Mstate) sqrt(.Machine$double.eps) else 1e-09, 
                numeriDeriv = "fd", eps.Hes = 1e-06, parscale = NULL, 
                step.max = 0.1, backtrackSteps = 2, knots = NULL, ObsTimes.knots = TRUE, 
                lng.in.kn = if (method == "piecewise-PH-GH") 6 else 5, 
                ord = 4, equal.strata.knots = TRUE, typeGH = if (ind.noadapt) "simple" else "adaptive", 
                GHk = if (ncol(Z) < 3 && nrow(Z) < 2000) 15 else 9, GKk = if (method == "piecewise-PH-GH" || length(Time) > nRisks * nT) 7 else 15, 
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
      warning("unknown names in 'control': ", paste(noNms, 
                                                    collapse = ", "))
    if (method == "Cox-PH-GH" && !con$only.EM) 
      stop("with method 'Cox-PH-GH' only the EM algorithm is used.\n")
    if (method == "Cox-PH-GH" && any(!is.na(match(c("iter.qN", 
                                                    "optimizer"), namc)))) 
      warning("method 'Cox-PH-GH' uses only the EM algorithm.\n")
    if (method %in% c("weibull-AFT-GH", "weibull-PH-GH", "spline-PH-GH", 
                      "spline-PH-Laplace")) {
      wk <- gaussKronrod(con$GKk)$wk
      sk <- gaussKronrod(con$GKk)$sk
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
      x <- c(x, list(P = P, st = c(t(st)), wk = wk, Ws.intF.vl = Ws.intF.vl, 
                     Ws.intF.sl = Ws.intF.sl))
      x <- switch(parameterization, 
                  value = c(x, list(Xs = Xs, Zs = Zs)),
                  slope = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv)),
                  both = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv, Xs = Xs, Zs = Zs)))
      
    #### method == "spline-PH-GH" || method == "spline-PH-Laplace" ####    
      if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
        strt <- if (is.null(survObject$strata)) 
          gl(1, length(Time)) # Generate factors by specifying the pattern of their levels.
        else survObject$strata
        nstrt <- length(levels(strt))
        split.Time <- split(Time, strt) # split divides the data in the vector x into the groups defined by f.
        # split(x, f, drop = FALSE, ...)
        # The replacement forms replace values corresponding to such a division. unsplit reverses the effect of split.
        
        ind.t <- if (LongFormat) { # ind.t vaut TRUE quand il y a événement pour l'ID, FALSE sinon
          if(!Mstate){
            unlist(tapply(idT, idT, function(x) c(rep(FALSE, length(x) - 1), TRUE)))
          }
          else { unlist(tapply(idT, idT, function(x) c(as.logical(data.Mstate[data.Mstate[ , ID.Mstate] %in% x, "status"])) )) }
        }
        else {
          rep(TRUE, length(Time))
        }
        kn <- if (con$equal.strata.knots) {
          kk <- if (is.null(con$knots)) {
            pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
            pp <- tail(head(pp, -1), -1)
            quantile(Time[ind.t], pp, names = FALSE) # quantile des temps d'événements (multi-états ou non)
          }
          else {
            con$knots
          }
          kk <- kk[kk < max(Time)]
          rr <- rep(list(sort(c(rep(range(Time, st), con$ord), kk))), nstrt)
          names(rr) <- names(split.Time)
          rr
        }
        else { # (i.e. si !con$equal.strata.knots) --> a été mis en place (le 21/11)
          spt <- if (length(Time) > nT & !CompRisk & !Mstate) 
            mapply(function(x, y) {
              x[unlist(tapply(y, y, function(z) c(rep(FALSE, length(z) - 1), TRUE)))]
            }, split.Time, split(idT, strt), SIMPLIFY = FALSE)
          else mapply(function(x, y) {x[y]}, split.Time, split(ind.t, strt)) # On a modifié cette ligne
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
    #### method == "piecewise-PH-GH" ####
    if (method == "piecewise-PH-GH") {
      wk <- gaussKronrod(con$GKk)$wk
      sk <- gaussKronrod(con$GKk)$sk
      nk <- length(sk)
      if (is.null(con$knots) || !is.numeric(con$knots)) {
        Q <- con$lng.in.kn + 1
        qs <- if (con$ObsTimes.knots) {
          unique(quantile(Time, seq(0, 1, len = Q + 1), 
                          names = FALSE)[-c(1, Q + 1)])
        }
        else {
          unique(quantile(Time[d == 1], seq(0, 1, len = Q - 
                                              1), names = FALSE))
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
        st[i, ] <- rep(P[i, ], each = nk) * skQ + rep(P1[i, 
                                                         ], each = nk)
      }
      y <- c(y, list(ind.D = ind))
      id.GK <- rep(seq_len(nY), rowSums(!is.na(st)))
      P <- c(t(P))
      data.id2 <- data.id[rep(seq_len(nY), each = nk * Q), 
                          ]
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
    #### method == "Cox-PH-GH" ####
    if (method == "Cox-PH-GH") {
      unqT <- sort(unique(Time[d == 1]))
      times <- lapply(Time, function(t) unqT[t >= unqT])
      ind.len <- sapply(times, length)
      indT <- rep(1:nrow(data.id), ind.len)
      data.id2 <- data.id[indT, ]
      data.id2[timeVar] <- pmax(unlist(times, use.names = FALSE) - 
                                  lag, 0)
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
    VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", 
                 lmeObject$sigma^2)[[1]]
    if (con$typeGH != "simple") {
      Vs <- vector("list", nY)
      inv.VC <- solve(VC)
      for (i in 1:nY) {
        Z.i <- Z[id == i, , drop = FALSE]
        Vs[[i]] <- solve(crossprod(Z.i)/lmeObject$sigma^2 + 
                           inv.VC)
      }
      con$inv.chol.VCs <- lapply(Vs, function(x) solve(chol(solve(x))))
      con$det.inv.chol.VCs <- sapply(con$inv.chol.VCs, det)
    }
    con$inv.chol.VC <- solve(chol(solve(VC)))
    con$det.inv.chol.VC <- det(con$inv.chol.VC)
    con$ranef <- b
    if (all(VC[upper.tri(VC)] == 0)) 
      VC <- diag(VC)
    init.surv <- initial.surv(Time, d, W, WintF.vl, WintF.sl, 
                              id, times = data[[timeVar]], method, parameterization, 
                              long = long, long.deriv = long.deriv, 
                              extra = list(W2 = x$W2, control = con, ii = idT, strata = survObject$strata), 
                              LongFormat = CompRisk | Mstate | length(Time) > nT)
    if (method == "Cox-PH-GH" && length(init.surv$lambda0) < 
          length(unqT)) 
      init.surv$lambda0 <- basehaz(survObject)$hazard
    initial.values <- c(list(betas = fixef(lmeObject), sigma = lmeObject$sigma, 
                             D = VC), init.surv) # valeurs initiales
    if (!is.null(init)) {
      nams1 <- names(init)
      nams2 <- names(initial.values)
      if (!is.list(init) || length(noNms <- nams1[!nams1 %in% 
                                                    nams2])) {
        warning("unknown names in 'init': ", paste(noNms, 
                                                   collapse = ", "))
      }
      else {
        initial.values[nams1] <- init
      }
    }
    rmObjs <- c(names(x), "y.long", "mfX", "mfZ", "data.id2")
    rm(list = rmObjs)
    gc()
    out <- switch(method, 
                  `Cox-PH-GH` = phGH.fit(x, y, id, initial.values, parameterization, derivForm, con),
                  `weibull-AFT-GH` = weibullAFTGH.fit(x, y, id, initial.values, scaleWB, parameterization, derivForm, con),
                  `weibull-PH-GH` = weibullPHGH.fit(x, y, id, initial.values, scaleWB, parameterization, derivForm, con),
                  `piecewise-PH-GH` = piecewisePHGH.fit(x, y, id, initial.values, parameterization, derivForm, con), 
                  `spline-PH-GH` = splinePHGH.fit(x, y, id, initial.values, parameterization, derivForm, con),
                  `ch-Laplace` = chLaplace.fit(x, y, id, initial.values, b, parameterization, derivForm, con))
    H <- out$Hessian
    if (any(is.na(H) | !is.finite(H))) {
      warning("infinite or missing values in Hessian at convergence.\n")
    }
    else {
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        warning("Hessian matrix at convergence is not positive definite.\n")
    }
    out$coefficients <- out$coefficients[!sapply(out$coefficients, 
                                                 is.null)]
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
    out$ID.Mstate <- data.Mstate[ ,ID.Mstate]
    out$LongFormat <- LongFormat
    out$assignY <- attr(lmeObject$fixDF, "assign")[-1]
    out$assignT <- survObject$assign
    out$call <- cl
    class(out) <- "jointModel"
    out
  }


#### initial.surv : OK ####
initial.surv <- function (Time, d, W, WintF.vl, WintF.sl, id, times, method, 
                          parameterization, long = NULL, long.deriv = NULL, extra = NULL, 
                          LongFormat) {
  old <- options(warn = (-1)) # "If warn is negative all warnings are ignored."
  on.exit(options(old)) # on.exit records the expression given as its argument as needing to be executed when the current function exits
  # (either naturally or as the result of an error). This is useful for resetting graphical parameters or performing other cleanup actions.

  if (!is.null(long)) {
    long.id <- tapply(long, id, tail, 1) # prend la dernière valeur de long pour chaque id
    if (parameterization == "value") 
      longD.id <- NULL
  }
  if (!is.null(long.deriv)) {
    longD.id <- tapply(long.deriv, id, tail, 1) # prend la dernière valeur de long.deriv pour chaque id
    if (parameterization == "slope") 
      long.id <- NULL
  }
  idT <- extra$ii # Dans jointModel, on a : extra = list(W2 = x$W2, control = con, ii = idT, strata = survObject$strata)
  
  
  WW <- if (!LongFormat) { # i.e. si on a pas de CompRisk, de Mstate, et pas de cov dpdtes du temps
    cbind(W, long.id, longD.id) # cbind(matrice des facteurs pronostiques, dernière valeur du niveau courant, de la pente du biomarqueur)
  }
  else {
    cbind(W, long.id[idT], longD.id[idT]) # Si on a LongFormat, on cbind facteurs pronostiques, et dernière valeur pour chaque ID1 du niveau courant et de la pente
  }
  
  if (method %in% c("Cox-PH-GH", "weibull-PH-GH", "piecewise-PH-GH", 
                    "spline-PH-GH", "spline-PH-Laplace")) {
    if (!LongFormat) {
      DD <- data.frame(id = id, Time = Time[id], d = d[id], 
                       times = times)
      if (!is.null(long)) {
        DD$long <- long * WintF.vl[id, , drop = FALSE]
        k <- ncol(DD$long)
      }
      if (!is.null(long.deriv)) {
        DD$longD <- long.deriv * WintF.sl[id, , drop = FALSE]
        l <- ncol(DD$longD)
      }
      dW <- as.data.frame(W[id, , drop = FALSE], row.names = row.names(DD))
      if (ncol(dW)) {
        names(dW) <- paste("W", seq_along(dW), sep = "")
        DD <- cbind(DD, dW)
      }
    }
    
    else { # i.e. if LongFormat
      DD <- data.frame(Time = Time, d = d) # tableau avec temps d'événements et indicateurs associés
      if (!is.null(long)) {
        DD$long <- as.vector(long.id[idT]) * WintF.vl # matrice avec en cov données ds interFact$value mulitpliées par la dernière valeur fittée du biomarqueur
        k <- ncol(DD$long) # Nb de cov ds interFact$value
      }
      if (!is.null(long.deriv)) {
        DD$longD <- as.vector(longD.id[idT]) * WintF.sl # Idem avec la pente
        l <- ncol(DD$longD) # Idem avec la pente
      }
      dW <- as.data.frame(W, row.names = row.names(DD)) #  Idem que W en fait ...
      if (ncol(dW)) {
        names(dW) <- paste("W", seq_along(dW), sep = "")
        DD <- cbind(DD, dW)
      }
      DD$strata <- extra$strata
    }# Fin de la construction de DD
    
    if (!LongFormat) {
      DD$start <- DD$times
      DD$stop <- unlist(lapply(split(DD[c("id", "start", "Time")], DD$id),
                                function(d) c(d$start[-1], d$Time[1])))
      DD$event <- ave(DD$d, DD$id, FUN = function(x) {
        if (length(x) == 1) {
          x
        }
        else {
          x[seq(length(x) - 1)] <- 0
          x
        }
      })
    }
    
    baseCovs <- if (ncol(dW)) {
      paste("+", paste(names(dW), collapse = " + "))
    }
    else NULL
    
    form <- if (!LongFormat) {
      switch(parameterization, value = paste("Surv(start, stop, event) ~", "long", baseCovs),
                               slope = paste("Surv(start, stop, event) ~", "longD", baseCovs),
                               both = paste("Surv(start, stop, event) ~", "long + longD", baseCovs))
    }
    
    else { # i.e. if LongFormat
      switch(parameterization, value = paste("Surv(Time, d) ~", "long", baseCovs),
                               slope = paste("Surv(Time, d) ~", "longD", baseCovs),
                               both = paste("Surv(Time, d) ~", "long + longD", baseCovs))
    }
    if (!is.null(DD$strata)) 
      form <- paste(form, "+ strata(strata)")
    form <- as.formula(form)
    cph <- coxph(form, data = DD)
    coefs <- cph$coefficients
    out <- switch(parameterization, value = list(alpha = coefs[1:k], gammas = coefs[-(1:k)]),
                                    slope = list(Dalpha = coefs[1:l], gammas = coefs[-(1:l)]),
                                    both = list(alpha = coefs[1:k], Dalpha = coefs[(k + 1):(k + l)], gammas = coefs[-(1:(k + l))]))
    if (method == "Cox-PH-GH") {
      out$lambda0 <- basehaz(cph, FALSE)$hazard
    }
    if (method == "weibull-PH-GH") {
      dat <- data.frame(Time = Time, d = d)
      init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
      coefs <- -init.fit$coef/init.fit$scale
      out$gammas <- c(coefs[1], out$gammas)
      out$sigma.t <- 1/init.fit$scale
    }
    if (method == "piecewise-PH-GH") {
      dat <- data.frame(Time = Time, d = d)
      cph. <- coxph(Surv(Time, d) ~ WW, data = dat, x = TRUE)
      init.fit <- piecewiseExp.ph(cph., knots = extra$control$knots)
      coefs <- init.fit$coef
      out$xi <- exp(coefs[grep("xi", names(coefs))])
    }
    if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
      if (is.null(extra$strata)) {
        dat <- data.frame(Time = Time, d = d, as.data.frame(WW))
        rn <- tapply(row.names(dat), idT, tail, 1)
        ind <- row.names(dat) %in% rn
        dat <- dat[ind, ]
        init.fit <- survreg(Surv(Time, d) ~ ., data = dat)
        coefs <- init.fit$coef
        xi <- 1/init.fit$scale
        phi <- exp(coefs[1])
        logh <- -log(phi * xi * dat$Time^(xi - 1))
        out$gammas.bs <- as.vector(lm.fit(extra$W2[ind, 
                                                   ], logh)$coefficients)
      }
      else {
        dat <- data.frame(Time = Time, d = d)
        dat <- cbind(dat, as.data.frame(WW))
        strata <- extra$strata
        split.dat <- split(dat, strata)
        gg <- NULL
        for (i in seq_along(split.dat)) {
          ii <- strata == levels(strata)[i]
          SpD.i <- split.dat[[i]]
          idT.i <- idT[ii]
          W2.i <- extra$W2[ii, ]
          rn <- tapply(row.names(SpD.i), idT.i, tail, 
                       1)
          ind <- row.names(SpD.i) %in% rn
          SpD.i <- SpD.i[ind, ]
          init.fit <- survreg(Surv(Time, d) ~ ., data = SpD.i)
          coefs <- init.fit$coef
          xi <- 1/init.fit$scale
          phi <- exp(coefs[1])
          logh <- -log(phi * xi * SpD.i$Time^(xi - 1))
          gg <- c(gg, as.vector(lm.fit(W2.i[ind, ], logh)$coefficients))
        }
        out$gammas.bs <- gg[!is.na(gg)]
      }
      out
    }
  }
  if (method == "weibull-AFT-GH") {
    dat <- data.frame(Time = Time, d = d)
    if (!is.null(long.id)) {
      long.id <- c(long.id) * WintF.vl
      k <- ncol(WintF.vl)
    }
    if (!is.null(longD.id)) {
      longD.id <- c(longD.id) * WintF.sl
      l <- ncol(WintF.sl)
    }
    WW <- cbind(W, long.id, longD.id)
    init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
    coefs <- -init.fit$coef
    nk <- if (is.null(W)) 
      1
    else ncol(W) + 1
    out <- switch(parameterization, 
                  value = list(gammas = coefs[1:nk], alpha = coefs[-(1:nk)], sigma.t = 1/init.fit$scale), 
                  slope = list(gammas = coefs[1:nk], Dalpha = coefs[-(1:nk)], sigma.t = 1/init.fit$scale),
                  both = list(gammas = coefs[1:nk], alpha = coefs[seq(nk + 1, nk + k)], Dalpha = coefs[-seq(1, nk + k)], sigma.t = 1/init.fit$scale))
  }
  if (method == "ch-Laplace") {
    dat <- data.frame(Time = Time, d = d)
    init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
    coefs <- -coef(init.fit)/init.fit$scale
    min.x <- min(logT)
    max.x <- max(logT)
    kn <- if (is.null(extra$control$knots)) {
      kk <- seq(0, 1, length.out = extra$control$lng.in.kn + 
                  2)[-c(1, extra$control$lng.in.kn + 2)]
      quantile(log(Time)[d == 1], kk, names = FALSE)
    }
    else {
      extra$control$knots
    }
    kn <- sort(c(rep(c(min.x, max.x), extra$control$ord), 
                 kn))
    W <- splineDesign(kn, log(Time), ord = extra$control$ord)
    nk <- ncol(W)
    nx <- NCOL(X)
    logH <- coefs[1] + logT/init.fit$scale
    coefs <- c(as.vector(lm.fit(W, logH)$coefficients), coefs[-1])
    out <- list(gammas = c(sort(coefs[1:nk]), if (nx > 1) coefs[seq(nk + 
                                                                      1, nk + nx - 1)] else NULL), alpha = coefs[nk + nx])
  }
  out
}

#### H.longSplinePH ####
H.longSplinePH <- 
  function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
      Ys <- as.vector(Xs %*% betas) + Zsb
      Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
      eta.s <- Ws.intF.vl.alph * Ys
    }
    if (parameterization %in% c("slope", "both")) {
      Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
        Zsb.deriv
      Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
      eta.s <- if (parameterization == "both") 
        eta.s + Ws.intF.sl.alph * Ys.deriv
      else Ws.intF.sl.alph * Ys.deriv
    }
    exp.eta.tw.P <- exp(eta.tw1) * P
    H1 <- XtX/sigma^2
    Int <- wk * exp(eta.ws + eta.s)
    H2 <- matrix(0, ncx, ncx)
    for (i in 1:ncx) {
      for (j in i:ncx) {
        XX <- if (parameterization == "value") {
          Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
        }
        else if (parameterization == "slope") {
          if (i %in% indFixed && j %in% indFixed) {
            ii <- match(i, indFixed)
            jj <- match(j, indFixed)
            Ws.intF.sl.alph^2 * Xs.deriv[, ii] * Xs.deriv[, 
                                                          jj]
          }
          else 0
        }
        else {
          if (i %in% indFixed && j %in% indFixed) {
            ii <- match(i, indFixed)
            jj <- match(j, indFixed)
            (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * 
               Xs.deriv[, ii]) * (Ws.intF.vl.alph * Xs[, 
                                                       j] + Ws.intF.sl.alph * Xs.deriv[, jj])
          }
          else if (i %in% indFixed && !j %in% indFixed) {
            ii <- match(i, indFixed)
            (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * 
               Xs.deriv[, ii]) * (Ws.intF.vl.alph * Xs[, 
                                                       j])
          }
          else if (!i %in% indFixed && j %in% indFixed) {
            jj <- match(j, indFixed)
            (Ws.intF.vl.alph * Xs[, i]) * (Ws.intF.vl.alph * 
                                             Xs[, j] + Ws.intF.sl.alph * Xs.deriv[, jj])
          }
          else {
            Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
          }
        }
        ki <- exp.eta.tw.P * rowsum(Int * XX, id.GK, reorder = FALSE)
        ki <- rowsum(ki, idT, reorder = FALSE)
        kii <- c((p.byt * ki) %*% wGH)
        H2[i, j] <- sum(kii, na.rm = TRUE)
      }
    }
    H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
    H1 + H2
  }

#### LogLik.splineGH ####
LogLik.splineGH <- 
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    gammas.bs <- thetas$gammas.bs
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    D <- thetas$D
    D <- if (diag.D) 
      exp(D)
    else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.tw1 <- if (!is.null(W1)) 
      as.vector(W1 %*% gammas)
    else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    if (parameterization %in% c("value", "both")) {
      Y <- as.vector(Xtime %*% betas) + Ztime.b
      Ys <- as.vector(Xs %*% betas) + Zsb
      eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * 
        Y
      eta.s <- c(Ws.intF.vl %*% alpha) * Ys
    }
    if (parameterization %in% c("slope", "both")) {
      Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + 
        Ztime.b.deriv
      Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
        Zsb.deriv
      eta.t <- if (parameterization == "both") 
        eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
      else eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * Y.deriv
      eta.s <- if (parameterization == "both") 
        eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv
      else c(Ws.intF.sl %*% Dalpha) * Ys.deriv
    }
    eta.ws <- as.vector(W2s %*% gammas.bs)
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE) # Contribution à la logvrais partie longit pour chaque ID
    log.p.yb <- rowsum(logNorm, id) # logvrais observée partie longit
    log.hazard <- eta.t # log(h) pour chaque ID, où h est le risque instantané
    log.survival <- -exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + 
                                                          eta.s), id.GK, reorder = FALSE)# log(S) pour chaque ID, où S est la survie
    dimnames(log.survival) <- NULL
    log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE) # logvrais observée partie survie
    log.p.b <- if (control$typeGH == "simple") { # Nous, on utilise control$typeGH == "adaptative"
      rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
    }
    else { # On utilise ceci :
      matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, 
                     TRUE), n, k, byrow = TRUE) # logvrais observée partie effets aléatoires
    }
    p.ytb <- exp(log.p.yb + log.p.tb + log.p.b) # vrais observée du modèle conjoint = exp(logvrais_longit|random + logvrais_surv|random + logvrais_random)
    if (control$typeGH != "simple") 
      p.ytb <- p.ytb * VCdets
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    log.p.yt <- log(p.yt)
    -sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
  }

#### Score.splineGH ####
Score.splineGH <- function (thetas) {
  thetas <- relist(thetas, skeleton = list.thetas)
  betas <- thetas$betas
  sigma <- exp(thetas$log.sigma)
  gammas <- thetas$gammas
  gammas.bs <- thetas$gammas.bs
  alpha <- thetas$alpha
  Dalpha <- thetas$Dalpha
  D <- thetas$D
  D <- if (diag.D) 
    exp(D)
  else chol.transf(D)
  eta.yx <- as.vector(X %*% betas)
  eta.tw1 <- if (!is.null(W1)) 
    as.vector(W1 %*% gammas)
  else rep(0, n)
  eta.tw2 <- as.vector(W2 %*% gammas.bs)
  if (parameterization %in% c("value", "both")) {
    Y <- as.vector(Xtime %*% betas) + Ztime.b
    Ys <- as.vector(Xs %*% betas) + Zsb
    WintF.vl.alph <- c(WintF.vl %*% alpha)
    Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
    eta.t <- eta.tw2 + eta.tw1 + WintF.vl.alph * Y
    eta.s <- Ws.intF.vl.alph * Ys
  }
  if (parameterization %in% c("slope", "both")) {
    Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + 
      Ztime.b.deriv
    Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
      Zsb.deriv
    WintF.sl.alph <- c(WintF.sl %*% Dalpha)
    Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
    eta.t <- if (parameterization == "both") 
      eta.t + WintF.sl.alph * Y.deriv
    else eta.tw2 + eta.tw1 + WintF.sl.alph * Y.deriv
    eta.s <- if (parameterization == "both") 
      eta.s + Ws.intF.sl.alph * Ys.deriv
    else Ws.intF.sl.alph * Ys.deriv
  }
  eta.ws <- as.vector(W2s %*% gammas.bs)
  exp.eta.tw.P <- exp(eta.tw1) * P
  mu.y <- eta.yx + Ztb
  logNorm <- dnorm(y, mu.y, sigma, TRUE)
  log.p.yb <- rowsum(logNorm, id)
  log.hazard <- eta.t
  Int <- wk * exp(eta.ws + eta.s)
  log.survival <- -exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
  dimnames(log.survival) <- NULL
  log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
  log.p.b <- if (control$typeGH == "simple") {
    rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
  }
  else {
    matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, 
                   TRUE), n, k, byrow = TRUE)
  }
  p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
  if (control$typeGH != "simple") 
    p.ytb <- p.ytb * VCdets
  dimnames(p.ytb) <- NULL
  p.yt <- c(p.ytb %*% wGH)
  p.byt <- p.ytb/p.yt
  post.b <- if (control$typeGH == "simple") {
    p.byt %*% (b * wGH)
  }
  else {
    sapply(seq_len(ncz), function(i) (p.byt * t(sapply(lis.b, 
                                                       "[", seq_len(k), i))) %*% wGH)
  }
  post.vb <- if (control$typeGH == "simple") {
    if (ncz == 1) {
      c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
    }
    else {
      (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function(x) x %o% 
                                         x))
    }
  }
  else {
    dd <- sapply(seq_len(ncz^2), function(i) (p.byt * t(sapply(lis.b2, 
                                                               "[", seq_len(k), i))) %*% wGH)
    bb <- apply(post.b, 1, function(x) x %o% x)
    dd - if (ncz == 1) 
      c(bb)
    else t(bb)
  }
  Zb <- if (ncz == 1) 
    post.b[id]
  else rowSums(Z * post.b[id, ], na.rm = TRUE)
  mu <- y - eta.yx
  tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
  sc1 <- -crossprod(X, y - eta.yx - Zb)/sigma^2
  sc2 <- numeric(ncx)
  for (i in 1:ncx) {
    ki <- exp.eta.tw.P * switch(parameterization, 
                                value = rowsum(Int * Ws.intF.vl.alph * Xs[, i], id.GK, reorder = FALSE), 
                                slope = {
                                  ii <- match(i, indFixed)
                                  if (is.na(ii)) 0 else rowsum(Int * Ws.intF.sl.alph * 
                                                                 Xs.deriv[, ii], id.GK, reorder = FALSE)
                                }, 
                                both = {
                                  ii <- match(i, indFixed)
                                  rowsum(Int * (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * 
                                                  if (is.na(ii)) 0 else Xs.deriv[, ii]), id.GK, 
                                         reorder = FALSE)
                                })
    ki <- c(rowsum(ki, idT, reorder = FALSE))
    kii <- c((p.byt * ki) %*% wGH)
    sc2[i] <- switch(parameterization, value = {
      ddd <- tapply(d * WintF.vl.alph * Xtime[, i], idT, sum)
      -sum(ddd - kii, na.rm = TRUE)
    }, slope = {
      ii <- match(i, indFixed)
      if (is.na(ii)) 0 else {
        ddd <- tapply(d * WintF.sl.alph * Xtime.deriv[, 
                                                      ii], idT, sum)
        -sum(ddd - kii, na.rm = TRUE)
      }
    }, both = {
      ii <- match(i, indFixed)
      ddd <- tapply(d * (WintF.vl.alph * Xtime[, i] + WintF.sl.alph * 
                           if (is.na(ii)) 0 else Xtime.deriv[, ii]), idT, 
                    sum)
      -sum(ddd - kii, na.rm = TRUE)
    })
  }
  score.y <- c(sc1 + sc2, -sigma * (-N/sigma + drop(crossprod(mu, 
                                                              mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb)/sigma^3))
  scgammas1 <- if (!is.null(W1)) {
    scg1 <- numeric(ncol(W1))
    for (jj in seq_along(scg1)) {
      tt <- rowsum(W1[, jj] * log.survival, idT, reorder = FALSE)
      scg1[jj] <- sum(c((p.byt * tt) %*% wGH), na.rm = TRUE)
    }
    -colSums(W1 * d, na.rm = TRUE) - scg1
  }
  else NULL
  scgammas2 <- numeric(nk)
  for (i in 1:nk) {
    kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
    kk <- rowsum(kk, idT, reorder = FALSE)
    scgammas2[i] <- -sum(W2[, i] * d) + sum(c((p.byt * kk) %*% 
                                                wGH))
  }
  scalpha <- if (parameterization %in% c("value", "both")) {
    rr <- numeric(ncol(WintF.vl))
    for (l in seq_along(rr)) {
      rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl[, l] * 
                                     Ys, id.GK, reorder = FALSE)
      rrr <- rowsum(rrr, idT, reorder = FALSE)
      rr[l] <- -sum((p.byt * (rowsum(d * WintF.vl[, l] * 
                                       Y, idT, reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
    }
    rr
  }
  else NULL
  scalpha.D <- if (parameterization %in% c("slope", "both")) {
    rr <- numeric(ncol(WintF.sl))
    for (l in seq_along(rr)) {
      rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.sl[, l] * 
                                     Ys.deriv, id.GK, reorder = FALSE)
      rrr <- rowsum(rrr, idT, reorder = FALSE)
      rr[l] <- -sum((p.byt * (rowsum(d * WintF.sl[, l] * 
                                       Y.deriv, idT, reorder = FALSE) - rrr)) %*% wGH, 
                    na.rm = TRUE)
    }
    rr
  }
  else NULL
  score.t <- c(scgammas1, scalpha, scalpha.D, scgammas2)
  score.b <- if (diag.D) {
    svD <- 1/D
    svD2 <- svD^2
    cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
    dim(cS.postVB) <- c(ncz, ncz)
    D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(post.b^2), 
                                                          na.rm = TRUE) * svD2)
  }
  else {
    svD <- solve(D)
    dD <- deriv.D(D)
    ndD <- length(dD)
    D1 <- sapply(dD, function(x) sum(svD * x))
    D2 <- t(sapply(dD, function(x) c(svD %*% x %*% svD)))
    cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
    out <- numeric(ndD)
    for (i in seq_along(dD)) {
      D.mat <- D2[i, ]
      dim(D.mat) <- c(ncz, ncz)
      out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + 
        sum((post.b %*% D.mat) * post.b, na.rm = TRUE)
    }
    J <- jacobian2(attr(D, "L"), ncz)
    drop(0.5 * (n * D1 - out) %*% J)
  }
  c(score.y, score.t, score.b)
}

#### chol.transf ####
chol.transf <- 
  function (x) {
    if (any(is.na(x) | !is.finite(x))) 
      stop("NA or infinite values in 'x'.\n")
    if (is.matrix(x)) {
      k <- nrow(x)
      U <- chol(x)
      U[cbind(1:k, 1:k)] <- log(U[cbind(1:k, 1:k)])
      U[upper.tri(U, TRUE)]
    }
    else {
      nx <- length(x)
      k <- round((-1 + sqrt(1 + 8 * nx))/2)
      mat <- matrix(0, k, k)
      mat[upper.tri(mat, TRUE)] <- x
      mat[cbind(1:k, 1:k)] <- exp(mat[cbind(1:k, 1:k)])
      res <- crossprod(mat)
      attr(res, "L") <- t(mat)[lower.tri(mat, TRUE)]
      res
    }
  }

#### deriv.D ####
deriv.D <- function (D) {
  ncz <- nrow(D)
  ind <- which(lower.tri(D, TRUE), arr.ind = TRUE)
  dimnames(ind) <- NULL
  nind <- nrow(ind)
  svD <- solve(D)
  lapply(1:nind, function(x, ind) {
    mat <- matrix(0, ncz, ncz)
    ii <- ind[x, , drop = FALSE]
    mat[ii[1], ii[2]] <- mat[ii[2], ii[1]] <- 1
    mat
  }, ind = ind[, 2:1])
}

#### dmvnorm ####
dmvnorm <- function (x, mu, Sigma, log = FALSE) {
  if (!is.matrix(x)) 
    x <- rbind(x)
  p <- length(mu)
  if (p == 1) {
    dnorm(x, mu, sqrt(Sigma), log = log)
  }
  else {
    t1 <- length(mu) == length(Sigma)
    t2 <- all(abs(Sigma[lower.tri(Sigma)]) < sqrt(.Machine$double.eps))
    if (t1 || t2) {
      if (!t1) 
        Sigma <- diag(Sigma)
      nx <- nrow(x)
      ff <- rowSums(dnorm(x, rep(mu, each = nx), sd = rep(sqrt(Sigma), 
                                                          each = nx), log = TRUE))
      if (log) 
        ff
      else exp(ff)
    }
    else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
      ss <- x - rep(mu, each = nrow(x))
      inv.Sigma <- evec %*% (t(evec)/ev)
      quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
      fact <- -0.5 * (p * log(2 * pi) + sum(log(ev)))
      if (log) 
        as.vector(fact - quad)
      else as.vector(exp(fact - quad))
    }
  }
}

#### dropAttr ####
dropAttr <- 
  function (mat) {
    d <- dim(mat)
    mat <- as.vector(mat)
    dim(mat) <- d
    mat
  }

#### fd.vec ####
fd.vec <- function (x, f, ..., eps = 1e-05) {
  n <- length(x)
  res <- matrix(0, n, n)
  ex <- pmax(abs(x), 1)
  f0 <- f(x, ...)
  for (i in 1:n) {
    x1 <- x
    x1[i] <- x[i] + eps * ex[i]
    diff.f <- c(f(x1, ...) - f0)
    diff.x <- x1[i] - x[i]
    res[, i] <- diff.f/diff.x
  }
  0.5 * (res + t(res))
}
  

#### gauher ####
gauher <- function (n) {
  m <- trunc((n + 1)/2)
  x <- w <- rep(-1, n)
  for (i in seq_len(m)) {
    z <- if (i == 1) {
      sqrt(2 * n + 1) - 1.85575 * (2 * n + 1)^(-0.16667)
    }
    else if (i == 2) {
      z - 1.14 * n^0.426/z
    }
    else if (i == 3) {
      1.86 * z - 0.86 * x[1]
    }
    else if (i == 4) {
      1.91 * z - 0.91 * x[2]
    }
    else {
      2 * z - x[i - 2]
    }
    for (its in seq_len(10)) {
      p1 <- 0.751125544464943
      p2 <- 0
      for (j in seq_len(n)) {
        p3 <- p2
        p2 <- p1
        p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * 
          p3
      }
      pp <- sqrt(2 * n) * p2
      z1 <- z
      z <- z1 - p1/pp
      if (abs(z - z1) <= 3e-14) 
        break
    }
    x[i] <- z
    x[n + 1 - i] <- -z
    w[i] <- 2/(pp * pp)
    w[n + 1 - i] <- w[i]
  }
  list(x = x, w = w)
}

#### gaussKronrod ####
gaussKronrod <- function (k = 15) {
  sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397, 
          0, 0.405845151377397, 0.741531185599394, 0.949107912342758, 
          -0.991455371120813, -0.864864423359769, -0.586087235467691, 
          -0.207784955007898, 0.207784955007898, 0.586087235467691, 
          0.864864423359769, 0.991455371120813)
  wk15 <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785, 
            0.209482141084728, 0.190350578064785, 0.140653259715526, 
            0.0630920926299786, 0.0229353220105292, 0.10479001032225, 
            0.169004726639268, 0.204432940075299, 0.204432940075299, 
            0.169004726639268, 0.10479001032225, 0.0229353220105292)
  wk7 <- c(0.12948496616887, 0.279705391489277, 0.381830050505119, 
           0.417959183673469, 0.381830050505119, 0.279705391489277, 
           0.12948496616887)
  if (k == 7) 
    list(sk = sk[1:7], wk = wk7)
  else list(sk = sk, wk = wk15)
}

#### gr.longSplinePH ####
gr.longSplinePH <- function (betas) {
  eta.yx <- as.vector(X %*% betas)
  if (parameterization %in% c("value", "both")) {
    Ys <- as.vector(Xs %*% betas) + Zsb
    WintF.vl.alph <- c(WintF.vl %*% alpha)
    Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
    eta.s <- Ws.intF.vl.alph * Ys
  }
  if (parameterization %in% c("slope", "both")) {
    Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
      Zsb.deriv
    WintF.sl.alph <- c(WintF.sl %*% Dalpha)
    Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
    eta.s <- if (parameterization == "both") 
      eta.s + Ws.intF.sl.alph * Ys.deriv
    else Ws.intF.sl.alph * Ys.deriv
  }
  exp.eta.tw.P <- exp(eta.tw1) * P
  sc1 <- -crossprod(X, y - eta.yx - Zb)/sigma^2
  Int <- wk * exp(eta.ws + eta.s)
  sc2 <- numeric(ncx)
  for (i in 1:ncx) {
    ki <- exp.eta.tw.P * switch(parameterization, 
                                value = rowsum(Int * Ws.intF.vl.alph * Xs[, i], id.GK, reorder = FALSE), 
                                slope = {
                                  ii <- match(i, indFixed)
                                  if (is.na(ii)) 0 else rowsum(Int * Ws.intF.sl.alph * 
                                                                 Xs.deriv[, ii], id.GK, reorder = FALSE)
                                }, 
                                both = {
                                  ii <- match(i, indFixed)
                                  rowsum(Int * (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * 
                                                  if (is.na(ii)) 0 else Xs.deriv[, ii]), id.GK, 
                                         reorder = FALSE)
                                })
    ki <- c(rowsum(ki, idT, reorder = FALSE))
    kii <- c((p.byt * ki) %*% wGH)
    sc2[i] <- switch(parameterization, value = {
      ddd <- tapply(d * WintF.vl.alph * Xtime[, i], idT, 
                    sum)
      -sum(ddd - kii, na.rm = TRUE)
    }, slope = {
      ii <- match(i, indFixed)
      if (is.na(ii)) 0 else {
        ddd <- tapply(d * WintF.sl.alph * Xtime.deriv[, 
                                                      ii], idT, sum)
        -sum(ddd - kii, na.rm = TRUE)
      }
    }, both = {
      ii <- match(i, indFixed)
      ddd <- tapply(d * (WintF.vl.alph * Xtime[, i] + WintF.sl.alph * 
                           if (is.na(ii)) 0 else Xtime.deriv[, ii]), idT, 
                    sum)
      -sum(ddd - kii, na.rm = TRUE)
    })
  }
  c(sc1 + sc2)
}

#### gr.survSplinePH ####
gr.survSplinePH <- function (thetas) {
  thetas <- relist(thetas, skeleton = list.thetas)
  gammas <- thetas$gammas
  alpha <- thetas$alpha
  Dalpha <- thetas$Dalpha
  gammas.bs <- thetas$gammas.bs
  eta.tw1 <- if (!is.null(W1)) 
    as.vector(W1 %*% gammas)
  else rep(0, n)
  eta.tw2 <- as.vector(W2 %*% gammas.bs)
  eta.t <- switch(parameterization, 
                  value = eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y,
                  slope = eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * Y.deriv,
                  both = eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y + c(WintF.sl %*% Dalpha) * Y.deriv)
  eta.s <- switch(parameterization,
                  value = c(Ws.intF.vl %*% alpha) * Ys,
                  slope = c(Ws.intF.sl %*% Dalpha) * Ys.deriv, 
                  both = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
  eta.ws <- as.vector(W2s %*% gammas.bs)
  exp.eta.tw.P <- exp(eta.tw1) * P
  Int <- wk * exp(eta.ws + eta.s)
  scgammas1 <- if (!is.null(W1)) {
    ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
    scg1 <- numeric(ncol(W1))
    for (jj in seq_along(scg1)) {
      tt <- rowsum(W1[, jj] * ki, idT, reorder = FALSE)
      scg1[jj] <- sum(c((p.byt * tt) %*% wGH), na.rm = TRUE)
    }
    -colSums(W1 * d, na.rm = TRUE) + scg1
  }
  else NULL
  scgammas2 <- numeric(nk)
  for (i in 1:nk) {
    kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
    kk <- rowsum(kk, idT, reorder = FALSE)
    scgammas2[i] <- -sum(W2[, i] * d) + sum(c((p.byt * kk) %*% 
                                                wGH))
  }
  scalpha <- if (parameterization %in% c("value", "both")) {
    rr <- numeric(ncol(WintF.vl))
    for (k in seq_along(rr)) {
      rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl[, k] * 
                                     Ys, id.GK, reorder = FALSE)
      rrr <- rowsum(rrr, idT, reorder = FALSE)
      rr[k] <- -sum((p.byt * (rowsum(d * WintF.vl[, k] * 
                                       Y, idT, reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
    }
    rr
  }
  else NULL
  scalpha.D <- if (parameterization %in% c("slope", "both")) {
    rr <- numeric(ncol(WintF.sl))
    for (k in seq_along(rr)) {
      rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.sl[, k] * 
                                     Ys.deriv, id.GK, reorder = FALSE)
      rrr <- rowsum(rrr, idT, reorder = FALSE)
      rr[k] <- -sum((p.byt * (rowsum(d * WintF.sl[, k] * 
                                       Y.deriv, idT, reorder = FALSE) - rrr)) %*% wGH, 
                    na.rm = TRUE)
    }
    rr
  }
  else NULL
  c(scgammas1, scalpha, scalpha.D, scgammas2)
}

#### jacobian2 ####
jacobian2 <- function (L, ncz) {
  ind <- which(lower.tri(matrix(0, ncz, ncz), TRUE), arr.ind = TRUE)
  dimnames(ind) <- NULL
  nind <- nrow(ind)
  id <- 1:nind
  rind <- which(ind[, 1] == ind[, 2])
  lind <- vector("list", length(rind))
  for (i in seq_along(rind)) {
    tt <- matrix(0, ncz - i + 1, ncz - i + 1)
    tt[lower.tri(tt, TRUE)] <- seq(rind[i], nind)
    tt <- tt + t(tt)
    diag(tt) <- diag(tt)/2
    lind[[i]] <- tt
  }
  out <- matrix(0, nind, nind)
  for (g in 1:ncz) {
    gind <- id[g == ind[, 2]]
    vals <- L[gind]
    for (j in gind) {
      k <- which(j == gind)
      out[cbind(lind[[g]][k, ], j)] <- if (j %in% rind) 
        vals[1] * vals
      else vals
    }
  }
  out[rind, ] <- 2 * out[rind, ]
  col.ind <- matrix(0, ncz, ncz)
  col.ind[lower.tri(col.ind, TRUE)] <- seq(1, length(L))
  col.ind <- t(col.ind)
  out[, col.ind[upper.tri(col.ind, TRUE)]]
}

#### nearPD ####
nearPD <- function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, 
                    maxits = 100) {
  if (!(is.numeric(M) && is.matrix(M) && identical(M, t(M)))) 
    stop("Input matrix M must be square and symmetric.\n")
  inorm <- function(x) max(rowSums(abs(x)))
  n <- ncol(M)
  U <- matrix(0, n, n)
  X <- M
  iter <- 0
  converged <- FALSE
  while (iter < maxits && !converged) {
    Y <- X
    T <- Y - U
    e <- eigen(Y, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    D <- if (length(d) > 1) 
      diag(d)
    else as.matrix(d)
    p <- (d > eig.tol * d[1])
    QQ <- Q[, p, drop = FALSE]
    X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
    U <- X - T
    X <- (X + t(X))/2
    conv <- inorm(Y - X)/inorm(Y)
    iter <- iter + 1
    converged <- conv <= conv.tol
  }
  X <- (X + t(X))/2
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
    d[d < Eps] <- Eps
    Q <- e$vectors
    o.diag <- diag(X)
    X <- Q %*% (d * t(Q))
    D <- sqrt(pmax(Eps, o.diag)/diag(X))
    X[] <- D * X * rep(D, each = n)
  }
  (X + t(X))/2
}

#### opt.longSplinePH ####
opt.longSplinePH <- 
  function (betas)   {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
      Y <- as.vector(Xtime %*% betas) + Ztime.b
      Ys <- as.vector(Xs %*% betas) + Zsb
      WintF.vl.alph <- c(WintF.vl %*% alpha)
      Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
      eta.t <- eta.tw2 + eta.tw1 + WintF.vl.alph * Y
      eta.s <- Ws.intF.vl.alph * Ys
    }
    if (parameterization %in% c("slope", "both")) {
      Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + 
        Ztime.b.deriv
      Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
        Zsb.deriv
      WintF.sl.alph <- c(WintF.sl %*% Dalpha)
      Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
      eta.t <- if (parameterization == "both") 
        eta.t + WintF.sl.alph * Y.deriv
      else eta.tw2 + eta.tw1 + WintF.sl.alph * Y.deriv
      eta.s <- if (parameterization == "both") 
        eta.s + Ws.intF.sl.alph * Ys.deriv
      else Ws.intF.sl.alph * Ys.deriv
    }
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    log.hazard <- eta.t
    log.survival <- -exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + 
                                                          eta.s), id.GK, reorder = FALSE)
    log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)
  }

#### opt.survSplinePH ####
opt.survSplinePH <- function (thetas) {
  thetas <- relist(thetas, skeleton = list.thetas)
  gammas <- thetas$gammas
  alpha <- thetas$alpha
  Dalpha <- thetas$Dalpha
  gammas.bs <- thetas$gammas.bs
  eta.tw1 <- if (!is.null(W1)) 
    as.vector(W1 %*% gammas)
  else rep(0, n)
  eta.tw2 <- as.vector(W2 %*% gammas.bs)
  eta.t <- switch(parameterization, 
                  value = eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y,
                  slope = eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * Y.deriv,
                  both = eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y + c(WintF.sl %*% Dalpha) * Y.deriv)
  eta.s <- switch(parameterization, 
                  value = c(Ws.intF.vl %*% alpha) * Ys,
                  slope = c(Ws.intF.sl %*% Dalpha) * Ys.deriv, 
                  both = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
  eta.ws <- as.vector(W2s %*% gammas.bs)
  log.hazard <- eta.t
  log.survival <- -exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + 
                                                        eta.s), id.GK, reorder = FALSE)
  dimnames(log.survival) <- NULL
  log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
  p.bytn <- p.byt * log.p.tb
  -sum(p.bytn %*% wGH, na.rm = TRUE)
}

#### splinePHGH.fit ####
splinePHGH.fit <- 
  function (x, y, id, initial.values, parameterization, derivForm, 
            control)   {
    logT <- as.vector(y$logT)
    Time <- exp(logT)
    d <- as.vector(y$d)
    strata <- y$strata
    y <- as.vector(y$y)
    X <- x$X
    Xtime <- x$Xtime
    Xs <- x$Xs
    Xtime.deriv <- x$Xtime.deriv
    Xs.deriv <- x$Xs.deriv
    Z <- x$Z
    Ztime <- x$Ztime
    Zs <- x$Zs
    Ztime.deriv <- x$Ztime.deriv
    Zs.deriv <- x$Zs.deriv
    W1 <- x$W
    if (!is.null(W1)) 
      rownames(W1) <- NULL
    W2 <- x$W2
    W2s <- x$W2s
    WW <- if (is.null(W1)) 
      W2
    else cbind(W2, W1)
    WintF.vl <- x$WintF.vl
    WintF.sl <- x$WintF.sl
    Ws.intF.vl <- x$Ws.intF.vl
    Ws.intF.sl <- x$Ws.intF.sl
    X <- dropAttr(X)
    Z <- dropAttr(Z)
    WW <- dropAttr(WW)
    WintF.vl <- dropAttr(WintF.vl)
    WintF.sl <- dropAttr(WintF.sl)
    Ws.intF.vl <- dropAttr(Ws.intF.vl)
    Ws.intF.sl <- dropAttr(Ws.intF.sl)
    if (parameterization == "value") {
      Xtime <- dropAttr(Xtime)
      Ztime <- dropAttr(Ztime)
      Xs <- dropAttr(Xs)
      Zs <- dropAttr(Zs)
    }
    else if (parameterization == "slope") {
      Xtime.deriv <- dropAttr(Xtime.deriv)
      Ztime.deriv <- dropAttr(Ztime.deriv)
      Xs.deriv <- dropAttr(Xs.deriv)
      Zs.deriv <- dropAttr(Zs.deriv)
    }
    else {
      Xtime <- dropAttr(Xtime)
      Ztime <- dropAttr(Ztime)
      Xs <- dropAttr(Xs)
      Zs <- dropAttr(Zs)
      Xtime.deriv <- dropAttr(Xtime.deriv)
      Ztime.deriv <- dropAttr(Ztime.deriv)
      Xs.deriv <- dropAttr(Xs.deriv)
      Zs.deriv <- dropAttr(Zs.deriv)
    }
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(WW)
    nk <- ncol(W2)
    N <- length(y)
    ni <- as.vector(tapply(id, id, length))
    n <- length(ni)
    nRisks <- x$nRisks
    CompRisk <- nRisks > 1
    idT <- x$idT
    XtX <- crossprod(X)
    ZtZ <- lapply(split(Z, id), function(x) crossprod(matrix(x, 
                                                             ncol = ncz)))
    names(ZtZ) <- NULL
    ZtZ <- matrix(unlist(ZtZ), n, ncz * ncz, TRUE)
    GH <- gauher(control$GHk)
    b <- as.matrix(expand.grid(rep(list(GH$x), ncz)))
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))
    wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * 
                                                           b))
    if (control$typeGH == "simple") {
      b <- sqrt(2) * t(control$inv.chol.VC %*% t(b))
      wGH <- wGH * control$det.inv.chol.VC
    }
    else {
      b <- sqrt(2) * b
      VCdets <- control$det.inv.chol.VCs
    }
    dimnames(b) <- NULL
    b2 <- if (ncz == 1) 
      b * b
    else t(apply(b, 1, function(x) x %o% x))
    Ztb <- Z %*% t(b)
    if (parameterization %in% c("value", "both")) {
      Ztime.b <- Ztime %*% t(b)
      Zsb <- Zs %*% t(b)
    }
    if (parameterization %in% c("slope", "both")) {
      if (length(indRandom) > 1 || indRandom) {
        Ztime.b.deriv <- Ztime.deriv %*% t(b[, indRandom, 
                                             drop = FALSE])
        Zsb.deriv <- Zs.deriv %*% t(b[, indRandom, drop = FALSE])
      }
      else {
        Ztime.b.deriv <- matrix(0, nrow(Ztime.deriv), k)
        Zsb.deriv <- matrix(0, nrow(Zs.deriv), k)
      }
    }
    wk <- rep(x$wk, length(logT))
    P <- as.vector(x$P)
    id.GK <- rep(seq_along(logT), each = control$GKk)
    if (control$typeGH != "simple") {
      lis.b <- vector("list", n)
      for (i in 1:n) {
        lis.b[[i]] <- t(control$inv.chol.VCs[[i]] %*% t(b)) + 
          rep(control$ranef[i, ], each = k)
        Ztb[id == i, ] <- Z[id == i, , drop = FALSE] %*% 
          t(lis.b[[i]])
      }
      lis.b2 <- lapply(lis.b, function(b) if (ncz == 1) 
        b * b
        else t(apply(b, 1, function(x) x %o% x)))
      for (i in seq_along(logT)) {
        if (parameterization %in% c("value", "both")) {
          bb <- t(lis.b[[idT[i]]])
          Ztime.b[i, ] <- Ztime[i, , drop = FALSE] %*% 
            bb
          Zsb[id.GK == i, ] <- Zs[id.GK == i, ] %*% bb
        }
        if (parameterization %in% c("slope", "both") && (length(indRandom) > 
                                                           1 || indRandom)) {
          bb <- t(lis.b[[idT[i]]][, indRandom, drop = FALSE])
          Ztime.b.deriv[i, ] <- Ztime.deriv[i, , drop = FALSE] %*% 
            bb
          Zsb.deriv[id.GK == i, ] <- Zs.deriv[id.GK == 
                                                i, ] %*% bb
        }
      }
    }
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- if (!is.null(W1)) 
      as.vector(initial.values$gammas)
    else NULL
    gammas.bs <- as.vector(initial.values$gammas.bs)
    alpha <- as.vector(initial.values$alpha)
    Dalpha <- as.vector(initial.values$Dalpha)
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) 
      dimnames(D) <- NULL
    else names(D) <- NULL
    environment(opt.survSplinePH) <- environment(gr.survSplinePH) <- environment()
    environment(opt.longSplinePH) <- environment(gr.longSplinePH) <- environment(H.longSplinePH) <- environment()
    environment(LogLik.splineGH) <- environment(Score.splineGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter + 1, ncx + 1)
    T.mat <- matrix(0, iter + 1, switch(parameterization, value = ncww + 
                                          ncol(WintF.vl), slope = ncww + ncol(WintF.sl), both = ncww + 
                                          ncol(WintF.vl) + ncol(WintF.sl)))
    B.mat <- if (diag.D) 
      matrix(0, iter + 1, ncz)
    else matrix(0, iter + 1, ncz * ncz)
    lgLik <- numeric(iter + 1)
    conv <- TRUE
    for (it in 1:iter) {
      Y.mat[it, ] <- c(betas, sigma)
      T.mat[it, ] <- switch(parameterization, value = c(gammas, 
                                                        alpha, gammas.bs), slope = c(gammas, Dalpha, gammas.bs), 
                            both = c(gammas, alpha, Dalpha, gammas.bs))
      B.mat[it, ] <- D
      eta.yx <- as.vector(X %*% betas)
      eta.tw1 <- if (!is.null(W1)) 
        as.vector(W1 %*% gammas)
      else rep(0, length(logT))
      eta.tw2 <- as.vector(W2 %*% gammas.bs)
      eta.ws <- as.vector(W2s %*% gammas.bs)
      if (parameterization %in% c("value", "both")) {
        Y <- as.vector(Xtime %*% betas) + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * 
          Y
        eta.s <- c(Ws.intF.vl %*% alpha) * Ys
      }
      if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + 
          Ztime.b.deriv
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
          Zsb.deriv
        eta.t <- if (parameterization == "both") 
          eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
        else eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * 
          Y.deriv
        eta.s <- if (parameterization == "both") 
          eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv
        else c(Ws.intF.sl %*% Dalpha) * Ys.deriv
      }
      mu.y <- eta.yx + Ztb
      logNorm <- dnorm(y, mu.y, sigma, TRUE)
      log.p.yb <- rowsum(logNorm, id, reorder = FALSE)
      dimnames(log.p.yb) <- NULL
      log.hazard <- eta.t
      log.survival <- -exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + 
                                                            eta.s), id.GK, reorder = FALSE)
      dimnames(log.survival) <- NULL
      log.p.tb <- rowsum(d * log.hazard + log.survival, idT, 
                         reorder = FALSE)
      log.p.b <- if (control$typeGH == "simple") {
        rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
      }
      else {
        matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), 
                       D, TRUE), n, k, byrow = TRUE)
      }
      p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
      if (control$typeGH != "simple") 
        p.ytb <- p.ytb * VCdets
      p.yt <- c(p.ytb %*% wGH)
      p.byt <- p.ytb/p.yt
      post.b <- if (control$typeGH == "simple") {
        p.byt %*% (b * wGH)
      }
      else {
        sapply(seq_len(ncz), function(i) (p.byt * t(sapply(lis.b, 
                                                           "[", seq_len(k), i))) %*% wGH)
      }
      post.vb <- if (control$typeGH == "simple") {
        if (ncz == 1) {
          c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
        }
        else {
          (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function(x) x %o% 
                                             x))
        }
      }
      else {
        dd <- sapply(seq_len(ncz^2), function(i) (p.byt * 
                                                    t(sapply(lis.b2, "[", seq_len(k), i))) %*% wGH)
        bb <- apply(post.b, 1, function(x) x %o% x)
        dd - if (ncz == 1) 
          c(bb)
        else t(bb)
      }
      log.p.yt <- log(p.yt)
      lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
      if (control$verbose) {
        cat("\n\niter:", it, "\n")
        cat("log-likelihood:", lgLik[it], "\n")
        cat("betas:", round(betas, 4), "\n")
        cat("sigma:", round(sigma, 4), "\n")
        if (!is.null(W1)) 
          cat("gammas:", round(gammas, 4), "\n")
        if (parameterization %in% c("value", "both")) 
          cat("alpha:", round(alpha, 4), "\n")
        if (parameterization %in% c("slope", "both")) 
          cat("Dalpha:", round(Dalpha, 4), "\n")
        cat("gammas.bs:", round(gammas.bs, 4), "\n")
        cat("D:", if (!diag.D) 
          round(D[lower.tri(D, TRUE)], 4)
          else round(D, 4), "\n")
      }
      if (it > 5 && lgLik[it] > lgLik[it - 1]) {
        thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 
                                                              1, ])
        thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
        check1 <- max(abs(thets2 - thets1)/(abs(thets1) + 
                                              control$tol1)) < control$tol2
        check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * 
          (abs(lgLik[it - 1]) + control$tol3)
        if (check1 || check2) {
          conv <- FALSE
          if (control$verbose) 
            cat("\n\nconverged!\ncalculating Hessian...\n")
          break
        }
      }
      if (iter == 0) 
        break
      Zb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
      mu <- y - eta.yx
      tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
      sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + 
                         tr.tZZvarb)/N)
      Dn <- if (control$typeGH == "simple") {
        matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), 
               ncz, ncz)
      }
      else {
        matrix(colMeans(dd, na.rm = TRUE), ncz, ncz)
      }
      Dn <- if (diag.D) 
        diag(Dn)
      else 0.5 * (Dn + t(Dn))
      Hbetas <- nearPD(H.longSplinePH(betas))
      scbetas <- gr.longSplinePH(betas)
      betasn <- betas - c(solve(Hbetas, scbetas))
      list.thetas <- list(gammas = gammas, alpha = alpha, Dalpha = Dalpha, 
                          gammas.bs = gammas.bs)
      list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
      thetas <- unlist(as.relistable(list.thetas))
      optz.surv <- optim(thetas, opt.survSplinePH, gr.survSplinePH, 
                         method = "BFGS", control = list(maxit = if (it < 
                                                                       5) 20 else 4, parscale = if (it < 5) rep(0.01, 
                                                                                                                length(thetas)) else rep(0.1, length(thetas))))
      betas <- betasn
      sigma <- sigman
      D <- Dn
      thetasn <- relist(optz.surv$par, skeleton = list.thetas)
      gammas <- if (!is.null(W1)) 
        thetasn$gammas
      else NULL
      alpha <- thetasn$alpha
      Dalpha <- thetasn$Dalpha
      gammas.bs <- thetasn$gammas.bs
    }
    list.thetas <- list(betas = betas, log.sigma = log(sigma), 
                        gammas = gammas, alpha = alpha, Dalpha = Dalpha, gammas.bs = gammas.bs, 
                        D = if (diag.D) log(D) else chol.transf(D))
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    lgLik <- -LogLik.splineGH(thetas)
    if (conv && !control$only.EM) {
      if (is.null(control$parscale)) 
        control$parscale <- rep(0.01, length(thetas))
      if (control$verbose) 
        cat("\n\nquasi-Newton iterations start.\n\n")
      out <- if (control$optimizer == "optim") {
        optim(thetas, LogLik.splineGH, Score.splineGH, method = "BFGS", 
              control = list(maxit = control$iter.qN, parscale = control$parscale, 
                             trace = 10 * control$verbose))
      }
      else {
        nlminb(thetas, LogLik.splineGH, Score.splineGH, scale = control$parscale, 
               control = list(iter.max = control$iter.qN, trace = 1 * 
                                control$verbose))
      }
      if ((conv <- out$convergence) == 0 || -out[[2]] > lgLik) {
        lgLik <- -out[[2]]
        thetas <- relist(out$par, skeleton = list.thetas)
        betas <- thetas$betas
        sigma <- exp(thetas$log.sigma)
        gammas <- if (!is.null(W1)) 
          thetas$gammas
        else NULL
        gammas.bs <- thetas$gammas.bs
        alpha <- thetas$alpha
        Dalpha <- thetas$Dalpha
        D <- thetas$D
        D <- if (diag.D) 
          exp(D)
        else chol.transf(D)
        it <- it + if (control$optimizer == "optim") 
          out$counts[1]
        else out$iterations
        eta.yx <- as.vector(X %*% betas)
        eta.tw1 <- if (!is.null(W1)) 
          as.vector(W1 %*% gammas)
        else rep(0, n)
        eta.tw2 <- as.vector(W2 %*% gammas.bs)
        exp.eta.tw <- exp(eta.tw1)
        eta.ws <- as.vector(W2s %*% gammas.bs)
        if (parameterization %in% c("value", "both")) {
          Y <- as.vector(Xtime %*% betas) + Ztime.b
          Ys <- as.vector(Xs %*% betas) + Zsb
          eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * 
            Y
          eta.s <- c(Ws.intF.vl %*% alpha) * Ys
        }
        if (parameterization %in% c("slope", "both")) {
          Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + 
            Ztime.b.deriv
          Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
            Zsb.deriv
          eta.t <- if (parameterization == "both") 
            eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
          else eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * 
            Y.deriv
          eta.s <- if (parameterization == "both") 
            eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv
          else c(Ws.intF.sl %*% Dalpha) * Ys.deriv
        }
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id)
        log.hazard <- eta.t
        log.survival <- -exp.eta.tw * P * rowsum(wk * exp(eta.ws + 
                                                            eta.s), id.GK, reorder = FALSE)
        dimnames(log.survival) <- NULL
        log.p.tb <- rowsum(d * log.hazard + log.survival, 
                           idT, reorder = FALSE)
        log.p.b <- if (control$typeGH == "simple") {
          rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
        }
        else {
          matrix(dmvnorm(do.call(rbind, lis.b), rep(0, 
                                                    ncz), D, TRUE), n, k, byrow = TRUE)
        }
        p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
        if (control$typeGH != "simple") 
          p.ytb <- p.ytb * VCdets
        p.yt <- c(p.ytb %*% wGH)
        p.byt <- p.ytb/p.yt
        post.b <- if (control$typeGH == "simple") {
          p.byt %*% (b * wGH)
        }
        else {
          sapply(seq_len(ncz), function(i) (p.byt * t(sapply(lis.b, 
                                                             "[", seq_len(k), i))) %*% wGH)
        }
        post.vb <- if (control$typeGH == "simple") {
          if (ncz == 1) {
            c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
          }
          else {
            (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, 
                                             function(x) x %o% x))
          }
        }
        else {
          dd <- sapply(seq_len(ncz^2), function(i) (p.byt * 
                                                      t(sapply(lis.b2, "[", seq_len(k), i))) %*% 
                         wGH)
          bb <- apply(post.b, 1, function(x) x %o% x)
          dd - if (ncz == 1) 
            c(bb)
          else t(bb)
        }
        Zb <- if (ncz == 1) 
          post.b[id]
        else rowSums(Z * post.b[id, ], na.rm = TRUE)
        if (control$verbose) 
          cat("\n\ncalculating Hessian...\n")
      }
    }
    Score <- Score.splineGH(unlist(thetas))
    Hessian <- if (control$numeriDeriv == "fd") {
      fd.vec(unlist(thetas), Score.splineGH, eps = control$eps.Hes)
    }
    else {
      cd.vec(unlist(thetas), Score.splineGH, eps = control$eps.Hes)
    }
    names(betas) <- names(initial.values$betas)
    if (!diag.D) 
      dimnames(D) <- dimnames(initial.values$D)
    else names(D) <- names(initial.values$D)
    names(gammas) <- colnames(W1)
    names(gammas.bs) <- if (length(levels(strata)) == 1) 
      paste("bs", 1:nk, sep = "")
    else {
      len.kn <- sapply(control$knots, length) - control$ord
      paste("bs", sapply(len.kn, seq_len), "(", rep(levels(strata), 
                                                    len.kn), ")", sep = "")
    }
    nm.alph <- colnames(x$WintF.vl)
    nm.alph <- if (!is.null(nm.alph)) {
      if (nm.alph[1] == "(Intercept)") 
        c("", nm.alph[-1])
      else nm.alph
    }
    else {
      "alpha"
    }
    nm.Dalph <- colnames(x$WintF.sl)
    nm.Dalph <- if (!is.null(nm.Dalph)) {
      if (nm.Dalph[1] == "(Intercept)") 
        c("", nm.Dalph[-1])
      else nm.Dalph
    }
    else {
      "alpha.s"
    }
    gg <- switch(parameterization, value = nm.alph, slope = nm.Dalph, 
                 both = c(nm.alph, nm.Dalph))
    if (parameterization %in% c("value", "both")) 
      names(alpha) <- nm.alph
    if (parameterization %in% c("slope", "both")) 
      names(Dalpha) <- nm.Dalph
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), 
              paste("T.", c(names(gammas), gg, names(gammas.bs)), sep = ""), 
              paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 
                                                                  1)/2), sep = "") else names(D), sep = ""))
    dimnames(Hessian) <- list(nams, nams)
    colnames(post.b) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, 
                             alpha = alpha, Dalpha = Dalpha, gammas.bs = gammas.bs, 
                             D = as.matrix(D)), Score = Score, Hessian = Hessian, 
         logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb, 
                                   Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb, 
                                   Ztimeb = if (parameterization %in% c("value", "both")) {
                                     rowSums(Ztime * post.b[idT, , drop = FALSE])
                                   } else NULL, Ztimeb.deriv = if (parameterization %in% 
                                                                     c("slope", "both")) {
                                     if (indRandom) {
                                       rowSums(Ztime.deriv * post.b[idT, indRandom, 
                                                                    drop = FALSE])
                                     } else rep(0, nrow(Ztime.deriv))
                                   } else NULL), iters = it, convergence = conv, n = n, 
         N = N, ni = ni, d = d, id = id)
  }


