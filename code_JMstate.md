---
title: "Example of joint multi-state model"
author: "Loïc Ferrer and Cécile Proust-Lima"
date: "Friday, February 27, 2015"
output: 
 html_document:
    fig_caption: yes
    fig_height: 4
    fig_width: 5
    highlight: tango
    theme: spacelab
    toc: yes
---

```{r, eval=TRUE, echo=FALSE, message=FALSE, warning=FALSE}
library(rmarkdown)
library(knitr)
```

Import two databases which contain longitudinal and survival data:
```{r, eval=TRUE}
load("data.RData")
ls()
```

```{r, eval=TRUE}
load("data.RData")
ls()
```

Load the packages and the function to estimate joint multi-state models:
```{r, message=FALSE}
library(mstate)
library(JM)
source("JMstateModel.R")
```

```{r, eval=TRUE, echo=FALSE}
# We take the 500 first ID
data_surv <- data_surv[data_surv$id <= 500, ]
data_long <- data_long[data_long$id <= 500, ]
```

```{r, eval=TRUE, echo=FALSE}
# We call X the covariate which is the same for the two processes
colnames(data_surv)[which(colnames(data_surv) == "X_s")] <- "X"
colnames(data_long)[which(colnames(data_long) == "X_l")] <- "X"
```


# Longitudinal sub-part
Plot the individual trajectories of the longitudinal responses:
```{r, eval=TRUE, warning=FALSE}
library(ggplot2)

plot_long <- (ggplot(data_long) +
                geom_line(aes(x=times, y=Y, group=id), color="grey30", alpha=0.8) +
                stat_smooth(aes(x=times, y=Y), method = "loess", size = 0.75) + 
                theme_bw() +
                xlab("Times") +
                ylab("Marker value"))
plot_long
```

Fit the longitudinal responses through a linear mixed model:
```{r, eval=TRUE}
lmeFit <- lme(fixed = Y ~ (1 + times) * X,
              data = data_long,
              random = ~ (1 + times) | id,      
              method = "REML",
              control = list(opt = "optim"))
summary(lmeFit)
```


# Multi-state sub-part
Construct the matrix of transitions:
```{r, eval=TRUE}
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, 3] <- 3
dimnames(tmat) <- list(from = c("State 0", "State 1", "State 2"),
                       to = c("State 0", "State 1", "State 2"))
tmat
```
The transition '0 -> 1' is called '1','0 -> 2' is called '2', '1 -> 2' is called '3'.

Define the covariate(s) in the multi-state sub-part:
```{r, eval=TRUE}
covs <- "X"
```
The *msprep()* function divides the survival database in order to have one line per transition at risk for each subject, with 'Tstart' the entry time in the current state, and 'Tstop' the time of transition or censorship; 'status' denotes if the transition has been performed:
```{r, eval=TRUE}
data_mstate <- msprep(time = c(NA, "time_of_Rec", "time_of_Death"),
                      status = c(NA, "Rec", "Death"),
                      data = data_surv,
                      trans = tmat,
                      keep = covs,
                      id = "id")
```

*expand.covs()* permits to define the covariates for each transition:
```{r, eval=TRUE}
data_mstate <- expand.covs(data_mstate, covs,
                           append = TRUE, longnames = FALSE)
head(data_mstate)
```

The *events()* function indicates the number of observed transitions and their percentages.
```{r, eval=TRUE}
events(data_mstate)
```

Multi-state model with proportional hazards:
```{r, eval=TRUE}
coxFit <- coxph(Surv(Tstart, Tstop, status) ~ X.1 + X.2 + X.3 + strata(trans),
                data = data_mstate, method = "breslow", x = TRUE, model = TRUE)
```

# Joint multi-state sub-part

Define the derived of the fixed and random parts in the mixed model, and indicates which covariates are kept:
```{r, eval=TRUE}
dForm <- list(fixed = ~ 1 + X,
              indFixed = c(2, 4),
              random = ~ 1,
              indRandom = 2)
```

Joint multi-state model with current level and current slope of the biomarker as dependence function:
```{r, eval=TRUE}
jointFit <- JMstateModel(lmeObject = lmeFit,
                         survObject = coxFit,
                         timeVar = "times",
                         parameterization = "both",
                         method = "spline-PH-aGH",
                         interFact = list(value = ~strata(trans) - 1,
                                          slope = ~strata(trans) - 1,
                                          data = data_mstate),
                         derivForm = dForm,
                         Mstate = TRUE,
                         data.Mstate = data_mstate,
                         ID.Mstate = "id",
                         control = list(GHk = 3, lng.in.kn = 3, iter.EM=0, only.EM = T))
summary(jointFit)
```