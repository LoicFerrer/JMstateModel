---
title: "Example of joint model for a longitudinal and a multi-state processes"
author: "Loïc Ferrer and Cécile Proust-Lima"
date: "April, Thuesday 07, 2015"
output: 
 html_document:
    fig_caption: yes
    fig_height: 4
    fig_width: 5
    highlight: tango
    theme: spacelab
    toc: yes
---



Import two databases which contain longitudinal and survival data:

```r
load("example_data.RData")
ls()
```

```
## [1] "data_long" "data_surv"
```

Load the packages and the function to estimate joint multi-state models:

```r
library(mstate)
library(JM)
source("JMstateModel.R")
```



# Longitudinal sub-part
Plot the first 500 individual trajectories of the longitudinal responses:

```r
library(ggplot2)
plot_long <- (ggplot(data_long[data_long$id <= 500, ]) +
                geom_line(aes(x = time, y = Y, group = id), color = "grey30", alpha = 0.8) +
                stat_smooth(aes(x = time, y = Y), method = "loess", size = 0.75) + 
                theme_bw() +
                xlab("Time") +
                ylab("Marker value"))
plot_long
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

Fit the longitudinal responses through a linear mixed model:

```r
lmeFit <- lme(fixed = Y ~ 1 + X + time + X:time,
              data = data_long,
              random = ~ (1 + time) | id,      
              method = "REML",
              control = list(opt = "optim"))
```


# Multi-state sub-part
Construct the matrix of transitions:

```r
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, 3] <- 3
dimnames(tmat) <- list(from = c("State_0", "State_1", "State_2"),
                       to = c("State_0", "State_1", "State_2"))
tmat
```

```
##          to
## from      State_0 State_1 State_2
##   State_0      NA       1       2
##   State_1      NA      NA       3
##   State_2      NA      NA      NA
```
The transition '0 -> 1' is called '1','0 -> 2' is called '2', '1 -> 2' is called '3'.

Define the covariate(s) in the multi-state sub-part:

```r
covs <- "X"
```
The *msprep()* function divides the survival database in order to have one line per transition at risk for each subject, with 'Tstart' the entry time in the current state, and 'Tstop' the time of transition or censorship; 'status' denotes if the transition has been performed:

```r
data_mstate <- msprep(time = c(NA, "time_of_State_1", "time_of_State_2"),
                      status = c(NA, "State_1", "State_2"),
                      data = data_surv,
                      trans = tmat,
                      keep = covs,
                      id = "id")
```

*expand.covs()* permits to define the covariates for each transition:

```r
data_mstate <- expand.covs(data_mstate, covs,
                           append = TRUE, longnames = FALSE)
head(data_mstate)
```

```
## An object of class 'msdata'
## 
## Data:
##   id from to trans   Tstart     Tstop      time status        X      X.1
## 1  1    1  2     1 0.000000  9.905201  9.905201      1 3.662981 3.662981
## 2  1    1  3     2 0.000000  9.905201  9.905201      0 3.662981 0.000000
## 3  1    2  3     3 9.905201 13.838360  3.933159      1 3.662981 0.000000
## 4  2    1  2     1 0.000000  5.139122  5.139122      0 1.192941 1.192941
## 5  2    1  3     2 0.000000  5.139122  5.139122      0 1.192941 0.000000
## 6  3    1  2     1 0.000000 10.048559 10.048559      1 1.549180 1.549180
##        X.2      X.3
## 1 0.000000 0.000000
## 2 3.662981 0.000000
## 3 0.000000 3.662981
## 4 0.000000 0.000000
## 5 1.192941 0.000000
## 6 0.000000 0.000000
```

The *events()* function indicates the number of observed transitions and their percentages.

```r
events(data_mstate)
```

```
## $Frequencies
##          to
## from      State_0 State_1 State_2 no event total entering
##   State_0       0     482     288      730           1500
##   State_1       0       0     311      171            482
##   State_2       0       0       0        0              0
## 
## $Proportions
##          to
## from        State_0   State_1   State_2  no event
##   State_0 0.0000000 0.3213333 0.1920000 0.4866667
##   State_1 0.0000000 0.0000000 0.6452282 0.3547718
##   State_2
```

Multi-state model with proportional intensities:

```r
coxFit <- coxph(Surv(Tstart, Tstop, status) ~ X.1 + X.2 + X.3 + strata(trans),
                data = data_mstate, method = "breslow", x = TRUE, model = TRUE)
```

# Joint multi-state sub-part

Define the derived of the fixed and random parts in the mixed model, and indicate which covariates are kept:

```r
dForm <- list(fixed = ~ 1 + X,
              indFixed = c(3, 4),
              random = ~ 1,
              indRandom = 2)
```

Joint multi-state model with current level and current slope of the biomarker as dependence function:

```r
jointFit_both <- JMstateModel(lmeObject = lmeFit,
                              survObject = coxFit,
                              timeVar = "time",
                              parameterization = "both",
                              method = "spline-PH-aGH",
                              interFact = list(value = ~ strata(trans) - 1,
                                               slope = ~ strata(trans) - 1,
                                               data = data_mstate),
                              derivForm = dForm,
                              Mstate = TRUE,
                              data.Mstate = data_mstate,
                              ID.Mstate = "id",
                              control = list(GHk = 9, lng.in.kn = 3))
summary(jointFit_both)
```


```
## 
## Call:
## JMstateModel(lmeObject = lmeFit, survObject = coxFit, timeVar = "time", 
##     parameterization = "both", method = "spline-PH-aGH", interFact = list(value = ~strata(trans) - 
##         1, slope = ~strata(trans) - 1, data = data_mstate), derivForm = dForm, 
##     control = list(GHk = 9, lng.in.kn = 3), Mstate = TRUE, data.Mstate = data_mstate, 
##     ID.Mstate = "id", verbose = TRUE)
## 
## Data Descriptives:
## Longitudinal Process		Event Process
## Number of Observations: 35977	Number of Events: 1081 (72.1%)
## Number of Groups: 1500
## 
## Joint Model Summary:
## Longitudinal Process: Linear mixed-effects model
## Event Process: Stratified relative risk model with spline-approximated
## 		baseline risk function
## Parameterization: Time-dependent + time-dependent slope 
## 
##    log.Lik      AIC      BIC
##  -32922.34 65920.68 66122.58
## 
## Variance Components:
##              StdDev    Corr
## (Intercept)  0.5835  (Intr)
## time         0.2461 -0.2699
## Residual     0.4790        
## 
## Coefficients:
## Longitudinal Process
##               Value Std.Err  z-value p-value
## (Intercept) -0.7792  0.0489 -15.9312 <0.0001
## X            0.5315  0.0227  23.3762 <0.0001
## time        -0.0701  0.0209  -3.3491  0.0008
## X:time       0.0172  0.0098   1.7499  0.0801
## 
## Event Process
##                                  Value Std.Err  z-value p-value
## X.1                             0.2760  0.0767   3.5989  0.0003
## X.2                            -0.0208  0.0943  -0.2207  0.8254
## X.3                            -0.2721  0.0897  -3.0339  0.0024
## Assoct:strata(trans)trans=1     0.9081  0.0725  12.5207 <0.0001
## Assoct:strata(trans)trans=2     0.2956  0.0689   4.2903 <0.0001
## Assoct:strata(trans)trans=3     0.1054  0.0702   1.5011  0.1333
## Assoct.s:strata(trans)trans=1   1.7036  0.4506   3.7808  0.0002
## Assoct.s:strata(trans)trans=2  -1.2282  0.6887  -1.7835  0.0745
## Assoct.s:strata(trans)trans=3  -0.0693  0.7814  -0.0887  0.9293
## bs1(trans=1)                   -7.5793  0.5568 -13.6111 <0.0001
## bs2(trans=1)                   -4.0473  0.4118  -9.8286 <0.0001
## bs3(trans=1)                   -4.5721  0.3413 -13.3960 <0.0001
## bs4(trans=1)                   -3.6716  0.2723 -13.4858 <0.0001
## bs5(trans=1)                   -3.4467  0.5853  -5.8892 <0.0001
## bs6(trans=1)                    0.9570  1.4974   0.6391  0.5227
## bs7(trans=1)                   -3.7318  4.7752  -0.7815  0.4345
## bs1(trans=2)                  -10.6728  2.0752  -5.1430 <0.0001
## bs2(trans=2)                   -4.6639  0.8913  -5.2326 <0.0001
## bs3(trans=2)                   -4.6580  0.5140  -9.0615 <0.0001
## bs4(trans=2)                   -3.0090  0.3158  -9.5279 <0.0001
## bs5(trans=2)                   -1.6966  0.4857  -3.4933  0.0005
## bs6(trans=2)                    0.0875  0.7518   0.1163  0.9074
## bs7(trans=2)                   -0.3098  1.2636  -0.2451  0.8064
## bs1(trans=3)                   -4.4036  4.8601  -0.9061  0.3649
## bs2(trans=3)                   -3.1485  1.3527  -2.3275  0.0199
## bs3(trans=3)                   -1.9357  0.5878  -3.2933  0.0010
## bs4(trans=3)                   -1.7222  0.3233  -5.3275 <0.0001
## bs5(trans=3)                   -0.7120  0.4510  -1.5788  0.1144
## bs6(trans=3)                    0.1252  0.6585   0.1901  0.8493
## bs7(trans=3)                    0.9051  0.8697   1.0407  0.2980
## 
## Integration:
## method: (pseudo) adaptive Gauss-Hermite
## quadrature points: 9 
## 
## Optimization:
## Convergence: 0
```
