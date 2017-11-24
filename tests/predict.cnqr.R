## Check that predict.cnqr() gives the same results as QYgX()
library(cnqr)
library(copsupp)
library(testthat)
tau <- space_taus(10)
data(egdat)
dat <- list(tr=egdat[1:500, ],
            val=egdat[501:750, ],
            test=egdat[751:1000, ])
xord <- 4:1
yvar <- 5
edges <- c(yvar, xord)

## Fit the cnqr object
xvine <- subset(egvine, xord)
fit <- cnqr(edges, dat=dat, sc=scorer(tau),
            basevine=xvine,
            families=c("bvncop", "frk", "gum"),
            verbose=TRUE)
summary(fit)
xyinfo <- xylink(fit)
uind <- fit$uind
uind_check <- lapply(dat, pcondseq, ord=xord, rv=xvine)
expect_identical(uind, uind_check)  # Passes the test.

## 1. Check each data type.
for (dt in names(dat)) {
    yhat_cnqr <- predict(fit, dat=dt)
    yhat_QYgX <- QYgX(tau, uind[[dt]], cops=xyinfo$cops, cpars=xyinfo$cpars, QY=identity)
    expect_identical(yhat_cnqr, yhat_QYgX)
}

## 2. Check computing all together
yhat_cnqr <- predict(fit, dat=names(dat))
yhat_QYgX <- QYgX(tau, do.call(rbind, uind), cops=xyinfo$cops, cpars=xyinfo$cpars, QY=identity)
expect_identical(yhat_cnqr, yhat_QYgX)

## 3. Try new data
set.seed(345)
datnew <- rrvine(100, egvine)
uindnew <- pcondseq(datnew, ord=xord, rv=xvine)
yhat_cnqr <- predict(fit, dat=datnew)
yhat_QYgX <- QYgX(tau, uindnew, cops=xyinfo$cops, cpars=xyinfo$cpars, QY=identity)
expect_identical(yhat_cnqr, yhat_QYgX)

## 4. Try different quantile levels
tau_new <- c(0.5, 0.8, 0.9)
yhat_cnqr <- predict(fit, dat=datnew, tau=tau_new)
yhat_QYgX <- QYgX(tau_new, uindnew, cops=xyinfo$cops, cpars=xyinfo$cpars, QY=identity)
expect_identical(yhat_cnqr, yhat_QYgX)
