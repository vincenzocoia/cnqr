## Test that character data inputs to xx.cnqr() functions work.

## 0. Setup
library(copsupp)
library(cnqr)
library(testthat)
set.seed(134)
d <- 3
ndat <- 3
dat <- lapply(1:ndat, function(i) matrix(runif(100*i*d), ncol=d))
names(dat) <- letters[1:length(dat)]
fit <- cnqr(1:d, dat, sc=scorer(space_taus(10)),
            families=c("indepcop", "bvncop",
                       "bvtcop", "mtcj", "gum", "frk", "joe", "bb1"))

## 1. predict.cnqr()
yhat1 <- predict(fit)
yhat2 <- predict(fit, dat=dat$a)
expect_identical(yhat1, yhat2)

yhat1 <- predict(fit, dat="c")
yhat2 <- predict(fit, dat=dat$c)
expect_identical(yhat1, yhat2)

yhat1 <- predict(fit, dat=c("tr", "c"))
yhat2 <- predict(fit, dat=rbind(dat$a, dat$c))
expect_identical(yhat1, yhat2)

## 2. calplot.cnqr()
cp1 <- calplot(fit)
cp2 <- calplot(fit, dat=dat$a)
expect_identical(cp1, cp2)

cp1 <- calplot(fit, dat="c")
cp2 <- calplot(fit, dat=dat$c)
expect_identical(cp1, cp2)

cp1 <- calplot(fit, dat=c("tr", "c"))
cp2 <- calplot(fit, dat=rbind(dat$a, dat$c))
expect_identical(cp1, cp2)

## 3. score.cnqr()
sc1 <- score(fit)
sc2 <- score(fit, dat=dat$a)
expect_identical(sc1, sc2)

sc1 <- score(fit, dat="c")
sc2 <- score(fit, dat=dat$c)
expect_identical(sc1, sc2)

sc1 <- score(fit, dat=c("tr", "c"))
sc2 <- score(fit, dat=rbind(dat$a, dat$c))
expect_identical(sc1, sc2)

