## Check that the score_eval() function works.
library(testthat)
library(cnqr)
n <- 10
set.seed(364)
y <- rnorm(n)
tau <- c(0.8, 0.9)
g <- exp
wtau <- function(p) 1-p
sc <- scorer(tau, g=g, wtau=wtau)
yhat <- matrix(qnorm(tau), nrow=n, ncol=length(tau), byrow=TRUE)

scmat <- score_eval(y, yhat, sc, cmb=FALSE)
scvec <- apply(scmat, 1, mean)

## 1. Basic mean
sc11 <- score_eval(y, yhat, sc)
sc12 <- mean(scvec)
expect_identical(sc11, sc12)

## 2. Weighted mean
sc21 <- score_eval(y, yhat, sc, w=1:n)
sc22 <- sum(1:n/sum(1:n)*scvec)
expect_identical(sc21, sc22)

## 3. Basic mean, se.
sc31 <- score_eval(y, yhat, sc, se=TRUE)[2]
names(sc31) <- NULL
sc32 <- sd(scvec)/sqrt(n)
expect_identical(sc31, sc32)

## 4. Weighted mean, se.
sc41 <- score_eval(y, yhat, sc, w=1:n, se=TRUE)[2]
names(sc41) <- NULL
sc42 <- sd(scvec)*sqrt(sum((1:n/sum(1:n))^2))
expect_identical(sc41, sc42)

## 5. Matrix remains unchanged if weights are indicated.
sc51 <- score_eval(y, yhat, sc, w=1:n, cmb=FALSE)
expect_identical(sc51, scmat)

## 6. First part of score matrix does not change if more quantiles are added:
sc <- scorer(c(tau, c(0.95, 0.99)), g=g, wtau=wtau)
scmat2 <- score_eval(y, cbind(yhat, yhat+1), sc, w=1:n, cmb=FALSE)
expect_identical(scmat, scmat2[, 1:length(tau)])
