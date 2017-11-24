## OBJECTIVE: Make a data matrix for use in examples in the package.
library(CopulaModel)
library(copsupp)
set.seed(555)

## 1. Specify a vine for the data.
G <- AtoG(Dvinearray(5))
copmat <- makevinemat("gum",
                      c("bvtcop", "frk"),
                      c("mtcj", "frk", "indepcop"),
                      c("bvncop", "joe", "mtcj", "frk"),
                      zerocol = TRUE)
cparmat <- makevinemat(3.1,
                       list(c(0.5, 4), 2.3),
                       list(4.2, 3.5, numeric(0)),
                       c(0.5, 2.2, 2.5, 1.6), zerocol = TRUE)
egvine <- rvine(G, copmat, cparmat)

## 2. Generate data
n <- 1000
egdat <- rrvine(n, egvine)

## Save
save(egdat, egvine, file="data/egdat.RData")
