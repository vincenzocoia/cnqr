#' #' Fit Bivariate Copula in PCBN
#' #'
#' #' Appends a fitted copula in the pairing order in a PCBN, via CNQR.
#' #'
#' #' @param v Vector of uniform-score observed outcomes.
#' #' @param ucond matrix of uniform conditional predictors.
#' #' @param prev_cops Names of previous copulas fitted in the pairing order,
#' #' matching the order of the columns in \code{ucond}.
#' #' Should contain the "flip qualifier" if it has one.
#' #' @param g Vectorized function to use for transformation of variables
#' #' before estimation. Fed into \code{\link{scoreq}}.
#' #' @param wtau Vectorized function indicating the weight over quantile indices,
#' #' used for estimation. Fed into \code{\link{scoreq}}. \code{NULL} for equal
#' #' weights (of 1).
#' #' @param prev_cpars List of previously fitted copula parameters.
#' #' @param cop Name of a copula family, without flip qualifier.
#' #' @param param0 Starting value for the copula parameter. If blank, I'll try
#' #' using VineCopula's
#' #' @param refit_all Logical; should all copula parameters in the chain
#' #' (i.e. in the pairing order) be re-estimated with the new addition?
#' #' @note \code{ucond} should have \code{length(prev_cops) + 1} columns at least.
#' #' Other columns besides those are ignored.
#' #'
#' #' The marginal distribution of the response, \code{QY} is assumed to be the
#' #' same for each observation.
#' #' @return A list containing
#' #'
#' #' \enumerate{
#' #'      \item \code{$cops}, a vector of copula models (with the most recently
#' #'      fitted shown last),
#' #'      \item \code{$cpars}, a list of copula parameters corresponding to
#' #'      the copulas.
#' #'      \item \code{$score}, the score of the resulting forecaster on the
#' #'      fitting set (i.e. the minimum CNQR objective function)
#' #' }
#' #' @export
#' append_copula <- function(v, ucond, tau, cop, QY = identity,
#'                           g = identity, wtau = NULL,
#'                           prev_cops = NULL, prev_cpars = NULL,
#'                           refit_all = FALSE) {
#'     if (is.null(wtau)) wtau <- function(tau) 1
#'     if (is.null(prev_cops)) prev_cops <- character(0)
#'     if (is.null(prev_cpars)) prev_cpars <- list()
#'     chainlen <- length(prev_cops) + 1
#'     ## Fit the new copula. Differentiate with the independence copula case,
#'     ##  where no fitting is needed.
#'     # if (cop == "indepcop" || cop == "indep") {
#'     #     ## marginal score
#'     #     yhat <- matrix(QY(tau), nrow=length(v), ncol=length(tau), byrow=TRUE)
#'     #     sc <- scoreq(QY(v), yhat, tau)
#'     #     ## Output
#'     #     return(list(cops = c(prev_cops, "indepcop"),
#'     #                 cpars = c(prev_cpars, list(numeric(0))),
#'     #                 score = sc))
#'     # }
#'     ## Get starting values and copula rotation
#'     if (length(prev_cops) == 0) {
#'         vcond <- v
#'     } else {
#'         vcond <- FYgX(v, ucond[, 1:length(prev_cops), drop=FALSE],
#'                       cops=prev_cops, cpars=prev_cpars, FY=identity)
#'     }
#'     initfit <- copsupp::fitbicop_lh(ucond[, chainlen], vcond, families=cop)
#'     ## Let's just use the copula family that comes out of fitbicop_lh(),
#'     ##  which may be different than the requested family (due to a
#'     ##  restriction of VineCopula's BiCopSelect()).
#'     cop <- initfit$cop
#'     copchain <- c(prev_cops, cop)
#'     ## Parameter starting value and parameter space:
#'     if (refit_all) {
#'         cparstart <- c(prev_cpars, initfit$cpar, recursive = TRUE)
#'         ## Parameter space
#'         paramspace <- copsupp::cparspace(c(prev_cops, cop))
#'         bnds <- copsupp::cparspace(c(prev_cops, cop), fn = FALSE)
#'     } else {
#'         cparstart <- initfit$cpar
#'         paramspace <- copsupp::cparspace(cop)
#'         bnds <- copsupp::cparspace(cop, fn = FALSE)
#'     }
#'     ncpar <- length(initfit$cpar)
#'     cparlen <- c(sapply(prev_cpars, length), ncpar, recursive=TRUE)
#'     ## Don't want starting value to be less than 0.01 units away from a
#'     ##  boundary. For example, with 'joe' copula, nlm() has problems
#'     ##  with cparstart=1.000001.
#'     cparstart <- pmax(bnds$lower + 0.01, cparstart)
#'     cparstart <- pmin(bnds$upper - 0.01, cparstart)
#'     ## Fit with CNQR:
#'     #### Get predictions as a function of parameter vector:
#'     if (refit_all) {
#'         yhat <- function(theta) {
#'             QYgX(tau, ucond[, 1:chainlen, drop=FALSE],
#'                  cops = c(prev_cops, cop),
#'                  cpars = cparvec2cpar(theta, cparlen),
#'                  QY = QY)
#'         }
#'     } else {
#'         yhat <- function(theta) {
#'             QYgX(tau, ucond[, 1:chainlen, drop=FALSE],
#'                  cops = c(prev_cops, cop),
#'                  cpars = c(prev_cpars, list(theta)),
#'                  QY = QY)
#'         }
#'     }
#'     #### Objective function:
#'     obj <- function(theta) scoreq(QY(v), yhat(theta), tau,
#'                                   g=g, wtau=wtau)
#'     #### Minimize objective (as long as there are parameters to estimate!)
#'     if (length(cparstart) != 0) {
#'         res <- rnlm(obj, cparstart, paramspace)
#'     } else {
#'         res <- list(estimate=numeric(0),
#'                     minimum=obj(numeric(0)))
#'     }
#'     ## Output.
#'     if (refit_all) {
#'         thetahat <- res$estimate
#'     } else {
#'         thetahat <- c(prev_cpars, res$estimate, recursive=TRUE)
#'     }
#'     return(list(cops = copchain,
#'                 cpars = copsupp::cparvec2cpar(thetahat, cparlen),
#'                 score = res$minimum))
#' }
#'
