#' #' Fit an Edge in CNQR
#' #'
#' #' Append an edge in the PCBN forecaster that fits a validation set best.
#' #'
#' #' @param vtr,vval Vector of uniform scores of the response. \code{vtr}
#' #' is the training data, and \code{vval} is the validation data.
#' #' @param ucondtr,ucondval Matrix of conditional predictor cdfs,
#' #' obtained by \code{\link{pcondseq}}.
#' #' @param tau Vector of quantile indices to score on.
#' #' @param QY Vectorized function of the marginal quantile function of the
#' #' response (which is assumed to be the same for all observations).
#' #' @param g Vectorized function to use for transformation of variables
#' #' before estimation. Fed into \code{\link{scoreq}}. If you wish to
#' #' differentially transform each observation, you can optionally make this
#' #' a list of two such functions -- the first function will be applied to
#' #' the training set, the second to the validation set.
#' #' @param wtau Vectorized function indicating the weight over quantile indices,
#' #' used for estimation. Fed into \code{\link{scoreq}}. \code{NULL} for equal
#' #' weights (of 1).
#' #' @param prev_cops Vector of copula names that link previous predictors
#' #' in the pairing order with the response. \code{NULL} if you're adding
#' #' the first predictor in a chain.
#' #' @param prev_cpars List of copula parameters corresponding to
#' #' \code{prev_cops}.
#' #' @param refit_all Should previous copula models be re-fit? Either \code{0},
#' #' \code{1}, or \code{2} -- see details.
#' #' @param families Vector of copula family names to try fitting in
#' #' the new edge.
#' #' @note The "base" forecaster for which a new predictor is being added
#' #' is determined by \code{prev_cops} and \code{prev_cpars}. That means the
#' #' \code{ucond}'s can be wider than the number of "base" copulas.
#' #' @details
#' #' Here are the codes for the \code{refit_all} argument.
#' #'
#' #' \itemize{
#' #'      \item \code{0}: Do not modify the previously fitted copulas.
#' #'      \item \code{1}: Once a copula model is chosen for the new edge, then
#' #'      re-estimate the previous copula parameters along with the new one.
#' #'      \item \code{2}: Re-estimate the previous copula parameters along with
#' #'      each candidate copula model for the new edge.
#' #' }
#' #' @return A list with the following entries:
#' #'
#' #' \enumerate{
#' #'      \item \code{$cops}: Vector of copula family names in the pairing
#' #'      order, up until (and including) the one just fitted.
#' #'      \item \code{$cpars}: List of fitted copula parameters corresponding
#' #'      to the copula families in \code{$cops}.
#' #'      \item \code{$scoretr}: The score on the training set.
#' #'      \item \code{$scoreval}: The score on the validation set.
#' #' }
#' #' @export
#' append_edge <- function(vtr, vval, ucondtr, ucondval, tau, QY = identity,
#'                         g = identity, wtau = NULL,
#'                         prev_cops = NULL, prev_cpars = NULL, refit_all = 0,
#'                         families = c("indepcop","bvncop","bvtcop","mtcj",
#'                                      "gum","frk","joe","bb1")) {
#'     if (is.null(wtau)) wtau <- function(tau) 1
#'     if (length(g) == 1){
#'         gtr <- g
#'         gval <- g
#'     } else {
#'         gtr <- g[[1]]
#'         gval <- g[[2]]
#'     }
#'     if (is.null(prev_cops)) prev_cops <- character(0)
#'     if (is.null(prev_cpars)) prev_cpars <- list()
#'     chainlen <- length(prev_cops) + 1
#'     ## Get fits for each copula family
#'     fits <- list()
#'     scores <- numeric(0)
#'     for (i in 1:length(families)) {
#'         ## Fit the copula
#'         fits[[i]] <- append_copula(vtr, ucondtr, tau=tau,
#'                                    cop=families[i],
#'                                    g=gtr, wtau=wtau,
#'                                    prev_cops=prev_cops,
#'                                    prev_cpars=prev_cpars,
#'                                    QY=QY,
#'                                    refit_all=refit_all==2)
#'         ## Score on the validation set
#'         yhat <- QYgX(tau, ucondval[, 1:chainlen, drop=FALSE],
#'                      cops=fits[[i]]$cops, cpars=fits[[i]]$cpars,
#'                      QY = QY)
#'         scores[i] <- scoreq(QY(vval), yhat, tau, g=gval, wtau=wtau)
#'     }
#'     ## Choose the fit with the smallest score. Some of the copula families might
#'     ##  be duplicates (because VineCopula's BiCopSelect() changes the inputted
#'     ##  copula sometimes, which is passed through to append_copula()), so we'll
#'     ##  take *one of the models* that has the minimum score.
#'     bestscore <- min(scores)
#'     bestfit <- fits[[which(scores == bestscore)[1]]]
#'     ## Now that the copula model has been chosen, do we need to re-fit all the
#'     ##  copulas in the chain together?
#'     if (refit_all == 1) {  #Yes -- re-fit all copulas.
#'         ## Get fitted copula
#'         fittedcop <- tail(bestfit$cops, 1)
#'         #### append_copula() doesn't accept flipped copulas, so get the base
#'         ####  copula family if necessary.
#'         trycop <- substr(fittedcop, 1, nchar(fittedcop)-1)
#'         if (exists(paste0("p", trycop))) fittedcop <- trycop
#'         ## Fit
#'         bestfit <- append_copula(vtr, ucondtr, tau=tau,
#'                                  cop=fittedcop,
#'                                  g=gtr, wtau=wtau,
#'                                  prev_cops=prev_cops,
#'                                  prev_cpars=prev_cpars,
#'                                  QY=QY,
#'                                  refit_all=TRUE)
#'         ## Score on the validation set
#'         yhat <- QYgX(tau, ucondval[, 1:chainlen, drop=FALSE],
#'                      cops=bestfit$cops, cpars=bestfit$cpars,
#'                      QY = QY)
#'         bestscore <- scoreq(QY(vval), yhat, tau, g=gval, wtau=wtau)
#'     }
#'     ## Return qualities of the best fit.
#'     return(list(cops = bestfit$cops,
#'                 cpars = bestfit$cpars,
#'                 scoretr = bestfit$score,
#'                 scoreval = bestscore))
#' }
