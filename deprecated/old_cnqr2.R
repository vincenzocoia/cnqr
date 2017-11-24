#' Composite Nonlinear Quantile Regression
#'
#' After fitting a vine on predictors, this function lengthens the
#' vine (or pair-copula Bayesian network
#' in general) by appending a response variable. Copulas that link
#' the response to the predictors are fit using CNQR.
#'
#' @param edges Integer vector; the pairing order of the predictors to the response.
#' Should start with the column number of the response variable, followed by the
#' column numbers of the predictors according to pairing order.
#' @param dattr,datval Data frame containing uniform scores of the data, each
#' variable having its own column. \code{dattr} contains training data to fit the
#' models, and \code{datval} contains validation data for model selection.
#' @param basevine Object of type \code{"rvine"} representing the R-vine fit
#' to the predictors. See \code{rvine()} and \code{fitrvine_basic()} in the
#' \code{copsupp} package for building such objects.
#' @param QY Vectorized function of the fitted quantile function of Y.
#' @param tauset Vector of quantile indices to fit the model to.
#' @param g Vectorized function to use for transformation of variables
#' before estimation. Fed into \code{\link{scoreq}}. If you wish to
#' differentially transform each observation, you can optionally make this
#' a list of two such functions -- the first function will be applied to
#' the training set, the second to the validation set.
#' @param wtau Vectorized function indicating the weight over quantile indices,
#' used for estimation. Fed into \code{\link{scoreq}}. \code{NULL} for equal
#' weights (of 1).
#' @param refit_all Should previous copula models be re-fit? Either \code{0},
#' \code{1}, or \code{2} -- see details.
#' @param cops List with vector entries of the copula families to try fitting
#' for each edge. \code{NULL} entries will be replaced with all copula families
#' in \code{families}. Or, leave the argument as \code{NULL} to fit all families in
#' \code{families} for all edges.
#' @param families Vector of copula families.
#' For those edges in \code{cop} that don't have copula families
#' specially selected (i.e. have \code{NULL} entries), these families are used.
#' @details
#' Here are the codes for the \code{refit_all} argument, regarding fitting
#' a new predictor/edge in the pairing order.
#'
#' \itemize{
#'      \item \code{0}: Do not modify the previously fitted copulas.
#'      \item \code{1}: Once a copula model is chosen for the new edge, then
#'      re-estimate the previous copula parameters along with the new one.
#'      \item \code{2}: Re-estimate the previous copula parameters along with
#'      each candidate copula model for the new edge.
#' }
#'
#' Fitting and scoring is always done on the scale of the response variable
#' (i.e. having quantile function \code{QY}), and not on the uniform scores.
#' @return
#' An object of class \code{"cnqr"}, which is a list of the following named
#' entries:
#'
#' \describe{
#'      \item{xord}{= \code{edges[-1]}, is the pairing order of the predictors.}
#'      \item{yvar}{= \code{edges[1]}, is the column number of the response.}
#'      \item{cops}{Character vector of fitted copula families when connecting
#'      predictors in the order of \code{$xord}.}
#'      \item{cpars}{List of numeric vector entries representing the fitted
#'      copula parameters for the families in \code{$cops}}
#'      \item{tau}{= the \code{tauset} argument, the quantile indices regressed on.}
#'      \item{basevine}{= the \code{basevine} argument, the vine fitted to the
#'      predictors.}
#'      \item{QY}{= the \code{QY} argument, the vectorized quantile function of the
#'      response.}
#'      \item{ucondtr,ucondval}{Matrix of conditional predictors,
#'      each having their own
#'      column, according to the pairing order. \code{$ucondtr} is
#'      calculated from the training set, and \code{$ucondval}
#'      is calculated from the validation set}
#'      \item{ytr,yval}{Vector of response values (training and validation).}
#'      \item{scorestr,scoresval}{Vector of scores of the forecasters
#'      starting from the marginal forecaster (which uses no predictors) to
#'      the full forecaster (which uses all indicated predictors).
#'      \code{scorestr} for the training set, and \code{scoresval} for
#'      the validation set.}
#' }
#' @examples
#' set.seed(123)
#' library(copsupp)
#' library(ggplot2)
#'
#' ## Get data
#' rv <- rvine(AtoG(CopulaModel::Dvinearray(5)), "frk", 3)
#' dattr <- rrvine(200, rv)
#' datval <- rrvine(200, rv)
#'
#' ## Fit CNQR with the order 5:1, by refitting none to refitting all copulas.
#' fit1 <- cnqr(5:1, dattr, datval, rv, QY=identity)
#' summary(fit1)
#' fit2 <- cnqr(5:1, dattr, datval, rv, QY=identity, refit_all=1)
#' summary(fit2)
#' fit3 <- cnqr(5:1, dattr, datval, rv, QY=identity, refit_all=2)
#' summary(fit3)
#'
#' ## Enforce a Frank copula on edge 1
#' fit4 <- cnqr(5:1, dattr, datval, rv, QY=identity,
#'              cops=list("frk", NULL, NULL, NULL))
#' summary(fit4)
#'
#' ## Try a different order
#' fit5 <- cnqr(c(5, 3, 2, 1, 4), dattr, datval, rv, QY=identity)
#' summary(fit5)
#' @import CopulaModel copsupp
#' @export
cnqr <- function(edges, dattr, datval, basevine, QY, tauset = space_taus(10),
                 g = identity, wtau = NULL,
                 refit_all = 0, cops = NULL,
                 families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                              "frk","joe","bb1")) {
    library(ggplot2)
    if (is.null(wtau)) wtau <- function(tau) 1
    if (length(g) == 1) {
        gtr <- g
        gval <- g
    } else {
        gtr <- g[[1]]
        gval <- g[[2]]
    }
    ## 1. Modify and get info from input
    vtr <- dattr[, edges[1]]
    Ytr <- QY(vtr)
    vval <- datval[, edges[1]]
    Yval <- QY(vval)
    ## Fill-in cops if not done already.
    nfam <- length(families)
    nedge <- length(edges) - 1  # Number of edges
    if (is.null(cops)) cops <- rep(list(NULL), nedge)
    cops <- as.list(lapply(cops, function(fams){
        if (is.null(fams)) families else fams
    }))
    ## Make cpars a list if it's not already.
    #     if (!is.list(cpars)) {
    #         cpars <- rep(list(NULL), nedge)
    #     }
    ## 2. Get sequential predictor cdfs
    ucondtr <- pcondseq(dattr, edges[-1], basevine)
    ucondval <- pcondseq(datval, edges[-1], basevine)
    ## 3. Get scores of the marginal forecaster
    yhat <- function(n) QYgX(tauset, matrix(ncol=0, nrow=n), QY=QY)
    scoretr <- scoreq(Ytr, yhat(length(Ytr)), tauset, g=gtr, wtau=wtau)
    scoreval <- scoreq(Yval, yhat(length(Yval)), tauset, g=gval, wtau=wtau)
    ## 4. Fit each edge, one at a time.
    recentfit <- list()
    for (i in 1:nedge) {
        recentfit <- append_edge(vtr, vval, ucondtr, ucondval, tau=tauset, QY=QY,
                                 g=g, wtau=wtau,
                                 prev_cops=recentfit$cops, prev_cpars=recentfit$cpars,
                                 refit_all=refit_all, families=cops[[i]])
        scoretr[i+1] <- recentfit$scoretr
        scoreval[i+1] <- recentfit$scoreval
    }
    ## Output
    res <- list(xord = edges[-1],
                yvar = edges[1],
                cops = recentfit$cops,
                cpars = recentfit$cpars,
                tau = tauset,
                basevine = basevine,
                QY = QY,
                ucondtr = ucondtr,
                ucondval = ucondval,
                ytr = Ytr,
                yval = Yval,
                scorestr = scoretr,
                scoresval = scoreval)
    class(res) <- "cnqr"
    return(res)
}
