
#' Fit a full vine model space using CNQR
#'
#' Uses CNQR to estimate/select a vine from a general model space. That is,
#' selects and estimates a copula family for each link between parameter
#' and response.
#'
#' @param edges Integer vector starting with the column number of the response
#' variable, followed by the column numbers of the predictors to use,
#' in the order that they are to be linked with the response.
#' @param dat Data frame or matrix containing the data (columns=variables), or
#' a list of such data frames. In the latter case, the first data frame is
#' used for parameter estimation (training data), the second data
#' frame is used for copula selection (validation data), and all other
#' data are not used in the fitting procedure, but are included in the
#' output. \code{dat} may
#' also be a single vector of response data, in the case you
#' have no predictors.
#' @param sc Scoring rule to use for the regression, as in the output
#' of \code{\link{scorer}}.
#' @param basevine Object of type \code{'rvine'} containing
#' the predictors, and not the response. If left blank, the predictors are
#' assumed to be independent.
#' @param cdf List of vectorized distribution functions of the data, where
#' the entries correspond respectively to the columns in the data frames in
#' \code{dat}. Or, if the distribution functions are all the same, \code{cdf}
#' can be that single function. You can ignore this argument if you don't
#' include any predictors in \code{edges}.
#' @param copspace List with vector entries of the copula families to try fitting
#' for each edge. \code{NULL} entries will be replaced with all copula families
#' in \code{families}. Or, leave the argument as \code{NULL} to fit all families in
#' \code{families} for all edges.
#' @param QY Quantile function of the response, which accepts a
#' vector of values (quantile levels) in (0,1). It should return
#' quantiles, either in the form of
#' a vector corresponding to the input,
#' or in the form of a matrix with columns corresponding to the inputted
#' quantile levels and rows corresponding to the observations
#' (thus allowing for each observation of the response to come from different
#' distributions).
#' @param refit After all the copula models have been selected and fit
#' (this is done sequentially), should the parameters of those copula
#' families be re-estimated, but this time altogether, using all the data?
#' \code{TRUE} if so.
#' @param verbose Logical; should messages be output to indicate what
#' \code{cnqr} is doing?
#' @param families Vector of copula families.
#' For those edges in \code{cop} that don't have copula families
#' specially selected (i.e. have \code{NULL} entries), these families are used.
#' @return
#' An object of class \code{'cnqr'}, which is also an object of type
#' \code{'rvine'} containing the final fitted vine.
#' In addition to the list entries present in \code{'rvine'}
#' objects, a \code{'cnqr'} object has the following named entries:
#'
#' \describe{
#'      \item{y}{List of response vectors.}
#'      \item{uind}{List of matrices of independent predictors.}
#'      \item{QY}{Either a list of quantile functions, or a vectorized
#'      quantile function. In the latter case, the marginal distribution
#'      is stationary -- that is, does not differ depending on the observation.}
#'      \item{cdf}{List of vectorized distribution functions, corresponding
#'      to the columns of the data frames/matrices inputted.}
#'      \item{yhat}{List of matrices of quantile predictions. Matrix columns
#'      correspond to the quantile levels, and rows correspond to observations.}
#'      \item{scorer}{The scoring rule used for the regression, as in the output
#'      of \code{\link{scorer}}.}
#'      \item{score}{A named vector of scores, one entry per data type.}
#' }
#'
#' Except for \code{$scorer} and possibly \code{$QY} (if response data
#' are stationary), each entry is a list/vector
#' corresponding to a different type of data set. The first one is always the
#' training data, named \code{"tr"}. If a second
#' data matrix is supplied in \code{dat}, the
#' second entry will correspond to this (the "validation
#' data"), and is named \code{"val"}.
#' Additional data supplied through \code{dat} will also appear in these lists,
#' and either be named according to their names in \code{dat}, or if
#' unnamed, will have the names \code{'dat3'}, \code{'dat4'}, etc.
#' The vector in \code{$score} also corresponds to these data sets.
#' @examples
#' data(egdat)
#' dat <- list(egdat[1:750, ], egdat[751:1000, ])
#' basevine <- subset(egvine, 1:4)
#' sc <- scorer(space_taus(10))
#'
#' summary(cnqr(5:1, dat, sc, basevine=basevine))
#' summary(cnqr(5, dat, sc))
#'
#' set.seed(123)
#' newdat <- rrvine(100, egvine)
#' dat <- c(dat, list(test=newdat))
#' fit <- cnqr(5:3, dat, sc, basevine=basevine,
#'             families=c("bvncop", "joe"))
#' summary(fit)
#' @import copsupp
#' @export
cnqr <- function(edges, dat, sc, basevine, cdf=identity,
                 QY=identity, copspace=NULL, refit=FALSE, verbose=FALSE,
                 families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                              "frk","joe","bb1", "bskewncop", "bskewncopp")) {
    ## --- Get appropriate basevine ---
    if (missing(basevine)) {
        ## basevine was not specified, which means that the predictors (if present)
        ##  are assumed to be independent.
        basevine <- rvine(matrix(edges[-1], nrow=1))
    }
    if (is.vector(dat) & !is.list(dat)) dat <- matrix(dat)
    if (is.data.frame(dat) | is.matrix(dat)) {
        # Only one data set has been input.
        novaldat <- TRUE  # There's no validation data.
        dat <- list(as.matrix(dat))
    } else {
        novaldat <- FALSE
        dat <- lapply(dat, as.matrix)
    }
    ## --- Manipulate data ---
    y <- lapply(dat, function(dat) dat[, edges[1]])
    ytr <- y[[1]]
    if (novaldat) yval <- ytr else yval <- y[[2]]
    if (length(cdf) == 1) cdf <- rep(list(cdf), ncol(dat[[1]]))
    dat <- lapply(dat, dat2udat, cdf=cdf)
    ## --- Learn about the Quantile Function ---
    ## Is the quantile function the same for each observation? Store the answer
    ##  in `stationary`. Either way, store the training and validation qf's
    ##  separately.
    if (is.list(QY)) {  # If QY is a list, then no.
        stationary <- FALSE
        QYtr <- QY[[1]]
        QYval <- QY[[2]]
    } else {
        QYtr <- QY
        QYval <- QY
        if (is.matrix(QY(0.5))) { # If QY outputs a matrix, then no.
            stationary <- FALSE
        } else {  # Otherwise, quantile function is the same.
            stationary <- TRUE
        }
    }
    ## --- Extract useful data quantities ---
    ## (a) Which are predictors (and what order)? What's the response?
    xlab <- edges[-1]
    ylab <- edges[1]
    ## (b) Extract the response; put marginals back in.
    ##  NOTE: Uniform responses are needed for getting starting values with MLE.
    if (novaldat) {  # Only one data set has been input.
        ## (c) Define the two training and validation matrices.
        dattr <- dat[[1]]
        datval <- dattr
        ## (d) Map predictors to independent uniform set.
        if (verbose) cat("Computing independent predictors for training data.\n")
        uindtr <- pcondseq(dattr, ord=xlab, rv=basevine)
        uindval <- uindtr
    } else {  # Both a training and validation set were input.
        ## (c) Define the two training and validation matrices.
        dattr <- dat[[1]]
        datval <- dat[[2]]
        ## (d) Map predictors to independent uniform set.
        if (verbose) cat("Computing independent predictors for training data.\n")
        uindtr <- pcondseq(dattr, ord=xlab, rv=basevine)
        if (verbose) cat("Computing independent predictors for validation data.\n")
        uindval <- pcondseq(datval, ord=xlab, rv=basevine)
    }
    vtr <- dattr[, ylab]
    vval <- datval[, ylab]
    ## --- Extract full model space ---
    res <- augment(basevine, a=ylab)
    d <- ncol(res$G)
    p <- length(xlab) # This is not necessarily the number of predictors in the vine (=d-1).
    ## Fill-in copspace if not done already.
    nfam <- length(families)
    if (is.null(copspace)) copspace <- rep(list(NULL), p)
    copspace <- as.list(lapply(copspace, function(fams){
        if (is.null(fams)) families else fams
    }))
    ## --- Fitting and selection procedure ---
    for (i in seq_len(p)) {
        if (verbose) cat(paste("--- Fitting edge", i, "of", p, "---\n"))
        ## --- Fit a copula to edge i ---
        ## Extract the required independent predictors
        this_uindtr <- uindtr[, 1:i, drop=FALSE]
        ## Get conditional PIT score of response and most recent predictor,
        ##  given the predictors that have already been fit.
        ucondtr <- uindtr[, i]
        fittedcops <- res$copmat[seq_len(i-1), d]
        if (length(fittedcops) == 0) {
            ## There are no copula families fit yet.
            vcondtr <- vtr
        } else {
            vcondtr <- FYgX(y=vtr,
                            ucond=this_uindtr[, -i, drop=F],  # rmv last ('active') col.
                            cops=fittedcops,
                            cpars=res$cparmat[seq_len(i-1), d],
                            FY=identity)
        }
        ## Loop through candidate copula families, fitting each one to the running vine.
        res_cand <- lapply(copspace[[i]], function(cop){
            if (verbose) cat(paste0("Fitting copula '", cop, "'. "))
            ## Get initial parameter estimates, and select copula reflection/permutation.
            init <- cpar_init(ucondtr, vcondtr, cop)
            ## Let's just use the copula family that comes out of cpar_init(),
            ##  which may be different than the requested family (due to a
            ##  restriction of VineCopula's BiCopSelect()).
            cop <- init$cop
            cpar <- list(init$cpar)
            ## Get parameter estimates using CNQR on this copula, with the
            ##  training data.
            cparhat <- cnqr_est(res, a=xlab[i], cop=cop, cpar_init=cpar, sc=sc,
                                y=ytr, uind=this_uindtr, QY=QYtr)
            if (verbose) cat(paste0("Parameter: (", paste(cparhat[[1]], collapse=", "), ")\n"))
            ## Augment running vine with this fit. The result is a candidate model.
            augment(res, a=xlab[i], cop=cop, cpar=cparhat, col=d)
        })
        ## Select the best candidate model on the validation set.
        if (verbose) cat("Selecting best copula family for this edge.\n")
        this_uindval <- uindval[, 1:i, drop=FALSE]
        res <- cnqr_sel(res_cand, sc=sc, y=yval, uind=this_uindval, QY=QYval)
        chosen_cop <- tail(xylink(res)$cops, 1)
        if (verbose) cat(paste0("Selected '", chosen_cop, "' copula.\n"))
    }
    ## --- Refit entire column --- (if asked)
    if (refit & p>0) {
        if (verbose) print("Refitting the selected copula families altogether.\n")
        ## Get copula families for each edge, and use their parameters as
        ##  starting values:
        cops <- res$copmat[seq_len(p), d]
        cpars <- res$cparmat[seq_len(p), d]
        ## Re-start `res` as having no links with the predictors.
        res <- augment(basevine, a=ylab)
        ## Get parameter estimates using *all* the data:
        if (novaldat) {
            ## There's no validation data. Just use training data.
            cparhat <- cnqr_est(res, a=xlab, cop=cops, cpar_init=cpars,
                                sc=sc, y=ytr, uind=uindtr, QY=QYtr)
        } else {
            ## There's separate validation data. Combine training and validation
            ##  data to use in the estimation.
            y <- c(ytr, yval)
            uind <- rbind(uindtr, uindval)
            if (stationary) {
                QYall <- QY
            } else {
                QYall <- function(tau) rbind(QYtr(tau), QYval(tau))
            }
            cparhat <- cnqr_est(res, a=xlab, cop=cops, cpar_init=cpars,
                                sc=sc, y=y, uind=uind, QY=QYall)
        }
        ## Bind parameter estimates to the vine.
        res <- augment(res, a=xlab, cop=cops, cpar=cparhat, col=d)
    }
    ## --- Convert vine to `cnqr` object ---
    if (verbose) cat("Computing predictions and scores.\n")
    ## (a) Get prediction matrices.
    yhattr <- QYgX(sc$tau, ucond=uindtr,
                   cops=res$copmat[seq_len(p), d],
                   cpars=res$cparmat[seq_len(p), d],
                   QY=QYtr)
    ## (b) Get scores
    scorestr <- score_eval(ytr, yhattr, sc)
    ## Repeat (a), (b) if there's validation data
    if (!novaldat) {
        yhatval <- QYgX(sc$tau, ucond=uindval,
                       cops=res$copmat[seq_len(p), d],
                       cpars=res$cparmat[seq_len(p), d],
                       QY=QYval)
        scoresval <- score_eval(yval, yhatval, sc)
        yhat <- list(tr=yhattr, val=yhatval)
        scores <- c(tr=scorestr, val=scoresval)
    } else {
        yhat <- list(tr=yhattr)
        scores <- c(tr=scorestr)
    }
    ## (c) Append items to the vine.
    if (novaldat) {
        res$y <- list(tr=ytr)
        res$uind <- list(tr=uindtr)
    } else {
        res$y <- list(tr=ytr, val=yval)
        res$uind <- list(tr=uindtr, val=uindval)
    }
    if (!stationary) {
        if (novaldat) {
            QY <- list(tr=QYtr)
        } else {
            QY <- list(tr=QYtr, val=QYval)
        }
    }
    res$QY <- QY
    res$cdf <- cdf
    res$yhat <- yhat
    res$scorer <- sc
    res$score <- scores
    class(res) <- c("cnqr", "rvine")
    if (length(dat) > 2) {
        if (stationary) {
            res <- adddat(res, dat[-(1:2)])
        } else {
            res <- adddat(res, dat[-(1:2)], QY=QY[-(1:2)])
        }
    }
    return(res)
}



#' @export
print.cnqr <- function(object) {
    d <- ncol(object$G)
    r <- nrow(object$G)
    yvar <- object$G[1, d]
    xord <- object$G[1+seq_len(r-1), d]
    xord <- xord[xord != 0]
    p <- length(xord)
    ## Output
    cat(paste0("CNQR fit for variable ", yvar, "|{",
               paste(xord, collapse=","), "}.\n\n"))
    cat("Quantile levels: ")
    cat(paste(signif(object$scorer$tau, 3), collapse=", "))
    invisible()
}

#' @export
summary.cnqr <- function(object) {
    d <- ncol(object$G)
    r <- nrow(object$G)
    yvar <- object$G[1, d]
    xord <- object$G[1+seq_len(r-1), d]
    xord <- xord[xord != 0]
    p <- length(xord)
    ## Output
    print(object)
    if (p > 0) {
        cat("\n\nCopulas:\n")
        cops <- object$copmat[seq_len(p), d]
        cpars <- object$cparmat[seq_len(p), d]
        cpars <- sapply(cpars, function(v) paste(signif(v, 3), collapse=", "))
        cops <- paste0(cops, "(", cpars, ")")
        namescops <- paste0("{", xord[1], ",", yvar, "}")
        for (i in 1+seq_len(p-1)) {
            namescops[i] <- paste0("{", xord[i], ",", yvar, "}|{",
                                   paste(xord[1:(i-1)], collapse = ","), "}")
        }
        names(cops) <- namescops
        print(cops)
    } else {
        cat("\n")
    }
    ## Score:
    cat("\nScores:\n")
    sc <- sapply(object$score, signif, digits=3)
    print(sc)
    invisible()
}
