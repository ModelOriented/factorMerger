#' Factor Merger - ...
#'
#' @importFrom data.table data.table
#'
#' @export
#'
merger <- function(response, factor, family = "gaussian",
                         subsequent = FALSE) {

    if (NROW(response) != NROW(factor)) {
        stop("Response and factor sizes do not match.")
    }

    factor <- as.factor(factor)
    # TODO: Make it insensitive to input types changes
    fm <- list(
        response = response,
        factor = factor,
        mergingList = list(`1` = list(groups = levels(factor),
                                factor = factor,
                                factorStats = list(),
                                modelStats = list(),
                                means = NA,
                                merged = NA))
    )

    class(fm) <- "factorMerger"

    if (NCOL(response) == 1) {
        fm$mergingList[[1]]$means <- calculateMeans(response, factor)
    }

    switch(family,
           "gaussian" = {
               class(fm) <- append(class(fm), "gaussianFactorMerger")
           },

           "survival" = {
               class(fm) <- append(class(fm), "survivalFactorMerger")
               stop("Survival analysis is not supported yet.")
           },
           "nonparametric" = {
               class(fm) <- append(class(fm), "nonparametricFactorMerger")
               stop("Non-parametric analysis is not supported yet.")
           },
           stop("Unknown family"))

    if (NCOL(factor) > 1) { # TODO: ...
        class(fm) <- append(class(fm), "multiClassFactorMerger")
        stop("Factor merging with multivariate factor is not supported yet.")
    }

    return(fm)
}

#' Show models statistics - ...
#'
#' @export
#'
stats <- function(object) {
    UseMethod("stats", object)
}

#' ---
stats.factorMerger <- function(factorMerger) {
    statsList <- lapply(factorMerger$mergingList, function(x) { x$modelStats })
    do.call(rbind, statsList)
}

#' Show levels means - ...
#'
#' @export
#'
means <- function(object) {
    UseMethod("means", object)
}

#' ---
#' @export
means.factorMerger <- function(factorMerger) {
    statsList <- lapply(factorMerger$mergingList,
                        function(x) { as.data.frame(x$means) })
    statsDf <- do.call(rbind, statsList) %>% unique()
    if (sum(complete.cases(statsDf)) == 0) {
        return(NULL)
    }
    rownames(statsDf) <- statsDf$level
    statsDf <- subset(statsDf, select = -level)
    return(statsDf)
}


#' Show merging history - ...
#'
#' @export
#'
mergingHistory <- function(object) {
    UseMethod("mergingHistory", object)
}

#' @export
mergingHistory.factorMerger <- function(factorMerger) {
    statsList <- sapply(factorMerger$mergingList,
                        function(x) { x$merged })
    do.call(rbind, statsList)
}

#' Factor Merger - ...
#'
#' @export
#'
print.factorMerger <- function(factorMerger) {
   stats <- round(stats(factorMerger), 4)
   mergList <- mergingHistory(factorMerger)
   mergList <- mergList[complete.cases(mergList),]
   mergList <- rbind(c("", ""), mergList)
   rownames(mergList) <- NULL
   df <- data.frame(mergList, stats)
   colnames(df)[1:2] <- c("groupA", "groupB")
   print(df)
}



node <- function(left, right = NULL, stat = NULL) {
    if (is.null(right)) {
        if (is.na(stat)) {
        return(list(stat = 1,
                   text = left))
        }
        else {
            return(list(stat = stat,
                        text = left))
        }
    }
    leftDiff <- left$stat - stat
    rightDiff <- right$stat - stat
    return(list(stat = stat,
               text = paste0("(", left$text, ": ", leftDiff,
                             ", ", right$text, ": ", rightDiff, ")")))
}


#' Merge factors - ...
#'
#' @export
#'
mergeFactors <- function(response, factor, family = "gaussian", subsequent = FALSE) {

    stopifnot(!is.null(response), !is.null(factor))
    if (NCOL(response) > 1 && subsequent) {
        warning("Subsequent merging with multivariate responseis not yet implemented. All-to-all merging run instead.")
        subsequent <- FALSE
    }

    if (is.data.frame(response)) {
        response <- as.matrix(response)
    }

    fm <- merger(response, factor, family)
    fm <- startMerging(fm, subsequent)
    while (canBeMerged(fm)) {
        fm <- mergePair(fm, subsequent)
    }
    return(fm)
}

