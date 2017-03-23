#' Factor Merger - ...
#'
#' @importFrom data.table data.table
#'
#' @export
#'
merger <- function(response, factor, family = "gaussian",
                         subsequent = FALSE) {

    stopifnot(NROW(response) == NROW(factor))

    factor <- as.factor(factor)
    map <- data.frame(`recoded` = paste0("(", abbreviate(levels(factor)), ")"),
                      `original` = levels(factor))
    rownames(map) <- NULL
    factor <- factor(factor, labels = map$recoded)

    fm <- list(
        response = response,
        factor = factor,
        map = map,
        mergingList = list(`1` = list(groups = levels(factor),
                                factor = factor,
                                factorStats = list(),
                                modelStats = list(),
                                groupStats = NA,
                                merged = NA))
    )

    class(fm) <- "factorMerger"

    switch(family,
           "gaussian" = {
               class(fm) <- append(class(fm), "gaussianFactorMerger")
           },

           "survival" = {
               stopifnot(!sum(response[, 1] < 0))
               stopifnot(length(unique(response[, 2])) == 2)
               class(fm) <- append(class(fm), "survivalFactorMerger")
           },

           "binomial" = {
               stopifnot(!sum(!response %in% c(0, 1)))
               class(fm) <- append(class(fm), "binomialFactorMerger")
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

#' Show levels statistic - ...
#'
#' @export
#'
groupsStats <- function(object) {
    UseMethod("groupsStats", object)
}

#' ---
#' @export
groupsStats.factorMerger <- function(factorMerger) {
    statsList <- lapply(factorMerger$mergingList,
                        function(x) { as.data.frame(x$groupStats) })
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
#' @importFrom dplyr rename
mergingHistory.factorMerger <- function(factorMerger) {
    statsList <- sapply(factorMerger$mergingList,
                        function(x) { x$merged })
    do.call(rbind, statsList) %>% as.data.frame(stringsAsFactors = FALSE) %>%
        rename(groupA = V1, groupB = V2)
}

#' Factor Merger - ...
#'
#' @export
#'
#' @importFrom knitr kable
#'
print.factorMerger <- function(factorMerger) {
   stats <- round(stats(factorMerger), 4)
   mergList <- mergingHistory(factorMerger)
   mergList <- mergList[complete.cases(mergList),]
   mergList <- rbind(c("", ""), mergList)
   rownames(mergList) <- NULL
   df <- data.frame(mergList, stats)
   colnames(df)[1:2] <- c("groupA", "groupB")
   cat("Factor levels were recoded as below:")
   cat(paste(c("", "", kable(factorMerger$map, output = FALSE)), collapse = "\n"))
   cat("\n\nFactor levels were merged in the following order:")
   cat(paste(c("", "", kable(df, output = FALSE)), collapse = "\n"))
   invisible(NULL)
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

