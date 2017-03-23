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

stats <- function(factorMerger) {
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
    statsDf <- do.call(rbind, statsList)
    statsDf <- subset(statsDf, !duplicated(level))

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
mergingHistory <- function(object, showStats = FALSE) {
    UseMethod("mergingHistory", object)
}

#' @export
#' @importFrom dplyr rename
mergingHistory.factorMerger <- function(factorMerger, showStats = FALSE) {
    mergingList <- sapply(factorMerger$mergingList,
                        function(x) { x$merged })
    mergingDf <- do.call(rbind, mergingList) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        rename(groupA = V1, groupB = V2)

    rownames(mergingDf) <- NULL
    if (showStats) {
        st <- round(stats(factorMerger), 4)
        mergingDf <- mergingDf[complete.cases(mergingDf), ]
        mergingDf <- rbind(c("", ""), mergingDf)
        mergingDf <- data.frame(mergingDf, st)
    }
    return(mergingDf)
}

call <- function(factorMerger) {
    return(
        paste0("Family: ", gsub('([[:upper:]])', ' \\1',
                                class(factorMerger)[length(class(factorMerger))]), ".")
        )
}

#' Factor Merger - ...
#'
#' @export
#'
#' @importFrom knitr kable
#'
print.factorMerger <- function(factorMerger) {
   df <- mergingHistory(factorMerger, TRUE)
   colnames(df)[1:2] <- c("groupA", "groupB")
   cat(call(factorMerger))
   cat("\nFactor levels were recoded as below:")
   cat(paste(c("", "", kable(factorMerger$map, output = FALSE)), collapse = "\n"))
   cat("\n\nFactor levels were merged in the following order:")
   cat(paste(c("", "", kable(df, output = FALSE)), collapse = "\n"))
   invisible(NULL)
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
    fmList <- startMerging(fm, subsequent)
    fm <- fmList$factorMerger
    while (canBeMerged(fm)) {

        fmList <- mergePair(fm, subsequent, fmList$factor, fmList$model)
        fm <- fmList$factorMerger
    }
    return(fmList$factorMerger)
}

