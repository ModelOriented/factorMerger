#' Factor Merger - ...
#'
#' @importFrom data.table data.table
#'
#' @export
#'
#'
#'

factorMerger <- function(response, factor, gaussian = TRUE,
                         subsequent = TRUE) {
    factor <- as.factor(factor)
    # TODO: Make it insensitive for changes in input types
    fm <- list(
        response = response,
        factor = factor,
        mergingList = list(`1` = list(groups = levels(factor),
                                factor = factor,
                                factorStats = list(),
                                modelStats = list()))
    )

    class(fm) <- "factorMerger"

    if (gaussian) { # Here t-test or Hotelling test is used
        class(fm) <- append(class(fm), "gaussianFactorMerger")
    } else { # TODO: Here we'll use Kruskal-Wallis test (equality of c.d.fs)
        class(fm) <- append(class(fm), "non-gaussianFactorMerger")
    }

    if (!is.null(dim(response))) { # MANOVA
        class(fm) <- append(class(fm), "multivariateFactorMerger")
        subsequent <- FALSE
    }

    if (subsequent) {
        class(fm) <- append(class(fm), "subsequentFactorMerger")
    } else {
        class(fm) <- append(class(fm), "allToAllFactorMerger")
    }

    if (!is.null(dim(class))) { # TODO: ...
        class(fm) <- append(class(fm), "multiClassFactorMerger")
    }

    return(fm)
}

print.factorMerger <- function(factorMerger) {
    lapply(factorMerger$mergingList, function(x) { print(x$groups) })
}

stats <- function(object) {
    UseMethod("stats", object)
}

stats.factorMerger <- function(factorMerger) {
    sapply(factorMerger$mergingList, function(x) { x$modelStats })
}

# ---

filterGroups <- function(response, factor, groupA, groupB) {
    xA <- response[factor == groupA]
    xB <- response[factor == groupB]
    return(list(xA, xB))
}

# ---

calculatePairStatistic <- function(factorMerger, factor,
                                   groupA, groupB) {
    UseMethod("calculatePairStatistic", factorMerger)
}

calculatePairStatistic.gaussianFactorMerger <- function(factorMerger, factor,
                                                        groupA, groupB) {
    groups <- filterGroups(factorMerger$response, factor, groupA, groupB)
    if (groupA == groupB) {
        return(-1)
    }
    return(t.test(groups[[1]], groups[[2]])$p.value)
}

# ---

calculateModelStatistic <- function(factorMerger) {
    UseMethod("calculateModelStatistic", factorMerger)
}

calculateModelStatistic.gaussianFactorMerger <- function(factorMerger) {
    y <- factorMerger$response
    x <- factorMerger$mergingList[[length(factorMerger$mergingList)]]$factor
    if (length(levels(x)) > 1) {
        return(logLik(lm(y ~ x))[1])
    } else {
        return(logLik(lm(y ~ 1))[1])
    }
}

# ---

startMerging <- function(factorMerger) {
    UseMethod("startMerging", factorMerger)
}

startMerging.subsequentFactorMerger <- function(factorMerger) {
    factorMerger$factor <- setIncreasingOrder(factorMerger$response,
                                              factorMerger$factor)
    factorMerger$mergingList[[1]]$factor <- factorMerger$factor
    groups <- levels(factorMerger$mergingList[[1]]$factor)
    factorMerger$mergingList[[1]]$groups <- groups
    noGroups <- length(groups)
    stats <- rep(NA, noGroups)
    names(stats) <- groups

    for (k in 1:max(1, (noGroups - 1))) {
        stats[k] <- calculatePairStatistic(factorMerger,
                                           factorMerger$factor,
                                           groups[k], groups[k + 1])
    }

    factorMerger$mergingList[[1]]$factorStats <- stats
    factorMerger$mergingList[[1]]$modelStats <- calculateModelStatistic(factorMerger)
    return(factorMerger)
}

startMerging.allToAllFactorMerger <- function(factorMerger) {
    groups <- factorMerger$mergingList[[1]]$groups
    noGroups <- length(groups)
    stats <- sapply(groups, function(y) {
        sapply(groups, function(x) {
            calculatePairStatistic(factorMerger, factorMerger$factor, x, y)
        })
    })
    colnames(stats) <- groups; rownames(stats) <- groups
    factorMerger$mergingList[[1]]$factorStats <- stats
    factorMerger$mergingList[[1]]$modelStats <- calculateModelStatistic(factorMerger)
    return(factorMerger)
}

# ---

canBeMerged <- function(factorMerger) {
    return(length(factorMerger$mergingList[[length(factorMerger$mergingList)]]$groups) > 1)
}

# ---

mergePair <- function(factorMerger) {
    UseMethod("mergePair", factorMerger)
}

mergePair.subsequentFactorMerger <- function(factorMerger) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    maxStat <- which.max(fs$factorStats)
    groupA <- names(maxStat)
    groupB <- names(fs$factorStats[maxStat + 1])
    groupAB <- paste0(groupA, groupB)
    groups <- fs$groups[-maxStat]
    groups[maxStat] <- groupAB
    factor <- mergeFactor(fs$factor, groupA, groupB, groups)
    factorStats <- fs$factorStats[-maxStat]
    names(factorStats) <- groups

    if (maxStat > 1) {
        factorStats[maxStat - 1] <- calculatePairStatistic(factorMerger, factor,
                                                           groups[maxStat],
                                                           groups[maxStat - 1])
    }

    if (maxStat < length(groups)) {
        factorStats[maxStat] <- calculatePairStatistic(factorMerger, factor,
                                                           groups[maxStat],
                                                           groups[maxStat + 1])
    }

    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")
    factorMerger$mergingList[["tmp"]] <- list(groups = groups,
                                              factor = factor,
                                              factorStats = factorStats,
                                              modelStats = NULL)

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        calculateModelStatistic(factorMerger)
    return(factorMerger)

}

mergePair.allToAllFactorMerger <- function(factorMerger) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    factorStats <- fs$factorStats
    maxInd <- which(factorStats == max(factorStats), arr.ind = TRUE)[1, ]
    maxStat <- factorStats[maxInd[1], maxInd[2]] # to można wrzucić do modelStats
    groups <- fs$groups
    groupA <- groups[maxInd[2]]; groupB <- groups[maxInd[1]]
    groups <- fs$groups[-maxInd]
    groups <- c(groups, groupAB <- paste0(groupA, groupB))
    factor <- mergeFactor(fs$factor, groupA, groupB, groups)
    colnames(factorStats)[maxInd[1]] <- groupAB; rownames(factorStats)[maxInd[1]] <- groupAB
    colnames(factorStats)[maxInd[2]] <- groupAB; rownames(factorStats)[maxInd[2]] <- groupAB
    for (i in 1:length(groups)) {
        factorStats[maxInd[1], i] <- calculatePairStatistic(factorMerger, factor,
                                                   groupAB, colnames(factorStats)[i])
        factorStats[i, maxInd[1]] <- calculatePairStatistic(factorMerger, factor,
                                                      groupAB, rownames(factorStats)[i])
    }

    factorStats <- factorStats[-maxInd[2], -maxInd[2]]
    factor <- factor(factor, levels = colnames(factorStats))
    groups <- colnames(factorStats)
    print(paste(step, groups))

    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")
    factorMerger$mergingList[["tmp"]] <- list(groups = groups,
                                              factor = factor,
                                              factorStats = factorStats,
                                              modelStats = NULL)

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        calculateModelStatistic(factorMerger)
    return(factorMerger)
    # zmienić "for" na apply
}

mergePair.multiClassFactorMerger <- function(factorMerger) {

}

# ---

mergeFactors <- function(response, factor, gaussian = TRUE, subsequent = TRUE) {
    fm <- factorMerger(response, factor, gaussian, subsequent)
    fm <- startMerging(fm)
    while (canBeMerged(fm)) {
        fm <- mergePair(fm)
    }
    return(fm)
}
