subsequentContrasts <- function(n) {
    matrix(rep(c(1, -1, rep(0, n - 1)), n - 1)[1:(n * (n - 1))],
           nrow = n,
           ncol = n - 1)
}

insertNAs <- function(vec, pos) {
    if (pos > length(vec)) {
        return(c(vec, NA))
    }
    diff <- length(vec) - pos + 2
    return(c(vec[0:(pos - 1)], rep(NA, diff)))
}


updateStatistics <- function(factorMerger, groups, factor) {
    UseMethod("updateStatistics", factorMerger)
}

updateStatistics.allToAllFactorMerger <- function(factorMerger, groups, factor) {
    noGroups <- length(groups)
    if (noGroups > 1) {
        statsTmp <- sapply(groups, function(y) {
            fac <- relevel(factor, ref = y)
            model <- calculateModel(factorMerger, fac)
            getPvals(model)
        })

        if (noGroups == 2) {
            statsTmp <- t(statsTmp)
        }
        stats <- matrix(NA, ncol = noGroups, nrow = noGroups)
        for (i in 1:noGroups) {
            stats[, i] <- insertNAs(statsTmp[, i], i)
        }
        colnames(stats) <- groups; rownames(stats) <- groups

    } else {
        stats <- NULL
    }
    return(stats)
}

startMerging <- function(factorMerger) {
    UseMethod("startMerging", factorMerger)
}

startMerging.subsequentFactorMerger <- function(factorMerger) {
    factorMerger$factor <- setIncreasingOrder(factorMerger$response,
                                              factorMerger$factor)
    contrasts(factorMerger$factor) <-
        subsequentContrasts(length(levels(factorMerger$factor)))
    factorMerger$mergingList[[1]]$factor <- factorMerger$factor
    groups <- levels(factorMerger$mergingList[[1]]$factor)
    factorMerger$mergingList[[1]]$groups <- groups
    model <- calculateModel(factorMerger, factorMerger$factor)
    stats <- c(getPvals(model), NA)
    names(stats) <- groups
    factorMerger$mergingList[[1]]$factorStats <- stats
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = calculateModelStatistic(model),
        pval = NA)
    return(factorMerger)
}

startMerging.allToAllFactorMerger <- function(factorMerger) {
    groups <- factorMerger$mergingList[[1]]$groups
    stats <- updateStatistics(factorMerger, groups, factorMerger$factor)
    factorMerger$mergingList[[1]]$factorStats <- stats
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = calculateModelStatistic(lm(factorMerger$response ~ factorMerger$factor)),
        pval = NA)
    return(factorMerger)
}

canBeMerged <- function(factorMerger) {
    ml <- factorMerger$mergingList
    return(length(ml[[length(ml)]]$groups) > 1)
}

mergePair <- function(factorMerger) {
    UseMethod("mergePair", factorMerger)
}

mergePair.subsequentFactorMerger <- function(factorMerger) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    maxInd <- which.max(fs$factorStats)
    maxStat <- fs$factorStats[maxInd]
    groupA <- names(maxInd)
    groupB <- names(fs$factorStats[maxInd + 1])
    factorMerger$mergingList[[step]]$merged <- c(groupA, groupB)
    groupAB <- paste0(groupA, groupB)
    groups <- fs$groups[-maxInd]
    groups[maxInd] <- groupAB
    factor <- mergeLevels(fs$factor, groupA, groupB, groupAB)

    if (length(unique(factor)) > 1) {
        contrasts(factor) <- subsequentContrasts(length(unique(factor)))
    }

    model <- calculateModel(factorMerger, factor)
    factorStats <- c(getPvals(model), NA)
    names(factorStats) <- groups

    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")
    factorMerger$mergingList[["tmp"]] <- list(groups = groups,
                                              factor = factor,
                                              factorStats = factorStats,
                                              modelStats = NULL,
                                              means = calculateMeans(
                                                  factorMerger$response, factor))

    names(factorMerger$mergingList)[step + 1] <- step + 1
    names(maxStat) <- NULL

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pval = maxStat)
    return(factorMerger)

}

mergePair.allToAllFactorMerger <- function(factorMerger) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    factorStats <- fs$factorStats
    maxInd <- which(factorStats == max(factorStats, na.rm = TRUE), arr.ind = TRUE)[1, ]
    maxStat <- factorStats[maxInd[1], maxInd[2]]
    groups <- fs$groups
    groupA <- groups[maxInd[1]]; groupB <- groups[maxInd[2]]
    groupAB <- paste0(groupA, groupB)
    factorMerger$mergingList[[step]]$merged <- c(groupA, groupB)
    factor <- mergeLevels(fs$factor, groupA, groupB, groupAB)
    groups <- levels(factor)
    noGroups <- length(groups)
    model <- calculateModel(factorMerger, factor)
    stats <- updateStatistics(factorMerger, groups, factor)
    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")
    factorMerger$mergingList[["tmp"]] <- list(groups = groups,
                                              factor = factor,
                                              factorStats = stats,
                                              modelStats = NULL,
                                              means = calculateMeans(
                                                  factorMerger$response, factor))

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pval = maxStat)
    return(factorMerger)
}

mergePair.multiClassFactorMerger <- function(factorMerger) {

}

getTreeWithEdgesLength <- function(factorMerger) {
    nodes <- list()
    initLevels <- levels(factorMerger$factor)
    for (level in initLevels) {
        nodes <- c(nodes, list(node(level)))
    }
    names(nodes) <- initLevels

    ml <- factorMerger$mergingList
    for (i in 1:(length(ml) - 1)) {
        merged <- ml[[i]]$merged
        nodes <- c(nodes, list(node(nodes[[merged[1]]], nodes[[merged[2]]], pval = ml[[i + 1]]$modelStats$pval)))
        names(nodes)[length(nodes)] <- paste0(merged[1], merged[2])
    }

    return(paste0(nodes[[length(nodes)]]$text, ":", pval, ";"))
}
