subsequentContrasts <- function(n) {
    matrix(rep(c(1, -1, rep(0, n - 1)), n - 1)[1:(n * (n - 1))],
           nrow = n,
           ncol = n - 1)
}
#' ----
startMerging <- function(factorMerger) {
    UseMethod("startMerging", factorMerger)
}

#' ----
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

#' ----
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
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = calculateModelStatistic(factorMerger),
        pval = NA)
    return(factorMerger)
}

#' ----
canBeMerged <- function(factorMerger) {
    ml <- factorMerger$mergingList
    return(length(ml[[length(ml)]]$groups) > 1)
}

#' ----
mergePair <- function(factorMerger) {
    UseMethod("mergePair", factorMerger)
}

#' ----
mergePair.subsequentFactorMerger <- function(factorMerger) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    maxInd <- which.max(fs$factorStats)
    maxStat <- fs$factorStats[maxInd]
    groupA <- names(maxInd)
    groupB <- names(fs$factorStats[maxInd + 1])
    factorMerger$mergingList[[step]]$merged <- c(groupA, groupB)
    # groupAB <- paste0("(", paste(groupA, groupB, sep = ","), ")")
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

#' ----
getNames <- function(object) {
    UseMethod("getNames", object)
}

#' ----
getNames.matrix <- function(vector) {
    return(colnames(vector))
}

#' ----
getNames.numeric <- function(numeric) {
    return(names(numeric))
}

#' ----
mergePair.allToAllFactorMerger <- function(factorMerger) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    factorStats <- fs$factorStats
    maxInd <- which(factorStats == max(factorStats), arr.ind = TRUE)[1, ]
    maxStat <- factorStats[maxInd[1], maxInd[2]] # to można wrzucić do modelStats
    groups <- fs$groups
    groupA <- groups[maxInd[2]]; groupB <- groups[maxInd[1]]
    # groupAB <- paste0("(", paste(groupA, groupB, sep = ","), ")")
    groupAB <- paste0(groupA, groupB)
    factorMerger$mergingList[[step]]$merged <- c(groupA, groupB) # trochę tu jest copy-paste...
    factor <- mergeLevels(fs$factor, groupA, groupB, groupAB)
    groups <- levels(factor)
    colnames(factorStats)[maxInd[1]] <- groupAB; rownames(factorStats)[maxInd[1]] <- groupAB
    colnames(factorStats)[maxInd[2]] <- groupAB; rownames(factorStats)[maxInd[2]] <- groupAB
    for (i in 1:length(groups)) {
        factorStats[maxInd[1], i] <- calculatePairStatistic(factorMerger, factor,
                                                            groupAB, colnames(factorStats)[i])
        factorStats[i, maxInd[1]] <- calculatePairStatistic(factorMerger, factor,
                                                            groupAB, rownames(factorStats)[i])
    }

    if (ncol(factorStats) > 2) {
        factorStats <- factorStats[-maxInd[2], -maxInd[2]]
    } else {
        factorStats <- factorStats[-maxInd[2], ][-maxInd[2]]
    }

    groups <- getNames(factorStats)
    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")
    factorMerger$mergingList[["tmp"]] <- list(groups = groups,
                                              factor = factor,
                                              factorStats = factorStats,
                                              modelStats = NULL,
                                              means = calculateMeans(
                                                  factorMerger$response, factor))

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(factorMerger),
                   pval = maxStat)
    return(factorMerger)
    # TODO: replace for with apply
}

#' ----
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
