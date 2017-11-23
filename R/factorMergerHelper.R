appendProjection <- function(factorMerger) {
    UseMethod("appendProjection", factorMerger)
}

appendProjection.default <- function(factorMerger) {
    return(factorMerger)
}

#' @importFrom dplyr filter left_join
reverseOrder <- function(factorMerger, isoMDSproj) {
    sim <- findSimilarities(factorMerger)
    colnames(isoMDSproj)[1] <- "proj"
    sim <- sim %>% filter(variable == levels(variable)[1]) %>%
        left_join(isoMDSproj, by = c("level" = "factor"))
    meansIncreasing <- sim[1, "mean"] < sim[nrow(sim), "mean"]
    projIncreasing <- sim[1, "proj"] < sim[nrow(sim), "proj"]
    if (xor(meansIncreasing, projIncreasing)) {
        return(-isoMDSproj$proj)
    }
    return(isoMDSproj$proj)
}

#' @importFrom dplyr left_join
appendProjection.gaussianFactorMerger <- function(factorMerger) {
    if (NCOL(factorMerger$response) > 1) {
        groupMeans <- calculateMeans(factorMerger$response,
                                     factorMerger$factor)
        tmpResponse <- MASS::isoMDS(dist(groupMeans[, -1]),
                                    k = 1, trace = FALSE)$points[, 1] %>%
            as.data.frame()
        tmpResponse$factor <- groupMeans$level
        tmpResponse[, 1] <- reverseOrder(factorMerger, tmpResponse)
        tmpResponse <- data.frame(factor = factorMerger$factor,
                                  stringsAsFactors = FALSE) %>%
            left_join(tmpResponse,
                      by = "factor")
        factorMerger$projectedResponse <- tmpResponse[, 2]
    }
    return(factorMerger)
}

convertToDistanceTensor <- function(modelsPvals, successive, labels) {
    if (successive) {
        modelsPvals <- data.frame(dist = modelsPvals,
                                  first = labels[-length(labels)],
                                  last = labels[-1],
                                  firstClusterLabel = labels[-length(labels)],
                                  lastClusterLabel = labels[-1],
                                  stringsAsFactors = FALSE)
        return(modelsPvals)
    }

    modelsPvals <- modelsPvals %>% matrix(ncol = length(labels))
    colnames(modelsPvals) <- labels
    rownames(modelsPvals) <- labels
    distances <- as.dist(modelsPvals)
    distances[distances == 0] <- max(distances) + 1
    return(distances)
}

#' @importFrom MASS isoMDS
startMerging <- function(factorMerger, successive, method) {

    factorMerger <- appendProjection(factorMerger)
    factorMerger$factor <- getIncreasingFactor(factorMerger)
    factor <- factorMerger$factor
    factorMerger$mergingList[[1]]$groupStats <-
        calculateGroupStatistic(factorMerger, factor)
    factorMerger$mergingList[[1]]$groups <- levels(factor)
    model <- calculateModel(factorMerger, factor)
    initStat <- calculateModelStatistic(model)
    factorMerger$initialModel <- model
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = initStat,
        pvalVsFull = 1,
        pvalVsPrevious = 1)

    if (method == "LRT") {
        return(
            list(
                factorMerger = factorMerger,
                factor = factor,
                model = model
            )
        )
    }

    # method "hclust"
    pairs <- getPairList(levels(factorMerger$factor), successive)
    modelsPvals <- sapply(pairs, function(x) {
        if (x[1] == x[2]) {
            return(1)
        }
        tmpFactor <- mergeLevels(factor, x[1], x[2])
        tmpModel <- calculateModel(factorMerger, tmpFactor)
        return(2 * initStat - 2 * calculateModelStatistic(tmpModel))
    })

    factorMerger$dist <- convertToDistanceTensor(modelsPvals,
                                                 successive,
                                                 levels(factorMerger$factor))

    return(factorMerger)
}

canBeMerged <- function(factorMerger) {
    ml <- factorMerger$mergingList
    return(length(ml[[length(ml)]]$groups) > 1)
}

getPairList <- function(groups, successive) {
    if (successive) {
        getsuccessivePairList(groups)
    } else {
        getAllPairList(groups)
    }
}

getsuccessivePairList <- function(groups) {
    noGroups <- length(groups)
    pairs <- t(cbind(groups[1:(noGroups - 1)], groups[2:(noGroups)]))
    return(split(pairs, rep(1:ncol(pairs), each = nrow(pairs))))
}

getAllPairList <- function(groups) {
    twoLevelList <- lapply(groups, function(x) {
        lapply(groups, function(y) c(x, y))
        })
    return(unlist(twoLevelList, recursive = FALSE))
}

clusterFactors <- function(dist) {
    return(hclust(d = dist, method = "complete"))
}

getLabel <- function(currentLabels, levels, num) {
    if (num < 1) {
        return(levels[-num])
    }
    return(paste(currentLabels[num, ], collapse = ""))
}

recodeClustering <- function(merge, levels, factor) {
    res <- matrix(NA, ncol = ncol(merge), nrow = nrow(merge))
    tmpLevels <- levels
    for (row in 1:nrow(merge)) {
        tmp1 <- getLabel(res, levels, merge[row, 1])
        tmp2 <- getLabel(res, levels, merge[row, 2])
        if (which(tmpLevels == tmp1) > which(tmpLevels == tmp2)) {
            t <- tmp1
            tmp1 <- tmp2
            tmp2 <- t
        }
        factor <- mergeLevels(factor, tmp1, tmp2)
        tmpLevels <- levels(factor)
        res[row, 1] <- tmp1
        res[row, 2] <- tmp2
    }
    return(res)
}

mergePairHClust <- function(factorMerger, factor) {
    step <- length(factorMerger$mergingList)
    merged <-  factorMerger$mergingHistory[step, ]
    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")

    factorMerger$mergingList[[step]]$merged <- merged
    prevModel <- calculateModel(factorMerger, factor)
    factor <- mergeLevels(factor, merged[1], merged[2])
    model <- calculateModel(factorMerger, factor)

    factorMerger$mergingList[["tmp"]] <-
        list(groups = levels(factor),
             modelStats = NULL,
             groupStats = calculateGroupStatistic(
                 factorMerger, factor))

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pvalVsFull = compareModels(factorMerger$initialModel, model),
                   pvalVsPrevious = compareModels(prevModel, model))

    return(
        list(factorMerger = factorMerger,
             factor = factor)
    )
}

mergePairLRT <- function(factorMerger, successive, factor, model) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    pairs <- getPairList(fs$groups, successive)
    modelsPvals <- sapply(pairs, function(x) {
        if (x[1] == x[2]) {
            return(-1)
        }
        tmpFactor <- mergeLevels(factor, x[1], x[2])
        tmpModel <- calculateModel(factorMerger, tmpFactor)
        return(compareModels(model, tmpModel))
    })

    whichMax <- which.max(modelsPvals)
    pval <- modelsPvals[whichMax]
    merged <- pairs[[whichMax]]
    factorMerger$mergingList[[step]]$merged <- merged
    factor <- mergeLevels(factor, merged[1], merged[2])
    model <- calculateModel(factorMerger, factor)

    factorMerger$mergingList[[step + 1]] <-
        list(groups = levels(factor),
             modelStats = NULL,
             groupStats = calculateGroupStatistic(
                 factorMerger, factor))

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pvalVsFull = compareModels(factorMerger$initialModel, model),
                   pvalVsPrevious = pval)

    return(
        list(factorMerger = factorMerger,
             factor = factor,
             model = model)
    )
}

# if reverse == TRUE, a group with the highest statistic will be the last one
getFinalOrder <- function(factorMerger, reverse = FALSE) {
    groups <- levels(factorMerger$factor)
    merging <- mergingHistory(factorMerger)
    noSteps <- nrow(merging)
    pos <- rep(1, length(groups))
    names(pos) <- groups
    for (step in 1:noSteps) {
        pos[names(pos) == merging[step, 2]] <-
            pos[names(pos) == merging[step, 2]] +
            max(pos[names(pos) == merging[step, 1]])
        names(pos)[names(pos) %in% merging[step, ]] <-
            paste(merging[step, ], collapse = "")
    }
    # Sort (approximately) according to group means
    stats <- groupsStats(factorMerger)
    stats$group <- rownames(stats)
    stats <- data.frame(group = groups, stringsAsFactors = F) %>%
        left_join(stats, by = "group")

    if (stats[1, 2] > stats[nrow(stats), 2] && reverse) {
        names(pos) <- groups[length(groups):1]
    } else {
        names(pos) <- groups
    }

    return(pos)
}

getStatNameInTable <- function(stat) {
    switch(stat,
           "loglikelihood" = return("model"),
           "p-value" = return("pvalVsFull"),
           "GIC" = return("GIC"))
}

#' @importFrom dplyr arrange
getFinalOrderVec <- function(factorMerger) {
    finalOrder <- data.frame(order = getFinalOrder(factorMerger))
    finalOrder$label <- rownames(finalOrder)
    finalOrder <- finalOrder %>% arrange(order)
    return(finalOrder$label %>% factor(levels = finalOrder$label))
}

#' Cut a Factor Merger Tree
#'
#' @description Splits factor levels into non-overlapping clusters
#' based on a \code{factorMerger} object.
#' If a \code{stat} is \code{"loglikelihood"} or {"p-value"}
#' then performs bottom-up search through models
#' on the merging path until spots a model scored worse
#' than the given threshold (\code{value}).
#' If \code{stat = "GIC"}, \code{value} is interpreted as
#' GIC penalty and optimal GIC model is returned..
#'
#' @param factorMerger object of a class \code{factorMerger}
#' @param stat statistic used in the bottom-up search. Available statistics are:
#' \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}.
#' @param value cut threshold or penalty (for GIC)
#'
#' @details By default, \code{cutree} returns factor partition
#' corresponding to the optimal GIC model (with the lowest GIC).
#'
#' @return Returns a factor vector - each observation is given a new cluster label.
#'
#' @export
cutTree <- function(factorMerger,
                    stat = "GIC",
                    value = 2) {
    stopifnot(!is.null(value) | stat == "GIC")
    stopifnot(stat != "GIC" | value > 0)
    mH <- mergingHistory(factorMerger, T)
    stopifnot(stat %in% c("loglikelihood", "p-value", "GIC"))
    if (stat == "GIC") {
        mH$GIC <- -2 * mH$model + as.numeric(value) * nrow(mH):1
        value <- min(mH$GIC)
    }
    statColname <- getStatNameInTable(stat)

    factor <- factorMerger$factor
    nMerges <- nrow(mH)
    if (nMerges < 2 | mH[1, statColname] == value) {
        return(factor)
    }

    for (i in 2:nMerges) {
        if (mH[i, statColname] >= value) {
            factor <- mergeLevels(factor, mH$groupA[i], mH$groupB[i])
        }
        else {
            return(factor)
        }
        if (stat == "GIC" && mH[i, stat] == value) {
            return(factor)
        }
    }
    return(factor)
}

#' Get optimal partition (clusters dictionary)
#'
#' @description Splits factor levels into non-overlapping
#' clusters based on a \code{factorMerger} object.
#' If a \code{stat} is \code{"loglikelihood"} or {"p-value"}
#' then performs bottom-up search through models
#' on the merging path until spots a model scored worse
#' than the given threshold (\code{value}).
#' If \code{stat = "GIC"}, \code{value} is interpreted as
#' GIC penalty and optimal GIC model is returned.
#'
#' @param factorMerger object of a class \code{factorMerger}
#' @param stat statistic used in the bottom-up search. Available statistics are:
#' \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}.
#' @param value cut threshold / GIC penalty
#'
#' @details By default, \code{cutree} returns factor partition
#' corresponding to the optimal GIC model (with the lowest GIC).
#'
#' @return Returns a dictionary in a data frame format.
#' Each row gives an original label of a factor level and its new (cluster) label.
#'
#' @export
getOptimalPartitionDf <- function(factorMerger,
                                  stat = "GIC",
                                  value = 2) {
    return(data.frame(orig = factorMerger$factor,
                      pred = cutTree(factorMerger, stat, value),
                      stringsAsFactors = FALSE) %>% unique())
}

#' Get optimal partition (final clusters names)
#'
#' @description Splits factor levels into non-overlapping
#' clusters based on a \code{factorMerger} object.
#' If a \code{stat} is \code{"loglikelihood"} or {"p-value"}
#' then performs bottom-up search through models
#' on the merging path until spots a model scored worse than
#' the given threshold (\code{value}).
#' If \code{stat = "GIC"}, \code{value} is interpreted as
#' GIC penalty and optimal GIC model is returned.
#'
#' @param factorMerger object of a class \code{factorMerger}
#' @param stat statistic used in the bottom-up search. Available statistics are:
#' \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}.
#' @param value cut threshold / GIC penalty
#'
#' @details By default, \code{cutree} returns factor partition
#' corresponding to the optimal GIC model (with the lowest GIC).
#'
#' @return Returns a vector with the final cluster names from the \code{factorMerger} object.
#'
#' @export
getOptimalPartition <- function(factorMerger,
                                stat = "GIC",
                                value = 2) {
    factor <- cutTree(factorMerger, stat, value)
    return(levels(factor))
}
