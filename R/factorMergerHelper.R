appendProjection <- function(factorMerger) {
    UseMethod("appendProjection", factorMerger)
}

appendProjection.factorMerger <- function(factorMerger) {
    return(factorMerger)
}

appendProjection.gaussianFactorMerger <- function(factorMerger) {
    if (NCOL(factorMerger$response) > 1) {
        tmpResponse <- MASS::isoMDS(dist(factorMerger$response), k = 1, trace = FALSE)$points[, 1]
        factorMerger$projectedResponse <- tmpResponse
    }
    return(factorMerger)
}

convertToDistanceMatrix <- function(modelsPvals, successive, labels) {
    if (successive) {
        m <- matrix(0, ncol = length(labels),
                    nrow = length(labels))
        tmp <- cbind(1:(nrow(m) - 1), 2:nrow(m))
        m[tmp] <- modelsPvals
        tmp <- cbind(2:nrow(m), 1:(nrow(m) - 1))
        m[tmp] <- modelsPvals
        modelsPvals <- m
    }

    modelsPvals <- modelsPvals %>% matrix(ncol = length(labels))
    colnames(modelsPvals) <- labels
    rownames(modelsPvals) <- labels
    distances <- as.dist(modelsPvals)
    distances[distances == 0] <- max(distances) + 1
    return(distances)
}

#' @importFrom MASS isoMDS
startMerging <- function(factorMerger, successive, method, penalty) {

    if (successive) {
        factorMerger$factor <- getIncreasingFactor(factorMerger)
    }
    factorMerger <- appendProjection(factorMerger)
    factor <- factorMerger$factor
    factorMerger$mergingList[[1]]$groupStats <-
        calculateGroupStatistic(factorMerger, factor)
    factorMerger$mergingList[[1]]$groups <- levels(factor)
    model <- calculateModel(factorMerger, factor)
    initStat <- calculateModelStatistic(model)
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = initStat,
        pval = 1,
        GIC = calculateGIC(model, length(levels(factor)), penalty))

    if (method == "LTR") {
        return(
            list(
                factorMerger = factorMerger,
                factor = factor,
                model = model
            )
        )
    }

    # method == "hclust"
    pairs <- getPairList(levels(factorMerger$factor), successive)
    modelsPvals <- sapply(pairs, function(x) {
        if (x[1] == x[2]) {
            return(1)
        }
        tmpFactor <- mergeLevels(factor, x[1], x[2])
        tmpModel <- calculateModel(factorMerger, tmpFactor)
        return(2 * initStat - 2 * calculateModelStatistic(tmpModel))
    })

    factorMerger$dist <- convertToDistanceMatrix(modelsPvals,
                                                 successive, levels(factorMerger$factor))

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
    twoLevelList <- lapply(groups, function(x) {lapply(groups, function(y) c(x, y))})
    return(unlist(twoLevelList, recursive = FALSE))
}

clusterFactors <- function(dist, successive) {
    if (successive) {
        return(hclust(d = dist, method = "single"))
    }
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

mergePairHClust <- function(factorMerger, factor, penalty) {
    step <- length(factorMerger$mergingList)
    merged <-  factorMerger$mergingHistory[step, ]
    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")

    factorMerger$mergingList[[step]]$merged <- merged
    factor <- mergeLevels(factor, merged[1], merged[2])
    model <- calculateModel(factorMerger, factor)

    factorMerger$mergingList[["tmp"]] <- list(groups = levels(factor),
                                              modelStats = NULL,
                                              groupStats = calculateGroupStatistic(
                                                  factorMerger, factor))

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pval = compareModels(calculateModel(factorMerger, factorMerger$factor), model),
                   GIC = calculateGIC(model, length(levels(factor)), penalty))

    return(
        list(factorMerger = factorMerger,
             factor = factor)
    )
}

mergePairLTR <- function(factorMerger, successive, factor, model, penalty) {
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

    factorMerger$mergingList[[step + 1]] <- list(groups = levels(factor),
                                              modelStats = NULL,
                                              groupStats = calculateGroupStatistic(
                                                  factorMerger, factor))

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pval = compareModels(calculateModel(factorMerger,
                                                       factorMerger$factor), model),
                   GIC = calculateGIC(model, length(levels(factor)), penalty))

    return(
        list(factorMerger = factorMerger,
             factor = factor,
             model = model)
    )
}

getFinalOrder <- function(factorMerger) {
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
    names(pos) <- groups
    return(pos)
}

#' @importFrom dplyr arrange
getFinalOrderVec <- function(factorMerger) {
    finalOrder <- data.frame(order = getFinalOrder(factorMerger))
    finalOrder$label <- rownames(finalOrder)
    finalOrder <- finalOrder %>% arrange(order)
    return(finalOrder$label %>% factor(levels = finalOrder$label))
}


getOptimalPartition <- function(factorMerger) {
    mH <- mergingHistory(factorMerger, T)
    nMerges <- which.min(mH$GIC) - 1
    factor <- factorMerger$factor
    if (nMerges == 0) {
        return(levels(factor))
    }
    for (i in 1:nMerges) {
        factor <- mergeLevels(factor, mH$groupA[i + 1], mH$groupB[i + 1])
    }
    return(levels(factor))
}
