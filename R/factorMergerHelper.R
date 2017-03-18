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

#' @importFrom MASS isoMDS
startMerging <- function(factorMerger, subsequent) {
    if (subsequent) {
        if (NCOL(factorMerger$response) > 1) {
            tmpResponse <- MASS::isoMDS(dist(factorMerger$response), k = 1, trace = FALSE)$points[, 1]
            factorMerger$projectedResponse <- tmpResponse
        } else {
            tmpResponse <- factorMerger$response
        }
            factorMerger$factor <- setIncreasingOrder(tmpResponse, factorMerger$factor)
    }

    factorMerger$mergingList[[1]]$groupStats <- calculateGroupStatistic(factorMerger, factorMerger$factor)
    factorMerger$mergingList[[1]]$factor <- factorMerger$factor
    factorMerger$mergingList[[1]]$groups <- levels(factorMerger$mergingList[[1]]$factor)
    model <- calculateModel(factorMerger, factorMerger$factor)
    factorMerger$mergingList[[1]]$model <- model
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = calculateModelStatistic(model),
        pval = 1)
    return(factorMerger)
}

canBeMerged <- function(factorMerger) {
    ml <- factorMerger$mergingList
    return(length(ml[[length(ml)]]$groups) > 1)
}

getPairList <- function(groups, subsequent) {
    if (subsequent) {
        getSubsequentPairList(groups)
    } else {
        getAllPairList(groups)
    }
}

getSubsequentPairList <- function(groups) {
    noGroups <- length(groups)
    pairs <- t(cbind(groups[1:(noGroups - 1)], groups[2:(noGroups)]))
    return(split(pairs, rep(1:ncol(pairs), each = nrow(pairs))))
}

getAllPairList <- function(groups) {
    twoLevelList <- lapply(groups, function(x) {lapply(groups, function(y) c(x, y))})
    return(unlist(twoLevelList, recursive = FALSE))
}

mergePair <- function(factorMerger, subsequent) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    pairs <- getPairList(fs$groups, subsequent)
    model <- fs$model
    modelsPvals <- sapply(pairs, function(x) {
        factor <- mergeLevels(fs$factor, x[1], x[2])
        tmpModel <- calculateModel(factorMerger, factor)
        return(compareModels(model, tmpModel))
    })

    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")

    whichMax <- which.max(modelsPvals)
    merged <- pairs[[whichMax]]
    factorMerger$mergingList[[step]]$merged <- merged
    factor <- mergeLevels(fs$factor, merged[1], merged[2])
    model <- calculateModel(factorMerger, factor)

    factorMerger$mergingList[["tmp"]] <- list(groups = levels(factor),
                                              factor = factor,
                                              modelStats = NULL,
                                              groupStats = calculateGroupStatistic(
                                                  factorMerger, factor),
                                              model = model)

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pval = modelsPvals[whichMax])

    return(factorMerger)

}

#' @export
getTreeWithEdgesLength <- function(factorMerger, stat) {
    nodes <- list()
    initLevels <- levels(factorMerger$factor)
    initStat <- factorMerger$mergingList$`1`$modelStats[stat]
    for (level in initLevels) {
        nodes <- c(nodes, list(node(level, stat = initStat)))
    }
    names(nodes) <- initLevels

    ml <- factorMerger$mergingList
    for (i in 1:(length(ml) - 1)) {
        merged <- ml[[i]]$merged
        nodes <- c(nodes, list(node(nodes[[merged[1]]], nodes[[merged[2]]],
                                    stat = ml[[i + 1]]$modelStats[stat])))
        names(nodes)[length(nodes)] <- paste0(merged[1], merged[2])
    }

    return(paste0(nodes[[length(nodes)]]$text, ":", nodes[[length(nodes)]]$stat, ";"))
}

getFinalOrder <- function(factorMerger) {
    groups <- levels(factorMerger$factor)
    merging <- data.frame(mergingHistory(factorMerger), stringsAsFactors = FALSE)
    noSteps <- nrow(merging)
    pos <- rep(1, length(groups))
    names(pos) <- groups
    for (step in 1:noSteps) {
        pos[names(pos) == merging[step, 2]] <-
            pos[names(pos) == merging[step, 2]] +
            max(pos[names(pos) == merging[step, 1]])
        names(pos)[names(pos) == merging[step, 1]] <-
            paste0(merging[step, 1], merging[step, 2])
        names(pos)[names(pos) == merging[step, 2]] <-
            paste0(merging[step, 1], merging[step, 2])
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
