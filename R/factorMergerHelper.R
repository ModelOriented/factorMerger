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

#' @importFrom MASS isoMDS
startMerging <- function(factorMerger, subsequent) {

    if (subsequent) {
        factorMerger$factor <- getIncreasingFactor(factorMerger)
    }

    factorMerger$mergingList[[1]]$groupStats <- calculateGroupStatistic(factorMerger, factorMerger$factor)
    factorMerger$mergingList[[1]]$groups <- levels(factorMerger$factor)
    model <- calculateModel(factorMerger, factorMerger$factor)
    factorMerger$mergingList[[1]]$modelStats <- data.frame(
        model = calculateModelStatistic(model),
        pval = 1,
        AIC = calculateAIC(model, length(levels(factorMerger$factor))))
    return(
        list(factorMerger = factorMerger,
             factor = factorMerger$factor,
             model = model)
        )
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

mergePair <- function(factorMerger, subsequent, factor, model) {
    step <- length(factorMerger$mergingList)
    fs <- factorMerger$mergingList[[step]]
    pairs <- getPairList(fs$groups, subsequent)
    # model <- fs$model
    modelsPvals <- sapply(pairs, function(x) {
        if (x[1] == x[2]) {
            return(-1)
        }
        tmpFactor <- mergeLevels(factor, x[1], x[2])
        tmpModel <- calculateModel(factorMerger, tmpFactor)
        return(compareModels(model, tmpModel))
    })

    factorMerger$mergingList <- c(factorMerger$mergingList,
                                  tmp = "tmp")

    whichMax <- which.max(modelsPvals)
    merged <- pairs[[whichMax]]
    factorMerger$mergingList[[step]]$merged <- merged
    factor <- mergeLevels(factor, merged[1], merged[2])
    model <- calculateModel(factorMerger, factor)

    factorMerger$mergingList[["tmp"]] <- list(groups = levels(factor),
                                              # factor = factor,
                                              modelStats = NULL,
                                              groupStats = calculateGroupStatistic(
                                                  factorMerger, factor))

    names(factorMerger$mergingList)[step + 1] <- step + 1

    factorMerger$mergingList[[step + 1]]$modelStats <-
        data.frame(model = calculateModelStatistic(model),
                   pval = modelsPvals[whichMax],
                   AIC = calculateAIC(model, length(levels(factor))))

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
