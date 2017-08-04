#' @importFrom magrittr "%>%"
#' @importFrom graphics text
#' @importFrom stats aggregate anova aov as.dist ave
#'           coef complete.cases glm hclust lm logLik
#'           predict proj qchisq quantile rbeta rbinom
#'           reshape residuals rexp rnorm runif step
#' @importFrom utils head tail
#' @importFrom proxy dist
#' @importFrom agricolae HSD.test
#' @importFrom dplyr select_ seqtequal
NULL

# Clean this below!
globalVariables(c("x1", "x2", "y1", "y2", "y0",
                  "pred", "group", "y",
                  "variable", "value", "stat", "significance",
                  "level", "hjust", "vjust", "label", "left", "n", "right",
                  "L1", "V0", "V05", "V1", "V2",
                  "xmin", "xpos", "y100", "y25", "xmax",
                  "y50", "y75", "ymax", "ymin", "ypos"))

cleanFactor <- function(factor) {
    factor <- as.factor(factor)
    levs <- levels(factor)
    factor <- droplevels(factor)
    if (!setequal(levels(factor), levs)) {
        warning("Dropped missing levels of the factor.")
    }
    return(factor)
}

merger <- function(response, factor,
                   family = "gaussian",
                   abbreviate) {

    stopifnot(NROW(response) == NROW(factor))

    factor <- cleanFactor(factor)

    if (abbreviate) {
        map <- data.frame(
            `recoded` = paste0("(", abbreviate(levels(factor)), ")"),
            `original` = levels(factor))
        rownames(map) <- NULL
        factor <- factor(factor, labels = map$recoded)

    }

    fm <- list(
        response = response,
        factor = factor,
        mergingList = list(`1` = list(groups = levels(factor),
                                modelStats = list(),
                                groupStats = NA,
                                merged = NA))
    )

    if (abbreviate) {
        fm[["map"]] <- map
    }

    class(fm) <- "factorMerger"

    switch(family,
           "gaussian" = {
               class(fm) <- append(class(fm), "gaussianFactorMerger")
           },

           "survival" = {
               stopifnot(class(response) == "Surv")
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

    if (NCOL(factor) > 1) {
        class(fm) <- append(class(fm), "multiClassFactorMerger")
        stop("Factor merging with multivariate factor is not supported yet.")
    }

    return(fm)
}

stats <- function(factorMerger) {
    statsList <- lapply(factorMerger$mergingList, function(x) x$modelStats)
    do.call(rbind, statsList)
}

#' Groups statistic
#'
#' @description Summary of statistics specific for a model for each group that
#' appeared in merging.
#'
#' @param factorMerger object of a class \code{factorMerger}
#'
#' @examples
#' randSample <- generateMultivariateSample(N = 100, k = 10, d = 3)
#' fm <- mergeFactors(randSample$response, randSample$factor)
#' groupsStats(fm)
#'
#' @export
groupsStats <- function(factorMerger) {
    statsList <- lapply(factorMerger$mergingList,
                        function(x) as.data.frame(x$groupStats))
    statsDf <- do.call(rbind, statsList)
    statsDf <- subset(statsDf, !duplicated(level))

    if (sum(complete.cases(statsDf)) == 0) {
        return(NULL)
    }

    rownames(statsDf) <- statsDf$level
    statsDf <- subset(statsDf, select = -level)
    return(statsDf)
}

#' Merging history
#'
#' @description Summarizes merging path by giving pairs of factor groups merged
#' in each iteration.
#' @export
#'
#' @param factorMerger Object of a class \code{factorMerger}
#' @param showStats If \code{TRUE} extends results with
#' the loglikelihood (column \code{model}),
#' p-value for the \code{LRT} tests against the full model (column \code{pval})
#' and Generalized Information Criterion value (column \code{GIC}).
#' By default \code{showStats} is set to \code{FALSE}.
#' @param round Logical. If \code{TRUE}, the default, statistics are rounded
#' @param penalty GIC penalty
#'
#' @examples
#' randSample <- generateMultivariateSample(N = 100, k = 10, d = 3)
#' fm <- mergeFactors(randSample$response, randSample$factor)
#' mergingHistory(fm, showStats = TRUE)
#'
#' @importFrom dplyr rename
mergingHistory <- function(factorMerger, showStats = FALSE,
                           penalty, round = TRUE) {
    mergingList <- sapply(factorMerger$mergingList,
                        function(x)  x$merged )
    mergingDf <- do.call(rbind, mergingList) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        rename(groupA = V1, groupB = V2)

    rownames(mergingDf) <- NULL
    if (showStats) {
        st <- stats(factorMerger)
        if (round) {
            st <- round(stats(factorMerger), 4)
        }
        mergingDf <- mergingDf[complete.cases(mergingDf), ]
        mergingDf <- rbind(c("", ""), mergingDf)
        mergingDf <- data.frame(mergingDf, st)
    }

    if (!missing(penalty)) {
        mergingDf$GIC <- -2 * mergingDf$model +
            penalty * nrow(mergingDf):1
    }

    return(mergingDf)
}

call <- function(factorMerger) {
    return(
        paste0("Family: ",
               gsub("([[:upper:]])", " \\1",
                    class(factorMerger)[length(class(factorMerger))]), ".")
        )
}


#' factorMerger
#'
#' @description \code{factorMerger} is the
#' base class of the factorMerger package. \code{factorMerger} stores information about
#' response, initial factor, its levels and their abbreviated names (field \code{map}).
#' \code{factorMerger} creates its own structure of inheritance connected with model family.
#'
#' When merging is applied, \code{factorMerger} shows which levels have been
#' merged together with
#' the matching summary statistics:
#' model loglikelihood, pvalue for the \code{LRT} test
#' against the full model and Generalized Information Criterion value.
#'
#' @export
#'
#' @param x object of a class \code{factorMerger}.
#' @param ... Other arguments
#'
#' @importFrom knitr kable
#'
print.factorMerger <- function(x, ...) {
   df <- mergingHistory(x, showStats = TRUE)
   colnames(df)[1:2] <- c("groupA", "groupB")
   cat(call(x))

   if ("map" %in% names(x)) {
       cat("\nFactor levels were recoded as below:")
       cat(paste(c("", "", kable(x$map, output = FALSE)), collapse = "\n"))
   }

   cat("\n\nFactor levels were merged in the following order:")
   cat(paste(c("", "", kable(df, output = FALSE)), collapse = "\n"))
   invisible(NULL)
}


#' Merge factors
#'
#' @description Performs step-wise merging of factor levels.
#'
#' @param response A response \code{vector/matrix} suitable for the model family.
#' @param factor A factor \code{vector}.
#' @param family Model family to be used in merging. Available models are: \code{"gaussian",}
#' \code{ "survival", "binomial"}.
#' By default \code{mergeFactors} uses \code{"gaussian"} model.
#' @param successive If \code{FALSE}, the default,
#' in each step of the merging procedure all possible pairs are compared.
#' Otherwise, factor levels are preliminarly sorted and only succesive pairs are compared.
#' @param method A string specifying method used during merging.
#' Two methods are availabel: \code{"hclust", "LRT"}. The default is \code{"LRT"}.
#' @param abbreviate Logical. If \code{TRUE}, the default, factor levels names
#' are abbreviated.
#'
#' @examples
#' randSample <- generateMultivariateSample(N = 100, k = 10, d = 3)
#' mergeFactors(randSample$response, randSample$factor)
#'
#' @export
#'
mergeFactors <- function(response, factor,
                         family = "gaussian",
                         successive = FALSE,
                         method = "LRT",
                         abbreviate = TRUE) {

    stopifnot(!is.null(response), !is.null(factor))


    if (is.data.frame(response)) {
        response <- as.matrix(response)
    }

    fm <- merger(response, factor, family, abbreviate)

    if (method == "LRT") {
        return(mergeLRT(fm, successive))
    }

    if (method == "hclust") {
        return(mergeHClust(fm, successive))
    }

    else {
        stop("Requested method of merging is not supported.")
    }
}

mergeLRT <- function(factorMerger, successive) {
    fmList <- startMerging(factorMerger, successive, "LRT")
    fm <- fmList$factorMerger

    while (canBeMerged(fm)) {
        fmList <- mergePairLRT(fm, successive, fmList$factor,
                               fmList$model)
        fm <- fmList$factorMerger
    }

    return(fmList$factorMerger)
}

getPairWithLowestDist <- function(distance, pos) {
    return(distance[pos, ] %>%
               select_("firstClusterLabel", "lastClusterLabel"))
}

replaceWithMergedPair <- function(distance, pos, toBeMerged) {
    mergedLabel <- paste0(toBeMerged, collapse = "")
    if (pos > 1) {
        distance[pos - 1, "lastClusterLabel"] <- mergedLabel
    }
    if (pos < nrow(distance)) {
        distance[pos + 1, "firstClusterLabel"] <- mergedLabel
    }
    return(distance)
}

updateRowAfterMerging <- function(distance, newPos, origPos,
                                  factorMerger, which) {
    initStat <- factorMerger$mergingList[[1]]$modelStats$model
    distance[newPos, which] <- distance[origPos, which]
    newPair <- distance[newPos, c("first", "last")]
    tmpFactor <- mergeLevels(factorMerger$factor,
                             newPair[1], newPair[2])
    tmpModel <- calculateModel(factorMerger, tmpFactor)
    distance[newPos, "dist"] <-
        2 * initStat - 2 * calculateModelStatistic(tmpModel)
    return(distance)
}

updateDistAfterMerging <- function(factorMerger, pos, toBeMerged) {
    distance <- factorMerger$dist
    if (pos > 1) {
        distance <- updateRowAfterMerging(distance, pos - 1, pos,
                                          factorMerger, "last")
    }
    if (pos < nrow(distance)) {
        distance <- updateRowAfterMerging(distance, pos + 1, pos,
                                          factorMerger, "first")
    }
    distance <- distance[-pos, ]
    return(distance)
}

addNewPairToMergingHistory <- function(factorMerger, toBeMerged) {
    toBeMerged <- unname(toBeMerged)
    if (length(factorMerger$mergingList) == 1) {
        return(as.matrix(toBeMerged))
    }
    mH <- factorMerger$mergingHistory
    mH <- rbind(mH, unname(as.matrix(toBeMerged)))
    return(mH)
}

mergeSuccessiveHClust <- function(factorMerger) {
    factor <- factorMerger$factor
    for (step in 1:(length(levels(factorMerger$factor)) - 1)) {
        # Take pair p of A and B with the lowest distance
        mergingPositionInTable <- which.min(factorMerger$dist$dist)
        toBeMerged <- getPairWithLowestDist(factorMerger$dist,
                                            mergingPositionInTable)
        # mergePairHClust with p (add to the merging list)
        factorMerger$mergingHistory <-
            addNewPairToMergingHistory(factorMerger, toBeMerged)
        fm <- mergePairHClust(factorMerger, factor)
        factorMerger <- fm$factorMerger
        factor <- fm$factor
        # Replace clusterLabels in the factorMerger$dist
        factorMerger$dist <- replaceWithMergedPair(factorMerger$dist,
                                                   mergingPositionInTable,
                                                   toBeMerged)

        factorMerger$dist <- updateDistAfterMerging(factorMerger,
                                                   mergingPositionInTable,
                                                   toBeMerged)
    }
    return(factorMerger)
}

mergeHClust <- function(factorMerger, successive) {
    factorMerger <- startMerging(factorMerger, successive, "hclust")
    if (successive) {
        return(mergeSuccessiveHClust(factorMerger))
    }

    # Original DMR
    clust <- clusterFactors(factorMerger$dist)
    factorMerger$mergingHistory <-
        recodeClustering(clust$merge,
                         clust$labels,
                         getIncreasingFactor(factorMerger))

    factor <- factorMerger$factor
    for (i in 1:nrow(factorMerger$mergingHistory)) {
        fm <- mergePairHClust(factorMerger, factor)
        factorMerger <- fm$factorMerger
        factor <- fm$factor
    }
    return(factorMerger)
}
