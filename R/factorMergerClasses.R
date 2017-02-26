#' Factor Merger - ...
#'
#' @importFrom data.table data.table
#'
#' @export
#'
merger <- function(response, factor, gaussian = TRUE,
                         subsequent = FALSE) {

    if (NROW(response) != NROW(factor)) {
        stop("Response and factor sizes do not match.")
    }

    factor <- as.factor(factor)
    # TODO: Make it insensitive to input types changes
    fm <- list(
        response = response,
        factor = factor,
        mergingList = list(`1` = list(groups = levels(factor),
                                factor = factor,
                                factorStats = list(),
                                modelStats = list(),
                                means = NA))
    )

    class(fm) <- "factorMerger"

    if (gaussian) { # Here t-test or Hotelling test is used

        if (!is.null(dim(response))) { # MANOVA
            class(fm) <- append(class(fm), "multivariateFactorMerger")
            if (subsequent) {
                warning("Subsequent merging is meaningless with multivariate response. All-to-all merging run instead.")
                subsequent <- FALSE
            }
        } else {
            fm$mergingList[[1]]$means <- calculateMeans(response, factor)
        }
        class(fm) <- append(class(fm), "gaussianFactorMerger")

    } else { # TODO: Here we'll use adonis{vegan}
        class(fm) <- append(class(fm), "nonparametricFactorMerger")
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

#' Factor Merger - ...
#'
#' @export
#'
print.factorMerger <- function(factorMerger) {
    lapply(factorMerger$mergingList, function(x) {
        cat("Merging p-value: ")
        cat(ifelse(is.na(x$modelStats$pval), " NULL", round(x$modelStats$pval, 3)))
        cat(", groups: ")
        cat(paste(x$groups, collapse = "|"))
        cat("\n")
        })
}

# ---

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

#' Show levels means - ...
#'
#' @export
#'
means <- function(object) {
    UseMethod("means", object)
}

#' ---
means.factorMerger <- function(factorMerger) {
    statsList <- lapply(factorMerger$mergingList,
                        function(x) { as.data.frame(x$means) })
    do.call(rbind, statsList) %>% unique()
}

#' Merge factors - ...
#'
#' @export
#'
mergeFactors <- function(response, factor, gaussian = TRUE, subsequent = FALSE) {

    if (is.null(response)) {
        stop('argument "response" is missing, with no default.')
    }
    if (is.null(factor)) {
        stop('argument "factor" is missing, with no default.')
    }

    fm <- merger(response, factor, gaussian, subsequent)
    fm <- startMerging(fm)
    while (canBeMerged(fm)) {
        fm <- mergePair(fm)
    }
    return(fm)
}
