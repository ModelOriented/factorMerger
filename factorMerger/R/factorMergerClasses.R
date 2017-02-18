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
    # TODO: Make it insensitive for changes in input types
    fm <- list(
        response = response,
        factor = factor,
        mergingList = list()
    )

    class(fm) <- "factorMerger"

    if (gaussian) { # Here t-test or Hotelling test is used
        class(fm) <- append(class(fm), "gaussianFactorMerger")
    } else { # TODO: Here we'll use ...
        class(fm) <- append(class(fm), "non-gaussianFactorMerger")
    }

    if (!is.vector(response)) { # MANOVA
        class(fm) <- append(class(fm), "multivariateFactorMerger")
        subsequent <- FALSE
    }

    if (subsequent) {
        class(fm) <- append(class(fm), "subsequentFactorMerger")
    } else {
        class(fm) <- append(class(fm), "allToAllFactorMerger")
    }

    if (!is.vector(class)) { # TODO: ...
        class(fm) <- append(class(fm), "multiClassFactorMerger")
    }

    return(fm)
}

getCurrentGroups <- function(factorMerger) {
    UseMethod("getCurrentGroups", factorMerger)
}

getCurrentGroups.subsequentFactorMerger <- function(factorMerger) {
    # returns a 2xk matrix of groups in a current state of merging

}

getCurrentGroups.allToAllFactorMerger <- function(factorMerger) {
    # returns a square matrix of groups in a current state of merging
}

size <- function(object) {
    UseMethod("size", object)
}

size.numeric <- function(object) {
    return(length(object))
}

size.matrix <- function(object) {
    return(ncol(object) * nrow(object))
}

calculateCurrentSimilarities <- function(factorMerger) {
    UseMethod("calculateCurrentSimilarities", factorMerger)
} # A może tu nie trzeba przeciążać? Zastanowić się!

mergePair <- function(factorMerger) {
    UseMethod("mergePair", factorMerger)
}

mergePair.subsequentFactorMerger <- function(factorMerger) {

}

mergePair.allToAllFactorMerger <- function(factorMerger) {

}

mergePair.multiClassFactorMerger <- function(factorMerger) {

}

mergeFactors <- function(response, factor) {
    fm <- factorMerger(response, factor)
    groups <- getCurrentGroups(fm)
    while (size(groups) > 1) {
        fm <- calculateCurrentSimilarities(fm)
        fm <- mergePair(fm)
        groups <- getCurrentGroups(fm)
    }
    return(fm)
}
