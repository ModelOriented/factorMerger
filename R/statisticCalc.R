#' Calculate pair statistic - ...
#'
calculatePairStatistic <- function(factorMerger, factor,
                                   groupA, groupB) {
    UseMethod("calculatePairStatistic", factorMerger)
}

#' Calculate pair statistic (gaussian case) - ...
#'
calculatePairStatistic.gaussianFactorMerger <- function(factorMerger, factor,
                                                        groupA, groupB) {
    groups <- filterGroups(factorMerger$response, factor, groupA, groupB)
    if (groupA == groupB) {
        return(-1)
    }
    return(t.test(groups[[1]], groups[[2]])$p.value)
}

#' Calculate pair statistic (multivariate response case) - ...
#'
#' @importFrom Hotelling hotel.test
#'

calculatePairStatistic.multivariateFactorMerger <- function(factorMerger, factor,
                                                            groupA, groupB) {
    groups <- filterGroups(factorMerger$response, factor, groupA, groupB)
    if (groupA == groupB) {
        return(-1)
    }
    return(hotel.test(groups[[1]], groups[[2]])$pval)
}

#' Calculate nonparametric pair statistic - ...
#'
#' @importFrom vegan adonis
#'

calculatePairStatistic.nonparametricFactorMerger <- function(factorMerger, factor,
                                                             groupA, groupB) {
    groups <- filterGroups(factorMerger$response, factor, groupA, groupB)
    if (groupA == groupB) {
        return(-1)
    }

    if (is.null(dim(groups[[1]]))) {
        X <- as.data.frame(c(groups[[1]], groups[[2]]))
    } else {
        X <- as.data.frame(rbind(groups[[1]], groups[[2]]))
    }

    g <- data.frame(group =
                        c(rep(groupA, NROW(groups[[1]])),
                          rep(groupB, NROW(groups[[2]])))
    )
    if (min(X) < 0) {
        X <- X - min(X)
    }

    return(adonis(X ~ group, data = g)$aov.tab$`Pr(>F)`[1])
}

# ---

#' Calculate model statistic - ...
#'
calculateModelStatistic <- function(factorMerger) {
    UseMethod("calculateModelStatistic", factorMerger)
}

#' Calculate model statistic (gaussian case) - ...
#'
calculateModelStatistic.gaussianFactorMerger <- function(factorMerger) {
    y <- factorMerger$response
    x <- factorMerger$mergingList[[length(factorMerger$mergingList)]]$factor
    if (length(levels(x)) > 1) {
        return(logLik(lm(y ~ x))[1])
    } else {
        return(logLik(lm(y ~ 1))[1])
    }
}

#' Calculate model statistic (multivariate response case) - ...
#'
calculateModelStatistic.multivariateFactorMerger <- function(factorMerger) {
    return(NA) # TODO
}

#' Calculate model statistic (nonparametric case) - ...
#'
calculateModelStatistic.nonparametricFactorMerger <- function(factorMerger) {
    return(NA) # TODO
}

# ---
