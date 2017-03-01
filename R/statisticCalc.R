# Hotelling test {car}: https://github.com/cran/car/blob/1425129f002cb91a38951ad8af641254368a4d96/R/Anova.R
HL <- function (eig, q, df.res) {
    test <- sum(eig)
    p <- length(eig)
    m <- 0.5 * (abs(p - q) - 1)
    n <- 0.5 * (df.res - p - 1)
    s <- min(p, q)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * (s * n + 1)
    return(c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2))
}

# pvalue for Hotelling test {car}:
# https://github.com/cran/car/blob/1425129f002cb91a38951ad8af641254368a4d96/R/Anova.R
HLPval <- function(linHyp.mlm) {
    SSPE.qr <- qr(linHyp.mlm$SSPE)
    eigs <- Re(eigen(qr.coef(SSPE.qr, linHyp.mlm$SSPH), symmetric = FALSE)$values)
    hotellVec <- HL(eigs, linHyp.mlm$df, linHyp.mlm$df.residual)
    return(stats::pf(hotellVec[2],
                     hotellVec[3],
                     hotellVec[4],
                     lower.tail = FALSE))
}

calculateModel <- function(factorMerger, factor) {
    UseMethod("calculateModel", factorMerger)
}

calculateModel.gaussianFactorMerger <- function(factorMerger, factor) {
    if (length(unique(factor)) > 1) {
        return(lm(factorMerger$response ~ factor))
    }
    return(lm(factorMerger$response ~ 1))
}

getPvals <- function(model) {
    UseMethod("getPvals", model)
}

getPvals.lm <- function(model) {
    return(summary(model)$coefficient[-1, 4])
}

#' Get hypothesis p-values
#'
#' @importFrom car linearHypothesis
#'
getPvals.mlm <- function(model) {
    contr <- contr.treatment(NROW(model$coefficients))
    return(apply(contr, 2, function(x) { # czy to na pewno jest dobrze?
        hyp <- linearHypothesis(model, x)
        HLPval(hyp)
    }))
}

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

#' Gaussian model loglikelihood
#'
#' http://r.789695.n4.nabble.com/Multiple-regression-information-criterion-td4689113.html
logLikLm <- function(obj) {
    E <- obj$residuals
    S <- cov(matrix(E, nrow = NROW(E)))
    Sinv <- solve(S)
    n <- NROW(E)
    p <- NCOL(E)
    return(- n * p / 2 * log(2 * pi) - n / 2 * log(det(S)) -
               1/2 * sum(diag(E %*% Sinv %*% t(E))))
}

#' Calculate model statistic - ...
#'
calculateModelStatistic <- function(model) {
    UseMethod("calculateModelStatistic", model)
}

#' Calculate model statistic (gaussian case) - ...
#'
calculateModelStatistic.lm <- function(model) {
    return(logLikLm(model))
}
