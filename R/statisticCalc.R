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

#' @importFrom survival Surv coxph
calculateModel.survivalFactorMerger <- function(factorMerger, factor) {
    survObject <- Surv(time = factorMerger$response[, 1],
                       event = factorMerger$response[, 2])
    if (length(unique(factor)) > 1) {
        return(coxph(survObject ~ factor))
    }
    return(coxph(survObject ~ 1))
}

calculateModel.binomialFactorMerger <- function(factorMerger, factor) {
    df <- data.frame(response = factorMerger$response,
                     factor = factor)
    if (length(unique(factor)) > 1) {
        mod <- glm(response ~ factor, family = "binomial", data = df)

    } else {
        mod <- glm(response ~ 1, family = "binomial", data = df)
    }

    class(mod) <- c("binomglm", class(mod))
    return(mod)
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

#' Gaussian model loglikelihood
#'
#' http://r.789695.n4.nabble.com/Multiple-regression-information-criterion-td4689113.html
logLik.mlm <- function(obj) {
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
calculateModelStatistic.default <- function(model) {
    return(logLik(model))
}

calculateModelStatistic.coxph <- function(model) {
    if (length(model$loglik) > 1) {
        return(model$loglik[2])
    }
    return(model$loglik)
}

compareModels <- function(model1, model2) {
    UseMethod("compareModels", model1)
}

compareModels.lm <- function(model1, model2) {
    return(anova(model1, model2)$`Pr(>F)`[2])
}

compareModels.coxph <- function(model1, model2) {
    return(anova(model1, model2)$`P(>|Chi|)`[2])
}

compareModels.binomglm <- function(model1, model2) {
    return(anova(model1, model2, test = "Chisq")$`Pr(>Chi)`[2])
}
