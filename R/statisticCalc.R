calculateGIC <- function(model, k, penatly = 2) {
    return(-2 * calculateModelStatistic(model) + penatly * k)
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

#' @importFrom survival Surv coxph coxph.control
calculateModel.survivalFactorMerger <- function(factorMerger, factor) {
    if (length(unique(factor)) > 1) {
        return(coxph(factorMerger$response ~ factor,
                     control = coxph.control(iter.max = 50)))
    }
    return(coxph(factorMerger$response ~ 1,
                 control = coxph.control(iter.max = 50)))
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

calculateModelStatistic <- function(model) {
    UseMethod("calculateModelStatistic", model)
}

calculateModelStatistic.default <- function(model) {
    return(logLik(model))
}

# http://r.789695.n4.nabble.com/Multiple-regression-information-criterion-td4689113.html
calculateModelStatistic.mlm <- function(obj) {
    E <- obj$residuals
    S <- cov(matrix(E, nrow = NROW(E)))
    Sinv <- solve(S)
    n <- NROW(E)
    p <- NCOL(E)
    return(- n * p / 2 * log(2 * pi) - n / 2 * log(det(S)) -
               1/2 * sum(diag(E %*% Sinv %*% t(E))))
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
