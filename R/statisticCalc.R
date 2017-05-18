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
        return(coxph(factorMerger$response ~ factor))
    }
    return(coxph(factorMerger$response ~ 1))
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

# https://rdrr.io/rforge/Atools/src/R/logLik.mlm.R
#' @importFrom mvtnorm dmvnorm
calculateModelStatistic.mlm <- function(obj) {

    resids <- residuals(obj)
    n <- nrow(resids)
    Sigma_ML <- crossprod(resids) /n
    ans <- sum(dmvnorm(resids, sigma = Sigma_ML, log = T))
    return(ans %>% as.numeric())
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
