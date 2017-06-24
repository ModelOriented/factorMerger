calculateGIC <- function(model, k, penatly = 2) {
    return(-2 * calculateModelStatistic(model) + penatly * k)
}

calculateModel <- function(factorMerger, factor) {
    UseMethod("calculateModel", factorMerger)
}

calculateModel.gaussianFactorMerger <- function(factorMerger, factor) {
    if (length(unique(factor)) > 1) {
        return(lm(factorMerger$response ~ factor - 1))
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
    Sigma_ML <- crossprod(resids) / n # sample covariance (https://en.wikipedia.org/wiki/Sample_mean_and_covariance)
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

calculateAnovaTable <- function(model) {
    UseMethod("calculateAnovaTable", model)
}

formatSingleRow <- function(num) {
    if (is.na(num)) {
        return("")
    }
    if (num == 0) {
        return("< 2.2e-16")
    }
    return(
        ifelse(num < 0.1, format(num, scientific = T), format(num, scientific = F))
    )
}

formatPvalue <- function(col) {
    formatCol <- Vectorize(formatSingleRow)
    return(
        formatCol(col)
    )
}

calculateAnovaTable.lm <- function(model) {
    anTable <- anova(model) %>% data.frame()
    anTable[, -ncol(anTable)] <- round(anTable[, -ncol(anTable)], 1)
    anTable[, ncol(anTable)] <- formatPvalue(anTable[, ncol(anTable)])
    anTable <- anTable[, c("Df", "F.value", "Pr..F.")]
    colnames(anTable)[2:3] <- c("F", "p-value")
    rownames(anTable)[2] <- "Res"
    return(anTable)
}

calculateAnovaTable.mlm <- function(model) {
    anTable <- anova(model) %>% data.frame()
    anTable[, -ncol(anTable)] <- round(anTable[, -ncol(anTable)], 1)
    anTable[, ncol(anTable)] <- formatPvalue(anTable[, ncol(anTable)])
    anTable <- anTable[, -(4:5)]
    colnames(anTable)[4] <- c("p-value")
    rownames(anTable)[2] <- "Res"
    return(anTable)
}

calculateAnovaTable.binomglm <- function(model) {
    anTable <- anova(model, test = "Chisq") %>% data.frame()
    anTable[, -ncol(anTable)] <- round(anTable[, -ncol(anTable)], 1)
    anTable[, ncol(anTable)] <- formatPvalue(anTable[, ncol(anTable)])
    anTable <- anTable[, -(2:3)]
    colnames(anTable)[2:3] <- c("ResDev", "p-value")
    return(anTable)
}

calculateAnovaTable.coxph <- function(model) {
    anTable <- anova(model, test = "Chisq") %>% data.frame()
    anTable[, -ncol(anTable)] <- round(anTable[, -ncol(anTable)], 1)
    anTable[, ncol(anTable)] <- formatPvalue(anTable[, ncol(anTable)])
    colnames(anTable)[4] <- c("p-value")
    return(anTable)
}

getTukeyGroups <- function(response, factor) {
    aovData <- aov(response ~ factor)
    hsd <- HSD.test(aovData, "factor")
    realGroupNames <- hsd$means %>% rownames() %>% as.character %>% sort()
    changedGroupNames <- hsd$groups$trt %>% as.character %>% sort()
    namesDict <- data.frame(real = realGroupNames,
                            changed = changedGroupNames,
                            stringsAsFactors = FALSE)
    hsd$groups$trt <- as.character(hsd$groups$trt)
    namesDict <- hsd$groups %>% left_join(namesDict, by = c("trt" = "changed"))

    groups <- data.frame(groups = hsd$groups$M)
    groupList <- apply(groups, 1, function(x) substring(x, 1:nchar(x), 1:nchar(x)))
    names(groupList) <- namesDict$real
    maxChar <- max(sapply(groupList, max))
    maxCharNum <- which(letters == maxChar)

    lettersUsed <- data.frame(letters[1:maxCharNum], stringsAsFactors = FALSE)
    lettersUsage <- apply(lettersUsed, 1, function(x) sapply(groupList, function(y) x %in% y)) %>%
        data.frame(stringsAsFactors = FALSE)

    colnames(lettersUsage) <- lettersUsed[,1]
    return(lettersUsage)
}
