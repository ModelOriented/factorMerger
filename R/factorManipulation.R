getIncreasingFactor <- function(factorMerger) {
    stats <- calculateGroupStatistic(factorMerger, factorMerger$factor)
    colnames(stats)[2] <- "stat"
    stats <- stats %>% arrange_("stat")
    return(factor(factorMerger$factor, levels = as.character(stats[, 1])))
}

setIncreasingOrder <- function(response, factor, family = "gaussian") {
    if (NROW(response) != NROW(factor)) {
        stop("Input data sizes do not match.")
    }

    if (NCOL(response) > 1) {
        tmpResponse <- isoMDS(dist(response), k = 1, trace = FALSE)$points
        return(setIncreasingOrder(tmpResponse, factor))
    }

    factor <- as.factor(as.character(factor))
    data <- data.frame(y = as.numeric(response), c = factor)
    newOrder <- aggregate(y ~ c, mean, data = data) %>%
        dplyr::arrange(y)
    factor(factor, levels = as.character(newOrder[, 1]))
}

mergeLevels <- function(factor, groupA, groupB, groupAB = NULL) {
    if (is.null(groupAB)) {
        groupAB <- paste0(groupA, groupB)
    }
    whichLevels <- which(levels(factor) %in% c(groupA, groupB))
    levels(factor)[whichLevels] <- groupAB
    factor
}

calculateGroupStatistic <- function(factorMerger, factor) {
    UseMethod("calculateGroupStatistic", factorMerger)
}

calculateGroupStatistic.default <- function(factorMerger, factor) {
    if (NCOL(factorMerger$response) == 1) {
        return(calculateMeans(factorMerger$response, factor))
    }
    else {
        if (!"projectedResponse" %in% names(factorMerger)) {
            factorMerger$projectedResponse <-
                MASS::isoMDS(dist(factorMerger$response),
                             k = 1, trace = FALSE)$points[, 1]
        }
        return(calculateMeans(factorMerger$projectedResponse, factor))
    }
}

#' @importFrom dplyr arrange left_join
calculateGroupStatistic.survivalFactorMerger <- function(factorMerger, factor) {
    model <- calculateModel(factorMerger, factor)
    coefs <- data.frame(level = factor,
                        coef = predict(model, type = "risk")) %>% unique()

    coefs <- coefs %>% arrange(coef)
    return(coefs)
}

#' @importFrom dplyr group_by summarize arrange
calculateMeans <- function(response, factor) {
    if (is.null(response)) {
        return(NA)
    }
    df <- data.frame(response, level = factor)
    df <- aggregate(. ~ level, function(x) mean(x, na.rm = TRUE), data = df)

    if (NCOL(response) == 1) {
        colnames(df)[2] <- "mean"
        df <- df %>% arrange(mean)
    }
    return(df)
}

#' @importFrom reshape2 melt
#' @importFrom dplyr rename
calculateMeansAndRanks <- function(response, factor) {
    means <- apply(as.data.frame(response), 2, function(x) {
        aggregate(x ~ level, function(x) mean(x, na.rm = T),
                  data = data.frame(x = x, level = factor))
    })

    means <- lapply(means, function(x) {
        df <- x %>% arrange(-x)
        df$rank <- ave(df$x, FUN = rank)
        df$rank <- nrow(df) - df$rank + 1
        df
    })

    attr(means, "varname") <- "variable"

    return(melt(means, id.vars = c("level", "rank"),
                value.name = "mean"))
}
