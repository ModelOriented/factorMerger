# Calculates points positions for the merging path plot
getTreeSegmentDf <- function(factorMerger, statisticColname, nodesPosition) {
    factor <- factorMerger$factor
    noGroups <- length(levels(factor))
    df <- nodesPosition[1:noGroups, ] %>% data.frame
    colnames(df) <- "y1"
    df$y2 <- df$y1
    df$label <- rownames(nodesPosition)[1:noGroups]
    df$x1 <- factorMerger$mergingList$`1`$modelStats[, statisticColname]
    df$x2 <- NA
    df$step <- 0
    labelsDf <- df
    pointsDf <- df
    pointsDf$significance <- ""
    merging <- mergingHistory(factorMerger)

    for (step in 1:nrow(merging)) {
        stepStats <- factorMerger$mergingList[[step + 1]]$modelStats
        statVal <- stepStats[, statisticColname]
        pair <- merging[step, ]
        whichDf <- which(df$label %in% pair)
        df[whichDf, "x2"] <- statVal
        pairLabel <- paste(pair, collapse = "")
        pairMean <- nodesPosition[rownames(nodesPosition) == pairLabel,]
        crosswise <- c(df$y1[whichDf], "", statVal, statVal, step)
        newVertex <- c(pairMean, pairMean, pairLabel, statVal, NA, step)
        df <- rbind(df, crosswise)
        df <- rbind(df, newVertex)
        pointsDf <- rbind(pointsDf, c(newVertex, ""))
        pointsDf[nrow(pointsDf), "significance"] <-
            getSignificanceStar(stepStats[, "pvalVsPrevious"])
    }

    lab <- df$label
    df <- df %>% subset(select = -label) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    df$label <- lab
    df <- df[complete.cases(df), ]

    pointsDf <- subset(pointsDf, select = c(x1, y1, significance, step))
    pointsDf[, -3] <- pointsDf[, -3] %>% apply(2, as.numeric) %>% as.data.frame()

    return(list(df = df,
                labelsDf = labelsDf,
                pointsDf = pointsDf))
}

getSignificanceStar <- function(pval) {
    if (pval > 0.1) {
        return("")
    }
    if (pval > 0.05) {
        return(".")
    }
    if (pval > 0.01) {
        return("*")
    }
    if (pval > 0.001) {
        return("**")
    }
    return("***")
}

optimalNumberOfMerges <- function(factorMerger, stat = "GIC", value = 2) {
    mH <- mergingHistory(factorMerger, T, F)
    if (stat == "GIC") {
        mH$GIC <- -2 * mH$model + as.numeric(value) * nrow(mH):1
        nSteps <- which.min(mH$GIC) - 1
    } else {
        statColname <- getStatNameInTable(stat)
        nSteps <- 1
        while (nSteps < nrow(mH) & mH[nSteps + 1, statColname] >= value) {
            nSteps <- nSteps + 1
        }
    }
    return(nSteps)
}

# Divides part of a segment plot into clusters
#' @importFrom dplyr arrange summarise group_by left_join
getClustersColors <- function(segment, factorMerger, clusterSplit, stat) {
    nSteps <- optimalNumberOfMerges(factorMerger,
                                    clusterSplit[[1]],
                                    clusterSplit[[2]])
    mH <- mergingHistory(factorMerger, T, F)
    bestModel <- mH[nSteps + 1, getStatNameInTable(stat)]
    segment <- lapply(segment, function(x)
        x %>% filter(step <= nSteps))
    if (length(segment$df[segment$df$x2 < bestModel, ]$x2) > 0) {
        segment$df[segment$df$x2 < bestModel, ]$x2 <- bestModel
    }
    if (length(segment$df[segment$df$x2 < bestModel, ]$x2) > 0) {
        segment$df[segment$df$x1 < bestModel, ]$x1 <- bestModel
    }
    # segment <- lapply(segment, function(x)x %>% filter(x1 >= bestModel))
    map <- getOptimalPartitionDf(factorMerger, clusterSplit[[1]], clusterSplit[[2]])
    map$pred <- as.character(map$pred)
    map$orig <- as.character(map$orig)
    segment$df <- segment$df %>% left_join(map, by = c("label" = "orig"))
    segment$labelsDf <- segment$labelsDf %>% left_join(map, by = c("label" = "orig"))
    segment$pointsDf$pred <- NA
    segment$pointsDf[1:nrow(segment$labelsDf), "pred"] <- segment$labelsDf$pred

    if (nSteps >= 1) {
        for (i in 1:nSteps) {
            pair <- mH[i + 1, 1:2]
            prediction <- (map$pred[map$orig %in% pair] %>% unique())
            segment$df[segment$df$step == i, "pred"] <- prediction
            segment$pointsDf[segment$pointsDf$step == i, "pred"] <- prediction
            map$orig[map$orig %in% pair] <- paste0(pair[1], pair[2])
        }
    }

    segment$labelsDf <- segment$labelsDf %>% arrange(y1)
    segment$labelsDf$pred <- factor(segment$labelsDf$pred,
                                    levels = segment$labelsDf$pred %>% unique())

    segment$pointsDf$pred <- factor(segment$pointsDf$pred,
                                    levels = segment$labelsDf$pred %>% unique())

    segment$df$pred <- factor(segment$df$pred,
                              levels = segment$labelsDf$pred %>% unique())

    return(segment)
}

getStatNameInTable <- function(stat) {
    switch(stat,
           "loglikelihood" = { return("model") },
           "p-value" = { return("pvalVsFull") },
           "GIC" = { return("GIC")} )
}

getLimits <- function(df, simplifiedNodes) {
    if (!simplifiedNodes) {
        return(c(min(df$y1) - 0.00001, max(df$y1) + 0.00001))
    } else {
        shift <- 0.6 / (nrow(df) - 1)
        return(c(min(df$y1) - shift, max(df$y1) + shift))
    }
}

getStatisticName <- function(factorMerger) {
    UseMethod("getStatisticName", factorMerger)
}

getStatisticName.gaussianFactorMerger <- function(factorMerger) {
    if (NCOL(factorMerger$response) > 1) {
        return ("isoMDS projection group means")
    }
    return("Group means")
}

getStatisticName.binomialFactorMerger <- function(factorMerger) {
    return("Group proportion of success")
}

getStatisticName.survivalFactorMerger <- function(factorMerger) {
    return("Initial survival model coefficient")
}
#' @importFrom dplyr left_join
getLabels <- function(labelsDf, factorMerger) {
    stats <- groupsStats(factorMerger)
    if (is.null(stats)) return(NULL)
    stats$label <- rownames(stats)
    labelsDf <- labelsDf %>% left_join(stats, by = "label")
    return(paste0(labelsDf$label, ": ", round(labelsDf[, ncol(labelsDf)], 2)))
}

nLabels <- 6

getChisqBreaks <- function(plotData, alpha) {
    right <- plotData$x1 %>% max()
    left <- plotData$x2 %>% min()
    breaks <- seq(left, right, qchisq(1 - alpha, 1))
    step <- ceiling(length(breaks) / nLabels)
    subs <- (1:nLabels) * step
    subs <- subs[subs <= length(breaks)]
    labels <- breaks %>% round()

    labels[setdiff(1:length(breaks), subs) %>% round()] <- ""
    return(
        list(
            breaks = breaks,
            labels = labels
        )
    )
}

calculateBoxPlotMoments <- function(df) {
    return(
        df %>% left_join(
            df %>% group_by(group) %>% summarize(y0 = min(y, na.rm = TRUE),
                                             y25 = quantile(y, 0.25, na.rm = TRUE),
                                             y50 = mean(y, na.rm = TRUE),
                                             y75 = quantile(y, 0.75, na.rm = TRUE),
                                             y100 = max(y, na.rm = TRUE)), by = "group")
    )

}

getMeansAndStds <- function(factorMerger, factor) {
    model <- lm(factorMerger$response ~ factor - 1)
    df <- data.frame(group = levels(factor))
    sumModel <- summary(model)
    df$mean <- sumModel$coefficients[, 1]
    df$left <- df$mean - 1.96 * sumModel$coefficients[, 2]
    df$right <- df$mean + 1.96 * sumModel$coefficients[, 2]
    df$group <- factor(df$group, levels = df$group)
    return(df)
}
