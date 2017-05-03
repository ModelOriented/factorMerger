#' Plot Factor Merger
#'
#' @param panel \code{c("all", "response", "GIC", "merging")}
#' @param show \code{c("loglikelihood", "p-value")}
#' @param nodesSpacing \code{c("equidistant", "effects", "modelSpecific")} # TODO: Change names
#' @param summary
#' @param color \code{c("none", "cluster", "summary")}
#' @param showDiagnostics Boolean
#' @param alpha significance interval
#' @param penalty GIC penalty
#' @param palette color palette
#'
#' @importFrom gridExtra grid.arrange
#'
#' @export
plot.factorMerger <- function(factorMerger, panel = "all",
                              show = "loglikelihood",
                              nodesSpacing = "equidistant",
                              summary = NULL, color = "none",
                              clusterSplit = list(stat = "GIC", value = 2),
                              markBestModel = TRUE, markStars = TRUE,
                              alpha = 0.05, penalty = 2,
                              mergingPalette = NULL,
                              responsePalette = NULL,
                              GICcolor = NULL) {

    stopifnot(panel %in% c("all", "response", "GIC", "merging"))
    stopifnot(show %in% c("loglikelihood", "p-value"))
    stopifnot(nodesSpacing %in% c("equidistant", "effects", "modelSpecific"))
    stopifnot(color %in% c("none", "cluster", "summary"))

    summary <- checkSummary(factorMerger, summary)
    responsePlot <- plotResponse(factorMerger, summary, color, responsePalette)
    if(!is.null(responsePalette)) {
        responsePlot <- ggpubr::ggpar(responsePlot, palette = responsePalette)
    }

    if (panel %in% c("all", "GIC") & show == "p-value") {
        warning("GIC plot is not supported with p-value yet.")
    }

    # Set colors
    if (color == "cluster") {
        colorsDf <- getOptimalPartitionDf(factorMerger,
                                          clusterSplit[[1]],
                                          clusterSplit[[2]])
    } else if (color == "summary") {
        colorsDf <- getColorsFromResponseDf(responsePlot)
    } else {
        colorsDf <- NULL
    }

    mergingPathPlot <- plotTree(factorMerger, show, nodesSpacing,
                                clusterSplit, markBestModel, markStars,
                                alpha, color, colorsDf, mergingPalette)

    if (!is.null(mergingPalette)) {
        mergingPathPlot <- ggpubr::ggpar(mergingPathPlot, palette = mergingPalette)
    }
    switch(panel,
           "merging" = {
               return(mergingPathPlot)
           },
           "all" = {
               return(grid.arrange(mergingPathPlot, responsePlot,
                                   plotGIC(factorMerger), ncol = 2,
                                   widths = c(7.5, 2.5), heights = c(7.5, 2.5)))
           },
           "response" = {
               return(grid.arrange(mergingPathPlot, responsePlot,
                                   ncol = 2,
                                   widths = c(7.5, 2.5)))
           },
           "GIC" = {
               return(grid.arrange(mergingPathPlot, plotGIC(factorMerger),
                                   ncol = 1,
                                   heights = c(7.5, 2.5)))
           })
}

# -------------------
# Merging path plot

#' @importFrom ggplot2 theme_classic theme element_line element_blank theme_minimal element_text
treeTheme <- function() {
    myTheme <- theme_minimal() +
        theme(panel.grid.major.x = element_line(color = "lightgrey",
                                                linetype = 2),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 12),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none")
    return(myTheme)
}

plotTree <- function(factorMerger, show, nodesSpacing,
                     clusterSplit, markBestModel, markStars,
                     alpha, color, colorsDf, palette) {
    stopifnot(show %in% c("loglikelihood", "p-value"))
    if (nodesSpacing == "equidistant") {
        return(
            plotSimpleTree(factorMerger, show, clusterSplit,
                           markBestModel, markStars,
                           alpha, color, colorsDf, palette)
        )
    }
    return(
        plotCustomizedTree(factorMerger, show, clusterSplit,
                           nodesSpacing, markBestModel, markStars,
                           alpha, color, colorsDf, palette, groupsStats(factorMerger))
    )
}

plotSimpleTree <- function(factorMerger, show, clusterSplit,
                           markBestModel, markStars,
                           alpha, color, colorsDf, palette) {
    nodesPosition <- getFinalOrder(factorMerger) %>% data.frame()
    mH <- mergingHistory(factorMerger)
    noStep <- nrow(mH)

    for (step in 1:noStep) {
        newLine <- mean(nodesPosition[rownames(nodesPosition) %in% mH[step, ],])
        nodesPosition <- rbind(nodesPosition, newLine)
        rownames(nodesPosition)[nrow(nodesPosition)] <-
            paste(mH[step, ], collapse = "")
    }
    return(plotCustomizedTree(factorMerger, show, clusterSplit,
                              "equidistant", markBestModel, markStars,
                              alpha, color, colorsDf, palette,
                              nodesPosition))

}

#' @importFrom dplyr mutate filter group_by count
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot geom_segment scale_x_log10 theme_bw scale_x_continuous
#' @importFrom ggplot2 coord_flip xlab ylab theme element_blank geom_vline geom_text
#' @importFrom ggplot2 geom_point aes geom_label scale_fill_manual scale_y_continuous labs
plotCustomizedTree <- function(factorMerger, show, clusterSplit,
                               nodesSpacing, markBestModel, markStars,
                               alpha, color, colorsDf, palette,
                               nodesPosition = NULL) {
    statisticColname <- getStatNameInTable(show)
    if (nodesSpacing == "modelSpecific") {
        nodesPosition <- applyModelTransformation(factorMerger, nodesPosition)
    }
    segment <- getTreeSegmentDf(factorMerger, statisticColname, nodesPosition)
    df <- segment$df
    pointsDf <- segment$pointsDf
    labelsDf <- segment$labelsDf
    labelsDf <- labelsDf %>% arrange(-y1)
    showY <- nodesSpacing != "equidistant"

    g <- df %>% ggplot() +
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.75) +
        geom_point(data = pointsDf, aes(x = x1, y = y1), size = 0.75)

    g <- g + scale_y_continuous(limits = getLimits(labelsDf, showY),
                                position = "right",
                                breaks = labelsDf$y1,
                                labels = getLabels(labelsDf, factorMerger)) +
        ylab(getStatisticName(factorMerger)) + xlab(show) +
        labs(title = "Merging path plot",
             subtitle = paste0("Optimal GIC partition: ",
                               paste(getOptimalPartition(factorMerger,
                                                         clusterSplit[[1]],
                                                         clusterSplit[[2]]),
                                     collapse = ":"))) + treeTheme()

    if (color == "cluster") {
        g <- colorCluster(g, segment, factorMerger, clusterSplit, show, palette)
    }

    g <- scaleAxis(g, show, alpha)

    if (markBestModel) {
        mark <- markOptimalModel(factorMerger, clusterSplit, show)
        g <- g + geom_vline(xintercept = mark$intercept, col = "mediumorchid3", linetype = "dotted") +
            geom_label(x = mark$labelIntercept, y = getLimits(labelsDf, showY)[1],
                       label = mark$label, alpha = 0.5, col = "mediumorchid3",
                       angle = 90, size = 3, fontface = "italic")
    }

    if (markStars) {
        g <- g + geom_text(data = pointsDf, aes(x = x1, y = y1,
                                                label = factor(significance)),
                           hjust = 1, vjust = 0.25, size = 5)
    }

    return(g)
}

colorCluster <- function(plot, segment, factorMerger, clusterSplit, show, palette) {
    segmentColoured <- getClustersColors(segment, factorMerger, clusterSplit, show)
    plot <- plot +
        geom_segment(data = segmentColoured$df,
                     aes(x = x1, y = y1, xend = x2, yend = y2, col = pred), size = 0.75) +
        geom_point(data = segmentColoured$pointsDf,
                   aes(x = x1, y = y1, col = pred), size = 0.75)
    nGroups <- length(unique(segmentColoured$labelsDf$pred))
    clusterColors <- factor(segmentColoured$labelsDf$pred, labels = hue_pal()(nGroups)) %>%
        as.character()
    plot <- plot + theme(axis.text.y = element_text(
        color = clusterColors[length(clusterColors):1],
        size = 12))
    return(plot)
}

scaleAxis <- function(plot, show, alpha) {
    if (show == "p-value") {
        plot <- plot + scale_x_log10()
    }

    if (show == "loglikelihood") {
        labBr <- getChisqBreaks(plot$data, alpha)
        plot <- plot +
            scale_x_continuous(breaks = labBr$breaks, labels = labBr$labels)
    }
    return(plot)
}

markOptimalModel <- function(factorMerger, clusterSplit, show) {
    optimalNSteps <- optimalNumberOfMerges(factorMerger,
                                           clusterSplit[[1]],
                                           clusterSplit[[2]])
    mH <- mergingHistory(factorMerger, T, F)
    intercept <- mH[optimalNSteps + 1, getStatNameInTable(show)]
    if (show == "p-value") {
        label <- paste0("alpha = ", intercept)
        labelIntercept <- log10(intercept)
    }
    if (show == "loglikelihood") {
        label <- paste0("loglikelihood = ", intercept)
        labelIntercept <- intercept
    }
    if (clusterSplit[[1]] == "GIC") {
        label <- paste0("min GIC")
    }
    return(list(
        label = label,
        intercept = intercept,
        labelIntercept = labelIntercept
    ))
}

applyModelTransformation <- function(object, nodesPosition) {
    UseMethod("applyModelTransformation", object)
}

applyModelTransformation.default <- function(object, nodesPosition) {
    warning("Model specific nodes spacing is not supported yet.")
    return(nodesPosition)
}

# -------------------
# Summary/response - top right plot

checkSummary <- function(object, summary) {
    UseMethod("checkSummary", object)
}

checkSummary.gaussianFactorMerger <- function(factorMerger, summary) {
    if (NCOL(factorMerger$response) > 1) {
        summarySet <-  c("heatmap", "profile")
    } else {
        summarySet <- c("means", "boxplot")
    }
    warnIfUnexpectedSummary(summarySet, summary)
}

checkSummary.binomialFactorMerger <- function(factorMerger, summary) {
    summarySet <- c("proportion")
    warnIfUnexpectedSummary(summarySet, summary)
}

checkSummary.survivalFactorMerger <- function(factorMerger, summary) {
    summarySet <- c("survival")
    warnIfUnexpectedSummary(summarySet, summary)
}

warnIfUnexpectedSummary <- function(summarySet, summary) {
    if (is.null(summary)) {
        return(summarySet[1])
    }
    if (!(summary %in% summarySet)) {
        warning(paste0("Summary '", summary, ", is not supported by supplied model family -- ", summarySet[1], " used insted."))
        return(summarySet[1])
    }
    return(summary)
}

plotResponse <- function(factorMerger, summary, color, palette) {
    switch(summary,
           "heatmap" = {
               return(plotHeatmap(factorMerger, color))
           },
           "profile" = {
               return(plotProfile(factorMerger, color))
           },
           "boxplot" = {
               return(plotBoxplot(factorMerger, color))
           },
           "means" = {
               return(plotMeansAndStds(factorMerger, color))
           },
           "survival" = {
               return(plotSurvival(factorMerger, color))
           },
           "proportion" = {
               return(plotProportion(factorMerger, color))
           })
}

#' @importFrom proxy dist
#' @importFrom dplyr arrange
findSimilarities <- function(factorMerger) {
    stats <- calculateMeansAndRanks(factorMerger$response,
                                    factorMerger$factor)
    varsToBePloted <- reshape(stats %>% subset(select = -mean),
                              idvar = "level",
                              timevar = "variable",
                              direction = "wide")
    distances <- dist(varsToBePloted[, -1] %>% t(), method = "manhattan")
    iso <- MASS::isoMDS(distances, k = 1, trace = FALSE)$points[, 1]
    iso <- data.frame(var = stats$variable %>% unique(), proj = iso) %>%
        arrange(proj)
    stats$variable <- factor(stats$variable,
                             levels = iso$var)
    return(stats)

}

#' Heatmap (multi-dimensional gaussian)
#'
#' @description Plots heatmap for each dimension of the response variable. Vector of means of factor levels for a given
#' dimension is scaled to have mean equal to zero and standard deviation equal to one.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_minimal scale_x_continuous theme
#' @importFrom ggplot2 coord_flip element_line element_blank scale_fill_distiller labs guides
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange
plotHeatmap <- function(factorMerger, color) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- findSimilarities(factorMerger)
    tmp <- tapply(df$mean, df$variable, function(x) scale(x) %>% as.numeric())
    df$mean <- (do.call(cbind, tmp) %>% reshape2::melt())$value
    df %>% ggplot() +
        geom_tile(aes(x = level, y = variable, fill = mean)) +
        coord_flip() + theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none") +
        scale_fill_distiller(palette = "Greys") +
        labs(title = "Heatmap", subtitle = "Group means by variables")
}

#' Profile plot (multi-dimensional gaussian)
#'
#' @description Plots rank plot - one series is a single factor level and one group
#' on the OX axis is a single dimension of the response.
#'
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_text theme_minimal theme scale_color_manual labs
#' @export
plotProfile <- function(factorMerger, color) {
    df <- findSimilarities(factorMerger)

    df$group <- factor(df$level,
                       levels = (df %>%
                                     filter(variable == levels(df$variable) %>%
                                                head(1)) %>% arrange(rank))$level
    )

    noLevels <- length(levels(df$level))
    df$rank <- factor(df$rank, levels = 1:noLevels)

    switch(color,
           "cluster" = {
               df <- df %>% left_join(getOptimalPartitionDf(factorMerger),
                                      by = c("level" = "orig"))
               g <- df %>% ggplot(aes(x = variable, y = rank,
                                      colour = pred, group = group, label = group))
           },
           "summary" = {
               g <- df %>% ggplot(aes(x = variable, y = rank,
                                      col = group, group = group, label = group))
           },
           "none" = {
               g <- df %>% ggplot(aes(x = variable, y = rank,
                                      col = group, group = group, label = group)) +
                   scale_color_manual(values = colorRamps::magenta2green(noLevels))
           })

    g <- g + geom_line(size = 1) +
        geom_text(data = subset(df,
                                variable == levels(df$variable) %>% tail(1)),
                  aes(x = variable),
                  size = 5.5, hjust = 0.8,  nudge_x = 0.5) +
        ylab("") + labs(title = "Profile plot", subtitle = "Variable means ranks") +
        theme_minimal() + theme(legend.position = "none",
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12))
    return(g)
}


