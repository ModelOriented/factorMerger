#' @importFrom ggplot2 ggplot aes aes_string coord_flip xlab ylab labs stat_summary guides position_dodge annotate unit ggtitle
#' @importFrom ggplot2 geom_point geom_step geom_bar geom_rect geom_ribbon geom_line geom_boxplot geom_errorbar geom_label geom_segment geom_vline geom_text geom_tile
#' @importFrom ggplot2 scale_y_discrete scale_y_continuous scale_x_discrete scale_color_manual scale_fill_brewer scale_fill_manual scale_fill_distiller scale_x_continuous scale_x_log10
#' @importFrom ggplot2 theme_classic theme theme_bw theme_minimal element_line element_blank element_text
#' @importFrom ggplot2 scale_color_brewer ggplot_build
NULL

#' @importFrom dplyr mutate filter group_by count arrange summarize left_join summarise
#' @importFrom dplyr group_by_ arrange_
NULL

#' @importFrom survival survexp
NULL

#' @importFrom reshape2 melt
NULL

#' @importFrom utils packageVersion
NULL

#' @importFrom ggpubr ggarrange
NULL

#' Plot Factor Merger
#'
#' @param x object of a class \code{factorMerger}.
#' @param panel Type of panels to be plot. Possible values are \code{c("all", "response", "GIC", "tree")}.
#' All types of plots include the Factor Merger Tree. Apart from the Factor Merger Tree there are
#' also two possible panels: the Response Plot (response summary, specific to the model family),
#' the GIC Plot (GIC vs. loglikelihood/p-value).
#' \itemize{
#' \item \code{"all"} plots all panels and a short summary of the full model,
#' \item \code{"response"} plots the Factor Merger Tree and the Response Plot,
#' \item \code{"GIC"} plots the Factor Merger Tree and the GIC Plot,
#' \item \code{"tree"} plots the Factor Merger Tree only.
#' }
#' @param statistic Statistic to be displayed on the OX axis of the Factor Merger Tree.
#' Possible values are \code{c("loglikelihood", "p-value")}.
#' If \code{"p-value"} is chosen p-value for the Likelihood Ratio Test against the full model is plot on the OX axis.
#' @param nodesSpacing Type of vertical nodes spacing in the Factor Merger Tree). May be chosen from
#'  \code{c("equidistant", "effects", "modelSpecific")}. \code{"effects"} arranges nodes according to
#'  the model coefficients estimatiors (e.g. in Gaussian case on the OY axis group means are plotted).
#'
#' # TODO: Implement "modelSpecific".
#'
#' @param responsePanel Response panel type -- accepts the following values dependent on the model family:
#' \itemize{
#' \item multi dimensional Gaussian: \code{c("heatmap", "profile")},
#' \item single dimensional Gaussian: \code{c("means", "boxplot","tukey")},
#' \item binomial: \code{c("proportion")},
#' \item survival: \code{c("survival")}
#' }
#'
#' @param colorClusters Boolean. If \code{TRUE}, the default, the Factor Merger Tree is colored according
#' to the optimal factor split (defined by \code{splitStatistic} and \code{splitThreshold} or
#' \code{splitStatistic} and \code{penalty}).
#' @param splitStatistic Statistic used in the optimal split definition. Possible values are:
#' \code{c("GIC", "loglikelihood", "p-value")}. If \code{"GIC"} is chosen, factor is split to minimize GIC
#' with the penalty \code{penalty}. Otherwise, chooses the very last partition whose corresponding statistic
#' (model loglikelihood or p-value for the LRT test) is not lower than \code{splitThreshold}.
#' @param splitThreshold Threshold used in the optimal split definition. Used only with
#' \code{splitStatistic = c("loglikelihood", "p-value")}.
#' @param penalty GIC penalty used for defining the optimal partition with \code{splitStatistic = "GIC"}.
#' The same penalty is used in the GIC plot.
#' @param showSplit Boolean. If \code{TRUE} plots vertical line crossing the optimal split.
#' @param showSignificance Boolean. If \code{TRUE}, the default, marks partitions
#' that are significantly worse than their predecessors on the Factor Merger Tree (uses the Likelihood Ratio Test).
#'
#' Significance codes are:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1.
#' @param chisqQuantile Significance level used if \code{panelGrid = TRUE}.
#'
#' @param palette Color palette used in the Factor Merger Tree and the Response Plot.
#' @param responsePanelPalette Additional color palette used in the Response Plot if
#' palettes for the Factor Merger Tree and the Response Plot are to be different.
#' @param gicPanelColor Color used in the GIC plot.
#' @param title Factor Merger Tree plot's title.
#' @param subtitle Factor Merger Tree plot's subtitle.
#' @param panelGrid Boolean. If \code{TRUE}, each interval on
#'  the OX axis of the Factor Merger Tree corresponds to the 1 - \code{chisqQuantile}
#'  quantile of chi-square distribution. Otherwise, panel is blank.
#' @param ... Other arguments
#'
#' @export
plot.factorMerger <- function(x, panel = "all",
                              statistic = "loglikelihood",
                              nodesSpacing = "equidistant",
                              colorClusters = TRUE,
                              splitStatistic = "GIC",
                              splitThreshold = NULL,
                              penalty = 2,
                              showSplit = FALSE,
                              showSignificance = TRUE,
                              title = "Factor Merger Tree",
                              subtitle = " ",
                              palette = NULL,
                              responsePanel = NULL,
                              responsePanelPalette = NULL,
                              gicPanelColor = NULL,
                              panelGrid = TRUE,
                              chisqQuantile = 0.05, ...) {

    stopifnot(panel %in% c("all", "response", "GIC", "tree"))
    stopifnot(statistic %in% c("loglikelihood", "p-value"))
    stopifnot(nodesSpacing %in% c("equidistant", "effects", "modelSpecific"))
    stopifnot(splitStatistic == "GIC" | !is.null(splitThreshold))

    clusterSplit <- getClusterSplit(splitStatistic, splitThreshold, penalty)

    responsePanel <- checkResponsePanel(x, responsePanel)

    if (is.null(responsePanelPalette) && !is.null(palette)) {
        responsePanelPalette <- palette
    }

    responsePlot <- plotResponse(x, responsePanel,
                                 colorClusters, clusterSplit,
                                 responsePanelPalette)

    # Set colors
    if (colorClusters) {
        colorsDf <- getOptimalPartitionDf(x,
                                          clusterSplit[[1]],
                                          clusterSplit[[2]])
    } else {
        colorsDf <- data.frame(
            orig = levels(x$factor),
            pred = levels(x$factor)
        )
    }

    mergingPathPlot <- plotTree(x, statistic, nodesSpacing,
                                clusterSplit, showSplit, showSignificance,
                                chisqQuantile, colorClusters, colorsDf, palette,
                                title, subtitle, panelGrid)

    switch(panel,
           "tree" = {
               return(mergingPathPlot)
           },
           "all" = {
               aovTable <- plotTable(calculateAnovaTable(x$initialModel))
               return(ggarrange(mergingPathPlot, responsePlot,
                                   plotGIC(x, gicPanelColor,
                                           penalty, statistic),
                                   aovTable,
                                ncol = 2, nrow = 2,  align = "hv",
                                widths = c(2, 1), heights = c(2, 1)))
           },
           "response" = {
               return(ggarrange(mergingPathPlot, responsePlot,
                                ncol = 2,  align = "hv",
                                widths = c(2, 1)))
           },
           "GIC" = {
               return(ggarrange(mergingPathPlot,
                                   plotGIC(x, gicPanelColor,
                                           penalty, statistic),
                                   nrow = 2, align = "hv",
                                   heights = c(2, 1)))
           })
}

getClusterSplit <- function(splitStatistic, splitThreshold, penalty) {
    stopifnot(splitStatistic %in% c("GIC", "loglikelihood", "p-value"))
    if (splitStatistic == "GIC") {
        return(
            list(stat = "GIC",
                 value = penalty)
        )
    }
    return(
        list(stat = splitStatistic,
             value = splitThreshold)
    )
}

# -------------------
# Factor Merger Tree

treeTheme <- function(panelGrid = TRUE) {
    myTheme <- theme_minimal() +
        theme(panel.grid.major.x = element_line(color = "lightgrey",
                                                linetype = 2),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.text.x = element_text(size = 12),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none")
    if (!panelGrid) {
        myTheme <- myTheme + theme(panel.grid.major.x = element_blank())
    }
    return(myTheme)
}

plotTree <- function(factorMerger, statistic, nodesSpacing,
                     clusterSplit, markBestModel, markStars,
                     alpha, color, colorsDf,
                     palette = NULL, title, subtitle, panelGrid) {
    stopifnot(statistic %in% c("loglikelihood", "p-value"))
    if (nodesSpacing == "equidistant") {
        return(
            plotSimpleTree(factorMerger, statistic, clusterSplit,
                           markBestModel, markStars,
                           alpha, color, colorsDf, palette,
                           title, subtitle, panelGrid)
        )
    }
    return(
        plotCustomizedTree(factorMerger, statistic, clusterSplit,
                           nodesSpacing, markBestModel, markStars,
                           alpha, color, colorsDf, palette,
                           groupsStats(factorMerger),
                           title, subtitle, panelGrid)
    )
}

plotSimpleTree <- function(factorMerger, statistic, clusterSplit,
                           markBestModel, markStars,
                           alpha, color, colorsDf, palette = NULL,
                           title, subtitle, panelGrid) {
    # We want to have reverse order of variables! TODO
    nodesPosition <- as.data.frame(getFinalOrder(factorMerger, TRUE))
    mH <- mergingHistory(factorMerger)
    noStep <- nrow(mH)

    for (step in 1:noStep) {
        which <- rownames(nodesPosition) %in% mH[step, ]
        newLine <- nodesPosition[which, ] %>%
            mean()
        nodesPosition <- rbind(nodesPosition, newLine)
        rownames(nodesPosition)[nrow(nodesPosition)] <-
            paste(mH[step, ], collapse = "")
    }

    return(plotCustomizedTree(factorMerger, statistic, clusterSplit,
                              "equidistant", markBestModel, markStars,
                              alpha, color, colorsDf, palette,
                              nodesPosition, title, subtitle, panelGrid))

}

#' @importFrom scales hue_pal
plotCustomizedTree <- function(factorMerger, statistic, clusterSplit,
                               nodesSpacing, markBestModel, markStars,
                               alpha, color, colorsDf, palette = NULL,
                               nodesPosition = NULL, title,
                               subtitle, panelGrid) {
    statisticColname <- getStatNameInTable(statistic)

    if (nodesSpacing == "modelSpecific") {
        nodesPosition <- applyModelTransformation(factorMerger, nodesPosition)
    }

    segment <- getTreeSegmentDf(factorMerger, statisticColname, nodesPosition)
    df <- segment$df
    pointsDf <- segment$pointsDf
    labelsDf <- segment$labelsDf %>%
        arrange(-y1)
    showY <- nodesSpacing != "equidistant"

    g <- ggplot(df) +
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.5) +
        geom_point(data = pointsDf, aes(x = x1, y = y1), size = 0.75) +
        scale_y_continuous(limits = getLimits(labelsDf, showY),
                           position = "right",
                           breaks = labelsDf$y1,
                           labels = getLabels(labelsDf,
                                              factorMerger),
                           expand = c(0, 0)) +
        ylab(getStatisticName(factorMerger)) + xlab(statistic) +
        labs(title = title,
             subtitle = subtitle) + treeTheme(panelGrid)
    if (color) {
        g <- addClustersColors(g, segment = segment, 
                               factorMerger = factorMerger,
                               clusterSplit = clusterSplit, 
                               statistic = statistic, 
                               palette = palette)
    } else {
        clusterColors <- getClustersColorsNames(palette,
                                                NROW(colorsDf),
                                                factor(labelsDf$label,
                                                       levels = labelsDf$label))
        g <- g + theme(axis.text.y = element_text(
            color = rev(clusterColors),
            size = 15))
    }

    g <- scaleAxis(g, statistic, alpha)

    if (markBestModel) {
        mark <- markOptimalModel(factorMerger, clusterSplit, statistic, alpha)
        g <- g + geom_vline(xintercept = mark$intercept,
                            col = "mediumorchid3", linetype = "dotted") +
            geom_label(x = mark$labelIntercept,
                       y = getLimits(labelsDf, showY)[1] +
                           0.1 * (getLimits(labelsDf, showY)[2] -
                                getLimits(labelsDf, showY)[1]),
                       label = mark$label,
                       alpha = 0.5,
                       col = "mediumorchid3",
                       angle = 90, size = 3,
                       fontface = "italic")
    }

    if (markStars) {
        g <- g + geom_text(data = pointsDf, aes(x = x1, y = y1,
                                                label = factor(significance)),
                           hjust = 1, vjust = 0.25, size = 5)
    }

    if (!is.null(palette)) {
        g <- g + scale_fill_brewer(palette = palette) +
            scale_color_brewer(palette = palette)
    }

    return(g)
}

getClustersColorsNames <- function(palette, nGroups, pred) {
    if (is.null(palette)) {
        return(factor(pred, labels = hue_pal()(nGroups)) %>%
                   as.character())
    }
    # Fake data
    df <- data.frame(x = 1:nGroups, y = 1, col = letters[1:nGroups])

    g <- ggplot(df, aes_string(x = "x", y = "y", color = "col")) +
        geom_point() +
        scale_color_brewer(palette = palette)
    colors <- ggplot_build(g)$data[[1]]$colour
    return(factor(pred, labels = colors) %>%
               as.character())
}

addClustersColors <- function(plot, segment,
                              factorMerger, clusterSplit,
                              statistic, palette = NULL) {
    segmentColoured <- getClustersColors(segment, factorMerger,
                                         clusterSplit, statistic)
    plot <- plot +
        geom_segment(data = segmentColoured$df,
                     aes(x = x1, y = y1,
                         xend = x2, yend = y2,
                         col = pred),
                     size = 0.75) +
        geom_point(data = segmentColoured$pointsDf,
                   aes(x = x1, y = y1, col = pred), size = 0.75)
    nGroups <- length(unique(segmentColoured$labelsDf$pred))
    clusterColors <- getClustersColorsNames(palette,
                                            nGroups,
                                            segmentColoured$labelsDf$pred)

    plot <- plot + theme(axis.text.y = element_text(
        color = rev(clusterColors),
        size = 15))

    return(plot)
}

scaleAxis <- function(plot, statistic, alpha) {
    if (statistic == "p-value") {
        plot <- plot + scale_x_log10()
    }

    if (statistic == "loglikelihood") {
        labBr <- getChisqBreaks(plot$data, alpha)
        plot <- plot +
            scale_x_continuous(breaks = labBr$breaks, labels = labBr$labels)
    }
    return(plot)
}

markOptimalModel <- function(factorMerger, clusterSplit, statistic, alpha) {
    optimalNSteps <- optimalNumberOfMerges(factorMerger,
                                           clusterSplit[[1]],
                                           clusterSplit[[2]])
    mH <- mergingHistory(factorMerger, T, F)
    intercept <- mH[optimalNSteps + 1, getStatNameInTable(statistic)]
    if (statistic == "p-value") {
        label <- paste0("alpha = ", round(intercept, 2))
        labelIntercept <- log10(intercept)
    }
    if (statistic == "loglikelihood") {
        label <- paste0("loglikelihood = ", round(intercept))
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
# Response - top right plot

checkResponsePanel <- function(object, responsePanel) {
    UseMethod("checkResponsePanel", object)
}

checkResponsePanel.gaussianFactorMerger <- function(factorMerger,
                                                    responsePanel) {
    if (NCOL(factorMerger$response) > 1) {
        responsePanelSet <-  c("profile", "heatmap", "frequency")
    } else {
        responsePanelSet <- c("means", "boxplot", "frequency", "tukey")
    }
    warnIfUnexpectedResponsePanel(responsePanelSet, responsePanel)
}

checkResponsePanel.binomialFactorMerger <- function(factorMerger,
                                                    responsePanel) {
    responsePanelSet <- c("proportion", "frequency")
    warnIfUnexpectedResponsePanel(responsePanelSet, responsePanel)
}

checkResponsePanel.survivalFactorMerger <- function(factorMerger,
                                                    responsePanel) {
    responsePanelSet <- c("survival", "frequency")
    warnIfUnexpectedResponsePanel(responsePanelSet, responsePanel)
}

warnIfUnexpectedResponsePanel <- function(responsePanelSet,
                                          responsePanel) {
    if (is.null(responsePanel)) {
        return(responsePanelSet[1])
    }
    if (!(responsePanel %in% responsePanelSet)) {
        warning(
            paste0(
                "ResponsePanel '", responsePanel,
                ", is not supported by supplied model family -- ",
                responsePanelSet[1], " used insted."
                )
            )
        return(responsePanelSet[1])
    }
    return(responsePanel)
}

plotResponse <- function(factorMerger, responsePanel,
                         colorClusters, clusterSplit,
                         responsePanelPalette) {
    switch(responsePanel,
           "heatmap" = {
               return(plotHeatmap(factorMerger,
                                  colorClusters, clusterSplit,
                                  responsePanelPalette))
           },
           "profile" = {
               return(plotProfile(factorMerger,
                                  colorClusters, clusterSplit,
                                  responsePanelPalette))
           },
           "boxplot" = {
               return(plotBoxplot(factorMerger,
                                  colorClusters, clusterSplit,
                                  responsePanelPalette))
           },
           "means" = {
               return(plotMeansAndConfInt(factorMerger,
                                          colorClusters, clusterSplit,
                                          responsePanelPalette))
           },
           "survival" = {
               return(plotSurvival(factorMerger,
                                   colorClusters, clusterSplit,
                                   responsePanelPalette))
           },
           "proportion" = {
               return(plotProportion(factorMerger,
                                     colorClusters, clusterSplit,
                                     responsePanelPalette))
           },
           "frequency" = {
               return(plotFrequency(factorMerger,
                                    FALSE, clusterSplit))
           },
           "tukey" = {
               return(plotTukey(factorMerger,
                                responsePanelPalette))
           })
}

findSimilarities <- function(factorMerger) {
    stats <- calculateMeansAndRanks(factorMerger$response,
                                    factorMerger$factor)
    varsToBePloted <- reshape(stats %>% subset(select = -mean),
                              idvar = "level",
                              timevar = "variable",
                              direction = "wide")
    distances <- dist(varsToBePloted[, -1] %>% t(), method = "manhattan")
    distances[distances == 0] <- 1e-4
    iso <- MASS::isoMDS(distances, k = 1, trace = FALSE)$points[, 1]
    iso <- data.frame(var = stats$variable %>% unique(), proj = iso) %>%
        arrange(proj)
    stats$variable <- factor(stats$variable,
                             levels = iso$var)
    return(stats)

}

#' Heatmap (multi-dimensional Gaussian)
#'
#' @description Plots heatmap for each dimension of the response variable. Vector of means of factor levels for a given
#' dimension is scaled to have mean equal to zero and standard deviation equal to one.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#' @param palette custom palette
#'
#' @export
plotHeatmap <- function(factorMerger, color, clusterSplit, palette = "Greys") {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- findSimilarities(factorMerger)
    tmp <- tapply(df$mean, df$variable, function(x) scale(x) %>% as.numeric())
    df$mean <- (do.call(cbind, tmp) %>% reshape2::melt())$value
    df %>% ggplot() +
        geom_tile(aes_string(x = "level", y = "variable", fill = "mean")) +
        coord_flip() + theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size = 12),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none") +
        xlab("") + scale_y_discrete(expand = c(0, 0)) +
        scale_fill_distiller(palette = palette) +
        labs(title = "Heatmap", subtitle = "Group means by variable")
}

#' Profile plot (multi-dimensional Gaussian)
#'
#' @description Plots rank plot - one series is a single factor level and one group
#' on the OX axis is a single dimension of the response.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#' @param palette custom palette
#'
#' @export
plotProfile <- function(factorMerger, color, clusterSplit, palette = NULL) {
    df <- findSimilarities(factorMerger)

    profileLevels <- df %>%
        filter(variable == levels(df$variable) %>%
                   head(1)) %>%
        arrange(rank) %>%
        subset(select = level)

    df$group <- factor(df$level, levels = profileLevels$level)

    noLevels <- length(levels(df$level))
    df$rank <- factor(df$rank, levels = noLevels:1)

    if (color) {
        df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                     clusterSplit[[1]],
                                                     clusterSplit[[2]]),
                               by = c("level" = "orig"))
        g <- df %>% ggplot(aes(x = variable, y = rank,
                               colour = pred, group = group, label = group))
    } else {
        g <- df %>% ggplot(aes(x = variable, y = rank,
                               col = group, group = group, label = group)) +
            scale_color_manual(values = colorRamps::magenta2green(noLevels))
    }

    g <- g + geom_line(size = 1) +
        geom_text(data =
                      subset(df,
                             variable == levels(df$variable) %>% tail(1)),
                  aes(x = variable),
                  size = 4, hjust = 0.8,  nudge_x = 0.5) +
        ylab("") + xlab("") +
        labs(title = "Profile plot", subtitle = "Variable means ranks") +
        theme_minimal() + theme(legend.position = "none",
              plot.title = element_text(size = 18),
              axis.text = element_text(size = 12),
              plot.subtitle = element_text(size = 12)) +
        scale_y_discrete(expand = c(1 / (2 * nrow(df)), 1 / (2 * nrow(df))))
    if (!is.null(palette) & color) {
        g <- g + scale_color_brewer(palette = palette) +
            scale_fill_brewer(palette = palette)
    }
    return(g)
}

#' Boxplot (single-dimensional Gaussian)
#'
#' @description Plots boxplot with mean as a summary statistic groupping observation by factor levels.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#' @param palette custom palette
#'
#' @export
plotBoxplot <- function(factorMerger, color, clusterSplit, palette = NULL) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- data.frame(group = factorMerger$factor, y = factorMerger$response)
    df <- calculateBoxPlotMoments(df)

    if (color) {
        df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                     clusterSplit[[1]],
                                                     clusterSplit[[2]]),
                               by = c("group" = "orig"))
        g <- df %>% ggplot(aes(y = y, x = group, group = group, fill = pred))
    } else {
        g <- df %>% ggplot(aes(y = y, x = group, group = group))
    }

    g <- g + geom_boxplot(aes(ymin = y0,
                         lower = y25,
                         middle = y50,
                         upper = y75,
                         ymax = y100), stat = "identity") +
        coord_flip() + treeTheme() + xlab("") + ylab("") +
        theme(axis.text.y = element_blank()) +
        labs(title = "Boxplot", subtitle = "Summary statistic: mean")

    if (!is.null(palette)) {
        g <- g +
            scale_color_brewer(palette = palette) +
            scale_fill_brewer(palette = palette)
    }

    return(g)
}

#' Means and standard deviation plot (single-dimensional Gaussian)
#'
#' @description For each factor level plots its mean and interval of the length equal to its standard deviation.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#' @param palette custom palette
#'
#' @export
plotMeansAndConfInt <- function(factorMerger, color,
                                clusterSplit, palette = NULL) {
    factor <- factor(factorMerger$factor,
                     levels = getFinalOrderVec(factorMerger))
    df <- getMeansAndStds(factorMerger, factor)

    if (color) {
        part <- getOptimalPartitionDf(factorMerger,
                                      clusterSplit[[1]],
                                      clusterSplit[[2]])
        part$orig <- factor(part$orig, levels = getFinalOrderVec(factorMerger))
        df <- df %>% left_join(part, by = c("group" = "orig"))
        g <- df %>% ggplot(aes(colour = pred, fill = pred,
                               x = as.factor(group),
                               y = mean, group = as.factor(group)))
    } else {
        g <- df %>% ggplot(aes(x = as.factor(group), colour = as.factor(group),
                               y = mean, group = as.factor(group)))
    }

    g <- g + geom_errorbar(aes(ymin = left, ymax = right),
                      width = .5,
                      position = position_dodge(.5)) + treeTheme() +
        geom_point(size = 2) + coord_flip() +
        theme(axis.title.x = element_text(),
              axis.text.y = element_blank()) +
        labs(title = "Group means",
             subtitle = "with 95% confidence intervals") +
        ylab("")

    if (!is.null(palette)) {
        g <- g +
            scale_color_brewer(palette = palette) +
            scale_fill_brewer(palette = palette)
    }

    return(g)
}

#' Proportion plot (binomial)
#'
#' @description Plots proportion of success for each factor level.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#' @param palette custom palette
#'
#' @export
plotProportion <- function(factorMerger, color, clusterSplit, palette = NULL) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- data.frame(group = factorMerger$factor, y = factorMerger$response)

    if (color) {
        df <- df %>% group_by(group) %>% summarise(mean = mean(y == 1))
        df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                     clusterSplit[[1]],
                                                     clusterSplit[[2]]),
                               by = c("group" = "orig"))
        g <- df %>% ggplot() +
            geom_bar(aes(x = group, y = mean, fill = pred), stat = "identity")
    } else {
        g <- df %>% ggplot() +
            geom_bar(aes(x = group, fill = as.factor(y)),
                     position = "fill")
    }

    g <- g + scale_y_continuous(labels = scales::percent, name = "") +
        scale_x_discrete(expand = c(0, 0)) +
        coord_flip() + treeTheme() +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
        labs(title = "Success ratio")

    if (!is.null(palette)) {
        g <- g +
            scale_fill_brewer(palette = palette)
    }

    return(g)
}

#' Survival plot (survival)
#'
#' @description Plots adjusted survival curves for coxph model for each group. Survival probabilities are calculated from \code{coxph} model.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#' @param palette custom palette
#'
#' @export
plotSurvival <- function(factorMerger, color, clusterSplit, palette = NULL) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor,
                                  levels = levels)
    df <- data.frame(response = factorMerger$response, factor = factorMerger$factor)

    model <- calculateModel(factorMerger, factorMerger$factor)
    pred <- survexp(~factor, data = df, ratetable = model)

    curve <- data.frame(time = rep(c(0,pred$time), length(levels)),
                        variable = factor(rep(levels, each=1+length(pred$time))),
                        surv = c(rbind(1, pred$surv)))

    if (color) {
      factorNew <- time <- surv <- NULL

      df$factorNew <- cutTree(factorMerger,
                              clusterSplit[[1]],
                              clusterSplit[[2]])
      keys <- unique(df[,c("factor", "factorNew")])
      both <- merge(curve, keys, by.x = "variable", by.y = "factor")
      curve <- summarise(
          group_by(
            both,
            variable = factorNew, time),
        surv = mean(surv, na.rm=TRUE))
    }

    g <- ggplot(curve, aes(x = time, y = surv, color=variable)) +
      geom_step() +
      treeTheme() +
      scale_y_continuous(limits = c(0,1)) +
      ylab("Survival rate") +
      labs(title = "Survival plot",
           subtitle = "")
    ggpubr::ggpar(g,  palette = palette)
}

getGICBreaks <- function(mH) {
    return(c(mH$GIC[1], min(mH$GIC), mH$GIC[length(mH$GIC)]) %>%
        unique() %>% round())
}

# -------------------
# GIC plot - bottom left plot

#' GIC plot
#'
#' @description Plots Generalized Information Criterion for models on the Factor Merger Tree.
#' @param factorMerger object of a class \code{factorMerger}
#' @param color GIC plot color.
#' @param statistic cluster split statistic
#' @param penalty GIC penalty
#' @importFrom grDevices hcl
#'
#' @export
plotGIC <- function(factorMerger, color, penalty = 2, statistic) {
    if (is.null(color)) {
        color <- "#762A83"
    }
    mH <- mergingHistory(factorMerger, T, F)
    mH$GIC <- -2 * mH$model + penalty * nrow(mH):1
    minGIC <- min(mH$GIC)
    yBreaks <- getGICBreaks(mH)
    minModel <- mH$model[which.min(mH$GIC)]
    g <- mH %>%
        ggplot(aes_string(x = getStatNameInTable(statistic),
                          y = "GIC")) +
        geom_line(col = color, size = 1) +
        geom_point(col = color, size = 1.5) +
        geom_point(x = minModel, y = minGIC,
                   col = color, size = 2.5) +
        geom_ribbon(aes_string(x = getStatNameInTable(statistic),
                               ymin = "minGIC", ymax = "GIC"),
                    fill = color, alpha = 0.2) +
        treeTheme() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(),
              panel.grid.major = element_blank()) +
        scale_y_continuous(position = "right",
                           breaks = yBreaks)

    annotation <- data.frame(
        xpos = max(mH[, getStatNameInTable(statistic)]),
        ypos = max(mH[, "GIC"]),
        text = paste0("GIC penalty = ", round(penalty, 1)),
        hjust = 1, vjust = 1)

    g <- g +
        geom_label(data = annotation,
                   aes(x = xpos, y = ypos, hjust = hjust,
                       vjust = vjust, label = text)) +
        ylab("")

    if (statistic == "p-value") {
        g <- g + scale_x_log10()
    }

    return(g)
}

ggDefaultPalette <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' TukeyHSD Plot
#'
#' @description TODO: Aga
#'
#' @param factorMerger object of a class \code{factorMerger}
#' @param palette RColorBrewer color palette
#'
#' @export
plotTukey <- function(factorMerger, palette = NULL) {

    response <- factorMerger$response
    factor <- factorMerger$factor
    tukeyGroups <- getTukeyGroups(response, factor)
    tukeyGroups$level <- rownames(tukeyGroups)
    tukeyGroups[tukeyGroups == FALSE] <- NA
    levelsOrder <- getFinalOrderVec(factorMerger)
    tukeyGroups <- tukeyGroups %>% melt(id.vars = "level") %>%
        filter(!is.na(value))
    
    cols <- scales::hue_pal()(length(levels(tukeyGroups$variable)))

    tukPlot <- tukeyGroups %>%
        ggplot(aes(variable, factor(level, levels = levelsOrder))) +
        geom_tile(aes(fill = paste(value, variable), color = "#000000"),
                  show.legend = FALSE, na.rm = TRUE) +
        scale_color_manual(values = "#000000") +
        scale_y_discrete(expand = c(0, 0)) +
        ggtitle(label = "Tukey HSD test", subtitle = " ") +
        xlab("") + treeTheme(FALSE) + theme(axis.text.y = element_blank())

    if (is.null(palette)) {
        return(tukPlot + 
                   scale_fill_manual(values = rev(cols)))
    }
    return(tukPlot + scale_fill_brewer(direction = -1, palette = palette))
}


# -------------------
# Frequency plot

#' Frequency plot
#'
#' @description Plots barplot with group frequencies.
#'
#' @param factorMerger object of a class \code{factorMerger}
#' @param color Boolean. If \code{TRUE}, the default, there is added aesthetic group corresponding
#' to the final cluster split.
#' @param clusterSplit final cluster split definition. A list with two fields:
#' \itemize{
#'     \item \code{stat} cluster statistic (available statistics are: \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}),
#'     \item \code{value} cut threshold / GIC penalty
#' }
#'
#' @export
plotFrequency <- function(factorMerger, color, clusterSplit) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- data.frame(group = factorMerger$factor)
    df <- df %>% count(group)

    if (color) {
        df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                     clusterSplit[[1]],
                                                     clusterSplit[[2]]),
                               by = c("group" = "orig"))
        g <- df %>% ggplot() +
            geom_bar(aes(x = group, y = n, fill = pred), stat = "identity")
    } else {
        g <- df %>% ggplot() +
            geom_bar(aes(x = group, y = n), fill = "navy", stat = "identity")
    }

    g + scale_y_continuous(name = "") +
        scale_x_discrete(expand = c(0, 0)) +
        coord_flip() +
        treeTheme() +
        theme(panel.grid.major.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank()) +
        xlab(" ") + ylab(" ") +
        labs(title = "Groups frequencies", subtitle = " ")
}

plotTable <- function(tab) {
    return(ggpubr::ggtexttable(tab, theme =
                                   ggpubr::ttheme(base_colour = "gray50",
                                                  base_size = 25,
                                                  "classic")))
}
