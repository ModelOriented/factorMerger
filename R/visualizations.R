#' Plot Factor Merger
#'
#' @param factorMerger object of a class \code{factorMerger}.
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
#' \item single dimensional Gaussian: \code{c("means", "boxplot")},
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
#' @param chisqQuantile Significance level. If \code{statistic = "likelihood"} each interval on
#'  the OX axis of the Factor Merger Tree corresponds to the 1 - \code{chisqQuantile}
#'  quantile of chi-square distribution.
#'
#' @param palette Color palette used in the Factor Merger Tree and the Response Plot.
#' @param responsePanelPalette Additional color palette used in the Response Plot if
#' palettes for the Factor Merger Tree and the Response Plot are to be different.
#' @param gicPanelColor Color used in the GIC plot.
#' @importFrom gridExtra grid.arrange
#'
#' @export
plot.factorMerger <- function(factorMerger, panel = "all",
                              statistic = "loglikelihood",
                              nodesSpacing = "equidistant",
                              colorClusters = TRUE,
                              splitStatistic = "GIC",
                              splitThreshold = NULL,
                              penalty = 2,
                              showSplit = FALSE,
                              showSignificance = TRUE,
                              chisqQuantile = 0.05,
                              palette = NULL,
                              responsePanel = NULL,
                              responsePanelPalette = NULL,
                              gicPanelColor = NULL) {

    stopifnot(panel %in% c("all", "response", "GIC", "tree"))
    stopifnot(statistic %in% c("loglikelihood", "p-value"))
    stopifnot(nodesSpacing %in% c("equidistant", "effects", "modelSpecific"))
    stopifnot(splitStatistic == "GIC" | !is.null(splitThreshold))

    clusterSplit <- getClusterSplit(splitStatistic, splitThreshold, penalty)

    responsePanel <- checkResponsePanel(factorMerger, responsePanel)
    responsePlot <- plotResponse(factorMerger, responsePanel, colorClusters, clusterSplit)

    if (is.null(responsePanelPalette) && !is.null(palette)) {
        responsePanelPalette <- palette
    }

    if(!is.null(responsePanelPalette)) {
        responsePlot <- ggpubr::ggpar(responsePlot, palette = responsePanelPalette)
    }

    # Set colors
    if (colorClusters) {
        colorsDf <- getOptimalPartitionDf(factorMerger,
                                          clusterSplit[[1]],
                                          clusterSplit[[2]])
    } else {
        colorsDf <- NULL
    }

    mergingPathPlot <- plotTree(factorMerger, statistic, nodesSpacing,
                                clusterSplit, showSplit, showSignificance,
                                chisqQuantile, colorClusters, colorsDf, palette)

    if (!is.null(palette)) {
        mergingPathPlot <- ggpubr::ggpar(mergingPathPlot, palette = palette)
    }
    switch(panel,
           "tree" = {
               return(mergingPathPlot)
           },
           "all" = {
               return(grid.arrange(mergingPathPlot, responsePlot,
                                   plotGIC(factorMerger, gicPanelColor, penalty, statistic), ncol = 2,
                                   widths = c(7.5, 2.5), heights = c(7.5, 2.5)))
           },
           "response" = {
               return(grid.arrange(mergingPathPlot, responsePlot,
                                   ncol = 2,
                                   widths = c(7.5, 2.5)))
           },
           "GIC" = {
               return(grid.arrange(mergingPathPlot, plotGIC(factorMerger, gicPanelColor, penalty, statistic),
                                   ncol = 1,
                                   heights = c(7.5, 2.5)))
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
              axis.text.x = element_text(size = 12),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none")
    return(myTheme)
}

plotTree <- function(factorMerger, statistic, nodesSpacing,
                     clusterSplit, markBestModel, markStars,
                     alpha, color, colorsDf, palette) {
    stopifnot(statistic %in% c("loglikelihood", "p-value"))
    if (nodesSpacing == "equidistant") {
        return(
            plotSimpleTree(factorMerger, statistic, clusterSplit,
                           markBestModel, markStars,
                           alpha, color, colorsDf, palette)
        )
    }
    return(
        plotCustomizedTree(factorMerger, statistic, clusterSplit,
                           nodesSpacing, markBestModel, markStars,
                           alpha, color, colorsDf, palette, groupsStats(factorMerger))
    )
}

plotSimpleTree <- function(factorMerger, statistic, clusterSplit,
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
    return(plotCustomizedTree(factorMerger, statistic, clusterSplit,
                              "equidistant", markBestModel, markStars,
                              alpha, color, colorsDf, palette,
                              nodesPosition))

}

#' @importFrom dplyr mutate filter group_by count
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot geom_segment scale_x_log10 theme_bw scale_x_continuous
#' @importFrom ggplot2 coord_flip xlab ylab theme element_blank geom_vline geom_text
#' @importFrom ggplot2 geom_point aes geom_label scale_fill_manual scale_y_continuous labs
plotCustomizedTree <- function(factorMerger, statistic, clusterSplit,
                               nodesSpacing, markBestModel, markStars,
                               alpha, color, colorsDf, palette,
                               nodesPosition = NULL) {
    statisticColname <- getStatNameInTable(statistic)
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
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.5) +
        geom_point(data = pointsDf, aes(x = x1, y = y1), size = 0.75)+
        scale_y_continuous(limits = getLimits(labelsDf, showY),
                                position = "right",
                                breaks = labelsDf$y1,
                                labels = getLabels(labelsDf, factorMerger)) +
        ylab(getStatisticName(factorMerger)) + xlab(statistic) +
        labs(title = "Factor Merger Tree",
             subtitle = paste0("Optimal GIC partition: ",
                               paste(getOptimalPartition(factorMerger,
                                                         clusterSplit[[1]],
                                                         clusterSplit[[2]]),
                                     collapse = ":"))) + treeTheme()

    if (color) {
        g <- addClustersColors(g, segment, factorMerger, clusterSplit, statistic, palette)
    }

    g <- scaleAxis(g, statistic, alpha)

    if (markBestModel) {
        mark <- markOptimalModel(factorMerger, clusterSplit, statistic, alpha)
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

addClustersColors <- function(plot, segment, factorMerger, clusterSplit, statistic, palette) {
    segmentColoured <- getClustersColors(segment, factorMerger, clusterSplit, statistic)
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

checkResponsePanel.gaussianFactorMerger <- function(factorMerger, responsePanel) {
    if (NCOL(factorMerger$response) > 1) {
        responsePanelSet <-  c("heatmap", "profile")
    } else {
        responsePanelSet <- c("means", "boxplot")
    }
    warnIfUnexpectedResponsePanel(responsePanelSet, responsePanel)
}

checkResponsePanel.binomialFactorMerger <- function(factorMerger, responsePanel) {
    responsePanelSet <- c("proportion")
    warnIfUnexpectedResponsePanel(responsePanelSet, responsePanel)
}

checkResponsePanel.survivalFactorMerger <- function(factorMerger, responsePanel) {
    responsePanelSet <- c("survival")
    warnIfUnexpectedResponsePanel(responsePanelSet, responsePanel)
}

warnIfUnexpectedResponsePanel <- function(responsePanelSet, responsePanel) {
    if (is.null(responsePanel)) {
        return(responsePanelSet[1])
    }
    if (!(responsePanel %in% responsePanelSet)) {
        warning(paste0("ResponsePanel '", responsePanel, ", is not supported by supplied model family -- ", responsePanelSet[1], " used insted."))
        return(responsePanelSet[1])
    }
    return(responsePanel)
}

plotResponse <- function(factorMerger, responsePanel, colorClusters, clusterSplit) {
    switch(responsePanel,
           "heatmap" = {
               return(plotHeatmap(factorMerger, colorClusters, clusterSplit))
           },
           "profile" = {
               return(plotProfile(factorMerger, colorClusters, clusterSplit))
           },
           "boxplot" = {
               return(plotBoxplot(factorMerger, colorClusters, clusterSplit))
           },
           "means" = {
               return(plotMeansAndStds(factorMerger, colorClusters, clusterSplit))
           },
           "survival" = {
               return(plotSurvival(factorMerger, colorClusters, clusterSplit))
           },
           "proportion" = {
               return(plotProportion(factorMerger, colorClusters, clusterSplit))
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

#' Heatmap (multi-dimensional Gaussian)
#'
#' @description Plots heatmap for each dimension of the response variable. Vector of means of factor levels for a given
#' dimension is scaled to have mean equal to zero and standard deviation equal to one.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_minimal scale_x_continuous theme
#' @importFrom ggplot2 coord_flip element_line element_blank scale_fill_distiller labs guides
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange
plotHeatmap <- function(factorMerger, color, clusterSplit) {
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
              axis.text.x = element_text(size = 12),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none") +
        xlab("") +
        scale_fill_distiller(palette = "Greys") +
        labs(title = "Heatmap", subtitle = "Group means by variables")
}

#' Profile plot (multi-dimensional Gaussian)
#'
#' @description Plots rank plot - one series is a single factor level and one group
#' on the OX axis is a single dimension of the response.
#'
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_text theme_minimal theme scale_color_manual labs
#' @export
plotProfile <- function(factorMerger, color, clusterSplit) {
    df <- findSimilarities(factorMerger)

    df$group <- factor(df$level,
                       levels = (df %>%
                                     filter(variable == levels(df$variable) %>%
                                                head(1)) %>% arrange(rank))$level
    )

    noLevels <- length(levels(df$level))
    df$rank <- factor(df$rank, levels = 1:noLevels)

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
        geom_text(data = subset(df,
                                variable == levels(df$variable) %>% tail(1)),
                  aes(x = variable),
                  size = 4, hjust = 0.8,  nudge_x = 0.5) +
        ylab("") + xlab("") +
        labs(title = "Profile plot", subtitle = "Variable means ranks") +
        theme_minimal() + theme(legend.position = "none",
              plot.title = element_text(size = 18),
              axis.text = element_text(size = 12),
              plot.subtitle = element_text(size = 12))
    return(g)
}

#' Boxplot (single-dimensional Gaussian)
#'
#' @description Plots boxplot with mean as a summary statistic groupping observation by factor levels.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_boxplot aes coord_flip labs ylab xlab
#' @importFrom dplyr group_by summarize left_join
plotBoxplot <- function(factorMerger, color, clusterSplit) {
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

    g + geom_boxplot(aes(ymin = y0,
                         lower = y25,
                         middle = y50,
                         upper = y75,
                         ymax = y100), stat = "identity") +
        coord_flip() + treeTheme() + xlab("") + ylab("") +
        theme(axis.text.y = element_blank()) +
        labs(title = "Boxplot", subtitle = "Summary statistic: mean")
}

#' Means and standard deviation plot (single-dimensional Gaussian)
#'
#' @description For each factor level plots its mean and interval of the length equal to its standard deviation.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_boxplot aes coord_flip labs geom_errorbar theme ylab position_dodge element_blank element_text
#' @importFrom dplyr group_by summarize left_join
plotMeansAndStds <- function(factorMerger, color, clusterSplit) {
    factor <- factor(factorMerger$factor, levels = getFinalOrderVec(factorMerger))
    df <- getMeansAndStds(factorMerger, factor)

    if (color) {
        df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                     clusterSplit[[1]],
                                                     clusterSplit[[2]]),
                               by = c("group" = "orig"))
        g <- df %>% ggplot(aes(colour = pred, fill = pred,
                               x = as.factor(group),
                               y = mean, group = as.factor(group)))
    } else {
        g <- df %>% ggplot(aes(x = as.factor(group), colour = as.factor(group),
                               y = mean, group = as.factor(group)))
    }

    g + geom_errorbar(aes(ymin = left, ymax = right),
                      width = .5,
                      position = position_dodge(.5)) + treeTheme() +
        geom_point(size = 2) + coord_flip() +
        theme(axis.title.x = element_text(),
              axis.text.y = element_blank()) +
        labs(title = "Summary statistics",
             subtitle = "Means and standard deviations of coefficients' estimators") +
        ylab("")
}

#' Proportion plot (binomial)
#'
#' @description Plots proportion of success for each factor level.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_bar aes coord_flip scale_fill_manual theme theme element_blank scale_y_continuous labs
#' @importFrom dplyr group_by summarize left_join summarise
plotProportion <- function(factorMerger, color, clusterSplit) {
    levels <- getFinalOrderVec(factorMerger)
    responseLevels <- factorMerger$response %>% as.factor() %>% levels()
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
        g <- df %>% ggplot() + geom_bar(aes(x = group, fill = as.factor(y)), position = "fill")
    }

    g + scale_y_continuous(label = scales::percent, name = "") +
        coord_flip() + treeTheme() +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
        labs(title = "Success ratio",
             subtitle = "")
}

#' Survival plot (survival)
#'
#' @description Plots \code{ggcoxadjustedcurves} from the \code{survminer} package.
#'
#' @importFrom ggplot2 labs
#'
#' @export
plotSurvival <- function(factorMerger, color, clusterSplit) {
    levels <- getFinalOrderVec(factorMerger)
    df <- data.frame(group = (factor(factorMerger$factor, levels = levels) %>%
                         as.character()))

    if (color) {
        df$group <- cutTree(factorMerger, clusterSplit[[1]], clusterSplit[[2]])
    }

    model <- calculateModel(factorMerger, df$group)

    g <- survminer::ggcoxadjustedcurves(model, data = data.frame(factorMerger$factor),
                                        individual.curves = TRUE,
                                        variable = df$group,
                                        curve.size = 1) +
        treeTheme() + labs(title = "Survival plot",
                           subtitle = "Adjusted survival curves for coxph model")
    return(g)
}

# -------------------
# GIC plot - bottom left plot

#' GIC plot
#'
#' @description Plots Generalized Information Criterion for models on the Factor Merger Tree.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line aes theme element_blank scale_y_continuous labs geom_point geom_ribbon
#'
#' @export
plotGIC <- function(factorMerger, color, penalty = 2, statistic) {
    if (is.null(color)) {
        color <- "#762A83"
    }
    mH <- mergingHistory(factorMerger, T, F)
    mH$GIC <- -2 * mH$model + penalty * nrow(mH):1
    minGIC <- min(mH$GIC)
    yBreaks <- c(mH$GIC[1], minGIC, mH$GIC[length(mH$GIC)]) %>%
        unique() %>% round()
    minModel <- mH$model[which.min(mH$GIC)]
    g <- mH %>% ggplot(aes_string(x = getStatNameInTable(statistic), y = "GIC")) +
        geom_line(col = color, size = 1) +
        geom_point(col = color, size = 1.5) +
        geom_point(x = minModel, y = minGIC, col = color, size = 2.5) +
        geom_ribbon(aes_string(x = getStatNameInTable(statistic),
                               ymin = "minGIC", ymax = "GIC"),
                    fill = color, alpha = 0.2) +
        treeTheme() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(),
              panel.grid.major = element_blank()) +
        scale_y_continuous(position = "right", breaks = yBreaks)
    if (statistic == "p-value") {
        g <- g + scale_x_log10()
    }
    return(g)
}


