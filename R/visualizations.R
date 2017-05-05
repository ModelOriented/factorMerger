#' Plot Factor Merger
#'
#' @param panel Type of panels to be plot. Possible values are \code{c("all", "response", "GIC", "merging")}.
#' All types of plots include the Merging Path Plot. Apart from the Merging Path Plot there are
#' also two possible panels: the Response Plot (response summary, specific to the model family),
#' the GIC Plot (GIC vs. loglikelihood/p-value).
#' \itemize{
#' \item \code{"merging"} plots the Merging Path Plot only,
#' \item \code{"response"} plots the Merging Path Plot and the Response Plot,
#' \item \code{"GIC"} plots the Merging Path Plot and the GIC Plot,
#' \item \code{"all"} plots all panels and a short summary of the full model.
#' }
#' @param show Statistic to be displayed on the OX axis. Possible values are \code{c("loglikelihood", "p-value")}.
#' If \code{"p-value"} is chosen p-value for the Likelihood Ratio Test against the full model is plot on the OX axis.
#' @param nodesSpacing Type of nodes vertical spacing. May be chosen from
#'  \code{c("equidistant", "effects", "modelSpecific")}. \code{"effects"} arranges nodes according to
#'  the model coefficients estimatiors (e.g. in Gaussian case on the OY axis group means are plotted).
#'
#' # TODO: Maybe different names? Implement "modelSpecific".
#' @param summary Response panel type -- accepts the following values dependent on the model family:
#' \itemize{
#' \item multi dimensional Gaussian: \code{c("heatmap", "profile")},
#' \item single dimensional Gaussian: \code{c("means", "boxplot")},
#' \item binomial: \code{c("proportion")},
#' \item survival: \code{c("survival")}
#' }
#' @param clusterSplit List of a length two:
#' \itemize{
#' \item \code{stat} - statistic used in the bottom-up search. Available statistics are:
#' \code{"loglikelihood"}, \code{"pvalue"}, \code{"GIC"}.
#' \item \code{value} cut threshold / GIC penalty
#' }
#' This parameter is considered only with \code{color = "cluster"}.
#' If \code{stat} is \code{"loglikelihood"} or \code{"pvalue"} performs bottom-up search through models
#' on the merging path until spots a model scored worse than the given \code{value}. If \code{stat = "GIC"}
#' treats \code{value} as GIC penalty and returns optimal GIC partition.
#' @param color Specifies how the Merging Path Plot should be colored. Possible values are \code{c("none", "cluster")}.
#' If \code{color = "cluster"} colors obtained tree according to the \code{clusterSplit}.
#' @param markBestModel Boolean. If \code{TRUE}, the default, plots vertical line crossing the optimal model
#' according to the \code{clusterSplit}.
#' @param markSignificanceStars Boolean. If \code{TRUE}, the default, marks models that are significantly
#' worse than their predecessors on the Merging Path Plot (uses the Likelihood Ratio Test).
#'
#' Significance codes are:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1.
#' @param alpha Significance level. If \code{show = "likelihood"} each interval on
#'  the OX axis of the Merging Path Plot corresponds to the 1 - \code{alpha} quantile of chi-square distribution.
#' @param penalty GIC penalty used in the GIC plot. If \code{clusterSplit$stat = "GIC"} it is
#' recommended to use \code{clusterSplit$value} as \code{penalty}.
#'
#' TODO: Is it a good idea?
#' @param mergingPalette Color palette used in the Merging Path Plot.
#' @param responsePalette Color palette used in the Response Plot.
#' @param GICcolor Color used in the GIC plot.
#' @importFrom gridExtra grid.arrange
#'
#' @export
plot.factorMerger <- function(factorMerger, panel = "all",
                              show = "loglikelihood",
                              nodesSpacing = "equidistant",
                              summary = NULL, color = "none",
                              clusterSplit = list(stat = "GIC", value = 2),
                              markBestModel = TRUE,
                              markSignificanceStars = TRUE,
                              alpha = 0.05, penalty = 2,
                              mergingPalette = NULL,
                              responsePalette = NULL,
                              GICcolor = NULL) {

    stopifnot(panel %in% c("all", "response", "GIC", "merging"))
    stopifnot(show %in% c("loglikelihood", "p-value"))
    stopifnot(nodesSpacing %in% c("equidistant", "effects", "modelSpecific"))
    stopifnot(color %in% c("none", "cluster"))

    summary <- checkSummary(factorMerger, summary)
    responsePlot <- plotResponse(factorMerger, summary, color, clusterSplit)
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
    } else {
        colorsDf <- NULL
    }

    mergingPathPlot <- plotTree(factorMerger, show, nodesSpacing,
                                clusterSplit, markBestModel, markSignificanceStars,
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
                                   plotGIC(factorMerger, GICcolor, penalty, show), ncol = 2,
                                   widths = c(7.5, 2.5), heights = c(7.5, 2.5)))
           },
           "response" = {
               return(grid.arrange(mergingPathPlot, responsePlot,
                                   ncol = 2,
                                   widths = c(7.5, 2.5)))
           },
           "GIC" = {
               return(grid.arrange(mergingPathPlot, plotGIC(factorMerger, GICcolor, penalty, show),
                                   ncol = 1,
                                   heights = c(7.5, 2.5)))
           })
}

getColorsFromResponseDf <- function(responsePlot) {

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
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.5) +
        geom_point(data = pointsDf, aes(x = x1, y = y1), size = 0.75)+
        scale_y_continuous(limits = getLimits(labelsDf, showY),
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
        mark <- markOptimalModel(factorMerger, clusterSplit, show, alpha)
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

markOptimalModel <- function(factorMerger, clusterSplit, show, alpha) {
    optimalNSteps <- optimalNumberOfMerges(factorMerger,
                                           clusterSplit[[1]],
                                           clusterSplit[[2]])
    mH <- mergingHistory(factorMerger, T, F)
    intercept <- mH[optimalNSteps + 1, getStatNameInTable(show)]
    if (show == "p-value") {
        label <- paste0("alpha = ", round(intercept, 2))
        labelIntercept <- log10(intercept)
    }
    if (show == "loglikelihood") {
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

plotResponse <- function(factorMerger, summary, color, clusterSplit) {
    switch(summary,
           "heatmap" = {
               return(plotHeatmap(factorMerger, color, clusterSplit))
           },
           "profile" = {
               return(plotProfile(factorMerger, color, clusterSplit))
           },
           "boxplot" = {
               return(plotBoxplot(factorMerger, color, clusterSplit))
           },
           "means" = {
               return(plotMeansAndStds(factorMerger, color, clusterSplit))
           },
           "survival" = {
               return(plotSurvival(factorMerger, color, clusterSplit))
           },
           "proportion" = {
               return(plotProportion(factorMerger, color, clusterSplit))
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
#' @importFrom magrittr %>%
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
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 12),
              legend.position = "none") +
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

    switch(color,
           "cluster" = {
               df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                            clusterSplit[[1]],
                                                            clusterSplit[[2]]),
                                      by = c("level" = "orig"))
               g <- df %>% ggplot(aes(x = variable, y = rank,
                                      colour = pred, group = group, label = group))
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

#' Boxplot (single-dimensional Gaussian)
#'
#' @description Plots boxplot with mean as a summary statistic groupping observation by factor levels.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_boxplot aes coord_flip labs ylab xlab
#' @importFrom dplyr group_by summarize left_join
plotBoxplot <- function(factorMerger, color = "none", clusterSplit) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- data.frame(group = factorMerger$factor, y = factorMerger$response)
    df <- calculateBoxPlotMoments(df)
    switch(color,
           "cluster" = {
               df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                            clusterSplit[[1]],
                                                            clusterSplit[[2]]),
                                      by = c("group" = "orig"))
               g <- df %>% ggplot(aes(y = y, x = group, group = group, fill = pred))
           },
           "none" = {
               g <- df %>% ggplot(aes(y = y, x = group, group = group))
           })

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
plotMeansAndStds <- function(factorMerger, color = "none", clusterSplit) {
    factor <- factor(factorMerger$factor, levels = getFinalOrderVec(factorMerger))
    df <- getMeansAndStds(factorMerger, factor)

    switch(color,
           "cluster" = {
               df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                            clusterSplit[[1]],
                                                            clusterSplit[[2]]),
                                      by = c("group" = "orig"))
               g <- df %>% ggplot(aes(colour = pred, fill = pred,
                                      x = as.factor(group),
                                      y = mean, group = as.factor(group)))
           },
           "none" = {
               g <- df %>% ggplot(aes(x = as.factor(group), colour = as.factor(group),
                                      y = mean, group = as.factor(group)))
           })

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
plotProportion <- function(factorMerger, color = "none", clusterSplit) {
    levels <- getFinalOrderVec(factorMerger)
    responseLevels <- factorMerger$response %>% as.factor() %>% levels()
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    df <- data.frame(group = factorMerger$factor, y = factorMerger$response)
    switch(color,
           "cluster" = {
               df <- df %>% group_by(group) %>% summarise(mean = mean(y == 1))
               df <- df %>% left_join(getOptimalPartitionDf(factorMerger,
                                                            clusterSplit[[1]],
                                                            clusterSplit[[2]]),
                                      by = c("group" = "orig"))
               g <- df %>% ggplot() +
                   geom_bar(aes(x = group, y = mean, fill = pred), stat = "identity")
           },
           "none" = {
               g <- df %>% ggplot() + geom_bar(aes(x = group, fill = as.factor(y)), position = "fill")
           }
    )
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

    if (color == "cluster") {
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
#' @description Plots Generalized Information Criterion for models on the Merging Path Plot.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line aes theme element_blank scale_y_continuous labs geom_point geom_ribbon
#'
#' @export
plotGIC <- function(factorMerger, color, penalty = 2, show) {
    if (is.null(color)) {
        color <- "#762A83"
    }
    mH <- mergingHistory(factorMerger, T, F)
    mH$GIC <- -2 * mH$model + penalty * nrow(mH):1
    minGIC <- min(mH$GIC)
    yBreaks <- c(mH$GIC[1], minGIC, mH$GIC[length(mH$GIC)]) %>% unique() %>% round()
    minModel <- mH$model[which.min(mH$GIC)]
    g <- mH %>% ggplot(aes_string(x = getStatNameInTable(show), y = "GIC")) +
        geom_line(col = color, size = 1) +
        geom_point(col = color, size = 1.5) +
        geom_point(x = minModel, y = minGIC, col = color, size = 2.5) +
        geom_ribbon(aes(x = model, ymin = minGIC, ymax = GIC),
                    fill = color, alpha = 0.2) +
        treeTheme() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(),
              panel.grid.major = element_blank()) +
        scale_y_continuous(position = "right", breaks = yBreaks)
    return(g)
}


