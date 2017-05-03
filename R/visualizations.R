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
                              show = "loglikelihood", nodesSpacing = "equidistant",
                              summary = NULL, color = "none",
                              clusterSplit = list(stat = "GIC", value = 2),
                              markBestModel = TRUE, markStars = TRUE,
                              alpha = 0.05, palette = NULL) {

    stopifnot(panel %in% c("all", "response", "GIC", "merging"))
    stopifnot(show %in% c("loglikelihood", "p-value"))
    stopifnot(nodesSpacing %in% c("equidistant", "effects", "modelSpecific"))
    stopifnot(color %in% c("none", "cluster", "summary"))

    summary <- checkSummary(factorMerger, summary)
    responsePlot <- plotResponse(factorMerger, summary, color)

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
                                alpha, color, colorsDf, palette)

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

checkSummary <- function(factorMerger, summary) {

}

plotResponse <- function(factorMerger, summary, color) {

}

#' @importFrom ggplot2 theme_classic theme element_line element_blank theme_minimal element_text
treeTheme <- function() {
    myTheme <- theme_minimal() +
        theme(panel.grid.major.x = element_line(color = "lightgrey",
                                                linetype = 2),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    return(myTheme)
}

plotTree <- function(factorMerger, show, nodesSpacing,
                     clusterSplit, markBestModel, markStars,
                     alpha, color, colorsDf, palette) {
    stopifnot(show %in% c("loglikelihood", "p-value"))
    if (nodesSpacing == "equidistance") {
        return(
            plotSimpleTree(factorMerger, show, clusterSplit,
                           markBestModel, markStars,
                           alpha, color, colorsDf, palette)
        )
    }
    return(
        plotCustomizedTree(factorMerger, show, clusterSplit,
                           nodesSpacing, markBestModel, markStars,
                           alpha, color, colorsDf, palette)
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
                              "equidistance", markBestModel, markStars,
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
    if (is.null(nodesPosition)) {
        nodesPosition <- groupsStats(factorMerger)
    }
    if (nodesSpacing == "modelSpecific") {
        nodesPosition <- applyModelTransformation(factorMerger, nodesPosition)
    }
    segment <- getTreeSegmentDf(factorMerger, statisticColname, nodesPosition)
    df <- segment$df
    pointsDf <- segment$pointsDf
    labelsDf <- segment$labelsDf
    labelsDf <- labelsDf %>% arrange(-y1)
    showY <- nodesSpacing != "equidistance"

    g <- df %>% ggplot() + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2)) +
        geom_point(data = pointsDf, aes(x = x1, y = y1), size = 0.75)

    g <- g + scale_y_continuous(limits = getLimits(labelsDf, showY),
                                position = "right",
                                breaks = labelsDf$y1,
                                labels = getLabels(labelsDf, factorMerger)) +
        ylab(getStatisticName(factorMerger))

    upperBreaks <- df$x1 %>% unique() %>% sort
    g <- g + xlab(show) +
        labs(title = "Merging path plot",
             subtitle = paste0("Optimal GIC partition: ",
                               paste(getOptimalPartition(factorMerger,
                                                         clusterSplit[[1]],
                                                         clusterSplit[[2]]),
                                     collapse = ":"))) + treeTheme()

    if (color == "cluster") {
        segmentColoured <- getClustersColors(segment, factorMerger, clusterSplit, show)
        g <- g +
            geom_segment(data = segmentColoured$df,
                         aes(x = x1, y = y1, xend = x2, yend = y2, col = pred)) +
            geom_point(data = segmentColoured$pointsDf,
                       aes(x = x1, y = y1, col = pred), size = 0.75)
        nGroups <- length(unique(segmentColoured$labelsDf$pred))
        clusterColors <- factor(segmentColoured$labelsDf$pred, labels = hue_pal()(nGroups)) %>%
            as.character()
        g <- g + theme(axis.text.y = element_text(color = clusterColors[length(clusterColors):1]))
    }

    if (show == "p-value") {
        g <- g + scale_x_log10()
    }

    if (show == "loglikelihood") {
        labBr <- getChisqBreaks(g$data, alpha)
        g <- g +
            scale_x_continuous(breaks = labBr$breaks, labels = labBr$labels)
    }

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
