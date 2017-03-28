customPalette <- "PRGn"
customPaletteValues <- c("#762A83", "#9970AB", "#C2A5CF", "#E7D4E8",
                         "#D9F0D3", '#A6DBA0', "#5AAE61", "#1B7837")

#' @importFrom ggplot2 theme_classic theme element_line element_blank theme_minimal element_text
treeTheme <- function(ticksColors) {
    myTheme <- theme_minimal() +
        theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.position = "none")
    if(!is.null(ticksColors)) {
        myTheme <- myTheme + theme(axis.text.y = element_text(color = ticksColors))
    } else {

    }
  return(myTheme)
}


#' @export
#'
plotTree <- function(factorMerger, stat = "model",
                     simplify = TRUE, alpha = 0.05,
                     showDiagnostics = TRUE) {
    stopifnot(alpha > 0 && alpha < 1)
    .plotTree(factorMerger, stat, NULL, simplify, alpha = alpha,
              showDiagnostics = showDiagnostics)
}

.plotTree <- function(factorMerger, stat,
                      levels, simplify,
                      alpha = 0.05,
                      showDiagnostics = TRUE) {
    UseMethod(".plotTree", factorMerger)
}

#' @export
.plotTree.factorMerger <- function(factorMerger, stat = "model",
                                   levels = NULL, simplify = TRUE,
                                   alpha = 0.05, showDiagnostics = TRUE) {
    stopifnot(stat %in% c("model", "pval"))

    if (simplify) {
        plotSimpleTree(factorMerger, stat, levels, alpha, showDiagnostics = showDiagnostics)
    }
    else {
        plotCustomizedTree(factorMerger, stat, groupsStats(factorMerger),
                           levels, alpha = alpha,
                           showDiagnostics = showDiagnostics)
    }
}

renameStat <- function(stat) {
    switch(stat,
           "pval" = {
               stat <- "log10(p-value)"
           },
           "model" = {
               stat <- "loglikelihood"
           },
           stat)
    return(stat)
}

getLimits <- function(df, showY) {
    if (showY) {
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
getLabels <- function(pointsDf, factorMerger) {
    stats <- groupsStats(factorMerger)
    if (is.null(stats)) return(NULL)
    stats$label <- rownames(stats)
    pointsDf <- pointsDf %>% left_join(stats, by = "label")
    return(paste0(pointsDf$label, ": ", round(pointsDf[, ncol(pointsDf)], 2)))
}

#' @importFrom dplyr left_join
getLabelsColors <- function(pointsDf, levels) {
    if (is.null(levels)) {
        return(NULL)
    }
    colors <- data.frame(
        label = levels,
        color = colorRamps::magenta2green(nrow(pointsDf)),
        stringsAsFactors = FALSE)
    return((pointsDf %>%
                left_join(colors, by = "label"))$color %>%
                       as.character())
}

getTreeSegmentDf <- function(factorMerger, stat, pos) {
    factor <- factorMerger$factor
    noGroups <- length(levels(factor))
    df <- pos[1:noGroups, ] %>% data.frame
    colnames(df) <- "y1"
    df$y2 <- df$y1
    df$label <- rownames(pos)[1:noGroups]
    df$x1 <- factorMerger$mergingList$`1`$modelStats[, stat]
    df$x2 <- NA
    pointsDf <- df
    merging <- mergingHistory(factorMerger)

    for (step in 1:nrow(merging)) {
        statVal <- factorMerger$mergingList[[step + 1]]$modelStats[, stat]
        pair <- merging[step, ]
        whichDf <- which(df$label %in% pair)
        df[whichDf, "x2"] <- statVal
        pairLabel <- paste(pair, collapse = "")
        pairMean <- pos[rownames(pos) == pairLabel,]
        crosswise <- c(df$y1[whichDf], "", statVal, statVal)
        newVertex <- c(pairMean, pairMean, pairLabel, statVal, NA)
        df <- rbind(df, crosswise)
        df <- rbind(df, newVertex)
    }

    df <- df %>% subset(select = -label) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    if (stat == "pval") {
        df$x1 <- log10(df$x1)
        df$x2 <- log10(df$x2)
        pointsDf$x1 <- log10(pointsDf$x1)
    }

    df$x2[is.na(df$x2)] <- min(df$x2, na.rm = TRUE) - ((max(df$x2, na.rm = TRUE) - min(df$x1)) / 20)

    return(list(df = df,
                pointsDf = pointsDf))
}

#' @importFrom dplyr mutate filter
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot geom_segment scale_x_log10 theme_bw scale_x_continuous
#' @importFrom ggplot2 coord_flip xlab ylab theme element_blank geom_vline geom_label
#' @importFrom ggplot2 geom_point aes geom_label scale_fill_manual scale_y_continuous
plotCustomizedTree <- function(factorMerger, stat = "model",
                               pos, levels = NULL,
                               showY = TRUE, alpha = 0.05,
                               showDiagnostics = TRUE) {

    segment <- getTreeSegmentDf(factorMerger, stat, pos)
    df <- segment$df
    pointsDf <- segment$pointsDf

    g <- df %>% ggplot() +
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2)) +
        geom_point(data = pointsDf, aes(x = x1, y = y1)) +
        scale_y_continuous(limits = getLimits(pointsDf, showY),
                           position = "right",
                           breaks = pointsDf$y1,
                           labels = getLabels(pointsDf, factorMerger)) +
        ylab(getStatisticName(factorMerger))

    upperBreaks <- df$x1 %>% unique() %>% sort
    g <- g + xlab(renameStat(stat)) + treeTheme(getLabelsColors(pointsDf, levels))

    if (showDiagnostics) {
        if (stat == "pval") {
            intercept <- log10(alpha)
            label <- paste0("alpha = ", alpha)
        }
        if (stat == "model") {
            aicMin <- mergingHistory(factorMerger, TRUE)[, c("model", "AIC")] %>%
                filter(AIC == min(AIC))
            intercept <- aicMin$model
            label <- paste0("min AIC")
            right <- g$data$x1 %>% max()
            left <- g$data$x2 %>% min()
            g <- g +
                scale_x_continuous(
                    breaks = seq(left, right, qchisq(1 - alpha, 1)),
                    labels = seq(left, right, qchisq(1 - alpha, 1)) %>% round())
        }
        y <- getLimits(pointsDf, showY)

        g <- g + geom_vline(xintercept = intercept, col = "mediumorchid3", linetype = "dotted") +
            geom_label(x = intercept, y = getLimits(pointsDf, showY)[1],
                      label = label, alpha = 0.5, col = "mediumorchid3",
                      angle = 90,
                      size = 3, fontface = "italic")
    }

    return(g)
}

plotSimpleTree <- function(factorMerger, stat = "model",
                           levels = NULL, alpha = 0.05,
                           showDiagnostics = TRUE) {
    pos <- getFinalOrder(factorMerger) %>% data.frame()
    merging <- mergingHistory(factorMerger)
    noStep <- nrow(merging)

    for (step in 1:noStep) {
        pos <- rbind(pos, mean(pos[rownames(pos) %in% merging[step, ],]))
        rownames(pos)[nrow(pos)] <- paste(merging[step, ], collapse = "")
    }
    return(plotCustomizedTree(factorMerger, stat, pos,
                              levels, showY = FALSE,
                              alpha, showDiagnostics = showDiagnostics))
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

#' @export
#' @importFrom gridExtra grid.arrange
appendToTree <- function(factorMerger, plot) {
    UseMethod("appendToTree", plot)
}

#' @export
appendToTree.default <- function(factorMerger, plot) {
    grid.arrange(.plotTree(factorMerger, showDiagnostics = FALSE, simplify = TRUE), plot, ncol = 2)
}

#' @export
appendToTree.profilePlot <- function(factorMerger, plot) {
    lev <- levels(plot$data$level)
    grid.arrange(.plotTree(factorMerger, levels = lev, showDiagnostics = FALSE, simplify = TRUE), plot, ncol = 2)
}

#' @importFrom ggplot2 ggplot aes geom_line geom_text theme_minimal theme scale_color_manual
#' @export
# TODO: dodać wybór między rank a średnią
plotProfile <- function(factorMerger) {
    stat <- findSimilarities(factorMerger)

    stat$level <- factor(stat$level,
                             levels = (stat %>%
                                           filter(variable == levels(stat$variable) %>%
                                                      head(1)) %>% arrange(rank))$level
        )

    noLevels <- length(levels(stat$level))
    stat$rank <- factor(stat$rank, levels = 1:noLevels)
    g <- stat %>% ggplot(aes(x = variable, y = rank, col = level, group = level, label = level)) +
        geom_line() +
        geom_text(data = subset(stat,
                                variable == levels(stat$variable) %>% tail(1)),
                  aes(x = variable),
                  size = 3.5, hjust = 0.8,  nudge_x = 0.1) +
        theme_minimal() + theme(legend.position = "none") +
        scale_color_manual(values = colorRamps::magenta2green(noLevels))
    class(g) <- append(class(g), "profilePlot")
    return(g)
}

#' @export
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_minimal scale_x_continuous theme
#' @importFrom ggplot2 coord_flip element_line element_blank scale_fill_distiller
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange
plotHeatmap <- function(factorMerger) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    findSimilarities(factorMerger) %>%
            ggplot() +
        geom_tile(aes(x = level, y = variable, fill = mean)) +
        coord_flip() + theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        scale_fill_distiller(palette = customPalette)
}

#' @export
#' @importFrom ggplot2 ggplot geom_boxplot aes coord_flip
#' @importFrom dplyr group_by summarize left_join
plotBoxplot <- function(factorMerger) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    data <- data.frame(x = factorMerger$factor, y = factorMerger$response)
    data %>% left_join(
        data %>% group_by(x) %>% summarize(y0 = min(y),
                                           y25 = quantile(y, 0.25),
                                           y50 = mean(y),
                                           y75 = quantile(y, 0.75),
                                           y100 = max(y)), by = "x") %>%
        ggplot(aes(y = y, x = x, group = x)) +
        geom_boxplot(aes(ymin = y0,
                         lower = y25,
                         middle = y50,
                         upper = y75,
                         ymax = y100), stat = "identity") +
        coord_flip() + treeTheme(NULL) +
        theme(axis.title = element_blank(), axis.text.y = element_blank())
}

#' @export
#' @importFrom ggplot2 ggplot geom_bar aes coord_flip scale_fill_manual theme theme element_blank scale_y_continuous
#' @importFrom dplyr group_by summarize left_join
plotProportion <- function(factorMerger) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    data <- data.frame(x = factorMerger$factor, y = factorMerger$response)
    data %>% ggplot() + geom_bar(aes(x = x, fill = as.factor(y)), position = "fill") +
        scale_y_continuous(label = scales::percent, name = "Success proportion") +
        coord_flip() + treeTheme(NULL) +
        scale_fill_manual(values = customPaletteValues[c(2, length(customPaletteValues) - 1)]) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank())
}

#' @export
plotSurvival <- function(factorMerger) {
    model <- calculateModel(factorMerger, factorMerger$factor)
    survminer::ggcoxadjustedcurves(model, data = data.frame(factorMerger$factor),
                        individual.curves = TRUE,
                        theme = treeTheme(NULL),
                        palette = "RdBu", curve.size = 1) +
        treeTheme(NULL)

    # variable = fdsfds
}
