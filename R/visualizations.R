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
              axis.title.y = element_blank(),
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
               stat <- "p-value"
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
getLabels <- function(labelsDf, factorMerger) {
    stats <- groupsStats(factorMerger)
    if (is.null(stats)) return(NULL)
    stats$label <- rownames(stats)
    labelsDf <- labelsDf %>% left_join(stats, by = "label")
    return(paste0(labelsDf$label, ": ", round(labelsDf[, ncol(labelsDf)], 2)))
}

#' @importFrom dplyr left_join
getLabelsColors <- function(labelsDf, levels) {
    if (is.null(levels)) {
        return(NULL)
    }
    colors <- data.frame(
        label = levels,
        color = colorRamps::magenta2green(nrow(labelsDf)),
        stringsAsFactors = FALSE)
    return((labelsDf %>%
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
    labelsDf <- df
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
        pointsDf <- rbind(pointsDf, newVertex)
    }

    df <- df %>% subset(select = -label) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    df <- df[complete.cases(df), ]

    pointsDf <- subset(pointsDf, select = c(x1, y1))
    pointsDf <- pointsDf %>% apply(2, as.numeric) %>% as.data.frame()

    return(list(df = df,
                labelsDf = labelsDf,
                pointsDf = pointsDf))
}

nLabels <- 5

getChisqBreaks <- function(plotData, alpha) {
    right <- plotData$x1 %>% max()
    left <- plotData$x2 %>% min()
    breaks <- seq(left, right, qchisq(1 - alpha, 1))
    labels <- seq(left, right, qchisq(1 - alpha, 1)) %>% round()
    labels[setdiff(1:length(breaks), seq(1, length(breaks), length.out = nLabels) %>% round())] <- ""
    return(
        list(
            breaks = breaks,
            labels = labels
        )
    )
}

#' @importFrom dplyr mutate filter
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot geom_segment scale_x_log10 theme_bw scale_x_continuous
#' @importFrom ggplot2 coord_flip xlab ylab theme element_blank geom_vline geom_label
#' @importFrom ggplot2 geom_point aes geom_label scale_fill_manual scale_y_continuous labs
plotCustomizedTree <- function(factorMerger, stat = "model",
                               pos, levels = NULL,
                               showY = TRUE, alpha = 0.05,
                               showDiagnostics = TRUE) {

    segment <- getTreeSegmentDf(factorMerger, stat, pos)
    df <- segment$df
    pointsDf <- segment$pointsDf
    labelsDf <- segment$labelsDf

    g <- df %>% ggplot() +
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2)) +
        geom_point(data = pointsDf, aes(x = x1, y = y1), size = 0.75) +
        scale_y_continuous(limits = getLimits(labelsDf, showY),
                           position = "right",
                           breaks = labelsDf$y1,
                           labels = getLabels(labelsDf, factorMerger)) +
        ylab(getStatisticName(factorMerger))

    upperBreaks <- df$x1 %>% unique() %>% sort
    g <- g + xlab(renameStat(stat)) + treeTheme(getLabelsColors(labelsDf, levels))

    if (showDiagnostics) {
        if (stat == "pval") {
            intercept <- alpha
            label <- paste0("alpha = ", alpha)
            g <- g + scale_x_log10()
        }
        if (stat == "model") {
            gicMin <- mergingHistory(factorMerger, TRUE)[, c("model", "GIC")] %>%
                filter(GIC == min(GIC))
            intercept <- gicMin$model
            label <- paste0("min GIC")
        }

        y <- getLimits(labelsDf, showY)

        g <- g + geom_vline(xintercept = intercept, col = "mediumorchid3", linetype = "dotted") +
            geom_label(x = intercept, y = getLimits(labelsDf, showY)[1],
                      label = label, alpha = 0.5, col = "mediumorchid3",
                      angle = 90,
                      size = 3, fontface = "italic")
    }

    labBr <- getChisqBreaks(g$data, alpha)
    g <- g +
        scale_x_continuous(breaks = labBr$breaks, labels = labBr$labels)

    g <- g + labs(title = "Merging path plot",
                  subtitle = paste0("Optimal GIC partition: ",
                                    paste(getOptimalPartition(factorMerger), collapse = ":")))

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
    grid.arrange(.plotTree(factorMerger,
                           levels = lev,
                           showDiagnostics = FALSE,
                           simplify = TRUE), plot, ncol = 2)
}

#' @export
appendToTree.survPlot <- function(factorMerger, plot) {
    lev <- levels(plot$data$variable)
    grid.arrange(.plotTree(factorMerger,
                           levels = lev,
                           showDiagnostics = FALSE,
                           simplify = TRUE), plot, ncol = 2)
}

#' @importFrom grid grid.draw grid.newpage
#' @importFrom ggplot2 ggplotGrob theme labs ylab
#'
#' @export
appendToTree.GICPlot <- function(factorMerger, plot) {
    grid.newpage()
    grid.draw(rbind(ggplotGrob(plot + ylab("") +
                                   labs(title = "Merging path plot (with GIC profile)",
                                        subtitle = paste0("Optimal GIC partition: ",
                                                          paste(getOptimalPartition(factorMerger), collapse = ":")))),
                    ggplotGrob(plotTree(factorMerger) +
                                   theme(title = element_blank())), size = "last"))
}


#' @importFrom ggplot2 ggplot aes geom_line geom_text theme_minimal theme scale_color_manual labs
#' @export
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
        scale_color_manual(values = colorRamps::magenta2green(noLevels)) +
        ylab("") +
        labs(title = "Profile plot", subtitle = "Variable means ranks")
    class(g) <- append(class(g), "profilePlot")
    return(g)
}

scaleStat <- function(df) {
    df <- split(df, df$variable)
    sapply
}

#' @export
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_minimal scale_x_continuous theme
#' @importFrom ggplot2 coord_flip element_line element_blank scale_fill_distiller labs guides
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange
plotHeatmap <- function(factorMerger) {
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
              legend.position = "none") +
        scale_fill_distiller(palette = customPalette) +
        labs(title = "Heatmap", subtitle = "Group means by variables")
}

#' @export
#' @importFrom ggplot2 ggplot geom_boxplot aes coord_flip labs
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
        theme(axis.title = element_blank(), axis.text.y = element_blank()) +
        labs(title = "Boxplot", subtitle = "Summary statistic: mean")
}

#' @export
#' @importFrom ggplot2 ggplot geom_boxplot aes coord_flip labs geom_errorbar theme ylab position_dodge element_blank element_text
#' @importFrom dplyr group_by summarize left_join
plotMeansAndStds <- function(factorMerger) {
    factor <- factor(factorMerger$factor, levels = getFinalOrderVec(factorMerger))
    model <- lm(factorMerger$response ~ factor - 1)
    df <- data.frame(group = levels(factor))
    sumModel <- summary(model)
    df$mean <- sumModel$coefficients[, 1]
    df$left <- df$mean - sumModel$coefficients[, 2]
    df$right <- df$mean + sumModel$coefficients[, 2]
    df$group <- factor(df$group, levels = df$group)

    ggplot(data = df, aes(x = as.factor(group), y = mean, group = as.factor(group))) +
        geom_errorbar(aes(ymin = left, ymax = right),
                      color = "black",
                      width = .5,
                      position = position_dodge(.5)) + treeTheme(NULL) +
        geom_point() + coord_flip() +
        theme(axis.title.x = element_text(), axis.text.y = element_blank()) +
        labs(title = "Summary statistics", subtitle = "Means and standard deviations of coefficients' estimators") +
        ylab("")

}

#' @export
#' @importFrom ggplot2 ggplot geom_bar aes coord_flip scale_fill_manual theme theme element_blank scale_y_continuous labs
#' @importFrom dplyr group_by summarize left_join
plotProportion <- function(factorMerger) {
    levels <- getFinalOrderVec(factorMerger)
    responseLevels <- factorMerger$response %>% as.factor() %>% levels()
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    data <- data.frame(x = factorMerger$factor, y = factorMerger$response)
    data %>% ggplot() + geom_bar(aes(x = x, fill = as.factor(y)), position = "fill") +
        scale_y_continuous(label = scales::percent, name = "") +
        coord_flip() + treeTheme(NULL) +
        scale_fill_manual(values = customPaletteValues[c(2, length(customPaletteValues) - 1)]) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
        labs(title = "Group success proportion",
             subtitle = paste0("Success: ", responseLevels[2],
                               " (green), failure: ", responseLevels[1], " (violet)"))
}

#' @importFrom ggplot2 labs
#'
#' @export
plotSurvival <- function(factorMerger) {
    model <- calculateModel(factorMerger, factorMerger$factor)
    g <- survminer::ggcoxadjustedcurves(model, data = data.frame(factorMerger$factor),
                        individual.curves = TRUE,
                        theme = treeTheme(NULL),
                        variable = factorMerger$factor,
                        palette = colorRamps::magenta2green(length(levels(factorMerger$factor))),
                        curve.size = 1) +
        treeTheme(NULL) + labs(title = "Survival plot", subtitle = "Adjusted survival curves for coxph model")
    class(g) <- append(class(g), "survPlot")
    return(g)
}

#' @importFrom ggplot2 ggplot geom_line aes theme element_blank scale_y_continuous labs geom_point geom_ribbon
#'
#' @export
plotGIC <- function(factorMerger) {
    mH <- mergingHistory(factorMerger, T)
    minGIC <- min(mH$GIC)
    minModel <- mH$model[which.min(mH$GIC)]
    g <- mH %>% ggplot(aes(x = model, y = GIC)) + geom_line(col = customPaletteValues[1], size = 1) +
        geom_point(x = minModel, y = minGIC, col = customPaletteValues[1], size = 2.5) +
        geom_ribbon(aes(x = model, ymin = minGIC, ymax = GIC), fill = customPaletteValues[1], alpha = 0.2) +
        treeTheme(NULL) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(),
              panel.grid.major = element_blank()) +
        scale_y_continuous(position = "right")
    class(g) <- append(class(g), "GICPlot")
    return(g)
}
