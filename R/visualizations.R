#' Generated sample visualization
#'
#' @param x generatesSample object
#'
#' @rdname plot.generatedSample
#' @importFrom ggplot2 ggplot geom_boxplot aes ylab xlab stat_summary labs
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize left_join
#' @importFrom data.table data.table
#' @export

plot.generatedSample <- function(data) {
    data <- data.frame(data)
    colnames(data) <- c("x", "y")
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

        xlab("Groups") + ylab("Response") +
        labs(title = "Generated sample - boxplot",
             subtitle = "Groups are sorted by means. Displayed summary statistic - mean in a group.")
}

#' @importFrom ape node.depth.edgelength
breaksAndLabelsCalc <- function(tr, shift, gridLength) {
    trHeight <- max(node.depth.edgelength(tr))
    br <- seq(0, trHeight, trHeight / gridLength)
    return(
        list(
            breaks = br,
            labels = as.character(round(br + as.numeric(shift), 2))
        )
    )
}

#' @importFrom ggplot2 theme_classic theme element_line element_blank theme_minimal
treeTheme <- function(showY) {

    if (showY) {
        return(theme_minimal() +
                   theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         legend.position = "none"))
    }
    else {
         return(theme_minimal() +
                     theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
                           panel.grid.major.y = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           panel.grid.minor.x = element_blank(),
                           legend.position = "none",
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank()))
    }
}


#' @export
#'
plotTree <- function(factorMerger, stat, color = FALSE) {
    UseMethod("plotTree", factorMerger)
}

#' @export
#' @importFrom magrittr %>%
plotTree.factorMerger <- function(factorMerger, stat = "model", color = FALSE) {
    stopifnot(stat %in% c("model", "pval"))
    means <- means(factorMerger)
    if (is.null(means)) {
        plotSimpleTree(factorMerger, stat, color)
    }
    else {
        plotCustomizedTree(factorMerger, stat, means, color)
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
        return(c(min(df$y1), max(df$y1)))
    } else {
        shift <- 0.6 / (nrow(df) - 1)
        return(c(min(df$y1) - shift, max(df$y1) + shift))
    }
}

#' @importFrom dplyr mutate filter
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot geom_segment scale_x_log10 theme_bw coord_flip xlab ylab theme element_blank geom_point aes geom_label scale_y_continuous
plotCustomizedTree <- function(factorMerger, stat = "model", pos, color = FALSE, showY = TRUE) {
    factor <- factorMerger$factor
    noGroups <- length(levels(factor))
    df <- pos[1:noGroups, ] %>%  data.frame
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
    df <- df %>% subset(select = -label) %>% filter(!is.na(x2)) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    if (stat == "pval") {
        df$x1 <- log10(df$x1)
        df$x2 <- log10(df$x2)
        pointsDf$x1 <- log10(pointsDf$x1)
    }

    g <- df %>% ggplot() +
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2)) +
        geom_point(data = pointsDf, aes(x = x1, y = y1)) +
        scale_y_continuous(limits = getLimits(pointsDf, showY), position = "right") +
        ylab("Group mean")

    stat <- renameStat(stat)
    g <- g + xlab(stat) + treeTheme(showY)

    if (color) {
         return(g + geom_label(data = pointsDf,
                                aes(x = x1, y = y1, label = pointsDf$label,
                                    fill = factor(pointsDf$label,
                                                  levels = getFinalOrderVec(factorMerger)))))
    } else {
        return(g + geom_label(data = pointsDf,
                              aes(x = x1, y = y1, label = pointsDf$label)))
    }

    # Jeżeli będziemy korzystać z ggrepel, to można zamienić na geom_label
}


plotSimpleTree <- function(factorMerger, stat = "model", color = FALSE) {
    pos <- getFinalOrder(factorMerger) %>% data.frame()
    merging <- mergingHistory(factorMerger)
    noStep <- nrow(merging)

    for (step in 1:noStep) {
        pos <- rbind(pos, mean(pos[rownames(pos) %in% merging[step, ],]))
        rownames(pos)[nrow(pos)] <- paste(merging[step, ], collapse = "")
    }
    return(plotCustomizedTree(factorMerger, stat, pos, color, showY = FALSE))
}


#' @export
bindPlots <- function(p1, p2) {
    grid.arrange(p1, p2, ncol = 2)
}

#' @importFrom ggplot2 ggplot aes geom_line geom_text theme_bw theme
#' @export
# TODO: dodać wybór między rank a średnią
plotProfile <- function(factorMerger) {
    stat <- calculateMeansByFactor(factorMerger$response,
                                   factorMerger$factor)

    stat$variable <- factor(stat$variable, levels = unique(stat$variable))
    stat$level <- factor(stat$level, levels = getFinalOrderVec(factorMerger))

    noLevels <- length(levels(stat$level))
    stat$rank <- factor(stat$rank, levels = 1:noLevels)
    stat %>% ggplot(aes(x = variable, y = rank, col = level, group = level, label = level)) +
        geom_line() +
        geom_text(data = subset(stat,
                                variable == levels(stat$variable) %>% tail(1)),
                  aes(x = variable),
                  size = 3.5, hjust = 0.8,  nudge_x = 0.1) +
        theme_minimal() + theme(legend.position = "none")
}

# TODO: Zrobić dendrogram na zmiennych, żeby zmienne podobn

#' @export
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_minimal scale_x_continuous theme
#' @importFrom ggplot2 coord_flip element_line element_blank scale_fill_distiller
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange
plotHeatmap <- function(factorMerger) {
    levels <- getFinalOrderVec(factorMerger)
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    data.frame(cbind(response = factorMerger$response), factor = factorMerger$factor) %>%
        melt(id.vars = "factor") %>% ggplot() +
        geom_tile(aes(x = factor, y = variable, fill = value)) +
        coord_flip() + theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        scale_fill_distiller(palette = "PRGn")
}


