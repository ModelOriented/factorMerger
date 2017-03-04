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

#' @export
#'
plotTree <- function(factorMerger, stat) {
    UseMethod("plotTree", factorMerger)
}

#' Plot
#' @importFrom ggtree read.tree ggtree gheatmap theme_tree2 geom_tiplab scale_x_ggtree
#' @export
#'
plotTree.multivariateFactorMerger <- function(factorMerger, stat = "model") {
    tr <- getTreeWithEdgesLength(factorMerger, stat)
    df <- data.frame(factorMerger$factor, factorMerger$response)
    df <- aggregate(df[, -1], list(df[, 1]), mean)
    dfRowNames <- df[, 1]
    df <- data.frame(df[, -1])
    rownames(df) <- dfRowNames
    shift <- factorMerger$mergingList[[length(factorMerger$mergingList)]]$modelStat[stat]
    tree <- read.tree(text = tr)
    brLabs <- breaksAndLabelsCalc(tree, shift, 5)
    off <- as.numeric(abs(shift)) / 1000

    ((ggtree(tree) + geom_tiplab(align = TRUE) +
            theme_tree2()) %>%
            gheatmap(df, offset = off, width = 0.5, colnames = FALSE)) %>%
        scale_x_ggtree(breaks = brLabs$breaks, labels = brLabs$labels) +
        scale_fill_distiller(palette = "PRGn") +
        xlab(paste0(stat, paste0(rep(" ", 50), collapse = ""))) # seriously - this paste ... needs rewriting
}

#' @export
#'
plotTree.default <- function(factorMerger, stat) {
    plotTree.multivariateFactorMerger(factorMerger, stat)
}

#' @export
#'
plot.factorMerger <- function(factorMerger) {
    plotTree(factorMerger, "model")
}

#'
#'
plotMyTree <- function(factorMerger, stat = "model") {
    means <- means(factorMerger)
    if (is.null(means)) {
        plotSimpleTree(factorMerger, stat)
    }
    else {
        plotCustomizedTree(factorMerger, stat, means)
    }
}

#' @importFrom dplyr mutate filter
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot geom_segment theme_bw coord_flip xlab ylab theme element_blank geom_point aes
plotCustomizedTree <- function(factorMerger, stat = "model", pos = NULL) {
    factor <- factorMerger$factor
    noGroups <- length(levels(factor))
    df <- pos[1:noGroups, ] %>%  data.frame
    colnames(df) <- "x1"
    df$x2 <- df$x1
    df$label <- rownames(pos)[1:noGroups]
    df$y1 <- factorMerger$mergingList$`1`$modelStats[, stat]
    df$y2 <- NA
    pointsDf <- df
    merging <- mergingHistory(factorMerger)
    for (step in 1:nrow(merging)) {
        statVal <- factorMerger$mergingList[[step + 1]]$modelStats[, stat]
        pair <- merging[step, ]
        whichDf <- which(df$label %in% pair)
        df[whichDf, "y2"] <- statVal
        pairLabel <- paste(pair, collapse = "")
        pairMean <- pos[rownames(pos) == pairLabel,]
        crosswise <- c(df$x1[whichDf], "", statVal, statVal)
        newVertex <- c(pairMean, pairMean, pairLabel, statVal, NA)
        df <- rbind(df, crosswise)
        df <- rbind(df, newVertex)
    }
    df %>% subset(select = -label) %>%
        filter(!is.na(y2)) %>% apply(2, as.numeric) %>% round(2) %>%
        as.data.frame() %>% ggplot() +
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2)) +
        theme_bw() + coord_flip() + ylab(stat) + xlab("") +
        theme(axis.text.y = element_blank()) +
        geom_point(data = pointsDf, aes(x = x1, y = y1)) +
        geom_label(data = pointsDf,
                   aes(x = x1, y = y1, label = pointsDf$label))
    # Jeżeli będziemy korzystać z ggrepel, to można zamienić na geom_label
}

plotSimpleTree <- function(factorMerger, stat = "model") {
    pos <- getFinalOrder(factorMerger) %>% data.frame()
    merging <- mergingHistory(factorMerger)
    noStep <- nrow(merging)

    for (step in 1:noSteps) {
        pos <- rbind(pos, mean(pos[rownames(pos) %in% merging[step, ],]))
        rownames(pos)[nrow(pos)] <- paste(merging[step, ], collapse = "")
    }
    return(plotCustomizedTree(factorMerger, stat, pos))
}

#' @importFrom ggplot2 ggplot aes geom_line geom_text theme_bw theme
# TODO: dodać wybór między rank a średnią
profilePlot <- function(factorMerger) {
    stat <- calculateMeansByFactor(factorMerger$response,
                                   factorMerger$factor)

    stat$variable <- factor(stat$variable, levels = unique(stat$variable))
    stat$level <- factor(stat$level,
                         levels = (stat %>%
                                       filter(variable == levels(stat$variable) %>%
                                                  tail(1)) %>% arrange(rank))$level
    )
    noLevels <- length(levels(stat$level))
    stat$rank <- factor(stat$rank, levels = 1:noLevels)
    stat %>% ggplot(aes(x = variable, y = rank, col = level, group = level, label = level)) +
        geom_line() +
        geom_text(data = subset(stat,
                                variable == levels(stat$variable) %>% tail(1)),
                  aes(x = variable),
                  size = 3.5, hjust = 0.8,  nudge_x = 0.1) +
        theme_bw() + theme(legend.position = "none")


}

