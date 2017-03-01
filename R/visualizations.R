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

