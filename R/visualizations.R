#' Generated sample visualization
#'
#' @param x generatesSample object
#'
#' @rdname plot.generatedSample
#' @importFrom ggplot2 ggplot geom_boxplot aes ylab xlab stat_summary labs
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize left_join
#'
#' @export

plot.generatedSample <- function(data) {
    colnames(data) <- c("y", "x")
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

#' Plot heatmap
#'
#' @param x generatesSample object
#'
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_bw scale_x_continuous theme
#' @importFrom ggplot2 coord_flip element_line element_blank scale_fill_brewer
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange
#'

plotHeatmap <- function(factorMerger) {
    levels <- (plot(factorMerger)$data %>% filter(!is.na(label)) %>%
        arrange(y))$label
    factorMerger$factor <- factor(factorMerger$factor, levels = levels)
    data.frame(cbind(response = factorMerger$response), factor = factorMerger$factor) %>%
        melt(id.vars = "factor") %>% ggplot() +
        geom_tile(aes(x = factor, y = variable, fill = value)) +
        coord_flip() + theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank()) +
        scale_fill_distiller(palette = "PRGn")
}

rootShifter <- function(x, shift) {
    lab <- round(x + as.numeric(shift), 2)
}

#' @importFrom ggtree read.tree ggtree gheatmap theme_tree2 geom_tiplab scale_x_ggtree
plot.multivariateFactorMerger <- function(factorMerger, stat = "model") {
    tr <- getTreeWithEdgesLength(factorMerger, stat)
    df <- data.frame(factorMerger$factor, factorMerger$response)
    df <- aggregate(df[, -1], list(df[, 1]), mean)
    rownames(df) <- df[, 1]
    df <- data.frame(df[, -1])
    shift <- factorMerger$mergingList[[length(factorMerger$mergingList)]]$modelStat[stat]
    (ggtree(read.tree(text = tr)) + geom_tiplab(align = TRUE) +
            theme_tree2() + xlab(stat) +
            scale_x_continuous(label = function(x) { rootShifter(x, shift) }))
}

plot.factorMerger <- function(fm) {
    plot.multivariateFactorMerger(fm)
}

plotTreeAndHeatmap <- function(factorMerger) {
    UseMethod("plotTreeAndHeatmap", factorMerger)
}

plotTreeAndHeatmap.factorMerger <- function(factorMerger) {
    plotTreeAndHeatmap.multivariateFactorMerger(factorMerger)
}

#' @importFrom gridExtra grid.arrange
plotTreeAndHeatmap.multivariateFactorMerger <- function(factorMerger) {
    g1 <- plot(factorMerger)
    g2 <- plotHeatmap(factorMerger) + xlab("")
    grid.arrange(g1, g2, ncol = 2)
}

