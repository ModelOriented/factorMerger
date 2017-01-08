#' Set increasing order
#'
#' Given a numeric vector and a factor vector changes levels order
#' in the factor. Means of the numeric variable in groups
#' appointed by the new order of the factor are increasing.
#'
#' @param numericVec numeric vector
#' @param factorVec factor vector (same length as numbericVec).
#'
#' @rdname setIncreasingOrder
#'

setIncreasingOrder <- function(numericVec, factorVec) {
    if (length(numericVec) != length(factorVec)) {
        stop("Vectors' lengths do not match.")
    }
    stopifnot(is.vector(numericVec))
    stopifnot(is.vector(numericVec))
    factorVec <- as.factor(as.character(factorVec))
    data <- data.table::data.table(y = numericVec, c = factorVec)
    newOrder <- aggregate(numericVec ~ c, mean, data = data) %>%
        dplyr::arrange(numericVec)
    factor(factorVec, levels = as.character(newOrder[, 1]))
}
