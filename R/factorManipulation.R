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
#' @export
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

#' Merge factor
#'
mergeFactor <- function(factor, groupA, groupB, groupAB) {
    # factor[factor == groupA] <- groupB
    whichLevels <- which(levels(factor) %in% c(groupA, groupB))
    # factor <- factor(factor, labels = groups)
    levels(factor)[whichLevels] <- groupAB
    factor
}

#' Calculate means by factor
#'
#' @rdname calculateMeans
#' @importFrom dplyr group_by summarize arrange
#'
calculateMeans <- function(numericVec, factorVec) {
    if (!is.null(dim(numericVec))) {
        return(NA)
    }
    data.frame(num = numericVec, level = factorVec) %>%
        group_by(level) %>% summarize(mean = mean(num)) %>%
        arrange(mean)
}

#' Filter groups - ...
#'
filterGroups <- function(response, factor, groupA, groupB) {
    response <- as.matrix(response)
    xA <- response[factor == groupA, ]
    xB <- response[factor == groupB, ]
    return(list(xA, xB))
}

