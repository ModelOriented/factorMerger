#' Generate sample
#'
#' Function generateSample() produces a random sample consisting
#' of k groups drawn from the same distribution with different
#' parameters.
#'
#' @param N sample size
#' @param k number of groups
#' @param distr type of distribution from \code{c("norm", "exp", "beta")}
#'
#' @examples
#' generateSample(100, 2)
#' generateSample(100, 2, "exp")
#'
#' @return
#' \code{data.table} with two columns: numeric variable
#' and factor variable.
#'
#' @rdname generateSample
#'
#' @importFrom data.table data.table
#'
#' @export
#'
#'
generateSample <- function(N, k, distr = "norm") {
    numericVec <- switch(distr,
        "norm" = rnorm(N),
        "exp" = rexp(N, 1),
        "beta" = rbeta(N, 1, 1),
        stop("Unknown distribution."))

    factorVec <- as.factor(sample(LETTERS[1:k],
                              size = N,
                              replace = TRUE))

    for (i in 1:k) {
        randomShift <- sample(seq(0, 1, 0.1), size = 1)
        let <- LETTERS[i]
        numericVec[factorVec == let] <- numericVec[factorVec == let] + randomShift
    }

    generatedSample <-
        data.table(numericVec = numericVec,
                   factorVec = setIncreasingOrder(numericVec, factorVec))
    class(generatedSample) <- append("generatedSample", class(generatedSample))
    generatedSample
}

generateMultivariateSample <- function(N, k, d = 2, distr = "norm") {
    res <- matrix(, nrow = N, ncol = d)
    tmp <- generateSample(N, k, distr)
    res[, 1] <- tmp$numericVec
    for (j in 2:d) {
        for (i in 1:k) {
            randomShift <- sample(seq(0, 1, 0.1), size = 1)
            normal <- rnorm(N)
            normal[tmp$factorVec == LETTERS[i]] <-
                normal[tmp$factorVec == LETTERS[i]] + randomShift
        }
        res[, j] <- normal
    }
    return(data.frame(res = res, factor = tmp$factorVec))
}
