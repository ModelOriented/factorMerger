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
#' \code{list} with two fields: numeric variable \code{response}
#' and factor variable \code{factor}.
#'
#' @rdname generateSample
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

    generatedSample <- list(
        factor = setIncreasingOrder(numericVec, factorVec),
        response = numericVec)

    class(generatedSample) <- append("generatedSample", class(generatedSample))
    return(generatedSample)
}

#' Generate multivariate sample
#'
#' @export
#'
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
    return(list(factor = tmp$factorVec, response = res))
}
