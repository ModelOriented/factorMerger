getLetters <- function(k) {
    reps <- rep(LETTERS, round(k / length(LETTERS) + 1))[1:k]
    prefix <- rep(c("", LETTERS), each = length(LETTERS))[1:k]
    return(paste0(prefix, reps))
}

#' Generate sample
#'
#' Produces a random sample of k groups drawn from the same distribution with different
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
generateSample <- function(N, k, distr = "gaussian") {
    numericVec <- switch(distr,
        "gaussian" = rnorm(N),
        "exp" = rexp(N, 1),
        "beta" = rbeta(N, 1, 1),
        "binomial" = rep(0, N),
        stop("Unknown distribution."))

    kLetters <- getLetters(k)
    factorVec <- as.factor(sample(kLetters,
                              size = N,
                              replace = TRUE))

    for (i in 1:k) {
        let <- kLetters[i]
        if (distr == "binomial") {
            numericVec[factorVec == let] <-
                rbinom(length(numericVec[factorVec == let]), 1, runif(1))
        } else {
            randomShift <- sample(seq(0, 1, 0.1), size = 1)
            numericVec[factorVec == let] <-
                numericVec[factorVec == let] + randomShift * 0.1
            randomShift <- sample(seq(0, 1, 0.1), size = 1)
            # numericVec[factorVec == let] <-
            #     numericVec[factorVec == let] * i
        }
    }

    generatedSample <- list(
        factor = setIncreasingOrder(numericVec, factorVec),
        response = numericVec)

    class(generatedSample) <- append("generatedSample", class(generatedSample))
    return(generatedSample)
}

#' Generate multivariate normal sample
#'
#' Produces a random sample of k groups and d dimensions drawn from the
#' normal distribution with different
#' parameters.
#'
#' @param N Sample size.
#' @param k Number of groups.
#' @param d Number of dimensions.
#'
#' @examples
#' generateMultivariateSample(N = 100, k = 10, d = 5)
#'
#' @return
#' \code{list} with two fields: matrix \code{response}
#' and factor variable \code{factor}.
#'
#' @export
#'
generateMultivariateSample <- function(N, k, d = 2) {
    tmp <- generateSample(N, k, "gaussian")
    if (d > 1) {
        res <- matrix(, nrow = N, ncol = d)
        res[, 1] <- tmp$response
        for (j in 2:d) {
            for (i in 1:k) {
                randomShift <- sample(seq(0, 1, 0.1), size = 1)
                normal <- rnorm(N)
                normal[tmp$factor == LETTERS[i]] <-
                    normal[tmp$factor == LETTERS[i]] + randomShift
            }
            res[, j] <- normal
        }
        return(list(factor = tmp$factor, response = res))
    } else {
        return(tmp)
    }

}
