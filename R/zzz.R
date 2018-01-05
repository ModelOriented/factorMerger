.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste0("Welcome to factorMerger ",
               packageVersion("factorMerger"), "!\n"))
}
