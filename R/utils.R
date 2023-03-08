#' Provide working directory for raer example files.
#'
#' @param path path to file
#'
#' @examples
#' raer_example("human.fasta")
#'
#' @return Character vector will path to an internal package file.
#' @export
raer_example <- function(path) {
    system.file("extdata", path, package = "raer", mustWork = TRUE)
}

# transpose a list
# https://stackoverflow.com/questions/30164803/fastest-way-to-transpose-a-list-in-r-rcpp
t_lst <- function(x) {
    split(unlist(x), sequence(lengths(x)))
}

check_tag <- function(x) {
    if (length(x) != 1 && nchar(x) != 2) {
        stop("supplied tag must by nchar of 2: ", x)
    }
}

chunk_vec <- function(x, n) {
    if (n == 1) {
        res <- list(`1` = x)
        return(res)
    }
    split(x, cut(seq_along(x), n, labels = FALSE))
}


# flatten top list, keeping names from inner list
unlist_w_names <- function(x) {
    nms <- unlist(lapply(x, names))
    res <- unlist(x, use.names = FALSE)
    names(res) <- nms
    res
}
