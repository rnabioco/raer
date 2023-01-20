
#' Generate a samtools style index for regions
#'
#' @param gr GRanges object
#'
#' @return Return an object of class RegionIndex which contains a pointer
#' to the region index.
#'
#' @importFrom methods as is new
#'
#' @rdname index_regions
#' @examples
#' library(rtracklayer)
#' bed_fn <- system.file("extdata", "regions.bed", package = "raer")
#' gr <- import(bed_fn)
#' gr <- tile(gr, width = 1) |> unlist()
#' strand(gr) <- sample(c("+", "-"), length(gr), replace = TRUE)
#' gr$ref <- sample(c("A", "T", "C", "G"), length(gr), replace = TRUE)
#' gr$alt <- sample(c("A", "T", "C", "G"), length(gr), replace = TRUE)
#' indexRegions(gr)
#' @export
setMethod(
  indexRegions, "GRanges",
  function(gr) {
    stopifnot(all(width(gr) == 1))

    nr <- length(gr);
    if(nr == 0) stop("No entries in GRanges")

    gr_nms <- c("ref", "alt")
    if(!all(gr_nms %in% names(mcols(gr)))){
      stop("GRanges must have a ref and alt columns")
    }

    if(any(strand(gr) == "*")){
      warning("missing strand not found in input, coercing strand to '+'")
      strand(gr) <- "+"
    }
    gr$idx <- seq(0, nr - 1) # is zero-based index

    obj <- .RegionIndex$new(
      gr = gr,
      open = FALSE
    )

    l <- list(as.character(seqnames(gr)),
              as.integer(start(gr)),
              as.integer(strand(gr)),
              as.character(mcols(gr)$ref),
              as.character(mcols(gr)$alt),
              as.integer(mcols(gr)$idx))

    tryCatch(
      {
        obj$.extptr <- .Call(".regidx_build", l)
        obj$open <- TRUE
      },
      error = function(err) {
        stop(conditionMessage(err), "\n  issue building index")
      }
    )

    return(obj)
  }
)

#' Close a RegionIndex index connection
#'
#' @param con RegionIndex class
#' @param ... present for consistency
#'
#' @rdname index_regions
#' @examples
#' library(rtracklayer)
#' bed_fn <- system.file("extdata", "regions.bed", package = "raer")
#' bed <- indexBed(bed_fn)
#' close(bed)
#' @export
close.RegionIndex <-
  function(con, ...) {
    stopifnot(!is.null(con$.extptr))

    if (!con$open) {
      return(con$open)
    }

    tryCatch(
      {
        con$.extptr <- .Call(".regidx_close", con$.extptr)
        con$open <- FALSE
      },
      error = function(err) {
        stop(conditionMessage(err), "\n  ", con)
      }
    )

    invisible(con)
  }

#' Print a RegionIndex index connection
#'
#' @param con RegionIndex class
#' @param ... present for consistency
#'
#' @rdname index_regions
#' @examples
#' library(rtracklayer)
#' bed_fn <- system.file("extdata", "regions.bed", package = "raer")
#' gr <- import(bed_fn)
#' i <- indexRegions(gr)
#' print(i)
#' @export
print.RegionIndex <-
  function(con, ...) {
    stopifnot(!is.null(con$.extptr))

    if (!con$open) {
      return(con$open)
    }

    tryCatch(
      {
        .Call(".print_regidx", con$.extptr)
      },
      error = function(err) {
        stop(conditionMessage(err), "\n  ", con)
      }
    )

    invisible(con)
  }
