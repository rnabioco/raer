
#' Generate a samtools style index for regions
#'
#' @param gr GRanges object
#'
#' @return Return an object of class RegionIndex which contains a pointer
#' to the region index.
#'
#' @importFrom methods as is new
#'
#' @rdname index_bed
#' @examples
#' bed_fn <- system.file("extdata", "regions.bed", package = "raer")
#' gr <- import(bed_fn)
#' indexRegions(gr)
#' @export
setMethod(
  indexRegions, "GRanges",
  function(gr) {
    stopifnot(all(width(gr) == 1))

    obj <- .RegionIndex$new(
      gr = gr,
      open = FALSE
    )
    l <- list(as.character(seqnames(gr)),
              as.integer(start(gr)),
              as.integer(seq_along(gr)))

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
#' @rdname index_bed
#' @examples
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
#' @rdname index_bed
#' @examples
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
