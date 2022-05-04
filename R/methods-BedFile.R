
#' Generate a samtools style index for a BedFile
#'
#' @param file path to bed file
#'
#' @returns Return an object of class BedFile which contains a pointer
#' to the bed file index.
#'
#' @importFrom methods as is new
#'
#' @rdname index_bed
#' @examples
#' bed_fn <- system.file("extdata", "regions.bed", package = "raer")
#' indexBed(bed_fn)
#' @export
setMethod(indexBed, "character",
          function(file)
          {
            stopifnot(file.exists(file))
            file <- normalizePath(path.expand(file))

            obj <- .BedFile$new(path = file,
                                open = FALSE)

            tryCatch({
              obj$.extptr <- .Call(".bedfile_open", file)
              obj$open <- TRUE
            }, error=function(err) {
              stop(conditionMessage(err), "\n  file: ", file)
            })

            return(obj)
})

#' Close a BedFile index connection
#'
#' @param con BedFile class
#' @param ... present for consistency
#' bed_fn <- system.file("extdata", "regions.bed", package = "raer")
#' bed <- indexBed(bed_fn)
#' close(bed)
#' @export
close.BedFile <-
  function(con, ...)
  {
    stopifnot(!is.null(con$.extptr))

    if(!con$open){
      return(con$open)
    }

    tryCatch({
      con$.extptr <- .Call(".bedfile_close", con$.extptr)
      con$open <- FALSE
    }, error=function(err) {
      stop(conditionMessage(err), "\n  file: ", con)
    })

    invisible(con)
  }

