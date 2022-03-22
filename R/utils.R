
is_null_extptr <- function(pointer) {
  stopifnot(is(pointer, "externalptr"))
  .Call(".isnull", pointer)
}
