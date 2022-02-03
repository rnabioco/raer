.onUnload <- function (libpath) {
  library.dynam.unload("ullr", libpath)
}