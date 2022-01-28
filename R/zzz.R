.onUnload <- function (libpath) {
  library.dynam.unload("eigvaldispRcpp", libpath)
}
