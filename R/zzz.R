.onLoad <- function(libname, pkgname) {
    message("Loading package: ", pkgname)
    message("From library: ", libname)
    message("Available functions: ", paste(ls(paste0("package:", pkgname)), collapse = ", "))
}
