.onAttach <- function (lib, pkg) {
  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
  ver <- as.character(ver)
  packageStartupMessage("\nvarianceExplainedPack ",
                        ver,
                        " loaded. \n",
                   #   "\nPlease cite as:\n   \n",
                        domain = NULL,
                        appendLF = TRUE)
}

.onLoad <- function(...) {
}

.onUnload <- function (libpath) {
}
