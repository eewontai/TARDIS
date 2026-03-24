#' @title Run T.A.R.D.I.S.
#' @description
#' Launches the GUI of T.A.R.D.I.S. that allows to input all parameters in an
#' intuitive way to perform targeted peak detection.
#'
#' @import shiny
#' @importFrom magrittr %>%

#' @export
runTardis <- function() {
    appDir <- system.file("tardis_app", package = "TARDIS")   # searches for tardis_app folder in TARDIS directory
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `TARDIS`.",
            call. = FALSE
        )
    }
    shiny::runApp(appDir, display.mode = "normal")   # run tardis_app
}
