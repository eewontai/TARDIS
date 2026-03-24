#' filter Spectra to single peak in single sample
#'
#' @param spectra a `Spectra` object
#' @param dataOrigin file that contains the spectra
#' @param rt_range `numeric(2)` rt range of compound
#' @param mz_range `numeric(2)` mz range of compound
#' @returns returns filtered `Spectra` object
filterSingle <- function(spectra, dataOrigin, rt_range, mz_range) {
    spectra <- spectra |>
        filterDataOrigin(dataOrigin) |>
        filterRt(rt_range) |>
        filterMzRange(mz_range)
    spectra
}

#' Sum Intensities of Spectra
#' @noRd
.sum_intensities <- function(x, ...) { # ... allows the function to accept any number of additional, arbitrary arguments.
    if (nrow(x)) {  # if data exist
        cbind(
            mz = NA_real_,  # mz does not matter anymore
            intensity = sum(x[, "intensity"], na.rm = TRUE)   # sum the intensity column
        )
    } else {
        cbind(mz = NA_real_, intensity = NA_real_)  # null
    }
}

#' @title Function to extract EIC from Spectra object
#'
#' @param spectra a `Spectra` object.
#'
#' @return two column matrix with rt and int
#'
#' @author Pablo Vangeenderhuysen
#'
#'
#' @export
extract_eic <- function(spectra) {
    sfs_agg <-
        addProcessing(spectra, .sum_intensities)  # add processing task, not run it yet
    # addProcessing extracts only the peaks matrix of the spectra, containg 2 columns; mz and int.
    # in 3d chromatogram, collapse the mz dimension, giving the eic (rt vs int).
    eic <-
        cbind(  # run it
            rtime(sfs_agg),  # retention time
            unlist(intensity(sfs_agg), use.names = FALSE)  # intensity
        ) # 2d matrix
    rownames(eic) <- NULL
    colnames(eic) <- c("rt", "int")
    eic
}
