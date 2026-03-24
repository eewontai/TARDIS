#' @title  Check if any samples are missing spectra
#' @description
#' Display warning if they do.
#' Currently warns user if a sample is detected that has less than 50% of the mean
#' of spectra in all samples.
#'
#' @param spectra `Spectra` object
#' Spectra are one-dimensional objects storing spectra, even from different files or samples, in a single list.
#' The data origin of each spectrum can be extracted with the dataOrigin() function. (.mzML file paths)
#' mzML: open source file format for raw mass spec data
#'
#' @importFrom ProtGenerics dataOrigin
#' @importFrom Spectra isEmpty
#' @export
#' @author Pablo Vangeenderhuysen
checkScans <- function(spectra) {
    scans_per_sample <- table(dataOrigin(spectra))  # table() makes a frequency table
    mean <- (mean(scans_per_sample))
    bad_runs <- which(scans_per_sample < 0.5 * mean)
    if (isEmpty(bad_runs) == FALSE) {
        names <- basename(names(bad_runs))  # basename(): only the file names, not full paths
        warning(paste("File", names, "contains less than 50% of the mean of scans
                   in the samples."))
    }
}
