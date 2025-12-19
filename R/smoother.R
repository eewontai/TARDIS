#' Simulate smoothing and integration for single chromatographic peak in single file
#'
#' @param file_path `character(1)` path to .mzML or .mzXML file
#' @param rt `numeric(1)` search retention time for target in minutes
#' @param mz `numeric(1)` search m/z of target compound
#' @param ppm `numeric(1)` allowed ppm deviance
#' @param rtdev `numeric(1)` RT window in seconds
#' @param smoother `character(1)` choose smoothing algortim, can be: "sgolay"
#'
#' @import MsExperiment
#' @importFrom Spectra MsBackendMzR
#' @importFrom Spectra filterMzRange
#' @importFrom Spectra filterEmptySpectra
#' @importFrom Spectra filterDataOrigin
#' @importFrom Spectra filterRt
#' @importFrom Spectra dataOrigin
#' @importFrom Spectra addProcessing
#' @importFrom signal sgolayfilt
#' @importFrom BiocParallel SnowParam
#' @importFrom xcms intensity
#' @importFrom pracma trapz
#' @importFrom manipulate manipulate
#' @importFrom manipulate slider
#' @importFrom manipulate picker
#'
#' @returns interactive plot of smoothing and integration
#'
#' @export


smoothingSim <- function(file_path, rt, mz, ppm, rtdev, smoother) {
    data <- MsExperiment::readMsExperiment(
        spectraFiles = file_path,
        backend = MsBackendMzR(),
        BPPARAM = SnowParam(workers = 1)
    )
    dbData <- data.frame(rt * 60, mz)
    colnames(dbData) <- c("tr", "m/z")
    range <- createRanges(data, dbData, ppm, rtdev)
    mzRange <- range[[1]]
    rtRange <- range[[2]]
    spectra <- data@spectra
    spectra <-
        filterRt(spectra, rtRange)
    spectra <-
        filterMzRange(spectra, mzRange)
    sfs_agg <-
        addProcessing(spectra, .sum_intensities)
    eic <-
        cbind(
            rtime(sfs_agg),
            unlist(intensity(sfs_agg), use.names = FALSE)
        )
    rt <- eic[, 1]
    int <- eic[, 2]
    manipulate(
        {
            if (smoother == "sgolay") {
                # smooth intensity
                int[which(is.na(int))] <- 0
                smoothed <- sgolayfilt(int, n = fl, p = p)
                smoothed[smoothed < 0] <- 0
            }
            # borders determined from smoothed data
            borders <- find_peak_points(rt, smoothed, dbData$tr)
            # choose data to integrate based on toggle
            if (smooth_on == "Yes") {
                y_use <- smoothed
                label <- "Smoothed"
                color <- "blue"
            } else {
                y_use <- int
                label <- "Raw"
                color <- "darkred"
            }
            # integrate selected data within smoothed borders
            x <- rt[borders["left"]:borders["right"]]
            y <- y_use[borders["left"]:borders["right"]]
            auc <- trapz(x, y)
            # plot raw + smoothed
            plot(
                rt, int,
                type = "b", col = "grey", pch = 16,
                ylim = c(0, max(int) * 1.1),
                xlab = "Retention time", ylab = "Intensity",
                main = paste0("Integration on ", label, " peak (AUC = ", round(auc, 2), ")")
            )
            lines(rt, smoothed, col = "blue", lwd = 2, lty = 2)
            lines(rt, y_use, col = color, lwd = 2)
            # show integrated region
            polygon(
                c(x[1], x, tail(x, 1)),
                c(0, y, 0),
                col = adjustcolor(color, alpha.f = 0.3),
                border = NA
            )
            # indicate peak maximum
            points(rt[borders["peak_index"]], smoothed[borders["peak_index"]],
                pch = 19, col = "red"
            )
        },
        fl = slider(1, length(int) - 1, step = 2, initial = 7, label = "Filter length (odd)"),
        p = slider(1, 7, step = 1, initial = 3, label = "Polynomial order"),
        smooth_on = picker("Yes", "No", label = "Use smoothing for integration?")
    )
}
