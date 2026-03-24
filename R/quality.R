#' @title qcoreCalculator
#' @description
#' Implementation of the work by William Kumler to calculate quality metrics:
#'  https://github.com/wkumler/MS_metrics
#'
#' @param rt `numeric` vector with retention times
#' @param int `numeric` vector with corresponding intensities
#' @author William Kumler
#'
#' @export
qscoreCalculator <- function(rt, int) {
    # Check for bogus EICs
    if (length(rt) < 5) {
        return(list(SNR = NA, peak_cor = NA))
    }
    # Calculate where each rt would fall on a beta dist (accounts for missed scans)
    # scale all rt to between a range of 0 and 1
    scaled_rts <- (rt - min(rt)) / (max(rt) - min(rt))
    # Create a couple different skews and test fit
    maybe_skews <- c(2.5, 3, 4, 5) # Add 7 to catch more multipeaks and more noise
    # Add 2 to catch very slopey peaks and more noise
    # which of the 4 maybe_skews has the max correlation with int?
    best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x) {
        cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), int)
    }))]
    # the best template
    perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
    # calculate final correlation score
    peak_cor <- cor(perf_peak, int)
    # Calculate the normalized residuals
    residuals <- int / max(int) - perf_peak / max(perf_peak)
    # Calculate the minimum SD, after normalizing for any shape discrepancy
    old_res_sd <- sd(residuals)
    # diff() calculates the difference between the elements (returns vector)
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(norm_residuals)
    while (new_res_sd < old_res_sd) {
        old_res_sd <- new_res_sd
        norm_residuals <- diff(residuals)
        new_res_sd <- sd(residuals)
    }
    # Calculate SNR
    SNR <- (max(int) - min(int)) / sd(norm_residuals * max(int))
    # Return the quality score
    return(list(SNR = SNR, peak_cor = peak_cor))
}




# https://link.springer.com/article/10.1186/s12859-023-05533-4
# We also calculated several novel metrics from the raw m/z/RT/intensity values by
# extracting the data points falling within each individual peak’s m/z and RT
# bounding box (values between the XCMS-reported min and max) separately for each
# file. The data points were then linearly scaled to fall within the 0–1 range by
# subtracting the minimum RT and dividing by the maximum RT, then each scaled RT
# was fit to a beta distribution with α values of 2.5, 3, 4, and 5, and a fixed β
# value of 5. This approach allowed us to approximate a bell curve with increasing
# degrees of right-skewness and the beta distribution was chosen because it is
# constrained between 0 and 1 and simple and speedy to generate in R. For each α
# value, Pearson’s correlation coefficient (r) was calculated between the beta
# distribution and the raw data, with the highest value returned as a metric for
# how peak-shaped the data were (Fig. 8). The beta distribution with the highest
# r was also then used to estimate the noise level within the peak by scaling both
# the beta distribution probability densities and the raw data intensity values as
# described above, then subtracting the scaled beta distribution from the scaled
# intensity values, producing the residuals of the fit (Fig. 8). The signal-to-noise
# ratio (SNR) was calculated by dividing the maximum original peak height by the
# standard deviation of the residuals multiplied by the maximum height of the
# original peak. This method of SNR calculation allowed us to rapidly estimate the
# noise within the peak itself rather than relying on background estimation using
# data points outside the peak, which may not exist or may be influenced by additional
# mass signals [2]. If there were fewer than 5 data points, a missing value was
# returned and dropped in subsequent summary calculations. Accessing the raw data
# values also allowed us to calculate the proportion of “missed” scans in a peak for
# which an RT exists at other masses in the same sample but for which no data was
# produced at the selected m/z ratio, divided by the total number of scans between
# the min and max RTs.
