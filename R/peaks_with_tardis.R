# TARDIS-GITHUB
#' @title TARDIS Peak Detection
#'
#' @description
#' Main function of the TARDIS package that is called in the Shiny app.
#' Given data files and a list of targeted compounds it returns the area of
#' those peaks, optional diagnostic plots and several other parameters.
#' See vignette for a detailed tutorial.
#'
#' @param file_path `character(1)` Path to the .mzML or .mzXML files
#'     containing LC-MS data.
#' @param lcmsData  `MsExperiment` MsExperiment containing the data to be
#'     preprocessed. Sampledata should at least include run type and should
#'     match the later provided "QC" or "sample" pattern.
#' @param dbData Output of [createTargetList()]
#' @param ppm `numeric(1)` Allowed deviance from given m/z of targets in ppm.
#' @param rtdev `numeric(1)` Allowed deviance from given retention time of
#'     compound, defines search window for the peak picking algorithm.
#' @param mass_range `numeric(2)` If the user uses data with overlapping mass
#'     windows, only one mass window at the time can be analyzed, specify this
#'     here.
#' @param polarity `character(1)` Ionisation mode to be considered, can be
#'     either "positive" or "negative"
#' @param output_directory `character(1)` Provide directory to store output
#' @param plots_samples `logical(1)` Create plots for all samples
#' @param plots_QC `logical(1)` Create plots for all QCs
#' @param diagnostic_plots `logical(1)` Create diagnostic plots of 5 QCs
#'     spread across the runs
#' @param batch_positions `list` Indicate start and end file of each batch,
#'     e.g. `list(c(1,20),c(21,40))`
#' @param QC_pattern `character(1)`  Pattern of QC files
#' @param sample_pattern `character(1)` Pattern of sample files
#' @param rt_alignment `logical(1)` Align retention time based on internal
#'     standard compounds in the QC samples.
#' @param int_std_id `character` Provide ID's of internal standard compounds for
#'     retention time alignment
#' @param screening_mode `logical(1)` Run the algorithm over 5 QCs to quickly
#'     check retention time shifts
#' @param smoothing `logical(1)` Smooth the peaks with [sgolayfilt()]
#' @param max_int_filter `numeric(1)` Disregard peaks with a max. int. lower
#'     than this value
#' @param num_cores `numeric(1)` Number of cores to use for parallelization
#'
#'
#' @import MsExperiment
#' @import xcms
#' @importFrom Spectra MsBackendMzR
#' @importFrom Spectra filterMzRange
#' @importFrom Spectra filterEmptySpectra
#' @importFrom Spectra filterDataOrigin
#' @importFrom Spectra filterRt
#' @importFrom Spectra filterPolarity
#' @importFrom Spectra dataOrigin
#' @importFrom Spectra addProcessing
#' @importFrom Spectra Spectra
#' @importFrom signal sgolayfilt
#' @importFrom xcms PeakGroupsParam
#' @importFrom xcms adjustRtime
#' @importFrom xcms applyAdjustedRtime
#' @importFrom xcms rtime
#' @importFrom xcms intensity
#' @importFrom pracma trapz
#' @importFrom BiocParallel SnowParam
#' @importFrom tidyr spread
#' @importFrom writexl write_xlsx
#' @importFrom dplyr summarise
#' @importFrom dplyr summarise_at
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @importFrom dplyr first
#' @importFrom dplyr mutate
#' @importFrom S4Vectors DataFrame
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' @importFrom stats loess
#' @import parallel
#'
#' @return returns `list` with auc table and feature table with summarized
#'     stats per compound. Outputs plots and other tables to output folder.
#'
#' @export
#'

# if these packages don't exist, install them
# and load them
if (!require("pbapply")) install.packages("pbapply")
library(pbapply)

if (!require("parallel")) install.packages("parallel")
library(parallel)

if (!require("tidyr")) install.packages("tidyr")
library(tidyr)

if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

if (!require("xcms")) install.packages("xcms")
library(xcms)


## the following functions are for repetitive code in tardisPeaks function

# function for retention time alignment
# need to fix!
rtAlignment <- function(minFraction, span, int_std, data_QC) {
  # The PeakGroupsParam class allows to specify all settings for the retention time adjustment based on house keeping peak groups present in most samples.
  # minFraction: numeric(1) between 0 and 1 defining the minimum required fraction of samples in which peaks for the peak group were identified.
  # span: numeric(1) defining the degree of smoothing (if smooth = "loess"). This parameter is passed to the internal call to loess.
  # peakGroupsMatrix: matrix of (raw) retention times for the peak groups on which the alignment should be performed.
  # peakGroupsMatrix: Each column represents a sample, each row a feature/peak group.
  param <- PeakGroupsParam(minFraction = minFraction,
                           span = span,
                           peakGroupsMatrix = int_std)
  # adjustRtime() performs the alignment (retention time correction) and stores the results alongside the raw data.
  data_QC <- adjustRtime(data_QC, param = param)
  # applyAdjustedRtime() replaces the original (raw) retention times with the newly calculated adjusted retention times.
  data_QC <- applyAdjustedRtime(data_QC)
  return (data_QC)
}

# function for Savitsky-Golay smoothing
smoothingSG <- function(p = 3,
                        tr,  		# ex. tr = dbData$tr[j]
                        all_files_i, 	# call function using all_files[i]
                        spectra_QC,
                        rt_input,
                        mz_input,
                        smoothing,
                        ratio = 0.5,
                        baseline_correction = FALSE) {
  # first, filter spectra by data and mz range
  filtered_spectra_data_mz <- spectra_QC |>
    filterDataOrigin(all_files_i) |>
    filterMzRange(mz_input)

  # to address peak cutoff in the final plots, increase rt window to the left or right
  for(j in 1:3){  # limit number of iterations to prevent eternal loop
    # next, filter the previous results with rt range
    filtered_spectra <- filtered_spectra_data_mz |>
      filterRt(rt_input)  # >= min, <= max (includes borders)

    eic <- extract_eic(filtered_spectra)
    # error check
    if (nrow(eic) == 0) break # Exit if no data found
    rt <- eic[, 1L]
    int <- eic[, 2L]

    # Baseline correction: change int values
    if (baseline_correction == TRUE && is.null(int) == FALSE){
      int <- int - rep(min(int), length(int))
    }

    max_val <- max(int, na.rm = TRUE)
    # If intensity at start or end is > 10% of the max, the peak is likely cut off
    is_cutoff_left  <- (!is.na(int[1])) && (int[1] > (max_val * 0.1))
    is_cutoff_right <- (!is.na(int[length(int)])) && (int[length(int)] > (max_val * 0.1))

    rt_left <- rt_input[1]
    rt_right <- rt_input[2]

    if (is_cutoff_left){
      rt_left <- rt_left - 10
    }
    if (is_cutoff_right){
      rt_right <- rt_right + 10
    }
    rt_input <- c(rt_left, rt_right)

    if(!is_cutoff_left && !is_cutoff_right){
      break
    }
  }
  length_error <- rt_input[2] - rt_input[1] + 1
  # EMERGENCY EXIT: If no data was found after filtering
  if (is.null(eic) || nrow(eic) == 0) {
    return(list(
      rt = rep(NA, length_error),
      int = rep(NA, length_error),
      border = rep(NA, 3)
    ))
  }

  # NA intensities are set to zero --> should change this so only NA's
  # at the edges get changed to zero, so the ones IN the peak will be
  # imputed
  int[which(is.na(int))] <- 0
  # if intensity length is under 7, lower filter length
  # to odd number <= intensity length
  if (length(int) < 7) {
    if (length(int) %% 2 == 0) {
      fl <- length(int) - 1
    } else {
      fl <- length(int)
    }
  } else {
    fl <- 7
  }

  # Savitzky-Golay requirement: n must be > p (else it will crash)
  if (fl <= p) fl <- p + (if (p %% 2 == 0) 1 else 2)

  smoothed <- sgolayfilt(int, p = p, n = fl)
  if (smoothing == TRUE) {
    int <- smoothed
    int[int < 0] <- 0
  }
  border <- find_peak_points(rt, smoothed, tr, .check = FALSE, ratio = ratio)  # edit - get ratio param from user
  # R cannot return multiple values, wrap them in a list is ok
  # return list of lists
  return(list(
    rt = if (is.null(rt)) rep(NA, length_error) else rt,
    int = if (is.null(int)) rep(NA, length_error) else int,
    border = if (is.null(border)) rep(NA, 3) else border
  ))
}

# a function to check if a peak is valid
checkValidPeak <- function(x, y, d, sample_name, int, rt, border) {

  # Force consistent scalar outputs
  if (length(unique(y)) > 1 && !any(is.na(border))) {

    auc <- pracma::trapz(x, y)
    pop <- length(x)
    qscore <- qscoreCalculator(x, y)

    MaxInt_val <- ifelse(length(int) >= border[3L], int[border[3L]], NA)
    foundRT_val <- ifelse(length(rt) >= border[3L], rt[border[3L]], NA)

    out <- data.frame(
      Component = as.character(d$ID),
      Sample = as.character(sample_name),
      AUC = as.numeric(auc),
      MaxInt = as.numeric(MaxInt_val),
      SNR = as.numeric(qscore[1]),
      peak_cor = as.numeric(qscore[2]),
      foundRT = as.numeric(foundRT_val),
      pop = as.numeric(pop),
      ID = as.character(d$ID),
      NAME = as.character(d$NAME),
      mz = as.numeric(d$`m/z`),
      tr = as.numeric(d$tr),
      stringsAsFactors = FALSE
    )

  } else {

    out <- data.frame(
      Component = as.character(d$ID),
      Sample = as.character(sample_name),
      AUC = NA_real_,
      MaxInt = NA_real_,
      SNR = NA_real_,
      peak_cor = NA_real_,
      foundRT = NA_real_,
      pop = NA_real_,
      ID = as.character(d$ID),
      NAME = as.character(d$NAME),
      mz = as.numeric(d$`m/z`),
      tr = as.numeric(d$tr),
      stringsAsFactors = FALSE
    )
  }

  colnames(out) <- c("Component", "Sample", "AUC", "MaxInt",
                     "SNR", "peak_cor", "foundRT", "pop", "ID", "NAME", "mz", "tr") # Force the names here

  return(out)
}

# function for data handling
dataHandling <- function(files, string, QC_pattern, polarity){
  data <- MsExperiment()
  experimentFiles(data) <-
    MsExperimentFiles(mzML = setNames(files, basename(
      tools::file_path_sans_ext(files)
    )))

  sampleData(data) <- DataFrame(sample_index = 1:length(files),
                                      spectraOrigin = files)
  if (string == "QC_files"){
    sampleData(data)$type <- "QC"
  }
  if (string == "files_batch"){
    # Define study and QC samples --> all not QC files are deemed study files
    sampleData(data)$type <- "study"
    sampleData(data)$type[grep(pattern = QC_pattern, files)] <- "QC"
  }
  sp <- Spectra(
    experimentFiles(data)[["mzML"]],
    backend = MsBackendMzR(),
    BPPARAM = SnowParam(workers = 1L)
  )
  if (polarity == "positive") {
    spectra(data) <- filterPolarity(sp, 1)
  } else if (polarity == "negative") {
    spectra(data) <- filterPolarity(sp, 0)
  }

  sampleData(data)$raw_file <- normalizePath(files)
  data <- linkSampleData(data, with = "sampleData.raw_file = spectra.dataOrigin")

  return (data)
}

# for cleaning dataframes, assigning colnames
standardize_results <- function(df) {

  required_cols <- c(
    "Component", "Sample", "AUC", "MaxInt",
    "SNR", "peak_cor", "foundRT", "pop",
    "ID", "NAME", "mz", "tr"
  )

  # Add missing columns
  missing <- setdiff(required_cols, colnames(df))
  for (m in missing) df[[m]] <- NA

  # Force atomic columns (CRUCIAL)
  df[] <- lapply(df, function(x) {
    if (is.list(x)) unlist(x) else x
  })

  # Reorder
  df <- df[, required_cols, drop = FALSE]

  colnames(df) <- c("Component", "Sample", "AUC", "MaxInt",
                     "SNR", "peak_cor", "foundRT", "pop", "ID", "NAME", "mz", "tr") # Force the names here
  return(df)
}

safe_bind <- function(x) {
  # 1. Remove NULLs and empty objects
  x <- Filter(function(df) is.data.frame(df) && nrow(df) > 0, x)
  if (length(x) == 0) return(NULL)

  # 2. Clean list-columns and ensure DF structure
  x <- lapply(x, function(df) {
    # Ensure it's a data frame first
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    # Clean columns without losing the DF structure
    for (col_name in colnames(df)) {
      if (is.list(df[[col_name]])) {
        # Replace list with first element or NA
        df[[col_name]] <- sapply(df[[col_name]], function(item) {
          if (length(item) == 0) return(NA)
          return(as.character(unlist(item)[1]))
        })
      }
    }
    return(df)
  })

  # 3. Combine first, THEN rename
  # bind_rows is smart—it will align columns by name automatically
  out <- dplyr::bind_rows(x)

  # 4. SAFETY CHECK: Only rename if column count matches
  expected_names <- c("Component", "Sample", "AUC", "MaxInt",
                      "SNR", "peak_cor", "foundRT", "pop", "ID", "NAME", "mz", "tr")

  if (ncol(out) == length(expected_names)) {
    colnames(out) <- expected_names
  } else {
    warning(paste("Column count mismatch! Expected 12, found", ncol(out)))
  }

  rownames(out) <- NULL
  return(out)
}

## jo: wouldn't it be better to call the function on a data object instead
## of a file path? The (advanced) user could eventually do some more quality
## checks on the data before?
## Pablo: Definitely something I should do, but I might do it a separate
## function, since for use with the GUI the file input is nice.

# If file_path is provided: The function goes to your hard drive and loads the raw files.
# If lcmsData is provided: The function skips the loading step and uses the data already sitting in your R memory.
tardisPeaks <-
  function(file_path = NULL,
           lcmsData = NULL,
           dbData,
           ppm = 5,
           rtdev = 18,
           mass_range = NULL,
           polarity = "positive",
           output_directory,
           plots_samples = FALSE,
           plots_QC = FALSE,
           diagnostic_plots = TRUE,
           batch_positions,
           QC_pattern = "QC",
           sample_pattern = "",
           rt_alignment = TRUE,
           int_std_id,
           screening_mode = FALSE,
           smoothing = TRUE,
           max_int_filter = NULL,
           num_cores = 1) { # edited GUI!
    # Setup the cluster (num_cores parameter)
    num_cores <- num_cores
    cl <- makeCluster(num_cores)
    # cluster is stopped if there is an error
    on.exit(stopCluster(cl), add = TRUE)

    if (is.null(file_path) == FALSE) {
      files <-
        list.files(file_path, full.names = T, pattern = "mzML|mzXML")
    }
    if (is.null(lcmsData) == FALSE) {
      if (polarity == "positive") {
        spectra(lcmsData) <- filterPolarity(spectra(lcmsData), 1)  # change the spectra in the lcmsData object
      } else if (polarity == "negative") {
        spectra(lcmsData) <- filterPolarity(spectra(lcmsData), 0)
      }
      suppressWarnings(
        lcmsData <- linkSampleData(lcmsData, with = "sampleData.raw_file = spectra.dataOrigin")
      )
    }
    if (is.null(mass_range) == FALSE) {
      dbData <-
        dbData[which(dbData$`m/z` < mass_range[2] &
                       dbData$`m/z` > mass_range[1]), ]
    }
    info_compounds <- dbData
    if (screening_mode == TRUE) {
      if (is.null(file_path) == FALSE) {
        QC_files <-
          files[grep(pattern = QC_pattern, files)]
        data_QC <- dataHandling(QC_files, "QC_files", QC_pattern, polarity)
      } else {
        data_QC <- lcmsData[which(sampleData(lcmsData)$type == QC_pattern)]
      }
      if (is.null(mass_range) == FALSE) {
        data_QC <- filterSpectra(data_QC, filterMzRange, mz = mass_range) |>
          filterSpectra(filterEmptySpectra)
      } else {
        data_QC <- data_QC
      }
      spectra_QC <- data_QC@spectra
      checkScans(spectra_QC)
      data_QC@spectra <- spectra_QC
      all_files <- unique(dataOrigin(spectra_QC))
      ## Create ranges for all compounds
      ranges <- createRanges(data_QC, dbData, ppm, rtdev)
      ## Get mz & rt ranges
      mzRanges <- ranges[[1L]]
      rtRanges <- ranges[[2L]]

      if (rt_alignment == TRUE) {
        ## Get the ranges for the internal standard compounds
        internal_standards_rt <-
          rtRanges[which(dbData$ID %in% int_std_id), ]
        internal_standards_mz <-
          mzRanges[which(dbData$ID %in% int_std_id), ]
        dbData_std <- dbData[which(dbData$ID %in% int_std_id), ]
        ## Get QC sample names
        sample_names <-
          lapply(data_QC@sampleData$spectraOrigin, basename)
        ## Initiate vectors
        int_std_foundrt <- c(length(sample_names))
        int_std <- c(dim(internal_standards_rt)[1] * length(sample_names))
        ## Retrieve foundRT of internal standards in QC's,
        ## loop over all samples and all internal standards
        # Parallelize the OUTER loop (j = internal standards)

        clusterEvalQ(cl, {
          library(MsExperiment)
          library(Spectra)
          library(signal)
          library(xcms)
          library(pracma)
          library(BiocParallel)
          library(tidyr)
          library(writexl)
          library(dplyr)
          library(S4Vectors)
          library(pbapply)
          library(parallel)
        })

        # Export everything the workers need to know
        # This includes variables AND the function smoothingSG
        # include ls("package:TARDIS"): to update worker processes
        clusterExport(cl, varlist = c(ls("package:TARDIS", all.names = TRUE), "sample_names", "dbData_std", "all_files",
                                      "spectra_QC", "internal_standards_rt",
                                      "internal_standards_mz", "smoothing",
                                      "int_std_foundrt", "smoothingSG", "batch_positions"), envir = environment())  # find variables/functions anywhere in the code
        results_list1 <- pblapply(1:nrow(internal_standards_rt), function(j) {
          local_found_rt <- numeric(length(sample_names))
          for (i in 1:length(sample_names)) {
            res <- smoothingSG(
              3,
              dbData_std$tr[j],
              all_files[i],
              spectra_QC,
              internal_standards_rt[j, ],
              internal_standards_mz[j, ],
              smoothing,
              0.5,
              TRUE
            )
            local_found_rt[i] <- res$rt[res$border[3L]]  # rt of peak
          }
          return(local_found_rt)
        }, cl = cl)   # results_list1: list of lists

        int_std <- do.call(rbind, results_list1)

        ## retention time adjustment
        data_QC <- rtAlignment(minFraction = 0.9,
                               span = 0.5,
                               int_std,  # no transpose!; peakGroupsMatrix: matrix with the retention times of the peak groups.
                               # peakGroupsMatrix: Each column represents a sample, each row a feature/peak group.
                               data_QC)
      }
      # Find all targets in x QC's

      spectra_QC <- data_QC@spectra
      all_files <- unique(dataOrigin(spectra_QC))
      sample_names <-
        lapply(data_QC@sampleData$spectraOrigin, basename)

      length_results_screening <- dim(rtRanges)[1] * length(sample_names)
      results_screening <-
        data.frame(
          Component = character(length_results_screening),
          Sample = character(length_results_screening),
          AUC = numeric(length_results_screening),
          SNR = numeric(length_results_screening),
          peak_cor = numeric(length_results_screening),
          foundRT = numeric(length_results_screening),
          pop = numeric(length_results_screening)
        )

      clusterExport(cl, varlist = c(ls("package:TARDIS", all.names = TRUE), "sample_names", "dbData", "all_files",
                                    "spectra_QC", "smoothing",
                                    "int_std_foundrt", "smoothingSG",
                                    "rtRanges", "mzRanges", "checkValidPeak",
                                    "plotDiagnostic", "diagnostic_plots",
                                    "output_directory"), envir = environment()) # find the variables in this function, not the global environment!

      # Parallelize the OUTER loop (j = internal standards)
      results_list2 <- pblapply(1:dim(rtRanges)[1], function(j) {  # for all target compounds
        compound_info <- dbData[j, ]  # id, name, mz, rt
        rt_list <- list()
        int_list <- list()
        x_list <- list()
        y_list <- list()
        results_screening_row <- vector("list", length(sample_names))

        compound_results <- list() # store each sample's results
        for (i in 1:length(sample_names)) {
          res <- smoothingSG(3,
                             dbData$tr[j],
                             all_files[i],
                             spectra_QC,
                             rtRanges[j, ],
                             mzRanges[j, ],
                             smoothing,
                             0.5,
                             TRUE)
          rt <- res$rt
          int <- res$int
          border <- res$border

          idx <- border[1L]:border[2L]
          x <- rt[idx]
          y <- int[idx]
          rt_list <- c(rt_list, list(rt))
          int_list <- c(int_list, list(int))
          x_list <- c(x_list, list(x))
          y_list <- c(y_list, list(y))

          results_screening_row[[i]] <- checkValidPeak(x,
                                                  y,
                                                  dbData[j, ],  # this is compound info
                                                  sample_names[i],
                                                  int,
                                                  rt,
                                                  border)
        }
        # Create and save the plot for the current component
        batchnr <- 1
        if (diagnostic_plots == TRUE) {
          plotDiagnostic(
            compound_info,
            output_directory,
            rt_list,
            int_list,
            x_list,
            y_list,
            batchnr,
            sample_names
          )
        }
        return(safe_bind(results_screening_row))
      }, cl=cl)

      # Combine and standardize result
      results_screening <- safe_bind(results_list2)
      results_screening <- standardize_results(results_screening)

      avg_metrics_table <- results_screening %>%
        group_by(Component) %>%
        summarise_at(vars(-Sample), list(~ if (is.numeric(.)) {
          mean(., na.rm = TRUE)
        } else {
          first(.)
        }))
      write.csv(avg_metrics_table,
                file = paste0(output_directory, "qc_screening.csv"))
      print("csv file saved!")

    } else {  # if screening mode is false
      ## Loop over the batches
      for (batchnr in 1:length(batch_positions)) {  # 1
        dbData <- info_compounds # need to reset? better to keep updated from last batch?
        if (is.null(file_path) == FALSE) {
          files_batch <-
            files[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
          data_batch <- dataHandling(files_batch, "files_batch", QC_pattern, polarity)
        } else {
          data_batch <- lcmsData[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
        }
        if (is.null(mass_range) == FALSE) {
          data_batch <- filterSpectra(data_batch, filterMzRange, mz = mass_range) |>
            filterSpectra(filterEmptySpectra)
        } else {
          data_batch <- data_batch
        }
        spectra_batch <- data_batch@spectra
        checkScans(spectra_batch)
        data_batch@spectra <- spectra_batch

        data_QC <-
          data_batch[which(sampleData(data_batch)$type == "QC")]  # among the data, select only QC samples
        if (is.null(mass_range == FALSE)) {
          spectra_QC <- data_QC@spectra |>
            filterMzRange(mass_range) |>
            filterEmptySpectra()
        } else {
          spectra_QC <- data_QC@spectra
        }
        all_files <- unique(dataOrigin(spectra_QC))
        ranges <- createRanges(data_QC, dbData, ppm, rtdev)
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        if (rt_alignment == TRUE) {
          ## Get the ranges for the internal standard compounds
          internal_standards_rt <-
            rtRanges[which(dbData$ID %in% int_std_id), ]
          internal_standards_mz <-
            mzRanges[which(dbData$ID %in% int_std_id), ]
          dbData_std <- dbData[which(dbData$ID %in% int_std_id), ]
          sample_names_batch <-
            lapply(data_batch@sampleData$spectraOrigin, basename)  # length of sample_names_batch: number of total samples
          sample_names_QC <-
            lapply(data_QC@sampleData$spectraOrigin, basename)  # length of sample_names_QC: number of QC samples
          int_std_foundrt <- c(length(sample_names_batch))
          int_std <- c(dim(internal_standards_rt)[1] * length(sample_names_batch))

          clusterExport(cl, varlist = c(ls("package:TARDIS", all.names = TRUE), "sample_names_batch", "sample_names_QC", "dbData_std", "all_files",
                                        "spectra_QC", "internal_standards_rt",
                                        "internal_standards_mz", "smoothing",
                                        "int_std_foundrt", "smoothingSG"), envir = environment()) # find the variables in this function & the global environment!

          # Parallelize the OUTER loop (j = internal standards)
          results_list3 <- pblapply(1:dim(internal_standards_rt)[1], function(j) {
            rt_list <- list()
            int_list <- list()
            x_list <- list()
            y_list <- list()
            local_found_rt <- numeric(length(sample_names_QC))
            for (i in 1:length(sample_names_QC)) {
              res <- smoothingSG(
                3,
                dbData_std$tr[j],
                all_files[i],
                spectra_QC,
                internal_standards_rt[j, ],
                internal_standards_mz[j, ],
                smoothing,
                0.5,
                TRUE
              )
              local_found_rt[i] <- res$rt[res$border[3L]]  # rt of peak
            }

            return(local_found_rt)
          }, cl=cl) # results_list3: list of lists
          int_std <- do.call(rbind, results_list3)

          results_screening_row <- NULL  # need this line for error handling (the computer keeps looking for results_screening_row and doesnt find it)
          data_QC <- rtAlignment(minFraction = 0.9,
                                    span = 0.5,
                                    int_std,
                                    data_QC)
        }
        ## Now, we try and find ALL compounds in the QC samples and save their
        ## foundRT to search the compounds at that RT in the sample files
        ## Skip this step if there aren't any QC's available.
        data_QC <-
          data_batch[which(sampleData(data_batch)$type == "QC")]
        if (length(data_QC) != 0) {
          sample_names_QC <-
            lapply(data_QC@sampleData$spectraOrigin, basename)
          if (is.null(mass_range) == FALSE) {
            spectra_QC <- data_QC@spectra |>
              filterMzRange(mass_range) |>
              filterEmptySpectra()
          } else {
            spectra_QC <- data_QC@spectra
          }
          all_files <- unique(dataOrigin(spectra_QC))

          length_results_QCs_batch <- dim(rtRanges)[1] * length(sample_names_QC)
          results_QCs_batch <-
            data.frame(
              Component = character(length_results_QCs_batch),
              Sample = character(length_results_QCs_batch),
              AUC = numeric(length_results_QCs_batch),
              SNR = numeric(length_results_QCs_batch),
              peak_cor = numeric(length_results_QCs_batch),
              foundRT = numeric(length_results_QCs_batch),
              pop = numeric(length_results_QCs_batch)
            )

          clusterExport(cl, varlist = c(ls("package:TARDIS", all.names = TRUE), "sample_names_batch", "sample_names_QC", "dbData", "all_files",
                                        "spectra_QC", "smoothing",
                                        "int_std_foundrt", "smoothingSG",
                                        "rtRanges", "mzRanges", "checkValidPeak",
                                        "plotDiagnostic", "plots_QC",
                                        "output_directory", "batch_positions"), envir = environment()) # find the variables in this function, not the global environment!

          # Parallelize the OUTER loop (j = internal standards)
          results_list4 <- pblapply(1:dim(rtRanges)[1], function(j) {
            compound_info <- dbData[j, ]
            rt_list <- list()
            int_list <- list()
            x_list <- list()
            y_list <- list()

            results_QCs_batch_row <- vector("list", length(sample_names_QC))

            compound_results <- list() # store each sample's results

            for (i in 1:length(sample_names_QC)) {
              res <- smoothingSG(3,
                                 dbData$tr[j],
                                 all_files[i],
                                 spectra_QC,
                                 rtRanges[j, ],
                                 mzRanges[j, ],
                                 smoothing,
                                 0.5,
                                 TRUE)
              rt <- res$rt
              int <- res$int
              border <- res$border

              idx <- border[1L]:border[2L]
              x <- rt[idx]
              y <- int[idx]
              rt_list <- c(rt_list, list(rt))
              int_list <- c(int_list, list(int))
              x_list <- c(x_list, list(x))
              y_list <- c(y_list, list(y))

              results_QCs_batch_row[[i]] <- checkValidPeak(x,
                                                      y,
                                                      dbData[j, ],
                                                      sample_names_QC[i],
                                                      int,
                                                      rt,
                                                      border)  # is a 1-row dataframe
            }
            # Create and save the plot for the current component
            batchnr <- 1
            if (plots_QC == TRUE) {
              plotQCs(
                compound_info,
                output_directory,
                rt_list,
                int_list,
                x_list,
                y_list,
                batchnr,
                sample_names_QC
              )
            }
            return(safe_bind(results_QCs_batch_row))
          }, cl=cl)

          # Combine and standardize result
          results_QCs_batch <- safe_bind(results_list4)
          results_QCs_batch <- standardize_results(results_QCs_batch)

          # Since the same metabolite (ID) is measured multiple times across the
          # batch, it calculates the average observed retention time (foundRT) for
          # each molecule. This creates a reference list of where the peaks actually
          # appeared in this specific batch.
          # need only ID and foundRT columns
          new_rt_avg <- results_QCs_batch %>%
            group_by(ID) %>%
            summarise(mean = mean(foundRT), na.rm = TRUE)
          dbData <- merge(dbData, new_rt_avg, by = "ID")
          dbData$trold <- dbData$tr
          dbData$tr <- new_rt_avg$mean
          ## If no RT is found, restore old RT
          for (k in 1:dim(dbData)[1]) {
            if (is.na(dbData$tr[k]) == TRUE) {
              dbData$tr[k] <- dbData$trold[k]
            }
          }
        }
        ## Next do the whole analysis for the samples in the same batch of the
        ## QC's to find ALL the compounds at the corrected RT. (SAMPLES + QC)
        ## Get sample data
        sample_names_batch <-
          lapply(data_batch@sampleData$spectraOrigin, basename)
        # Create ranges around new RT
        ranges <- createRanges(data_batch, dbData, ppm, rtdev)
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        if (is.null(mass_range) == FALSE) {
          spectra <- data_batch@spectra |>
            filterMzRange(mz = mass_range) |>
            filterEmptySpectra()
        } else {
          spectra <- data_batch@spectra
        }

        length_results_samples <- dim(rtRanges)[1] * length(sample_names_batch)
        results_samples <-
          data.frame(
            Component = character(length_results_samples),
            Sample = character(length_results_samples),
            AUC = numeric(length_results_samples),
            SNR = numeric(length_results_samples),
            peak_cor = numeric(length_results_samples),
            foundRT = numeric(length_results_samples),
            pop = numeric(length_results_samples)
          )

        all_files <- unique(dataOrigin(spectra))

        clusterExport(cl, varlist = c(ls("package:TARDIS", all.names = TRUE), "sample_names_batch", "sample_names_QC", "dbData", "all_files",
                                      "spectra", "smoothing",
                                      "smoothingSG",
                                      "rtRanges", "mzRanges", "checkValidPeak",
                                      "plotSamples", "plotDiagnostic",
                                      "plots_samples", "diagnostic_plots",
                                      "output_directory", "batch_positions"), envir = environment()) # find the variables in this function, not the global environment!

        # Parallelize the OUTER loop (j = internal standards)
        results_list5 <- pblapply(1:dim(rtRanges)[1], function(j){
          compound_info <- dbData[j, ]
          rt_list <- list()
          int_list <- list()
          x_list <- list()
          y_list <- list()

          results_samples_row <- vector("list", length(sample_names_batch))
          for (i in 1:length(sample_names_batch)) {
            res <- smoothingSG(3,
                               dbData$tr[j],
                               all_files[i],
                               spectra,
                               rtRanges[j, ],
                               mzRanges[j, ],
                               smoothing,
                               0.5,
                               TRUE)
            rt <- res$rt
            int <- res$int
            border <- res$border

            idx <- border[1L]:border[2L]
            x <- rt[idx]
            y <- int[idx]
            rt_list <- c(rt_list, list(rt))
            int_list <- c(int_list, list(int))
            x_list <- c(x_list, list(x))
            y_list <- c(y_list, list(y))

            results_samples_row[[i]] <- checkValidPeak(x,
                                                  y,
                                                  dbData[j, ],
                                                  sample_names_batch[i],
                                                  int,
                                                  rt,
                                                  border)

          }
          if (plots_samples == TRUE) {
            plotSamples(
              compound_info,
              output_directory,
              rt_list,
              int_list,
              x_list,
              y_list,
              batchnr,
              sample_names_batch
            )
          }
          if (diagnostic_plots == TRUE) {
            plotDiagnostic(
              compound_info,
              output_directory,
              rt_list,
              int_list,
              x_list,
              y_list,
              batchnr,
              sample_names_batch
            )
          }
          return (safe_bind(results_samples_row))
        }, cl=cl)

        # Combine dataframes
        results_samples <- safe_bind(results_list5)
        results_samples <- standardize_results(results_samples)
      }
      results <- results_samples  # 1 batch, same copy
      results <- standardize_results(results)

      # if (is.null(max_int_filter) == FALSE &&
      #     max_int_filter != 0) {
      #   results <- results[which(results$MaxInt >= max_int_filter), ]
      # }
      # Use is.numeric and !is.na to ensure the value is actually a number
      if (is.numeric(max_int_filter) && !is.na(max_int_filter) && max_int_filter > 0) {
        results <- results[which(results$MaxInt >= max_int_filter), ]
      }

      # This prevents: duplicated (Component, Sample) pairs,
      # pivot chaos, column explosions
      results <- results %>%
        group_by(Component, Sample) %>%
        summarise(
          AUC = mean(AUC, na.rm = TRUE),
          MaxInt = mean(MaxInt, na.rm = TRUE),
          SNR = mean(SNR, na.rm = TRUE),
          peak_cor = mean(peak_cor, na.rm = TRUE),
          foundRT = mean(foundRT, na.rm = TRUE),
          pop = mean(pop, na.rm = TRUE),
          ID = first(ID),
          NAME = first(NAME),
          mz = first(mz),
          tr = first(tr),
          .groups = "drop"
        )

      auc_table <- results %>%
        dplyr::select(Component, Sample, AUC) %>%
        tidyr::pivot_wider(
          names_from = Sample,
          values_from = AUC,
          values_fill = NA   # prevents missing structure issues
        )
      write.csv(auc_table, file = paste0(output_directory, "auc_table.csv"))

      pop_table <- results %>%
        dplyr::select(Component, Sample, pop) %>%
        tidyr::pivot_wider(
          names_from = Sample,
          values_from = pop,
          values_fill = NA
        )
      write.csv(pop_table, file = paste0(output_directory, "pop_table.csv"))

      SNR_table <- results %>%
        dplyr::select(Component, Sample, SNR) %>%
        tidyr::pivot_wider(
          names_from = Sample,
          values_from = SNR,
          values_fill = NA
        )
      write.csv(SNR_table, file = paste0(output_directory, "snr_table.csv"))

      int_table <- results %>%
        dplyr::select(Component, Sample, MaxInt) %>%
        tidyr::pivot_wider(
          names_from = Sample,
          values_from = MaxInt,
          values_fill = NA
        )
      write.csv(int_table, file = paste0(output_directory, "int_table.csv"))

      peakcor_table <- results %>%
        dplyr::select(Component, Sample, peak_cor) %>%
        tidyr::pivot_wider(
          names_from = Sample,
          values_from = peak_cor,
          values_fill = NA
        )
      write.csv(peakcor_table,
                file = paste0(output_directory, "peakcor_table.csv"))


      # summarize feature table based on QC's
      avg_metrics_table <- NULL
      if (length(data_QC) != 0) {
        QC_results <- results[grep("QC", results$Sample), ]
        avg_metrics_table <- QC_results %>%
          group_by(Component) %>%
          summarise_at(vars(-Sample), list(~ if (is.numeric(.)) {
            mean(., na.rm = TRUE)
          } else {
            first(.)
          }))
        avg_metrics_table[] <- lapply(avg_metrics_table, function(x) {
          if (is.list(x)) unlist(x) else x
        })
        write_xlsx(avg_metrics_table,
                   paste0(output_directory, "feat_table.xlsx"))
      }

      # save input parameters to .csv

      input_params <- data.frame(
        "ppm" = .collapse_safe(ppm),
        "rtdev" = .collapse_safe(rtdev),
        "mass_range_low" = .collapse_safe(mass_range[1]),
        "mass_range_high" = .collapse_safe(mass_range[2]),
        "polarity" = .collapse_safe(polarity),
        "batch_positions" = .collapse_safe(batch_positions),
        "QC_pattern" = .collapse_safe(QC_pattern),
        "sample_pattern" = .collapse_safe(sample_pattern),
        "int_std_id" = .collapse_safe(int_std_id),
        "screening_mode" = .collapse_safe(screening_mode),
        "rt_alignment" = .collapse_safe(rt_alignment),
        "plots_samples" = .collapse_safe(plots_samples),
        "plots_QC" = .collapse_safe(plots_QC),
        "diagnostic_plots" = .collapse_safe(diagnostic_plots),
        "max_int_filter" = .collapse_safe(max_int_filter),
        "smoothing" = .collapse_safe(smoothing),
        stringsAsFactors = FALSE
      )

      write.csv(
        t(input_params),
        file = paste0(output_directory, "input_params.csv"),
        row.names = TRUE
      )

      return(list(auc_table, avg_metrics_table))
    }
  stopCluster(cl)  # stop parallel processes
  }
