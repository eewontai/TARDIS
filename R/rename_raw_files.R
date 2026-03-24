# this function syncs your raw data files with your experimental notes.
# It reads a Thermo XCalibur sequence list (the .csv) and renames the physical
# .raw files on your computer by appending the "Sample Name" to the original "File Name."


#' Rename Thermo .raw files
#' Script to read XCalibur sequence list and rename files according to run type
#' (sample, QC, blank, etc...)
#' sequence list: XCalibur instructions for the machine
#' @param data_path_raw_files path to raw files. Important: all runs from the
#'     sequence list need to have a corresponding .raw file in your input folder.
#' @param data_path_list Path to exported sequence list in .csv format,
#'     all columns may be exported
#'
#' @importFrom stringr str_replace_all
#'
#' @export
#'
#' @author Pablo Vangeenderhuysen

renameRawFiles <- function(data_path_raw_files, data_path_list) {
    file_names_original <- list.files(data_path_raw_files)
    # gsub() searches for a regular expression in a string and replaces it
    file_names <- gsub(".raw", "", list.files(data_path_raw_files))

    # read.csv: reads a file in table format and creates data frame
    # arguments: file to read, skip 1 line of data file before beggining to read data,
    # do not check the names of variables in data frame to ensure that they are syntactically valid variable names
    sequence <- read.csv(data_path_list, skip = 1, check.names = FALSE) # read sequence list
    sequence[, "save_name"] <- NA
    # fill in 'save_name' column
    for (k in 1:length(file_names)) {
        if (file_names[k] == sequence$`File Name`[k]) {
            sequence$save_name[k] <- paste0(
                file_names[k], "_",
                sequence$`Sample Name`[k], ".raw"
            )
        }
    }
    sequence$save_name <- str_replace_all(sequence$save_name, "/", "")  # "/" is illegal character - errors
    new_names <- sequence$save_name
    file.rename(
        paste0(data_path_raw_files, file_names_original),  # paste0: concatenate strings
        paste0(data_path_raw_files, new_names)  # rename - change to this
    )

}
