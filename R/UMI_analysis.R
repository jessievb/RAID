#' # Function for UMI analysis from iSeq split data table.

# Each row comes from a sequence from the FASTQ file that contained a
# common
#'  sequence. The first question is: how many PCR duplicates can be found?
#'  A PCR duplicate will result a duplicate UMI (and barcodes) sequences.
#'  Practically this means if a row is a similar to another one, this is a PCR
#'  duplication.
#'
#'  # Count the duplicate UMI in the whole file & filter duplicates.

#' Three functions are defined:
#' 1. import the split reads into the R environment
#' 2. count the occurence of a specific row (which contains 1 read)
#' 3. count the frequency of the occurences (direct usable for UMI frequency
#'    plot)
#'
#' Function 2 and 3 each save a TSV file in the output/data/ folder .

#' @export
raid_import_split_data <- function(
    split_data_file_path) {
  data.table::fread(
      input = split_data_file_path, sep = "\t", stringsAsFactors = FALSE,
      colClasses = c("character", "character"),
      header = FALSE,
      col.names = c("UMI", "Barcode_1"))
}

## ~12 minutes for ~50M observations
## If iSeq_split_data data.frame is not imported in environment set
## import to TRUE and specify split_data_file_path
#' @export
raid_umi_count <- function(split_data_file_path,
    raid_barcode_tbl = raid_split_data, output_filename = "umi_count_",
    output_dir = "output/data/") {
  # using dplyr package to easily summarize the counts for each row.
  raid_umi_count <- dplyr::tbl_df(raid_barcode_tbl) %>%
      dplyr::group_by(UMI, Barcode_1, add = TRUE) %>%
      dplyr::summarize(row_occurence_count = n())

    sampleID <- basename(
      tools::file_path_sans_ext(file_path_sans_ext(split_data_file_path)))

  output_filename <- utils::write.table(
      raid_umi_count,
      file = file.path(output_dir, paste0(output_filename, sampleID, ".tsv")), sep = "\t",
      row.names = FALSE, col.names = TRUE)

  raid_umi_count
}

#' @export
raid_umi_count_frequency <- function(split_data_file_path,
    UMI_count_input = UMI_count,
    output_filename = "umi_count_frequency_",
    output_dir = "output/data/") {
  UMI_count_frequency <- dplyr::tbl_df(UMI_count_input) %>%
      dplyr::group_by(row_occurence_count) %>%
      dplyr::summarize(count_frequency = n())

  sampleID <- basename(
      tools::file_path_sans_ext(file_path_sans_ext(split_data_file_path)))
  
  # TODO: use data.table::fwrite
  # (https://github.com/Rdatatable/data.table/issues/580)
  utils::write.table(
      UMI_count_frequency, file = file.path(output_dir, paste0(output_filename, sampleID, ".tsv")), sep = "\t",
      row.names = FALSE, col.names = TRUE)

  UMI_count_frequency
}
