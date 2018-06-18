#' Count the number of unique UMI per Barcode_1 and Barcode_2 columns.

#' The data table that should be used is the iSeq_UMI_count table, which
#' contains all the unique (collapsed) UMI,
#' barcode_1 and barcode_2 information. The column row_occurence_count can be
#' ignored.
#' @export
raid_barcode_count <- function(split_data_file_path,
    UMI_count_tbl = UMI_count, output_filename = "output/data/barcode_count_") {
  barcode_count <- UMI_count_tbl %>%
      dplyr::group_by(Barcode_1) %>%
      dplyr::summarize(antibody_count = n())


  sampleID <- basename(
      tools::file_path_sans_ext(file_path_sans_ext(split_data_file_path)))

  utils::write.table(barcode_count, file = paste0(output_filename, sampleID, ".tsv"),
                     sep = "\t", row.names = FALSE, col.names = TRUE)

  barcode_count
}

#' Add antibody and sample specific information to the table with barcode
#' counts.
#'
#' The function makes use of the dplyr package 'left_join' option. If there is
#' a direct match between shared columns of data table and the input tables
#' antibody_barcode_index and well_barcode_index.
#' @export
raid_barcode_match <- function(split_data_file_path,
    barcode_count_tbl = barcode_count,
    output_filename = "output/data/barcode_count_matched_") {
 
  antibody_barcode_tbl <- tbl_df(
      read.table("config/antibody_barcode_index.txt", header = TRUE,
                 stringsAsFactors = FALSE))
  
  barcode_count_tbl$Barcode_1 <- as.character(barcode_count_tbl$Barcode_1)
  

  barcode_count_matched <- barcode_count_tbl %>%
      dplyr::left_join(antibody_barcode_tbl, copy = TRUE) 


  sampleID <- basename(
      tools::file_path_sans_ext(file_path_sans_ext(split_data_file_path)))

  utils::write.table(
      barcode_count_matched, file = paste0(output_filename, sampleID, ".tsv"), sep = "\t",
      row.names = FALSE, col.names = TRUE)

  barcode_count_matched

}

#' @export
raid_barcode_match_na <- function(
    barcode_count_matched_tbl = barcode_count_matched) {
  barcode_count_NA <- barcode_count_matched_tbl %>%
      dplyr::filter(is.na(Ab_barcode_nr))
}

#' @export
raid_barcode_match_filtered <- function(
    barcode_count_matched_tbl = barcode_count_matched) {
  barcode_count_matched_filtered <- barcode_count_matched_tbl %>%
      dplyr::filter(!is.na(Ab_barcode_nr))
}
