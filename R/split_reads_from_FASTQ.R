## 1) Load reads from separate small FASTQ.gz files
## 2) The reads will be split in UMI and BARCODE1 strings
## 3) also a column with sample-ID will be added before 
## 4) export a .tsv file with columns: - UMI  - BARCODE1 - sample-ID

## First function uses Shortread package to extract all reads from the FASTQ files. 
## The fastq files can be zipped (fastq.gz), if doesn`t work have a look at 
## if there is a change in the shortread package and file extensions that can be used.
RAID_load_reads <- function(fastq_file) {
  FASTQ_data <- ShortRead::readFastq(fastq_file)
  reads <- ShortRead::sread(FASTQ_data)
  reads <- as.character(reads)
}

## The ANCHOR defines the common sequence of all antibody barcodes
## The regular expression defines the 10 bp barcode directly before the ANCHOR
## and the 15 bp UMI sequence before the barcode sequence.
## Note: 
## In comparison with the iSeq data processing pipeline: the anchor is shorter
## and there is no Barcode 2 defined behind the ANCHOR.
ANCHOR <- "(ATCAGTCAACAG)"
REGULAR_EXPRESSION_PATTERN <- paste0(
  "[ACTGN]+([ACTGN]{15,15})([ACTGN]{10,10})", ANCHOR,
  "[ACTGN]+")



#' Split reads from FASTQ files
#'
#' First, the iseq_split_reads function uses shortread package function "readFastq" to load all reads from
#' a .fastq or .fastq.gz file. 
#' Then the reads containing an anchor sequence are extracted 
#' and these reads are split to retreive - UMI - Barcode_1 sequences (and stored in a table with 2 columns)
#' Finally a column with the sample_ID is added.
#' @param fastq_file fastq formatted files in a directory path. Default settings or more information from ShortRead::readFastq(fastq_file)
#' @return a tsv file with three columns: UMI , Barcode_1 , sample_ID
#' @seealso \code{\link{ShortRead::readFastq}}
#' @export
#' @examples
#' RAID_split_reads("experiment/data/sample-ID.fastq.gz")

raid_split_reads <- function(fastq_file) {
  ## Create log files in working directory 
  futile.logger::flog.appender(
      futile.logger::appender.file("output/exp_log/experimental_info.log"),
      name = "1")
  futile.logger::flog.layout(futile.logger::layout.tracearg, name = "1")
  futile.logger::flog.appender(
      futile.logger::appender.file("output/exp_log/Run_info.log"), name = "2")
  futile.logger::flog.layout(futile.logger::layout.tracearg, name = "2")

  futile.logger::flog.info("Loading reads from '%s' ...", fastq_file,
                           name = "2")

  reads <- RAID_load_reads(fastq_file)
  ## option to test modifications in script and only use 10 reads for example
  # reads <- reads[1:10]

  futile.logger::flog.info("Split reads from '%s' ...", fastq_file, name = "2")

  anchor_matches <- stringr::str_match(reads, REGULAR_EXPRESSION_PATTERN)


  ## Determine which row numbers contain `NA` (not 100% matched).
  row_index_for_aregexec <- which(is.na(anchor_matches[, 1]))

  # Use the list with row numbers to perform approximate matching on the
  # left sequences. with aregexec the 'rownumbers' an index is create
  # where there is approximate match with max.distance of 2 With
  # max.distance the number of mismatches that are allowed can be set. In
  # our system we allow 10% mismatch.
  approximate_matches_index <- utils::aregexec(
      REGULAR_EXPRESSION_PATTERN, reads[row_index_for_aregexec],
      max.distance = 1)

  # Creates a matrix with 5 columns | Read | UMI | Barcode_1 | anchor |
  # Barcode_2 | of all 'approximated matches' via `regmatches())` that
  # only uses the reads not having 100% match
  # (`reads[row_index_for_aregexec]`) and only matches with an allowed
  # distance of 2 via list of `approximate_matches_index`.
  approximate_matches <- matrix(unlist(regmatches(reads[row_index_for_aregexec],
                                                  approximate_matches_index)),
                               byrow = TRUE, ncol = 4)

    anchor_matches <- rbind(anchor_matches, approximate_matches)

  colnames(
      anchor_matches) <- c("Read", "UMI", "Barcode_1", "anchor")

 

  ## Creates 'clean' table (usable by dplyr) from matched sequences that
  ## contains 3 columns UMI, Barcode_1 (= Antibody) , sample-ID
  ## Use folder above fastq-file as sample-ID (for iSeq)
  #sampleID <- basename(dirname(fastq_file))
  ## use fastq file name as sample-ID
  sampleID <- basename(
      tools::file_path_sans_ext(file_path_sans_ext(fastq_file)))

  data <- dplyr::data_frame(
              anchor_matches[, 2], anchor_matches[, 3]) %>% 
      na.omit()
      
  colnames(data) <- c("UMI", "Barcode_1")

 

  ## use sampleID to save new .tsv file
  save_split_data_name <- paste0("split_", sampleID, ".tsv")
  location_dir <- dirname(fastq_file)

  utils::write.table(data, file = file.path(location_dir, save_split_data_name),
                     sep = "\t", row.names = FALSE, col.names = FALSE)
  ## check if nrow(reads) works properly
  futile.logger::flog.info(
      "Total reads from: %s = %d", fastq_file, length(reads), name = "1")
  futile.logger::flog.info(
      "Total split reads from: %s = %d", fastq_file, nrow(data), name = "1")
  futile.logger::flog.info(
      "Saved data file in %s ...", location_dir, name = "2")
  
remove(reads, anchor_matches)

  data
  

}
