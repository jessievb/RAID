# General Description

The R package *RAID* contains functions that are used to process
sequencing reads from the RAID technology. The RAID package requires FASTQ files from
individual cells as input.

The package allows to:

1. Import reads from FASTQ files and extracts the following
nucleotidesequences from the
reads:
    - Unique Molecular Identifier (UMI)
    - Antibody barcode (BARCODE1)
2. Create a data table with two columns of DNA sequences (UMI, BARCODE1)
3. Report the frequency of duplicate (or more) UMI DNA sequences.
4. Create count table (Unique UMI) per antibody per cell.
5. Log the run and experiment information.

# Install and load RAID package in R

## Quick and easy install directly from GitHub

The 'devtools' package allows installing the RAID package from GitHub
directly.

```r
# First start R and install and load the devtools package
install.packages('devtools')
library(devtools)

# Then install the RAID package from the repository jessievb/RAID
devtools::install_github("jessievb/RAID")

# Finally load the package
library(RAID)
```
 ## Optional: install package after manual download

 #### Download package from GitHub

- Windows: [Direct
download](https://github.com/jessievb/RAID/archive/master.zip)
from github

- Linux: Download via command line and git. If git is not installed,
download it
[here](https://git-scm.com/download/linux). Then run the following
command: `git
clone git@github.com:jessievb/RAID.git`

```
# move into the folder with R packages (any folder you like)
cd ~/my-R-packages/

# download the RAID package using git
git clone git@github.com:jessievb/RAID.git
```


#### Install/deploy RAID package in R

'Devtools' package also allows you to install a package from a local
folder. Extra information on devtools package you can find [here](devtoolsinfo).

- Start R
- Make sure working directory = the package directory
- Run documentation and installation of the RAID package (only once needed)

```
setwd(~/my-R-packages/RAID/) # Make sure working directory = folder with
R-package
devtools::document() # to create documentation
devtools::install() # install package
```

```r
# Finally load the package
library(RAID)
```

# Functions description

### (I) raid_split_reads()

RAID_split_reads function uses ShortRead package to extract all reads from
(zipped)
FASTQ files.

1. Load reads from .FASTQ files
2. (Approximate) matches reads to anchor sequence (see note)
3. Then, the UMI sequence and Barcode 1 are extracted from the reads via a
regular
expression
4. Export: "splitreads.tsv": table with - UMI  - Barcode 1
5. Run information is logged to output/exp_log/Run_info.log


> NOTE: The ANCHOR defines the common sequence of all antibody barcodes.
The regular
expression defines the 10 bp barcode directly before the ANCHOR and the 15
bp UMI
sequence before the barcode sequence. In comparison with
[ID-seq](https://github.com/jessievb/IDseq): the anchor is shorter and
there is no Barcode 2 defined behind the ANCHOR.  
ANCHOR <- "(ATCAGTCAACAG)"  
REGULAR_EXPRESSION_PATTERN <- paste0( 
"[ACTGN]+([ACTGN]{15,15})([ACTGN]{10,10})",
ANCHOR, "[ACTGN]+")   



The file 'split_reads_from_FASTQ.R' can be modified according to the
design of
the barcode.


### (II) raid_analysis_splitreads()

This functions combines the following functions in order:
 1. raid_import_split_data()
 2. raid_umi_count()
 3. raid_umi_count_frequency()
 4. raid_barcode_count()
Also, it adds run and experiment information to the .log files in the
output/exp_log/ folder. Finally it creates table with counts, and a table
with
UMI-duplicate rates in the output/data/ folder

# Example workflow

## 1. Create the following experiment directory (via command line in linux)
```sh
mkdir -p /home/Experiment_ID/{data,config,output/{data,exp_log,figures}}
```

## 2. Split reads


### In R (or R-studio)

```r
library(RAID)
setwd(~/Experiment_ID/)
raid_split_reads(fastq_file="data/sample_1/sample_name.fastq.gz")
```


### From command line
Start any number of processes, depending on the number of fastq.gz files
to process:

```sh
# check if the command works. Should print all FASTQ filenames found in the
indicated. Indicate behind -P how many cores should be used.
find Path_to_folders_with_fastqfiles/*/sample_name*__R1_0*.fastq.gz -name
"*.fastq.gz" | xargs -P 4 -i -- echo "'{}'"

find Path_to_folders_with_fastqfiles/*/sample_name*__R1_0*.fastq.gz -name
"*.fastq.gz" | xargs -P 4 -i -- R -e 'library(RAID);
setwd("~/Experiment_ID");
system.time(raid_split_reads(fastq_file="'{}'")); quit(save="no")'

```

## 3. Analyse split reads (count and add antibody information to counts)


> Note: Make sure the config directory contains an antibody_barcode_index.txt
The antibody_barcode_index.txt file should contain the antibody specific
nucleotide sequences a barcode number and the name of the Antibody.  


Example:  

| Barcode_1  | Ab_barcode_nr | Ab_name        |
|------------|---------------|----------------|
| GTATCGTCGT | 1             | Name_Antibody1 |
| GTGTATGCGT | 2             | Name_Antibody2 |
| …          |               |                |          

From the command line multiple files can be analysed in parallel:

```sh

find Path_to_folders_with.splitfiles/split*.tsv -name "*.tsv" | nice -15
xargs -P 5
-i -- echo "'{}'"

find Path_to_folders_with.splitfiles/split*.tsv -name "*.tsv" | nice -15
xargs -P 12
-i R -e 'library(RAID); setwd("~/Experiment_ID");
system.time(RAID::raid_analysis_splitreads(experiment_dir="~/Experiment_ID",
data="'{}'")); quit(save="no")'

```


## 4. Create annotated single count matrix with antibody barcode counts per cell

As for each cell 1 file was generated, the data needs to be merged. One can use the following commands in R to obtain such count table:

```r

library(dplyr)
# replace "folder/" with directory where barcode_count_matched_split_files
are
located (twice in this script)

# get all file names starting with "barcode_count_matched_split_"
file_names_single_cells <- list.files(path =
"~/Experiment_ID/output/data/", pattern
="^barcode_count_matched_split_")

# create empty data frame with column of all Ab_barcode_numbers. Replace
c(1:100,"x") with a list of other numbers if needed
antibody_counts_all_cells <- data.frame(Ab_barcode_nr =c(1:100,"x"))

# For loop to read file, adjust column names and appende to
antibody_counts_all_cell
dataframe
for (single_cell_number in c(1:length(file_names_single_cells))){

  # import table and add correct colnames
  single_cell_table <-
read.table(paste0("~/Experiment_ID/output/data/",file_names_single_cells[single_cell_number]))
  colnames(single_cell_table) <- c("Barcode_1", "antibody_count",
"Ab_barcode_nr","Ab_name")

  # get string with single-cell id from the filename
  single_cell_ID <- gsub(pattern = "barcode_count_matched_split_",
replacement =
"",x = file_names_single_cells[single_cell_number])
  single_cell_ID <-gsub(pattern = ".tsv", replacement = "", x =
single_cell_ID)

  # select columns from imported table
  single_cell_table <- single_cell_table %>%
select(Ab_barcode_nr,antibody_count)

  # replace "antibody_count" column name with single-cell ID name
  colnames(single_cell_table)  <- c("Ab_barcode_nr",single_cell_ID )

  # join the antibody_counts_all_cells with the importedtable
  antibody_counts_all_cells <- antibody_counts_all_cells %>%
left_join(single_cell_table, by = "Ab_barcode_nr")

}

# save csv file
write.csv(antibody_counts_all_cells, file =
"~/Experiment_ID/data/expID_counts.csv",
col.names = TRUE)


q()

```
