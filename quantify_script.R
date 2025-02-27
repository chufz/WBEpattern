################ QUANTIFICATION ###################################
## cc Carolin Huber, Inga Haalck
## R E M A R K S ## 
# this script is used to build the calibration curves and perform the quantification 
# input: mzmine output, either already narrowed down to "annotated_only" or do that in line 18-21
# remember: "quantify_function.R" needs to be loaded when running this script 
# spiking levels of calibration curve need to be adjusted manually 
###################################################################
# Variables to define
mzmine_input <- "mzmineclean_neg.csv" # csv output from mzmine
cal_pattern <- "Cal" # pattern that distinguishes calibration samples
is_pattern <- "^IS" # pattern that distinguishes internal standards
blank_pattern <- "blank" # pattern that distinguishes blank samples
sample_pattern <- "_24hWW_batch[1-3]_s|24hcomp|Cal|nomatrix" # pattern that distinguishes samples that should be quantified
sample_is_pattern <- "_24hWW_batch[1-3]_s|24hcomp" # pattern that distinguishes samples that should be used to calculate RSD of IS
###################################################################
blank_factor <- 1.5 # factor for blank correction
IS_threshold <- 80 # Detection rate of IS necessary
sample_threshold <- 3 # how often a compound needs to be found to be quantified
IS_rsd_threshold <- 45 # RSD allowed for a IS to be used for quantification
enr_factor <- 1 # enrichment factor that should be calculated in
###################################################################
#Output names to define
quant_output_pdf <- "quantcurves_pos.pdf" # pdf file with calibration curves
quant_csv <- "quantified_pos.csv" # csv file with quantified values
r_square_output <- "r_square_pos.csv" # output table with Rquares
###################################################################
# levels can be used from filenames
# levels <- as.numeric(gsub("ngmL.mzML.height", "", unlist(lapply(strsplit(colnames(cal), "_"), "[[", 6))))
# or defined labels directly in this script:
levels[1] <- 0.1
levels[2] <- 10
levels[3] <- 0.2
levels[4] <- 0.05
levels[5] <- 0.5
levels[6] <- 1
levels[7] <- 2
levels[8] <- 5
###################################################################
# load packages & quantify_function
library(tidyverse)
source("quantify_function.R")

# get the new annotated lists
mzmine <- read.csv(mzmine_input) # mzmine output

# only annotated rows
mzmine <- mzmine[which(mzmine$compound_db_identity.compound_name != ""  ),]
metadata <- data.frame(name=mzmine$compound_db_identity.compound_name, rt=mzmine$rt)

# only peak height table
int_IS <- mzmine[grep(is_pattern,metadata$name),grep("mzML.height", colnames(mzmine))]
metadata_IS <- metadata[grep(is_pattern,metadata$name),]
int <- mzmine[-grep(is_pattern,metadata$name),grep("mzML.height", colnames(mzmine))]
metadata <- metadata[-grep(is_pattern,metadata$name),]

#get the tables ready for quantificaion
cal <- int[,grep(cal_pattern, colnames(int))]
cal_IS <- int_IS[,grep(cal_pattern, colnames(int_IS))]
sample_frame <- int[, grep(sample_pattern, colnames(int))]
blank_frame <- int[, grep(blank_pattern, colnames(int))]
sample_is_frame <- int_IS[, grep(sample_pattern, colnames(int_IS))]

# blank correction 
sample_frame <- blankCorrection(t(sample_frame), t(blank_frame), factor=blank_factor)
sample_frame <- t(sample_frame)

# IS with a DR of > 80%
DR_IS <- DR(sample_is_frame)
sample_frame <- sample_frame[which(DR(cal)> IS_threshold),]
blank_frame <- blank_frame[which(DR(cal)> IS_threshold),]
metadata <- metadata[which(DR(cal)> IS_threshold),]
cal <- cal[which(DR(cal)> IS_threshold),]

# compounds that are found in at least 3 samples
metadata <- metadata[which(DR(sample_frame)> sample_threshold),]
cal <- cal[which(DR(sample_frame)> sample_threshold),]
blank_frame <- blank_frame[which(DR(sample_frame)> sample_threshold),]
sample_frame <- sample_frame[which(DR(sample_frame)> sample_threshold),]

# IS deviation 
matrix_is_frame <- int_IS[, grep(sample_is_pattern, colnames(int_IS))]
RSD_IS <- apply(matrix_is_frame, 1, RSD)

# only is with a deviation of RSD < 45%
sample_is_frame <- sample_is_frame[which(RSD_IS < IS_rsd_threshold),]
cal_IS <- cal_IS[which(RSD_IS < IS_rsd_threshold),]
metadata_IS <- metadata_IS[which(RSD_IS < IS_rsd_threshold),]

# find best internal standard for each target
IS_assignment  <- findBestIS(t(sample_frame), t(sample_is_frame), metadata$rt, metadata_IS$rt)
################################################################################
# get list of IS used 
metadata_IS_reset <- metadata_IS
rownames(metadata_IS_reset) <- 1:nrow(metadata_IS_reset)

IS_assigned <- cbind(metadata$name, IS_assignment)

# Add a column "IS_number" that goes from 1 to 24 in `metadata_IS_reset`
metadata_IS_reset$IS_number <- 1:nrow(metadata_IS_reset)

IS_assigned_names <- merge(IS_assigned, metadata_IS_reset, 
                                by.x = "IS_assignment", by.y = "IS_number", 
                                all.x = TRUE)

IS_assigned_names <- IS_assigned_names %>%
  select(V1, name)
################################################################################
# quantify
quant <- quantify(t(sample_frame), t(sample_is_frame), t(cal), levels, t(cal_IS), IS_assignment, enr_factor=enr_factor,compounds=metadata$name,
                  output_pdf=quant_output_pdf)
q <- data.frame(quant)
colnames(q) <- metadata$name
q <- data.frame(cbind(rownames(quant), q))
write_csv(q, quant_csv)

# run the linear function to get summary & r square values 
summary <- linear(t(sample_frame), t(sample_is_frame), t(cal), levels, t(cal_IS), IS_assignment, enr_factor=1,compounds=metadata$name)
r_square <- data.frame(summary)

# get compound names, remove the header and add the compounds as a new column 
col_names <- colnames(q)
col_names <- col_names[-1]
r_square$compound <- col_names

# write csv file for further evaluation 
write.csv(r_square, file = r_square_output)








