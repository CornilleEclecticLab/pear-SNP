suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
  library(data.table)
  # library(parallel)
})

rm(list=ls())

# Set working directory to the location of the script (if running interactively)
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
}

# Set parameters ----------------------------------------------------------
tajimad_dir <- "../Tajima_D_link.NOGIT./s02.split_pixy_by_pop"
omega_dir <- "../output/s06.mask_results"

tajimad_cutoff_dir   <- "../../run_pixy_TajimaD/output/s04.catch_outliers/"
omegaplus_cutoff_dir <- "../../run_OmegaPlus/output/s04.catch_outliers/"

tajimad_cutoff_dir2   <- "../../run_pixy_TajimaD_refer_on_communis/output/s04.catch_outliers/"
omegaplus_cutoff_dir2 <- "../../run_OmegaPlus_refer_on_communis/output/s04.catch_outliers/"


# Check if the data directory exists
if (!dir.exists(tajimad_dir)) {
  stop(paste("Tajima_D data directory", tajimad_dir, "not found."))
}
if (!dir.exists(omega_dir)) {
  stop(paste("OmegaPlus data directory", omega_dir, "not found."))
}

# Initialize parallel processing
num_cores <- min(parallel::detectCores() - 1, 4)
options(mc.cores = num_cores)

# Convert Chr ID ----------------------------------------------------------
# Read the chromosome to ID mapping file
chr_map_file <- "../input/chromosome_accession_ID_map.txt"
chr_map <- read.table(chr_map_file, header = FALSE, col.names = c("Accession", "ChrID"))

convert_chr_to_id <- function(data, chr_map) {
  # # merge chr ID map with data
  # data_with_id <- merge(data, chr_map, by.x = "Chr", by.y = "Accession", all.x = TRUE)
  
  # # If no mapped value, use original
  # data_with_id$Chr <- ifelse(is.na(data_with_id$ChrID), data_with_id$Chr, data_with_id$ChrID)
  
  # # rm ChrID column
  # data_with_id <- data_with_id[, !names(data_with_id) %in% c("ChrID")]
  #  return(data_with_id)
  data <- merge(data, chr_map, by.x = "Chr", by.y = "Accession", all.x = TRUE)
  data[, Chr := ifelse(is.na(ChrID), Chr, ChrID)]
  data[, ChrID := NULL]  # Ensure this column is removed after use
  return(data)
}



# Load data ---------------------------------------------------------------
#### Load Tajima_D data ####
tajimad_file_list <- list.files(path = tajimad_dir, pattern = "^Pixy\\..*\\.1k_tajima_d\\.txt", full.names = TRUE)

tajimad_data <- tajimad_file_list %>%
  lapply(function(file_path) {
    # Read the data
    print(paste0("Reading file: ", file_path))
    data <- fread(file_path, header = TRUE, sep = "\t", fill = TRUE)
    # Header:
      # pop     chromosome      window_pos_1    window_pos_2    tajima_d 
      # no_sites        raw_pi  raw_watterson_theta     tajima_d_stdev

    # Check the number of columns and set column names accordingly
    if (ncol(data) == 9) {
      colnames(data) <- c("Population", "Chr", "Start", "End", "Score", "no_sites", "raw_pi", "raw_watterson_theta", "tajima_d_stdev")
      } else {
      stop(paste("Unexpected number of columns in file:", file_path))
    }
    data$Pos <- as.integer(0.5 * (data$Start + data$End))  # Calculate the midpoint of the window
    data$Pos <- as.integer(data$Pos)
    
    # Filter out lines that are empty or start with "//"
    data <- data[!is.na(data$Score), ]  # Remove rows with NA in tajima_d

    # data$Population <- strsplit(basename(file_path), "\\.")[[1]][2]
    # data$Chr <- strsplit(basename(file_path), "\\.")[[1]][3]
    
    data <- convert_chr_to_id(data, chr_map)

    # Set the chromosome level
    data$Chr <- factor(data$Chr, levels = paste0("Chr", 1:17))
    
    # Set Method
    data$Method <- "Pixy_TajimaD"
    
    return(data)
  }) %>%
  bind_rows()

tajimad_data <- tajimad_data %>%
  drop_na(Score) %>%
  select(Chr, Start, End, Pos, Method, Population, Score)



#### Load Omega Data ####
omega_file_list <- list.files(path = omega_dir, pattern = "^OmegaPlus_Report.", full.names = TRUE)

omega_data <- omega_file_list %>%
  lapply(function(file_path) {
    
    # Read the data
    print(paste0("Reading file: ", file_path))
    data <- fread(file_path, header = FALSE, sep = "\t", fill = TRUE)
    
    # Check the number of columns and set column names accordingly
    if (ncol(data) == 5) {
      colnames(data) <- c("Pos", "Score", "Start", "End", "Valid")
    } else if (ncol(data) == 2) {
      colnames(data) <- c("Pos", "Score")
    } else {
      stop(paste("Unexpected number of columns in file:", file_path))
    }
    
    # Filter out lines that are empty or start with "//"
    data <- data[!(data$Pos == "" | grepl("^//", data$Pos)), ]
    data$Pos <- as.integer(as.numeric(data$Pos))
    
    data$Population <- strsplit(basename(file_path), "\\.")[[1]][2]
    data$Chr <- strsplit(basename(file_path), "\\.")[[1]][3]
    
    data <- convert_chr_to_id(data, chr_map)
    
    # Set the chromosome level
    data$Chr <- factor(data$Chr, levels = paste0("Chr", 1:17))
    
    # Set Method
    data$Method <- "OmegaPlus"
    
    return(data)
  }) %>%
  bind_rows()

# Remove duplicated ranges
omega_data <- omega_data %>%
  mutate(score_start_end = paste0(Score, "-", Start, "-", End)) %>%
  filter(!duplicated(score_start_end))

# Purge data, remove invalid score
if ("Valid" %in% colnames(omega_data)) {
  omega_data <- omega_data %>%
    filter(Valid != 0, !is.na(Valid))
}

omega_data <- omega_data %>%
  drop_na(Score) %>%
  select(Chr, Start, End, Pos, Method, Population, Score)

combined_data <-bind_rows(tajimad_data, omega_data)

head(combined_data)
dim(combined_data)


##########  TEST  #########################################################
## /!\ Randoly select 5000 observations
# set.seed(123)
# omega_data_sample <- omega_data %>%
#   sample_n(5000)

# tajimad_data_sample <- tajimad_data %>%
#   sample_n(5000)
###########################################################################



#### Cumulate position ####

combined_data <- combined_data %>%
  # Get chromosome size
  group_by(Chr) %>%
  summarise(chr_len = max(Pos)) %>%
  
  # Calulate cumulative position
  mutate(Tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  
  # Add this information to data
  left_join(combined_data, ., by = c("Chr" = "Chr")) %>%
  arrange(Chr, Pos) %>%
  mutate(Cum_pos = Pos+ Tot)


#### Calculate boundry and midpoints of chromosome for labeling ####
chr_label <- combined_data %>%
  group_by(Chr) %>%
  summarize(
    end = max(Cum_pos)  # The maximum cumulative position for this chromosome
  ) %>%
  mutate(
    start = lag(end, default = 1),  # Start is the last max Cum_pos (default 1 for the first chromosome)
    midpoints = (start + end) / 2   # Calculate the midpoint for labeling
  )


#### Prune data ####
prune_data <- function(df, method="Pixy_TajimaD", window_size=10) { # Keep one mark for every window_size (3) marks
  df <- as.data.table(df)
  df_method<- df[df$Method == method, ]
  df_other <- df[df$Method != method, ]
  
  if (nrow(df_method) == 0) {return(df)}

  df_method[, window_id := ceiling(seq_len(.N) / window_size), by = Population]
  
  pruned_df <- df_method[, .SD[which.max(Score)], by = .(Population, window_id)]
  
  final_df <- rbind(pruned_df[, !"window_id"], df_other, fill = TRUE)
  
  return(final_df)
}

## Pruning
print(paste0("Pruning data from  ", nrow(combined_data), "..."))
# combined_data <- prune_data(combined_data, "Pixy_TajimaD", 20)
print(paste0("Done, purning data, left ", nrow(combined_data), "..." ))


populations <- na.omit(unique(combined_data$Population))

# Initialize
cutoff_df <- data.frame(Population = character(),  tajimad = numeric(),  omegaplus = numeric(),  stringsAsFactors = FALSE)

# Loop through populations and read cutoff values
for (pop in populations) {
  print(paste0("Reading cutoff for ", pop))
  # Generate file paths for each method
  tajimad_file <- file.path(tajimad_cutoff_dir, paste0(pop, ".Pixy_TajimaD.0_0040.cutoff.txt"))
  omegaplus_file <- paste0(omegaplus_cutoff_dir, pop, ".OmegaPlus.0_0055.cutoff.txt")
  tajimad_file2 <- paste0(tajimad_cutoff_dir2, pop, ".Pixy_TajimaD.0_0040.cutoff.txt")
  omegaplus_file2 <- paste0(omegaplus_cutoff_dir2, pop, ".OmegaPlus.0_0055.cutoff.txt")
  
  # Initialize cutoff values to NA
  tajimad_cutoff <- NA
  omegaplus_cutoff <- NA
  
  # Try to read Tajima_D cutoff value
  if (file.exists(tajimad_file)) {
    tajimad_cutoff <- read.table(tajimad_file, header = FALSE, skip = 1)[1, 1]
  } else if (file.exists(tajimad_file2)) {
    tajimad_cutoff <- read.table(tajimad_file2, header = FALSE, skip = 1)[1, 1]
  }
  
  # Try to read OmegaPlus cutoff value
  if (file.exists(omegaplus_file)) {
    omegaplus_cutoff <- read.table(omegaplus_file, header = FALSE, skip = 1)[1, 1]
  } else if (file.exists(omegaplus_file2)) {
    omegaplus_cutoff <- read.table(omegaplus_file2, header = FALSE, skip = 1)[1, 1]
  }
  
  # Add to the dataframe
  cutoff_df <- rbind(
    cutoff_df,
    data.frame(
      Population = pop,
      tajimad = tajimad_cutoff,
      omegaplus = omegaplus_cutoff
    )
  )
}


# Add cutoff values to the combined data
combined_data <- combined_data %>%
  left_join(cutoff_df, by = "Population") %>%
  mutate(
    Cutoff = case_when(
      Method == "Pixy_TajimaD" ~ tajimad,
      Method == "OmegaPlus" ~ omegaplus
    ),
    Color_group = case_when(
      Method == "Pixy_TajimaD"  & Score < Cutoff ~ Population,
      Method == "OmegaPlus" & Score > Cutoff ~ Population,
      TRUE ~ "Below_Cutoff"
    ),
    Color_group = factor(Color_group, levels = c(unique(Population), "Below_Cutoff"))
  )



# Save data ---------------------------------------------------------------
binary_cache_file <- "combined_data_Omega_TajimaD_binary.rds"

save(combined_data, cutoff_df, file = binary_cache_file)

print("Data saved to binary file")
q()
