# Convert pixy output to long format
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()

  for (i in 1:length(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])

    df <- read_delim(pixy_files[i], delim = "\t")

    if (stat_file_type == "pi") {
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
    } else {
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
    }

    pixy_df[[i]] <- df
  }

  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

pixy_folder <- "../output/s02.run_pixy_by_chromosome.mask_input"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)