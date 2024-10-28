pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  
  for (i in 1:length(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])  # Extract file type
    
    df <- read_delim(pixy_files[i], delim = "\t")
    
    if (stat_file_type == "pi") {
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)  # Add pop2 as NA for pi files
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
