library(VennDiagram)

# Set working directory to the location of the script and clean up
if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list=ls())
}



# Define file paths
file1 <- "../output/s05.parse_MCScanX/comm_pyri.parsed.collinearity.txt"
file2 <- "../output/s05.parse_MCScanX/pyri_comm.parsed.collinearity.txt"


# Read files as data frames
df1 <- read.table(file1, header = FALSE, stringsAsFactors = FALSE)
df2 <- read.table(file2, header = FALSE, stringsAsFactors = FALSE)


# Ensure correct structure
colnames(df1) <- c("Gene1", "Gene2")
colnames(df2) <- c("Gene1", "Gene2")

# Create unique gene pair identifiers
pairs1 <- apply(df1, 1, function(x) paste(x[1], x[2], sep = "-"))
pairs2 <- apply(df2, 1, function(x) paste(x[1], x[2], sep = "-"))


# Compute set sizes and intersection
venn.plot <- draw.pairwise.venn(
  area1 = length(pairs1),
  area2 = length(pairs2),
  cross.area = length(intersect(pairs1, pairs2)),
  category = c(expression("Target" ~ italic("P. pyrifolia")), expression("Target" ~ italic("P. communis"))),
  fill = c("blue", "red"),
  alpha = 0.5,
  cat.cex = 0.8
)

# Display the Venn diagram
grid.newpage()
grid.draw(venn.plot)
