library(VennDiagram)
library(dplyr)

# Clean up
rm(list = ls())


# Load data
df <- read.table("../input/s06.plot_Venn/500outliers.comm_Dessert.merged.bed.genes.txt", 
                 header = F)
names(df) <- c("GeneID", "APP")

# Example dataframe
# df <- data.frame(GeneID = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6"),
#                  APP = c("APP1", "APP2", "APP1", "APP3", "APP2", "APP3"))

# Group by APP and create a list of GeneIDs for each APP
app_gene_list <- df %>% 
  group_by(APP) %>% 
  summarise(genes = list(GeneID)) %>% 
  pull(genes)

# Assign names to the list based on APP names
names(app_gene_list) <- unique(df$APP)


# Check the number of APPs
num_apps <- length(app_gene_list)

# Draw the Venn diagram based on the number of APPs
if (num_apps == 2) {
  venn.plot <- venn.diagram(
    x = app_gene_list,
    category.names = names(app_gene_list),
    filename = NULL,
    fill = c("blue", "red")
  )
} else if (num_apps == 3) {
  venn.plot <- venn.diagram(
    x = app_gene_list,
    category.names = names(app_gene_list),
    filename = NULL,
    fill = c("blue", "red", "green")
  )
} else if (num_apps == 4) {
  venn.plot <- venn.diagram(
    x = app_gene_list,
    category.names = names(app_gene_list),
    filename = NULL,
    fill = c("blue", "red", "green", "yellow")
  )
} else {
  stop("Too many APP groups for Venn Diagram. Please limit to 2-4 groups.")
}

# Plot the Venn diagram
grid.draw(venn.plot)


grid.text(
  "Venn Diagram of Genes Matched by Top 500 Outliers comm_Dessert Positive Selection",  # Title label
  x = 0.5,  # Position in x-axis (0.5 is the center)
  y = 0.95,  # Position in y-axis (near the top)
  gp = gpar(fontsize = 12, fontface = "bold")  # Text properties like fontsize and boldness
)
