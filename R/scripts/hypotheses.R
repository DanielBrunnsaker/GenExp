
library(dplyr)
library(corrplot)
library(gplots)
library(UpSetR)
library(ggplot2)
library(grid)
library(gridExtra)
library(patchwork)  # or use cowplot if you prefer

# Set the directory containing the CSV files
folder_path <- "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/results/coefficients/20241114"
# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
# Initialize an empty list to store the means
average_list <- list()
# Loop through each file
for (file in csv_files) {
  # Read the CSV file
  data <- read.csv(file, row.names = 1)
  # Calculate the average of each column
  column_means <- rowMeans(data, na.rm = TRUE)
  # Use the file name (without path and extension) as the column name
  col_name <- tools::file_path_sans_ext(basename(file))
  # Store the column means in the list with the file name as the key
  average_list[[strsplit(col_name, '_')[[1]][1]]] <- column_means
}
# Combine all column means into a single data frame
average_data <- as.data.frame(average_list)
average_data <- subset(average_data, select = -c(tryptophan, glycine,phenylalanine))

coefs = average_data
# Let nTargets = ncol(coefs)
nTargets <- ncol(coefs)
non_zero <- coefs != 0

overlap_info <- data.frame(
  i = integer(),
  j = integer(),
  n_i = integer(),
  n_j = integer(),
  overlap = integer(),
  negneg = integer(),
  pospos = integer(),
  diffSign = integer(),
  stringsAsFactors = FALSE
)

# 1. Suppose your real target order is simply colnames(coefs):
ordered_targets <- colnames(coefs)

# 2. Convert target_i and target_j into factors with the desired ordering:
overlap_info$target_i <- factor(overlap_info$target_i, levels = ordered_targets)
overlap_info$target_j <- factor(overlap_info$target_j, levels = ordered_targets)


for(i in seq_len(nTargets)){
  for(j in seq_len(nTargets)){
    nz_i <- non_zero[, i]
    nz_j <- non_zero[, j]
    
    overlap_ij <- sum(nz_i & nz_j)
    both_neg   <- sum(coefs[, i] < 0 & coefs[, j] < 0)
    both_pos   <- sum(coefs[, i] > 0 & coefs[, j] > 0)
    diff_sign  <- sum((coefs[, i] > 0 & coefs[, j] < 0) | (coefs[, i] < 0 & coefs[, j] > 0))
    
    overlap_info <- rbind(
      overlap_info,
      data.frame(
        i        = i,
        j        = j,
        n_i      = sum(nz_i),
        n_j      = sum(nz_j),
        overlap  = overlap_ij,
        negneg   = both_neg,
        pospos   = both_pos,
        diffSign = diff_sign
      )
    )
  }
}

# Attach names of the targets for clearer facet labels
overlap_info$target_i <- colnames(coefs)[overlap_info$i]
overlap_info$target_j <- colnames(coefs)[overlap_info$j]

df_single <- data.frame()

for(k in seq_len(nrow(overlap_info))){
  
  rowk <- overlap_info[k, ]
  
  name_i   <- rowk$target_i
  name_j   <- rowk$target_j
  
  n_i      <- rowk$n_i
  overlap  <- rowk$overlap
  negneg   <- rowk$negneg
  pospos   <- rowk$pospos
  diffSign <- rowk$diffSign
  
  # If target i selected zero features, the entire donut is conceptually notShared.
  # We can just store "notShared = 1" (which becomes invisible).
  if(n_i == 0){
    df_single <- rbind(
      df_single,
      data.frame(target_i=name_i, target_j=name_j, slice="notShared", amount=1)
    )
    next
  }
  
  # fraction of i's features that are shared
  frac_shared <- overlap / n_i
  
  # Among overlap, subdivide by sign
  total_sign <- pospos + negneg + diffSign
  if(total_sign == 0){
    # no sign overlap => no pospos/negneg/diffSign
    pospos_frac   <- 0
    negneg_frac   <- 0
    diffSign_frac <- 0
  } else {
    pospos_frac   <- frac_shared * (pospos   / total_sign)
    negneg_frac   <- frac_shared * (negneg   / total_sign)
    diffSign_frac <- frac_shared * (diffSign / total_sign)
  }
  
  not_shared_frac <- 1 - (pospos_frac + negneg_frac + diffSign_frac)
  
  # Add one row per slice
  df_single <- rbind(
    df_single,
    data.frame(target_i=name_i, target_j=name_j, slice="pospos",   amount=pospos_frac),
    data.frame(target_i=name_i, target_j=name_j, slice="negneg",   amount=negneg_frac),
    data.frame(target_i=name_i, target_j=name_j, slice="diffSign", amount=diffSign_frac),
    data.frame(target_i=name_i, target_j=name_j, slice="notShared",amount=not_shared_frac)
  )
}




p_main <- ggplot(df_single, aes(x0=0, y0=0, r0=0.3, r=1.0, amount=amount, fill=slice)) +
  geom_arc_bar(stat="pie", color=NA) +
  # Switch the row names to the left side:
  facet_grid(target_i ~ target_j, switch="y") +
  coord_fixed() +
  theme_void() +
  # Adjust facet labels to avoid overlapping:
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    strip.text.x = element_text(size=5, angle=45, hjust=1),  # Increase bottom margin of facet labels
    strip.placement = "outside",  # Keep facet labels outside the plotting area
    strip.text.y = element_text(size=5),
    #plot.margin = margin(t = 50, r = 10, b = 10, l = 10),  # Add more space at the top of the plot
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  scale_fill_manual(
    values = c(
      "pospos"    = "#1b9e77",
      "negneg"    = "#d95f02",
      "diffSign"  = "#7570b3",
      "notShared" = "#00000000"  # invisible
    )
  )

df_bar <- overlap_info %>%
  select(target_i, i, n_i) %>%
  distinct()  # one row per target_i

p_bar <- ggplot(df_bar, aes(y=rev(target_i), x=n_i)) +
  geom_col(fill="darkgrey", color = 'black') +
  theme_minimal() +
  #labs(x="# Non-Zero Features", y=NULL) +
  # Flip coords so that bars go horizontally, and the row order goes top->down
  #coord_flip() +
  # Make sure the row order matches the facet grid: "Target_1" on top:
  scale_y_discrete(limits=rev(levels(df_bar$target_i))) +
  theme(
    # Remove extra grid lines if you want less clutter:
    panel.grid.major.y = element_blank(),
    # Or keep them if you find them helpful:
    # panel.grid.major.x = element_line(color="grey70"),
    axis.text.y = element_blank(),   # we won't label each bar, since the row label is in the main plot
    axis.ticks.y = element_blank()
  )

p_combined <- p_main + p_bar + 
  plot_layout(ncol=2, widths = c(4, 1))

p_combined


ggsave("/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/R/plots/plot.pdf", plot = p_combined, device = "pdf", width = 8, height = 6, units = "in")




