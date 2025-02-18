


gdata = read.csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102/results/data/lactic_acid_smoothed.txt', sep = '\t')
layout = read.csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102/protocol/plate_layout/layout.tsv', sep = '\t')
layout$well <- gsub("([A-Z])0+(\\d+)", "\\1\\2", layout$well)

merged_df <- merge(layout[c('well','Summary')],gdata, by.x = "well", by.y = 'Well')


df = merged_df

# Filter out rows where the Summary column is empty or missing
df <- df %>%
  filter(!is.na(Summary) & Summary != "")


# Ensure columns are selected properly by specifying numeric indices
averaged_data <- df %>%
  group_by(Summary) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))

# Reshape data for plotting
averaged_data_long <- averaged_data %>%
  pivot_longer(cols = -Summary, names_to = "Timepoint", values_to = "Average") %>%
  mutate(Timepoint = as.numeric(gsub("X", "", Timepoint)))  # Convert timepoints to numeric

# Plot the time series of averaged replicates
ggplot(averaged_data_long, aes(x = Timepoint, y = Average, color = Summary, group = Summary)) +
  geom_line(size = 1) +
  labs(
    title = "Time Series of Averaged Replicates",
    x = "Timepoint",
    y = "Average Value",
    color = "Group"
  ) +
  theme_minimal()



# Generate 3x3 versions of the dataset with added noise
set.seed(123)  # Ensure reproducibility
replicates <- bind_rows(
  lapply(1:9, function(i) {
    df %>%
      mutate(across(where(is.numeric), ~ . + rnorm(n(), 0, 0.01))) %>%
      mutate(Replicate = paste0("Replicate ", i))
  })
)

# Filter out rows where the Summary column is empty or missing
replicates <- replicates %>%
  filter(!is.na(Summary) & Summary != "")

# Group by Summary and Replicate, then calculate means
averaged_data <- replicates %>%
  group_by(Summary, Replicate) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Reshape data for plotting
averaged_data_long <- averaged_data %>%
  pivot_longer(cols = -c(Summary, Replicate), names_to = "Timepoint", values_to = "Average") %>%
  mutate(Timepoint = as.numeric(gsub("X", "", Timepoint)))  # Convert timepoints to numeric

# Plot the 3x3 grid with the legend underneath
ggplot(averaged_data_long, aes(x = Timepoint, y = Average, color = Summary, group = Summary)) +
  geom_line(size = 1) +
  facet_wrap(~ Replicate, ncol = 3) +  # Create a 3x3 grid
  labs(
    title = "Time Series of Averaged Replicates",
    x = "Timepoint",
    y = "Average Value",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Customize facet labels
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",  # Move the legend to the bottom
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

















# Simulate bar plot data for each replicate
simulate_bar_data <- function() {
  data.frame(
    Group = c("Control", "Low", "High"),
    Value = c(runif(1, 40, 60), runif(1, 60, 80), runif(1, 80, 100)),  # Random values for each replicate
    p_value = c(runif(1, 0.01, 0.1), runif(1, 0.001, 0.05), runif(1, 0.0001, 0.01))  # Simulated p-values
  ) %>%
    mutate(Significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ ""
    ))
}

# Function to create bar plot inset for each replicate
create_bar_plot <- function(bar_data) {
  ggplot(bar_data, aes(x = Group, y = Value, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = Significance, y = Value + 5), size = 4, vjust = 0) +  # Add stars above bars
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = "black"),  # White background with black border
      plot.margin = margin(2, 2, 2, 2)  # Adjust spacing around inset
    )
}

# Simulate replicates for time series with added noise
set.seed(123)
replicates <- bind_rows(
  lapply(1:9, function(i) {
    df %>%
      mutate(across(where(is.numeric), ~ . + rnorm(n(), 0, 0.01))) %>%
      mutate(Replicate = paste0("Replicate ", i))
  })
)

# Filter out rows where the Summary column is empty or missing
replicates <- replicates %>%
  filter(!is.na(Summary) & Summary != "")

# Group by Summary and Replicate, then calculate means
averaged_data <- replicates %>%
  group_by(Summary, Replicate) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Reshape data for plotting
averaged_data_long <- averaged_data %>%
  pivot_longer(cols = -c(Summary, Replicate), names_to = "Timepoint", values_to = "Average") %>%
  mutate(Timepoint = as.numeric(gsub("X", "", Timepoint)))

# Generate a list of individual plots with unique insets for each replicate
replicate_plots <- lapply(unique(averaged_data_long$Replicate), function(replicate) {
  # Subset data for the current replicate
  replicate_data <- averaged_data_long %>%
    filter(Replicate == replicate)
  
  # Simulate bar plot data for this replicate
  bar_data <- simulate_bar_data()
  
  # Create the time series plot
  time_series_plot <- ggplot(replicate_data, aes(x = Timepoint, y = Average, color = Summary, group = Summary)) +
    geom_line(size = 1) +
    labs(x = NULL, y = NULL, title = replicate, color = "Group") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, hjust = 0.5)
    )
  
  # Combine the time series plot with the inset bar plot in the upper-left corner
  time_series_plot + 
    inset_element(create_bar_plot(bar_data), left = 0.05, bottom = 0.6, right = 0.5, top = 0.95)  # Taller box
})

# Combine all replicate plots into a 3x3 grid
final_plot <- wrap_plots(replicate_plots, ncol = 3) +
  plot_annotation(
    title = "Time Series with Bar Plot Insets and Simulated Data",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

# Display the final plot
print(final_plot)






