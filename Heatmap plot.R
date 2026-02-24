# Get frequencies in order 
freq_levels <- unique(final_results$Level)
freq_levels <- freq_levels[freq_levels != "Average"]  
freq_levels <- c(freq_levels, "Average")            

# Convert Level to factor with desired order
final_results$Level <- factor(final_results$Level, levels = freq_levels)
smape_wins$Level    <- factor(smape_wins$Level, levels = freq_levels)
mase_wins$Level     <- factor(mase_wins$Level, levels = freq_levels)


# Heatmap by country and frequency
ggplot(final_results,
       aes(x = Level, y = country, fill = Best_SMAPE)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "Base"  = "red",
      "BU"    = "blue",
      "Struc" = "green",
      "OLS"   = "orange",
      "SAM" = "purple"
    )
  ) +
  labs(
    title = "Best Method by SMAPE",
    x = "Frequency",
    y = "Country",
    fill = "Winning method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(final_results,
       aes(x = Level, y = country, fill = Best_MASE)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "Base"  = "red",
      "BU"    = "blue",
      "Struc" = "green",
      "OLS"   = "orange",
      "SAM" = "purple"
    )
  ) +
  labs(
    title = "Best Method by MASE",
    x = "Frequency",
    y = "Country",
    fill = "Winning method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# plot the Winner SMAPE and MASE
ggplot(smape_wins,
       aes(x = Level, y = wins, fill = Best_SMAPE)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(
      "Base"  = "red",
      "BU"    = "blue",
      "Struc" = "green",
      "OLS"   = "orange",
      "SAM" = "purple"
    )
  ) +
  labs(
    title = "SMAPE: win counts",
    x = "Frequency",
    y = "Number of countries",
    fill = "Method"
  ) +
  theme_minimal()

ggplot(mase_wins,
       aes(x = Level, y = wins, fill = Best_MASE)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(
      "Base"  = "red",
      "BU"    = "blue",
      "Struc" = "green",
      "OLS"   = "orange",
      "SAM" = "purple"
    )
  ) +
  labs(
    title = "MASE: win counts",
    x = "Frequency",
    y = "Number of countries",
    fill = "Method"
  ) +
  theme_minimal()
