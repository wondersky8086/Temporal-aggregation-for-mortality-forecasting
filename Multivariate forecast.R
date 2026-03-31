# ===============
# Load libraries
# ===============

require(readr)
require(ftsa)
require(thief)
require(forecast)
require(ISOweek)
require(tidyr)
require(ggplot2)
require(dplyr)
require(viridis)
require(purrr)
require(rainbow)
require(patchwork)
require(grid)

# ======================
# Load weekly STMF data
# ======================

stmf_data <- read_csv("stmf.csv", skip = 1, col_names = TRUE)

age_cols <- c("D0_14","D15_64","D65_74","D75_84","D85p")

scale_levels <- c(
  week  = 1,
  two_w = 2,
  four_w = 4,
  qtr   = 13,
  biann = 26,
  ann   = 52
)

results_all <- data.frame()

# ===================
# Accuracy functions
# ===================

smape <- function(actual, forecast){
  mean(2 * abs(actual - forecast) /
         (abs(actual) + abs(forecast)), na.rm=TRUE)
}

mase <- function(actual, forecast, train){
  scale_factor <- mean(abs(diff(train)), na.rm=TRUE)
  mean(abs(actual - forecast), na.rm=TRUE) / scale_factor
}

# ==================
# Main Calculations
# ==================

countries <- unique(stmf_data$CountryCode)
# countries <- c("AUS","CAN","SWE")
results_all <- data.frame()

for(country_code in countries){
  
  cat("Processing:", country_code, "\n")
  
  tryCatch({
  
  # Prepare country data
  stmf_cty <- stmf_data %>%
    filter(CountryCode == country_code,
           Sex == "b",
           Week <= 52) %>%
    select(Year, Week, all_of(age_cols)) %>%
    arrange(Year, Week)
  
  train_data <- stmf_cty %>% filter(Year <= 2023)
  test_data  <- stmf_cty %>% filter(Year == 2024)
  
  if(nrow(train_data) == 0 | nrow(test_data) == 0){
    warning("Insufficient data for ", country_code, "- skip")
    next
  }
  
  base_fc_scales <- list()
  test_scales <- list()
  train_scales <- list()
  
  # Forecast each scale
  for(scale_name in names(scale_levels)){
    
    f <- scale_levels[scale_name]
    
    train_agg <- train_data %>%
      mutate(group = rep(1:ceiling(n()/f), each=f, length.out=n())) %>%
      group_by(group) %>%
      summarise(across(all_of(age_cols), sum)) %>%
      ungroup()
    
    test_agg <- test_data %>%
      mutate(group = rep(1:ceiling(n()/f), each=f, length.out=n())) %>%
      group_by(group) %>%
      summarise(across(all_of(age_cols), sum)) %>%
      ungroup()
    
    train_scales[[scale_name]] <- train_agg
    test_scales[[scale_name]]  <- test_agg
    
  # Fit FTSM
    train_matrix <- t(as.matrix(train_agg[, age_cols]))
    
    rownames(train_matrix) <- age_cols
    colnames(train_matrix) <- seq_len(ncol(train_matrix))
    
    age_grid <- seq_len(length(age_cols))
    
    h <- nrow(test_agg)
    
    fc_fts <- forecast(ftsm(fts(age_grid, train_matrix), order = 3), h = h)
    
    base_fc_scales[[scale_name]] <- fc_fts$mean$y
    
  }
  
  # Reconcile across scales per age band
  for(i in seq_along(age_cols)){
    
    age <- age_cols[i]
    
    age_scale_fc <- list()
    
    for(scale_name in names(scale_levels)){
      
      f <- scale_levels[scale_name]
      
      fc_vec <- base_fc_scales[[scale_name]][i, ]
      
      age_scale_fc[[scale_name]] <- ts(fc_vec, frequency = 52 / f)
      
    }

    rec_bu  <- reconcilethief(age_scale_fc, comb="bu")
    rec_ols <- reconcilethief(age_scale_fc, comb="ols")
    rec_str <- reconcilethief(age_scale_fc, comb="struc")
    
  # Accuracy evaluation
    scale_idx <- 1
    
    for(scale_name in names(scale_levels)){
      
      f <- scale_levels[scale_name]
      
      actual <- ts(test_scales[[scale_name]][[age]], frequency = 52 / f)
      
      base_mean <- age_scale_fc[[scale_name]]
      bu_mean   <- rec_bu[[scale_idx]]
      ols_mean  <- rec_ols[[scale_idx]]
      str_mean  <- rec_str[[scale_idx]]
      
      scale_idx <- scale_idx + 1
      
      train_vec <- train_scales[[scale_name]][[age]]
      
      results_all <- rbind(
        results_all,
        data.frame(
          Country = country_code,
          AgeBand = age,
          Scale   = scale_name,
          Method  = c("Base","BU","OLS","STRUC"),
          SMAPE   = c(
            smape(actual, base_mean),
            smape(actual, bu_mean),
            smape(actual, ols_mean),
            smape(actual, str_mean)
          ),
          MASE    = c(
            mase(actual, base_mean, train_vec),
            mase(actual, bu_mean, train_vec),
            mase(actual, ols_mean, train_vec),
            mase(actual, str_mean, train_vec)
          )
        )
      )
      
    }
    
  }
  
  })
  
}

# Compute best SMAPE and MASE per AgeBand and Scale
best_results <- results_all %>%
  group_by(Country, AgeBand, Scale) %>%
  summarise(
    Best_SMAPE = Method[which.min(SMAPE)],
    Min_SMAPE  = min(SMAPE),
    Best_MASE  = Method[which.min(MASE)],
    Min_MASE   = min(MASE),
    .groups = "drop"
  )

results_all <- results_all %>%
  left_join(best_results, by = c("Country","AgeBand","Scale"))

write.csv(results_all,
          "ftsm_temporal_reconciliation_results.csv",
          row.names = FALSE)

# ======
# Plots
# ======

# Prepare ordered factors and data for plot
scale_order <- c("week", "two_w", "four_w", "qtr", "biann", "ann")
results_all$Scale <- factor(results_all$Scale, levels = scale_order)

age_order <- age_cols
results_all$AgeBand <- factor(results_all$AgeBand, levels = age_order)

method_colors <- c(
  "Base"  = "red",
  "BU"    = "blue",
  "STRUC" = "green",
  "OLS"   = "orange"
)

# Function - a plot for one scale
plot_best_method_grid <- function(df, measure, 
                        show_x = FALSE, show_y = FALSE, show_legend = FALSE) { 
                        ggplot(df, aes(x = AgeBand, y = Country, fill = BestMethod)) + 
                        geom_tile(color = "white") + 
                        scale_fill_manual(values = method_colors) + 
                        coord_fixed(ratio = 0.3)+ 
                        labs( title = paste0(unique(df$Scale)), x = if(show_x) "Age Band" else NULL, 
                              y = if(show_y) "Country" else NULL, fill = "Best Method" ) + 
                        theme_minimal() + 
                        theme( axis.text.x = if(show_x) element_text(angle = 0, hjust = 0.5, size = 7) 
                               else element_blank(), 
                               axis.text.y = if(show_y) element_text(size = 7) 
                               else element_blank(), 
                               axis.ticks = element_blank(), 
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(), 
                               legend.position = if(show_legend) "bottom" else "none", 
                               legend.title = element_text(size = 6), 
                               legend.text = element_text(size = 5), 
                               legend.key.size = unit(0.3, "cm"), 
                               legend.spacing.x = unit(0.1, "cm"), 
                               legend.spacing.y = unit(0.1, "cm"), 
                               legend.box.margin = margin(0, 0, 0, 0), 
                               legend.margin = margin(0, 0, 0, 0), 
                               plot.title = element_text(hjust = 0.5, size = 10) ) } 
plot_best_method_all_scales <- function(results_all, measure) {
  prepare_best_method_data <- function(results_all, measure) {
    results_all %>%
      group_by(Country, AgeBand, Scale) %>%
      summarise(
        BestMethod = Method[which.min(.data[[measure]])],
        .groups = "drop"
      )
  }
  plot_data <- prepare_best_method_data(results_all, measure) 
  split_data <- split(plot_data, plot_data$Scale) 
  plots <- map2(split_data, seq_along(split_data), function(df, idx) {
    n <- length(split_data)
    ncol <- 3
    show_x <- idx > (n - ncol)
    show_y <- (idx - 1) %% ncol == 0
    show_legend <- idx == n
    plot_best_method_grid(df, measure = measure, 
                           show_x = show_x, show_y = show_y, 
                           show_legend = show_legend ) }) 
  wrap_plots(plots, ncol = 3, nrow = 2) + 
    plot_annotation( title = paste0("Best Method by ", measure, " Across Scales"), 
                     theme = theme( plot.title = 
                                    element_text( hjust = 0.5, size = 14, face = "bold", 
                                                  margin = margin(b = 10) ) ) ) }

smape_grid <- plot_best_method_all_scales(results_all, "SMAPE")
mase_grid  <- plot_best_method_all_scales(results_all, "MASE")

# Save as png
savefig("smape_grid", height = 38, width = 25, toplines = 0.5, type = "png")
print(smape_grid)
dev.off()

savefig("mase_grid", height = 38, width = 25, toplines = 0.5, type = "png")
print(mase_grid)
dev.off()

# Method performance across temporal aggregation levels

plot_smape_data <- results_all %>% 
  group_by(Country, AgeBand, Scale) %>% 
  slice_min(SMAPE, with_ties = FALSE) %>% 
  ungroup()

savepdf("Overall Performance by Method", width = 12, height = 10, toplines = 0.8)
plot_smape_data %>%
  count(Scale, Method) %>%
  ggplot(aes(Method, n, fill = Method)) +
  geom_col() + 
  scale_fill_manual(values = method_colors, drop = FALSE) +
  facet_wrap(~Scale) +
  labs(title = "Overall Performance by Method", x = NULL, y = NULL) +
  theme(
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

dev.off()

# Shows method dominance across ages
  win_table <- plot_smape_data %>%
    count(AgeBand, Best_SMAPE)
  savepdf("Dominant method by age band", width = 12, height = 10, toplines = 0.8)
  plot(ggplot(win_table,
       aes(AgeBand, Best_SMAPE, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c()+
    labs(
      title = "Dominant method by age band",
      x = NULL,
      y = NULL
    ) +
    theme(
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)))
  dev.off()
  
# Performance pattern across ages
  savepdf("Average SMAPE per age band by frequency", width = 12, height = 10, toplines = 0.8)
  plot(results_all %>%
         group_by(AgeBand, Method, Scale) %>%
         summarise(SMAPE = mean(SMAPE), .groups = "drop") %>%
         ggplot(aes(AgeBand, SMAPE, color = Method, group = Method)) +
         geom_line() +
         geom_point() +
         scale_color_manual(values = method_colors) +
         facet_wrap(~Scale, ncol = 3) +
         labs(
           title = "Average SMAPE per age band by frequency",
           x = NULL,
           y = NULL
         ) +
         theme(
           axis.ticks.x = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
           legend.title = element_blank(),
           legend.position = "bottom",
           plot.title = element_text(hjust = 0.5)
         ))
  dev.off()
  
  # Performance pattern across ages - MASE
  savepdf("Average MASE per age band by frequency", width = 12, height = 10, toplines = 0.8)
  plot(results_all %>%
         group_by(AgeBand, Method, Scale) %>%
         summarise(MASE = mean(MASE), .groups = "drop") %>%
         ggplot(aes(AgeBand, MASE, color = Method, group = Method)) +
         geom_line() +
         geom_point() +
         scale_color_manual(values = method_colors) +
         facet_wrap(~Scale, ncol = 3) +
         labs(
           title = "Average MASE per age band by frequency",
           x = NULL,
           y = NULL
         ) +
         theme(
           axis.ticks.x = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
           legend.title = element_blank(),
           legend.position = "bottom",
           plot.title = element_text(hjust = 0.5)
         ))
  dev.off()
