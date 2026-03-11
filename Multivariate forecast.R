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
  
  }, 
  error = function(e){
    warning("Processing failed for ", country_code, ": ", e$message)
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

# Prepare ordered factors
scale_order <- c("week", "two_w", "four_w", "qtr", "biann", "ann")
results_all$Scale <- factor(results_all$Scale, levels = scale_order)

age_order <- age_cols
results_all$AgeBand <- factor(results_all$AgeBand, levels = age_order)

# Function to create a plot for one scale
plot_best_method <- function(df, measure = "SMAPE") {
  ggplot(df, aes(x = AgeBand, y = Country, fill = BestMethod)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c(
        "Base"  = "red",
        "BU"    = "blue",
        "STRUC" = "green",
        "OLS"   = "orange"
      )) +
    theme_minimal() +
    labs(
      title = paste0("Best Method by ", measure, " - Scale: ", unique(df$Scale)),
      x = "Age Band",
      y = "Country",
      fill = "Best Method"
    )
}

# Prepare data for SMAPE
plot_smape_data <- results_all %>%
  group_by(Country, AgeBand, Scale) %>%
  summarise(BestMethod = Method[which.min(SMAPE)], .groups = "drop")

# Prepare data for MASE
plot_mase_data <- results_all %>%
  group_by(Country, AgeBand, Scale) %>%
  summarise(BestMethod = Method[which.min(MASE)], .groups = "drop")

# Split by scale and generate plots
smape_plots <- plot_smape_data %>%
  split(.$Scale) %>%
  map(~ plot_best_method(.x, measure = "SMAPE"))

mase_plots <- plot_mase_data %>%
  split(.$Scale) %>%
  map(~ plot_best_method(.x, measure = "MASE"))


for(i in seq_along(scale_order)) {
  
  p <- smape_plots[[i]]
  
  savepdf(paste0("SMAPE_", scale_order[i]), width = 12, height = 10, toplines = 0.8)
  plot(p)  
  dev.off()
}

for(j in seq_along(scale_order)) {
  
  p <- mase_plots[[j]]
  
  savepdf(paste0("MASE_", scale_order[j]), width = 12, height = 10, toplines = 0.8)
  plot(p)  
  dev.off()
}