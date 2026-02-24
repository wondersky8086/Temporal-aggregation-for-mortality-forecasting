run_country_analysis <- function(country_code, stmf_data) {
  message("Running analysis for ", country_code)
  # ---------------------------
  # Filter data
  # ---------------------------
  stmf_cty <- stmf_data %>%
    filter(CountryCode == country_code, Sex == "b", Week <= 52) %>%
    select(Year, Week, RTotal) %>%
    mutate(
      iso_week   = sprintf("%04d-W%02d", Year, Week),
      week_start = ISOweek::ISOweek2date(paste0(iso_week, "-1"))
    ) %>%
    arrange(week_start)
  
  # Training data: up to end of 2023
  train_data <- stmf_cty %>%
    filter(Year >= 2015 & Year <= 2022)
  
  # Test data: year 2024 only
  test_data <- stmf_cty %>%
    filter(Year == 2023 | Year == 2024)
  
  # create 2-week rolling sums
  test_2w <- test_data %>%
    mutate(two_week_group = rep(1:52, each = 2)) %>%
    group_by(two_week_group) %>%
    summarise(RTotal = sum(RTotal)) %>%
    pull(RTotal)
  
  # create 4-week rolling sums
  test_4w <- test_data %>%
    mutate(four_week_group = rep(1:26, each = 4)) %>%
    group_by(four_week_group) %>%
    summarise(RTotal = sum(RTotal)) %>%
    pull(RTotal)
  
  # create quarter rolling sums
  test_qtr <- test_data %>%
    mutate(qtr_group = rep(1:8, each = 13)) %>%
    group_by(qtr_group) %>%
    summarise(RTotal = sum(RTotal)) %>%
    pull(RTotal)
  
  # create Biannual rolling sums
  test_biann <- test_data %>%
    mutate(biann_group = rep(1:4, each = 26)) %>%
    group_by(biann_group) %>%
    summarise(RTotal = sum(RTotal)) %>%
    pull(RTotal)
  
  # create annual rolling sums
  test_ann <- test_data %>%
    mutate(ann_group = rep(1:2, each = 52)) %>%
    group_by(ann_group) %>%
    summarise(RTotal = sum(RTotal)) %>%
    pull(RTotal)
  
  # Convert to time series
  train_ts <- ts(train_data$RTotal,
                 frequency = 52,
                 start = c(min(train_data$Year), 1))
  
  test_ts <- ts(test_data$RTotal,
                frequency = 52,
                start = c(min(test_data$Year), 1))
  
  test_ts_2w <- ts(test_2w,
                   frequency = 26,
                   start = c(min(test_data$Year), 1))
  
  test_ts_4w <- ts(test_4w,
                   frequency = 13,
                   start = c(min(test_data$Year), 1))
  
  test_ts_qtr <- ts(test_qtr,
                    frequency = 4,
                    start = c(min(test_data$Year), 1))
  
  test_ts_biann <- ts(test_biann,
                      frequency = 2,
                      start = c(min(test_data$Year), 1))
  
  test_ts_ann <- ts(test_ann,
                    frequency = 1,
                    start = c(min(test_data$Year), 1))
  
  # Combine all test sthetaf into a named list
  test_ts_list <- list(
    week  = test_ts,
    two_w = test_ts_2w,
    four_w = test_ts_4w,
    qtr   = test_ts_qtr,
    biann = test_ts_biann,
    ann   = test_ts_ann
  )
  
  
  # ------------------------------
  # Aggregation and Reconciliation
  # ------------------------------
  
  # Create all aggregate levels
  totalagg <- tsaggregates(train_ts)
  
  # Base forecasts (with 95% CI)
  base_theta <- list()
  for (i in seq_along(totalagg)) {
    train_end_time <- time(totalagg[[i]])[length(totalagg[[i]])]
    test_end_time  <- max(time(test_ts_list[[i]]))
    
    h_steps <- ceiling( (test_end_time - train_end_time) * frequency(totalagg[[i]]) )
    base_theta[[i]] <- thetaf(totalagg[[i]],
                              h = h_steps,
                              level = 95)
    
  }
  
  # Reconcile forecasts
  rec_bu_theta <- reconcilethief(base_theta,comb="bu")
  rec_struc_theta <- reconcilethief(base_theta,comb="struc")
  rec_ols_theta <- reconcilethief(base_theta,comb="ols")
  rec_sam_theta <- reconcilethief(base_theta,comb="sam")
  
  # ------------------------
  # Plot comparison with CI
  # ------------------------
  
  savepdf(country_code, width = 12, height = 10, toplines = 0.8)
  par(mfrow = c(2, 3), mar = c(3, 3, 2, 1), oma = c(0, 0, 3, 0))
  
  for (i in length(totalagg):1) {
    # Base
    plot(base_theta[[i]], main = names(totalagg)[i],
         xlim = c(2015.0, 2025), flwd = 1,
         shadecols = adjustcolor("red", alpha.f = 0.1))
    lines(base_theta[[i]]$mean, col = "red", lwd = 1.5)
    lines(test_ts_list[[i]], col = "black", lwd = 1.5)
    # Reconciled
    lines(rec_bu_theta[[i]]$mean, col = "blue", lwd = 2)
    lines(rec_struc_theta[[i]]$mean, col = "green", lwd = 2)
    lines(rec_ols_theta[[i]]$mean, col = "orange", lwd = 2)
    
    legend("topleft", legend = c("Base", "Rec_bu_theta","Rec_struc_theta","Rec_ols_theta","actual"),
           col = c("red", "blue","green","orange","black"), lty = 1, lwd = 2, bty = "n", cex = 0.8)
  }
  mtext(paste0(country_code, " Results comparison - theta method"),
        side = 3, line = 0.5, outer = TRUE, cex = 1.2)
  
  dev.off()
  
  # -------------------------------
  # Accuracy metrics: MASE and RMAE
  # -------------------------------
  
  results_theta <- data.frame(
    Level = names(test_ts_list),
    SMAPE_Base = NA, SMAPE_BU = NA, SMAPE_Struc = NA, SMAPE_OLS = NA, SMAPE_SAM = NA, 
    MASE_Base = NA, MASE_BU = NA, MASE_Struc = NA, MASE_OLS = NA, MASE_SAM = NA, 
    Best_SMAPE = NA,
    Best_MASE  = NA,
    Freq = NA,
    stringsAsFactors = FALSE
  )
  
  # Loop over all frequencies
  for (i in seq_along(test_ts_list)) {
    t_i <- test_ts_list[[i]]
    f_i <- frequency(t_i)
    
    # Find matching index in base_theta by frequency
    match_idx <- which(sapply(base_theta, function(x) frequency(x$x)) == f_i)
    if (length(match_idx) == 0) {
      warning(paste("No matching model found for frequency", f_i))
      next
    }
    idx <- match_idx[1]
    
    # Errors
    err_base  <- as.numeric(t_i) - as.numeric(base_theta[[idx]]$mean)
    err_bu    <- as.numeric(t_i) - as.numeric(rec_bu_theta[[idx]]$mean)
    err_struc <- as.numeric(t_i) - as.numeric(rec_struc_theta[[idx]]$mean)
    err_ols   <- as.numeric(t_i) - as.numeric(rec_ols_theta[[idx]]$mean)
    err_sam   <- as.numeric(t_i) - as.numeric(rec_sam_theta[[idx]]$mean)
    
    # MAE
    mae_base  <- mean(abs(err_base), na.rm = TRUE)
    mae_bu    <- mean(abs(err_bu), na.rm = TRUE)
    mae_struc <- mean(abs(err_struc), na.rm = TRUE)
    mae_ols   <- mean(abs(err_ols), na.rm = TRUE)
    mae_sam   <- mean(abs(err_sam), na.rm = TRUE)
    
    # MASE
    train_ts_i <- aggregate(train_ts, nfrequency = f_i)
    naive_errors <- diff(train_ts_i, lag = 1)
    scale_factor <- mean(abs(naive_errors), na.rm = TRUE)
    
    mase_base  <- mae_base / scale_factor
    mase_bu    <- mae_bu / scale_factor
    mase_struc <- mae_struc / scale_factor
    mase_ols   <- mae_ols / scale_factor
    mase_sam   <- mae_sam / scale_factor
    
    # sMAPE (symmetric mean absolute percentage error)
    smape_base  <- mean( 2 * abs(err_base) / 
                           (abs(as.numeric(t_i)) + abs(as.numeric(base_theta[[idx]]$mean))), na.rm = TRUE)
    smape_bu    <- mean( 2 * abs(err_bu) / 
                           (abs(as.numeric(t_i)) + abs(as.numeric(rec_bu_theta[[idx]]$mean))), na.rm = TRUE)
    smape_struc <- mean( 2 * abs(err_struc) / 
                           (abs(as.numeric(t_i)) + abs(as.numeric(rec_struc_theta[[idx]]$mean))), na.rm = TRUE)
    smape_ols   <- mean( 2 * abs(err_ols) / 
                           (abs(as.numeric(t_i)) + abs(as.numeric(rec_ols_theta[[idx]]$mean))), na.rm = TRUE)
    smape_sam   <- mean( 2 * abs(err_sam) / 
                           (abs(as.numeric(t_i)) + abs(as.numeric(rec_sam_theta[[idx]]$mean))), na.rm = TRUE)
    
    # Collect metrics into named vectors
    smape_vec <- c(
      Base  = smape_base,
      BU    = smape_bu,
      Struc = smape_struc,
      OLS   = smape_ols,
      SAM   = smape_sam
    )
    
    mase_vec <- c(
      Base  = mase_base,
      BU    = mase_bu,
      Struc = mase_struc,
      OLS   = mase_ols,
      SAM   = mase_sam
    )
    
    best_smape_method <- names(smape_vec)[which.min(smape_vec)]
    best_mase_method  <- names(mase_vec)[which.min(mase_vec)]
    
    # Store results
    results_theta[i, 2:14] <- c(
      smape_base, smape_bu, smape_struc, smape_ols, smape_sam,
      mase_base, mase_bu, mase_struc, mase_ols, mase_sam,
      best_smape_method,
      best_mase_method,
      f_i
    )
    
  }
  
  results_theta <- results_theta[order(-as.numeric(results_theta$Freq)), ]
  results_theta$country <- country_code
  
  # -------------------------------
  # Add average across frequencies
  # -------------------------------
  
  # Columns for SMAPE and MASE
  smape_cols <- grep("^SMAPE_", names(results_theta), value = TRUE)
  mase_cols  <- grep("^MASE_",  names(results_theta), value = TRUE)
  
  # Convert columns to numeric 
  results_theta[smape_cols] <- lapply(results_theta[smape_cols], as.numeric)
  results_theta[mase_cols]  <- lapply(results_theta[mase_cols], as.numeric)
  
  # Compute average values
  avg_smape <- colMeans(results_theta[, smape_cols], na.rm = TRUE)
  avg_mase  <- colMeans(results_theta[, mase_cols],  na.rm = TRUE)

  avg_values <- c(avg_smape, avg_mase)
  
  # Copy first row for structure
  avg_row <- results_theta[1, , drop = FALSE]
  avg_row[,] <- NA
  
  # Fill in average metrics
  avg_row$Level <- "Average"
  avg_row[smape_cols] <- avg_smape
  avg_row[mase_cols]  <- avg_mase
  
  # Best methods based on averages
  best_smape_method <- names(avg_smape)[which.min(avg_smape)]
  best_mase_method  <- names(avg_mase)[which.min(avg_mase)]
  avg_row$Best_SMAPE <- sub("SMAPE_", "", best_smape_method)
  avg_row$Best_MASE  <- sub("MASE_",  "", best_mase_method)
  
  # Other columns
  avg_row$Freq <- "MIXED"
  avg_row$country <- country_code
  
  # Append to results
  results_theta <- rbind(results_theta, avg_row)

  write.csv(results_theta,
           paste0("results_theta_", country_code, ".csv"),
           row.names = FALSE)
  
  return(results_theta)
}
