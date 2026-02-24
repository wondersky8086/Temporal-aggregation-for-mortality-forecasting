# ---------------------------
# Packages
# ---------------------------

require(readr)
require(dplyr)
require(lubridate)
require(ISOweek)
require(forecast)
require(thief)
require(ggplot2)
require(patchwork)
require(xtable)

# ---------------------------
# Load weekly STMF data
# ---------------------------

stmf_data <- read_csv("stmf.csv", skip = 1, col_names = TRUE)

# ---------------------------
# Select Countries
# ---------------------------

# countries <- c("AUS", "CAN", "SWE")
 
countries <- stmf_data %>%
  distinct(CountryCode) %>%
  pull(CountryCode) %>%
  sort()


all_results <- lapply(countries, function(cty) {
  tryCatch(
    run_country_analysis(cty, stmf_data = stmf_data),
    error = function(e) {
      message("Failed for country: ", cty)
      NULL
    })
})

names(all_results) <- countries

final_results <- bind_rows(all_results) %>%
  filter(!is.na(country)) %>%
  select(-contains("MAE"))


# SMAPE table for all levels
smape_table <- final_results %>%
  select(country, Level, contains("SMAPE"))

# MASE table for all levels
mase_table <- final_results %>%
  select(country, Level, contains("MASE"))

smape_wins <- final_results %>%
  count(Level, Best_SMAPE, name = "wins")

mase_wins <- final_results %>%
  count(Level, Best_MASE, name = "wins")

# SMAPE LaTex table

smape_table_print <- smape_table %>%
  mutate(across(contains("SMAPE_"), ~ as.numeric(.))) %>%  
  mutate(across(contains("SMAPE_"), ~ round(., 3)))       

xt_smape <- xtable(smape_table_print,
                   caption = "Forecast SMAPE by Country and Level",
                   label = "tab:smape")


digits(xt_smape) <- c(0, 0, 0, rep(3,  ncol(smape_table) - 3),0)

print(xt_smape, include.rownames = FALSE, booktabs = TRUE)


# MASE LaTex table

mase_table_print <- mase_table %>%
  mutate(across(contains("MASE_"), ~ as.numeric(.))) %>%  
  mutate(across(contains("MASE_"), ~ round(., 3)))       

xt_mase <- xtable(mase_table_print,
                   caption = "Forecast MASE by Country and Level",
                   label = "tab:mase")


digits(xt_mase) <- c(0, 0, 0, rep(3,  ncol(mase_table) - 3),0)

print(xt_mase, include.rownames = FALSE, booktabs = TRUE)

