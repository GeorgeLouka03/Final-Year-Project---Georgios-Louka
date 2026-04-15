############################################################
# FINAL YEAR PROJECT - U.S. HURRICANE CAT MODEL
# Georgios Louka
############################################################

# 0. Setup and package loading
# 1. Load and clean data
# 2. Build hurricane-only sample
# 3. Frequency model
# 4. Severity model
# 5. Wind-loss regression and imputation
# 6. Monte Carlo simulation of annual losses
# 7. XoL pricing
# 8. Parametric CAT pricing
# 9. Sensitivity analysis
# 10. Export tables / plots

############################################################
# FINAL YEAR PROJECT - SETUP
############################################################

# Run once if needed to install the new packages:
# install.packages(c(
#   "readxl", "janitor", "dplyr", "ggplot2", "broom",
#   "fitdistrplus", "extRemes", "lubridate", "scales", "patchwork"
# ))

library(readxl)
library(janitor)
library(dplyr)
library(ggplot2)
library(broom)
library(fitdistrplus)
library(extRemes)
library(lubridate)   # Needed for arrival-time (seasonality) analysis
library(scales)      # Needed for the log-scale survival plot
library(patchwork)   # Needed for the side-by-side sensitivity plots

############################################################
# 1. LOAD AND CLEAN DATA
############################################################

path <- "Final Year Project Data.xlsx"

# Check sheet names first
excel_sheets(path)

# Read first sheet
df_raw <- read_excel(path, sheet = 2)

# Standardise names
df <- df_raw %>%
  clean_names()

# Quick inspection
names(df)
glimpse(df)

############################################################
# 2. BUILD ANALYSIS SAMPLE
############################################################

unique(df$max_winds_kt)

df_clean <- df %>%
  mutate(
    year = as.integer(year),
    max_winds_kt = as.numeric(max_winds_kt),
    min_pressure_mb = as.numeric(min_pressure_mb),
    u_s_damage_million = as.numeric(u_s_damage_million),
    landfall_yes_no = trimws(landfall_yes_no),
    typea = trimws(typea)
  ) %>%
  filter(
    year >= 1991,
    year <= 2024,
    landfall_yes_no == "Yes"
  )

# Hurricanes only: HU + MH (Aligns with >= 64 knots methodology)
df_hu <- df_clean %>%
  filter(typea %in% c("HU", "MH"))

# Basic checks
summary_hu <- df_hu %>%
  summarise(
    n_rows = n(),
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    min_wind = min(max_winds_kt, na.rm = TRUE),
    max_wind = max(max_winds_kt, na.rm = TRUE),
    missing_wind = sum(is.na(max_winds_kt)),
    missing_loss = sum(is.na(u_s_damage_million))
  )

print(summary_hu)

########################################################################################################################

############################################################
# 3. FREQUENCY MODEL
# Covers Methodology Section 2.1 and Results Section 4.1
############################################################

t0 <- 1991

############################################################
# 3.1 Annual hurricane counts N_t
# Supports Section 3.1 (Data) and 2.1.1 (NHPP Framework)
############################################################

annual_counts <- df_hu %>%
  count(year, name = "N") %>%
  mutate(t = year - t0)

print(annual_counts)

############################################################
# 3.2 Descriptive statistics
# Supports Section 3.2 (Validating the Poisson Model)
############################################################

freq_stats <- annual_counts %>%
  summarise(
    mean_N = mean(N),
    var_N = var(N),
    dispersion = var_N / mean_N
  )

print(freq_stats)

############################################################
# 3.3 Fit HPP and NHPP models
# Corresponds to Section 2.1.1 (Intensity Functions)
############################################################

# Homogeneous Poisson model
m0 <- glm(N ~ 1, family = poisson(link = "log"), data = annual_counts)

# Non-homogeneous Poisson model with log-linear trend
m1 <- glm(N ~ t, family = poisson(link = "log"), data = annual_counts)

############################################################
# 3.4 Parameter estimation
# Corresponds to Section 2.1.2 (Parameter Estimation via MLE)
############################################################

tidy_m0 <- tidy(m0)
tidy_m1 <- tidy(m1)

print(tidy_m0)
print(tidy_m1)

# Extract estimated baseline HPP intensity
mu_hat <- as.numeric(exp(coef(m0)["(Intercept)"]))
print(mu_hat)

# Extract NHPP coefficients
alpha_hat <- coef(m1)[1]
beta_hat  <- coef(m1)[2]

cat("Estimated alpha_hat =", alpha_hat, "\n")
cat("Estimated beta_hat  =", beta_hat, "\n")

############################################################
# 3.5 Model comparison and testing
# Corresponds to Section 2.1.3 (Model Selection and Validation)
############################################################

# Likelihood ratio test: H0: beta = 0
lrt <- anova(m0, m1, test = "Chisq")
print(lrt)

# AIC comparison
aic_tab <- AIC(m0, m1)
print(aic_tab)

# Pearson dispersion ratio for NHPP fit
dispersion_ratio <- sum(residuals(m1, type = "pearson")^2) / m1$df.residual
print(dispersion_ratio)

############################################################
# 3.6 Fitted values and plot
# Supports Section 4.1 (Frequency Model Calibration Results)
############################################################

annual_counts <- annual_counts %>%
  mutate(
    mu_hat_hpp  = fitted(m0),
    mu_hat_nhpp = fitted(m1)
  )

ggplot(annual_counts, aes(x = year)) +
  geom_point(aes(y = N), size = 2) +
  geom_line(aes(y = mu_hat_hpp, linetype = "HPP")) +
  geom_line(aes(y = mu_hat_nhpp, linetype = "NHPP")) +
  labs(
    x = "Year",
    y = "Annual landfalling hurricane count",
    title = "Observed hurricane counts and fitted Poisson models",
    linetype = "Model"
  ) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

############################################################
# 3.7 Frequency simulation inputs
# Corresponds to Section 2.1.4 (Simulation of Hurricane Frequencies)
############################################################

# Baseline simulation uses HPP because beta is not statistically significant
cat("Baseline frequency model uses HPP with constant intensity mu_hat =", mu_hat, "\n")

# Projected NHPP intensities for sensitivity analysis only
future_years <- data.frame(
  year = 2025:2030,
  t = (2025:2030) - t0
)

future_years <- future_years %>%
  mutate(
    mu_hat_nhpp = exp(alpha_hat + beta_hat * t)
  )

print(future_years)

############################################################
# 3.8 Intra-annual arrival-time analysis (exploratory)
# Corresponds to Results Section 4.1 (Appendix A)
############################################################

get_midpoint <- function(date_range, year_val) {
  parts <- strsplit(date_range, "–")[[1]]
  
  start <- paste(parts[1], year_val)
  end   <- paste(parts[2], year_val)
  
  start_date <- dmy(start)
  end_date   <- dmy(end)
  
  mid_date <- start_date + (end_date - start_date) / 2
  
  yday(mid_date) / 365
}

df_hu <- df_hu %>%
  rowwise() %>%
  mutate(arrival_time = get_midpoint(datesb_utc, year)) %>%
  ungroup()

# Histogram of within-year arrival times (Upgraded for Appendix)
ggplot(df_hu, aes(x = arrival_time)) +
  geom_histogram(bins = 20, fill = "#4A90E2", color = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Hurricane Arrival Times Within the Year",
    subtitle = "Peak activity is heavily concentrated between July (0.5) and October (0.8)",
    x = "Time Within Year (Fraction: 0 = Jan 1, 1 = Dec 31)",
    y = "Frequency (Number of Landfalls)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(face = "italic", color = "grey30", margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank() # Removes the smaller background grid lines for a cleaner look
  )

# Kernel density estimate of within-year intensity
ggplot(df_hu, aes(x = arrival_time)) +
  geom_density(size = 1.2) +
  labs(
    title = "Estimated intensity function within year",
    x = "Time (fraction of year)",
    y = "Density"
  ) +
  theme_minimal()


# Residual diagnostics (NHPP model)
par(mfrow = c(1,1))

plot(residuals(m1, type="pearson"),
     main="Pearson Residuals (NHPP Model)",
     ylab="Residuals",
     xlab="Observation index")

abline(h=0, col="red")

########################################################################################################################

############################################################
# 4. SEVERITY MODEL (WIND SPEED)
# Covers Methodology Sections 2.2 - 2.2.2 and Results Section 4.2
#
# IMPORTANT:
# This section models WIND severity (W), not insured loss severity (L).
# Therefore we use all hurricanes with valid positive wind speeds,
# regardless of whether insured loss data is missing or zero.
############################################################

library(dplyr)
library(ggplot2)
library(fitdistrplus)
library(extRemes)

############################################################
# 4.1 Prepare wind-speed data
# Supports Methodology Section 3.3 (Severity Data)
############################################################

wind_data <- df_hu %>%
  dplyr::filter(!is.na(max_winds_kt), max_winds_kt > 0) %>%
  dplyr::select(year, storm_name, max_winds_kt)

wind <- wind_data$max_winds_kt

# Basic descriptive statistics
wind_stats <- wind_data %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_wind = mean(max_winds_kt),
    sd_wind = sd(max_winds_kt),
    min_wind = min(max_winds_kt),
    q25 = quantile(max_winds_kt, 0.25),
    median_wind = median(max_winds_kt),
    q75 = quantile(max_winds_kt, 0.75),
    max_wind = max(max_winds_kt),
    skewness_proxy = mean((max_winds_kt - mean(max_winds_kt))^3) / sd(max_winds_kt)^3
  )

print(wind_stats)

############################################################
# 4.2 Exploratory plots for wind severity
############################################################

# Histogram of wind speeds
ggplot(wind_data, aes(x = max_winds_kt)) +
  geom_histogram(bins = 20, fill = "grey70", color = "black") +
  labs(
    title = "Distribution of landfall wind speeds",
    x = "Maximum sustained wind speed (kt)",
    y = "Frequency"
  ) +
  theme_minimal()

# Boxplot
ggplot(wind_data, aes(y = max_winds_kt)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of landfall wind speeds",
    y = "Maximum sustained wind speed (kt)"
  ) +
  theme_minimal()

############################################################
# 4.3 Lognormal fit (baseline model)
# Corresponds to Section 2.2.1: W ~ Lognormal(mu_W, sigma_W^2)
############################################################

fit_lnorm <- fitdist(wind, "lnorm")

print(summary(fit_lnorm))

# Extract MLEs
mu_w_hat <- as.numeric(fit_lnorm$estimate["meanlog"])
sigma_w_hat <- as.numeric(fit_lnorm$estimate["sdlog"])

cat("Estimated mu_W (meanlog) =", mu_w_hat, "\n")
cat("Estimated sigma_W (sdlog) =", sigma_w_hat, "\n")

# AIC
aic_lnorm <- fit_lnorm$aic
cat("Lognormal AIC =", aic_lnorm, "\n")

############################################################
# 4.4 Lognormal diagnostic plots
############################################################

# Standard fitdistrplus diagnostics
plot(fit_lnorm)

# QQ plot manually (cleaner for report)
lnorm_theoretical <- qlnorm(
  ppoints(length(wind)),
  meanlog = mu_w_hat,
  sdlog = sigma_w_hat
)

wind_sorted <- sort(wind)

qq_df <- data.frame(
  theoretical = lnorm_theoretical,
  empirical = wind_sorted
)

ggplot(qq_df, aes(x = theoretical, y = empirical)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "QQ plot: empirical wind speeds vs fitted Lognormal",
    x = "Theoretical quantiles (Lognormal)",
    y = "Empirical quantiles"
  ) +
  theme_minimal()

############################################################
# 4.5 GPD tail fit (robustness check)
# Corresponds to Section 2.2.2: Diagnostic Tail Severity Validation
############################################################

# Candidate thresholds in knots
thresholds <- c(120, 125, 130, 135, 140)

gpd_results <- list()
gpd_mle_table <- data.frame()

for (u in thresholds) {
  
  n_exc <- sum(wind > u)
  
  fit <- fevd(
    x = wind,
    threshold = u,
    type = "GP",
    method = "MLE"
  )
  
  sigma_u_hat <- fit$results$par["scale"]
  xi_hat      <- fit$results$par["shape"]
  
  gpd_results[[paste0("u_", u)]] <- fit
  
  gpd_mle_table <- rbind(
    gpd_mle_table,
    data.frame(
      Threshold = u,
      Num_Exceedances = n_exc,
      Scale = as.numeric(sigma_u_hat),
      Shape = as.numeric(xi_hat)
    )
  )
}

print(gpd_mle_table)

############################################################
# 4.6 Diagnostics for each threshold
############################################################

par(mfrow = c(1,1))

for (u in thresholds) {
  
  fit <- gpd_results[[paste0("u_", u)]]
  
  plot(fit, type = "qq", main = paste("Q-Q Plot (u =", u, ")"))
  plot(fit, type = "probprob", main = paste("P-P Plot (u =", u, ")"))
}

par(mfrow = c(1, 1))

############################################################
# 4.7 Choose a reporting threshold
# Example: use 130 or 135 if estimates look most stable
############################################################

u_main <- 130
gpd_fit <- gpd_results[[paste0("u_", u_main)]]

print(summary(gpd_fit))

gpd_params <- gpd_fit$results$par
gpd_scale_hat <- as.numeric(gpd_params["scale"])
gpd_shape_hat <- as.numeric(gpd_params["shape"])

cat("Chosen threshold u_main =", u_main, "\n")
cat("Estimated GPD scale =", gpd_scale_hat, "\n")
cat("Estimated GPD shape =", gpd_shape_hat, "\n")

############################################################
# 4.8 Tail comparison plot
# Corresponds to Results Section 4.2 (Appendix B)
############################################################

# Empirical survival function of wind speeds 
# (Requires 'scales' package loaded in Setup)

tail_df <- data.frame(wind = wind) %>%
  dplyr::arrange(wind) %>%
  dplyr::mutate(
    emp_cdf = (dplyr::row_number() - 0.5) / dplyr::n()
  )

ggplot(tail_df, aes(x = wind, y = 1 - emp_cdf)) +
  geom_point(color = "#4A90E2", alpha = 0.8, size = 2.5) +
  geom_vline(xintercept = 130, linetype = "dashed", color = "#D0021B", linewidth = 1) +
  scale_y_log10(labels = scales::label_number()) + 
  labs(
    title = "Empirical Survival Function of Wind Speeds",
    subtitle = "The dashed red line indicates the 130-knot threshold used for GPD tail fitting.",
    x = "Maximum Sustained Wind Speed (knots)",
    y = "Survival Probability (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(face = "italic", color = "grey30", margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank()
  )

############################################################
# 4.9 Summary table for Results Section 4.2
############################################################

severity_summary <- data.frame(
  Model = c("Lognormal", "GPD"),
  Parameter_1 = c(mu_w_hat, gpd_scale_hat),
  Parameter_2 = c(sigma_w_hat, gpd_shape_hat),
  AIC = c(aic_lnorm, NA)
)

print(severity_summary)

############################################################
# 4.10 Interpretation flags for report writing
############################################################

cat("\nINTERPRETATION NOTES:\n")
cat("- Wind severity is modelled using all valid wind observations.\n")
cat("- Insured loss missingness does not affect this section.\n")
cat("- Lognormal is the baseline model.\n")
cat("- GPD is a tail robustness check only.\n")

########################################################################################################################


############################################################
# 5. WIND -> INSURED LOSS REGRESSION
# Covers Methodology Section 2.2.3 and Results Section 4.3
############################################################

############################################################
# 5.1 Prepare regression dataset
# Use only storms with valid positive reported insured losses
############################################################

loss_data <- df_hu %>%
  dplyr::filter(
    !is.na(u_s_damage_million),
    u_s_damage_million > 0,
    !is.na(max_winds_kt),
    max_winds_kt > 0
  ) %>%
  dplyr::select(
    year, storm_name, max_winds_kt,
    u_s_damage_million, tropical_cyclone_report_status
  ) %>%
  dplyr::mutate(
    log_loss = log(u_s_damage_million),
    log_wind = log(max_winds_kt)
  )

############################################################
# 5.2 Fit log-log regression
############################################################

lm_model <- lm(log_loss ~ log_wind, data = loss_data)

print(summary(lm_model))

gamma0_hat <- as.numeric(coef(lm_model)[1])
gamma1_hat <- as.numeric(coef(lm_model)[2])
sigma_eps  <- summary(lm_model)$sigma

cat("gamma0_hat =", gamma0_hat, "\n")
cat("gamma1_hat =", gamma1_hat, "\n")
cat("Residual SD (sigma_eps) =", sigma_eps, "\n")

############################################################
# 5.3 Diagnostics
############################################################

ggplot(loss_data, aes(x = log_wind, y = log_loss)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Log-Log Relationship between Wind Speed and Insured Loss",
    x = "log(Wind speed)",
    y = "log(Insured Loss)"
  ) +
  theme_minimal()

ggplot(loss_data, aes(x = log_wind, y = resid(lm_model))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs log(Wind)",
    x = "log(Wind speed)",
    y = "Residuals"
  ) +
  theme_minimal()

############################################################
# 5.4 DETERMINISTIC Loss Imputation (For Baseline Comparison Only)
#
# Note: This deterministic imputation (no epsilon) is used 
# strictly to calculate the baseline "deterministic" average 
# discussed in Results Section 4.4, proving the necessity 
# of the full stochastic Monte Carlo loop later.
############################################################

df_hu <- df_hu %>%
  dplyr::mutate(
    predicted_log_loss = gamma0_hat + gamma1_hat * log(max_winds_kt),
    predicted_loss = exp(predicted_log_loss), # Deterministic (no epsilon)
    
    report_final = !is.na(tropical_cyclone_report_status) &
      tropical_cyclone_report_status == "Final",
    
    loss_final = dplyr::case_when(
      !is.na(u_s_damage_million) & u_s_damage_million > 0 ~ u_s_damage_million,
      !is.na(u_s_damage_million) & u_s_damage_million == 0 & report_final ~ 0,
      TRUE ~ predicted_loss
    )
  )

############################################################
# 5.5 Check classification counts
############################################################

loss_check <- df_hu %>%
  dplyr::mutate(
    loss_source = dplyr::case_when(
      !is.na(u_s_damage_million) & u_s_damage_million > 0 ~ "Observed positive loss",
      !is.na(u_s_damage_million) & u_s_damage_million == 0 &
        tropical_cyclone_report_status == "Final" ~ "Observed final zero loss",
      TRUE ~ "Imputed loss"
    )
  ) %>%
  dplyr::count(loss_source)

print(loss_check)

summary(df_hu$loss_final)

############################################################
# 5.6 DETERMINISTIC AGGREGATE LOSS MODELLING 
# (Compound Poisson Framework Analytical Moments)
############################################################

# Severity moments from completed historical storm-level losses
# This provides the analytical comparison for Section 4.4
loss_mean_hat <- mean(df_hu$loss_final, na.rm = TRUE)
loss_var_hat  <- var(df_hu$loss_final, na.rm = TRUE)

cat("Estimated mu_hat (Frequency) =", mu_hat, "\n")
cat("Estimated E[L] =", loss_mean_hat, "\n")
cat("Estimated Var(L) =", loss_var_hat, "\n")

# Compound Poisson moments
# E[S] = mu * E[L]
agg_mean_hat <- mu_hat * loss_mean_hat

# Var(S) = mu * Var(L) + mu * (E[L])^2
agg_var_hat <- mu_hat * loss_var_hat + mu_hat * (loss_mean_hat^2)

cat("Analytical E[S] (Deterministic) =", agg_mean_hat, "\n")
cat("Analytical Var(S) =", agg_var_hat, "\n")
cat("Analytical SD(S) =", sqrt(agg_var_hat), "\n")

compound_summary <- data.frame(
  mu_hat = mu_hat,
  severity_mean = loss_mean_hat,
  severity_var = loss_var_hat,
  aggregate_mean = agg_mean_hat,
  aggregate_var = agg_var_hat,
  aggregate_sd = sqrt(agg_var_hat)
)

print(compound_summary)

########################################################################################################################

############################################################
# 6. BASELINE MONTE CARLO SIMULATION OF ANNUAL AGGREGATE LOSSES
# Covers Methodology Section 2.3.2 and Results Section 4.4
############################################################

set.seed(123)

M <- 10000
annual_loss_sim_base <- numeric(M)
storm_count_sim <- numeric(M)
max_wind_sim <- numeric(M)

for (s in 1:M) {
  
  # Step 1: simulate annual hurricane count
  N_s <- rpois(1, mu_hat)
  storm_count_sim[s] <- N_s
  
  if (N_s == 0) {
    annual_loss_sim_base[s] <- 0
    max_wind_sim[s] <- 0
  } else {
    
    # Step 2: simulate wind speeds from fitted Lognormal model
    W_s <- rlnorm(N_s, meanlog = mu_w_hat, sdlog = sigma_w_hat)
    max_wind_sim[s] <- max(W_s)
    
    # Step 4: STOCHASTIC wind-to-loss mapping (Insured Losses)
    epsilon_s <- rnorm(N_s, mean = 0, sd = sigma_eps)
    logL_s <- gamma0_hat + gamma1_hat * log(W_s) + epsilon_s
    L_s <- exp(logL_s)
    
    # Step 5: aggregate annual insured loss
    annual_loss_sim_base[s] <- sum(L_s)
  }
}

############################################################
# 6.1 OUTPUTS: BASELINE SIMULATED ANNUAL LOSS DISTRIBUTION
############################################################

sim_summary_base <- data.frame(
  Mean = mean(annual_loss_sim_base),
  Median = median(annual_loss_sim_base),
  SD = sd(annual_loss_sim_base),
  Q90 = as.numeric(quantile(annual_loss_sim_base, 0.90)),
  Q95 = as.numeric(quantile(annual_loss_sim_base, 0.95)),
  Q99 = as.numeric(quantile(annual_loss_sim_base, 0.99))
)

print(sim_summary_base)

# Compare analytical vs baseline simulated moments
moment_compare_base <- data.frame(
  Measure = c("Mean", "Variance", "Standard Deviation"),
  Analytical = c(
    agg_mean_hat,
    agg_var_hat,
    sqrt(agg_var_hat)
  ),
  Simulated = c(
    mean(annual_loss_sim_base),
    var(annual_loss_sim_base),
    sd(annual_loss_sim_base)
  )
)

print(moment_compare_base)
print(max_wind_sim)

############################################################
# 7. CAT BOND PRICING RESULTS (RISK-ADJUSTED Q-MEASURE)
# Covers Methodology Sections 2.4 - 2.6 and Results Section 4.5
############################################################

# Contract Parameters
d <- as.numeric(quantile(annual_loss_sim_base, 0.90)) # 95th Percentile Attachment
K <- as.numeric(quantile(max_wind_sim, 0.90))         # 95th Percentile Parametric Trigger
L_limit <- 20000                                      # $20 Billion Limit
r <- 0.05                                             # Risk-free rate (5%)
lambda_risk <- 0.50                                   # 50% Safety Loading for tail risk

cat("Calibrated Indemnity Attachment (d) = $", d/1000, "Billion\n")
cat("Calibrated Parametric Trigger (K) = ", K, "Knots\n")

# 1. Payoffs
indemnity_payoff <- pmin(pmax(annual_loss_sim_base - d, 0), L_limit)
parametric_payoff <- ifelse(max_wind_sim > K, L_limit, 0)
basis_risk <- parametric_payoff - indemnity_payoff

# 2. Pure expected loss (under P measure)
exp_payoff_ind <- mean(indemnity_payoff)
exp_payoff_par <- mean(parametric_payoff)

# 3. Risk-Adjusted Discounted Fair Price (Q-measure)
price_ind <- exp(-r) * exp_payoff_ind * (1 + lambda_risk)
price_par <- exp(-r) * exp_payoff_par * (1 + lambda_risk)

# 4. Trigger probabilities 
p_trigger_ind <- mean(annual_loss_sim_base > d)
p_trigger_par <- mean(max_wind_sim > K)

# 5. Efficiency Ratio (Section 2.6.3)
eff_ratio_ind <- price_ind / exp_payoff_ind
eff_ratio_par <- price_par / exp_payoff_par

# Summary Table
pricing_summary <- data.frame(
  Structure = c("Indemnity XoL", "Parametric Wind"),
  Trigger_Prob = c(p_trigger_ind, p_trigger_par),
  Expected_Payout_USDm = c(exp_payoff_ind, exp_payoff_par),
  Fair_Market_Price_USDm = c(price_ind, price_par),
  Efficiency_Ratio = c(eff_ratio_ind, eff_ratio_par),
  Mean_Basis_Risk = c(NA, mean(basis_risk))
)

print(pricing_summary)

############################################################
# Plot 1: Histogram of Baseline Simulated Losses
############################################################

sim_df_base <- data.frame(annual_loss = annual_loss_sim_base)

ggplot(sim_df_base, aes(x = annual_loss / 1000)) + 
  geom_histogram(bins = 50, fill = "#4A90E2", color = "white", alpha = 0.8) +
  geom_vline(xintercept = d / 1000, color = "#D0021B", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = (d / 1000) + 5, y = 800, 
           label = paste0("Attachment Point: $", round(d/1000, 1), "B"), 
           color = "#D0021B", hjust = 0, fontface = "bold") +
  coord_cartesian(xlim = c(0, quantile(sim_df_base$annual_loss, 0.99) / 1000)) +
  labs(
    title = "Simulated Annual Aggregate Insured Hurricane Losses (Baseline)",
    subtitle = "Showing the 95th percentile Indemnity Attachment Point",
    x = "Annual Aggregate Insured Loss (Billion USD)",
    y = "Frequency (Simulated Years)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.subtitle = element_text(face = "plain") 
  )

############################################################
# Plot 2: Density of baseline simulated annual losses
############################################################

ggplot(sim_df_base, aes(x = annual_loss)) +
  geom_density(linewidth = 1.2) +
  coord_cartesian(xlim = c(0, quantile(sim_df_base$annual_loss, 0.99))) +
  labs(
    title = "Density of Baseline Simulated Annual Aggregate Insured Losses",
    x = "Annual aggregate insured loss (million USD)",
    y = "Density"
  ) +
  theme_minimal()

############################################################
# Plot 3: ECDF of baseline simulated annual losses
############################################################

ggplot(sim_df_base, aes(x = annual_loss)) +
  stat_ecdf(geom = "step") +
  coord_cartesian(xlim = c(0, quantile(sim_df_base$annual_loss, 0.99))) +
  labs(
    title = "Empirical CDF of Baseline Simulated Annual Aggregate Insured Losses",
    x = "Annual aggregate insured loss (million USD)",
    y = "Cumulative probability"
  ) +
  theme_minimal()

############################################################
# Plot 4: COLORED Basis Risk Distribution
############################################################

br_df <- data.frame(Basis_Risk = basis_risk) %>%
  mutate(Risk_Type = case_when(
    Basis_Risk > 0 ~ "Parametric Over-pays",
    Basis_Risk < 0 ~ "Parametric Under-pays",
    TRUE ~ "Perfect Match (Zero Payouts)"
  ))

ggplot(br_df, aes(x = Basis_Risk / 1000, fill = Risk_Type)) +
  geom_histogram(bins = 40, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c(
    "Parametric Over-pays" = "#F5A623",      # Orange (Inefficiency)
    "Parametric Under-pays" = "#D0021B",     # Red (Dangerous Under-payment)
    "Perfect Match (Zero Payouts)" = "#7ED321" # Green (Aligned Outcomes)
  )) +
  labs(
    title = "Distribution of Basis Risk (Parametric minus Indemnity)",
    subtitle = "Positive values indicate the Parametric bond triggered while Insured Losses were lower",
    x = "Basis Risk (Billion USD)",
    y = "Frequency (Simulated Years)",
    fill = "Outcome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic", color = "grey30"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black")
  )
########################################################################################################################


########################################################################################################################

############################################################
# 8. SENSITIVITY ANALYSIS (CLIMATE IMPACTS)
# Covers Methodology Section 2.6 and Results Section 4.6
############################################################
lambda_risk <- 0.50  # The 50% safety loading from your methodology

# 1. Create a reusable simulation function (Strictly Aligned)
simulate_pricing_year <- function(sim_mu, sim_mu_w, sim_sigma_w, M = 10000) {
  
  sim_losses <- numeric(M)
  sim_max_wind <- numeric(M)
  
  for (s in 1:M) {
    N_s <- rpois(1, sim_mu)
    
    if (N_s == 0) {
      sim_losses[s] <- 0
      sim_max_wind[s] <- 0
    } else {
      # No arbitrary cap (160) - allow natural extreme tail events
      W_s <- rlnorm(N_s, meanlog = sim_mu_w, sdlog = sim_sigma_w)
      sim_max_wind[s] <- max(W_s)
      
      # Pure Stochastic Regression (Insured Losses)
      epsilon_s <- rnorm(N_s, mean = 0, sd = sigma_eps)
      logL_s <- gamma0_hat + gamma1_hat * log(W_s) + epsilon_s
      L_s <- exp(logL_s)
      
      sim_losses[s] <- sum(L_s)
    }
  }
  
  # Calculate expected payoffs using original 90th percentile triggers (d and K)
  ind_payoff <- pmin(pmax(sim_losses - as.numeric(d), 0), L_limit)
  par_payoff <- ifelse(sim_max_wind > K, L_limit, 0)
  
  return(list(
    Exp_Ind = mean(ind_payoff),
    Exp_Par = mean(par_payoff),
    Prob_Ind = mean(sim_losses > as.numeric(d)),
    Prob_Par = mean(sim_max_wind > K)
  ))
}

# 2. Define Climate Scenarios
# Scenario A: Baseline (using your strictly calibrated mu_hat)
base_res <- simulate_pricing_year(mu_hat, mu_w_hat, sigma_w_hat)

# Scenario B: +20% Frequency (e.g., warmer oceans leading to more storms)
freq_res <- simulate_pricing_year(mu_hat * 1.20, mu_w_hat, sigma_w_hat)

# Scenario C: +5% Severity (e.g., more energy leading to stronger peak winds)
sev_res <- simulate_pricing_year(mu_hat, mu_w_hat * 1.05, sigma_w_hat)

# 3. Compile Results Table with RISK LOADING (1 + lambda_risk)
sensitivity_table <- data.frame(
  Scenario = c("Baseline", "+20% Frequency", "+5% Severity (Wind)"),
  Indemnity_Prob = c(base_res$Prob_Ind, freq_res$Prob_Ind, sev_res$Prob_Ind),
  Parametric_Prob = c(base_res$Prob_Par, freq_res$Prob_Par, sev_res$Prob_Par),
  
  Indemnity_Price_USDm = exp(-r) * (1 + lambda_risk) * c(base_res$Exp_Ind, freq_res$Exp_Ind, sev_res$Exp_Ind),
  Parametric_Price_USDm = exp(-r) * (1 + lambda_risk) * c(base_res$Exp_Par, freq_res$Exp_Par, sev_res$Exp_Par)
)

print(sensitivity_table)

############################################################
# Plot 5 & 6: COMBINED SENSITIVITY PLOTS (Side-by-Side)
############################################################

# 1. Create the Probabilities Plot (Saved as 'p5')
plot_data_prob <- data.frame(
  Scenario = rep(sensitivity_table$Scenario, 2),
  Structure = rep(c("Indemnity XoL", "Parametric Bond"), each = 3),
  Probability = c(sensitivity_table$Indemnity_Prob, sensitivity_table$Parametric_Prob)
)
plot_data_prob$Scenario <- factor(plot_data_prob$Scenario, levels = c("Baseline", "+20% Frequency", "+5% Severity (Wind)"))

p5 <- ggplot(plot_data_prob, aes(x = Scenario, y = Probability, fill = Structure)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.9, color = "black", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(Probability * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Indemnity XoL" = "#4A90E2", "Parametric Bond" = "#F5A623")) +
  scale_y_continuous(limits = c(0, max(plot_data_prob$Probability) * 1.2)) +
  labs(
    title = "A) Trigger Probabilities",
    x = "Climate Scenario",
    y = "Probability of Triggering"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black")
  )

# 2. Create the Prices Plot (Saved as 'p6')
plot_data_price <- data.frame(
  Scenario = rep(sensitivity_table$Scenario, 2),
  Structure = rep(c("Indemnity XoL", "Parametric Bond"), each = 3),
  Price_Billion = c(sensitivity_table$Indemnity_Price_USDm, sensitivity_table$Parametric_Price_USDm) / 1000
)
plot_data_price$Scenario <- factor(plot_data_price$Scenario, levels = c("Baseline", "+20% Frequency", "+5% Severity (Wind)"))

p6 <- ggplot(plot_data_price, aes(x = Scenario, y = Price_Billion, fill = Structure)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.9, color = "black", linewidth = 0.8) +
  geom_text(aes(label = paste0("$", round(Price_Billion, 2), "B")), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Indemnity XoL" = "#4A90E2", "Parametric Bond" = "#F5A623")) +
  scale_y_continuous(limits = c(0, max(plot_data_price$Price_Billion) * 1.2)) + 
  labs(
    title = "B) Expected Fair Prices",
    x = "Climate Scenario",
    y = "Discounted Price (Billion USD)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black")
  )

# 3. Combine them side-by-side using Patchwork
combined_plot <- p5 + p6 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") 

# Display the final combined plot
print(combined_plot)


############################################################
# APPENDIX D: DIFFERENTIAL RISK LOADING
# Real-world pricing scenario (different lambdas)
############################################################

# Define the real-world lambdas
lambda_ind <- 0.60  # 60% loading for Indemnity (moral hazard / trapped capital)
lambda_par <- 0.40  # 40% loading for Parametric (transparency / rapid settlement)

# Recalculate Risk-Adjusted Discounted Fair Price (Q-measure)
price_ind_real <- exp(-r) * exp_payoff_ind * (1 + lambda_ind)
price_par_real <- exp(-r) * exp_payoff_par * (1 + lambda_par)

# Recalculate Efficiency Ratios
eff_ratio_ind_real <- price_ind_real / exp_payoff_ind
eff_ratio_par_real <- price_par_real / exp_payoff_par

# Create Summary Table for Appendix D
appendix_d_table <- data.frame(
  Structure = c("Indemnity XoL", "Parametric Wind"),
  Expected_Payout_USDm = c(exp_payoff_ind, exp_payoff_par),
  Risk_Loading = c("60%", "40%"),
  Fair_Market_Price_USDm = c(price_ind_real, price_par_real),
  Efficiency_Ratio = c(eff_ratio_ind_real, eff_ratio_par_real)
)

cat("\n--- APPENDIX D: DIFFERENTIAL RISK LOADING RESULTS ---\n")
print(appendix_d_table)


