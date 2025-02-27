library(speff2trial)  
library(dplyr)        
library(ggplot2)      


# Dataset
data("ACTG175", package = "speff2trial")

# Relevant columns and handling missing values
actg_data = ACTG175 %>%
  select(cd420, age, wtkg, karnof, cd40) %>%
  na.omit()

head(actg_data)

# Noise levels
noise_levels = c(0.1, 0.5, 1.0, 2.0, 5.0)

# True Model
true_model = lm(cd420 ~ age + wtkg + karnof + cd40, data = actg_data)

# Coefficients, R-squared, and standard errors
true_coefs = coef(true_model)
r_squared = summary(true_model)$r.squared
SE = summary(true_model)$coefficients[, "Std. Error"]

print(true_coefs)
print(r_squared)
print(SE)

set.seed(123)

# Adding noise to predictors
add_noise = function(data, noise_sd) {
  data %>%
    mutate(
      age_noisy = age + rnorm(n(), mean = 0, sd = noise_sd),
      wtkg_noisy = wtkg + rnorm(n(), mean = 0, sd = noise_sd),
      karnof_noisy = karnof + rnorm(n(), mean = 0, sd = noise_sd),
      cd40_noisy = cd40 + rnorm(n(), mean = 0, sd = noise_sd)
    )
}

# Linear model with noisy predictors 
fit_model = function(data) {
  lm(cd420 ~ age_noisy + wtkg_noisy + karnof_noisy + cd40_noisy, data = data)
}


# Result with intercept
results_intercept = lapply(seq_along(noise_levels), function(i) {
  noise_sd = noise_levels[i]
  noisy_data = add_noise(actg_data, noise_sd)
  model = fit_model(noisy_data)
  
  data.frame(
    noise_sd = noise_sd,
    intercept = coef(model)["(Intercept)"],
    age_coef = coef(model)["age_noisy"],
    wtkg_coef = coef(model)["wtkg_noisy"],
    karnof_coef = coef(model)["karnof_noisy"],
    cd40_coef = coef(model)["cd40_noisy"],
    r_squared = summary(model)$r.squared
  )
}) %>%
  bind_rows()

# Re-assign simple numeric row names
rownames(results_intercept) = 1:nrow(results_intercept)
print(results_intercept)

# Results with no intercept(This is what we will be using from now due to high values of the intercept)
results_no_intercept = lapply(seq_along(noise_levels), function(i) {
  noise_sd = noise_levels[i]
  noisy_data = add_noise(actg_data, noise_sd)
  model = fit_model(noisy_data)
  
  data.frame(
    noise_sd = noise_sd,
    age_coef = coef(model)["age_noisy"],
    wtkg_coef = coef(model)["wtkg_noisy"],
    karnof_coef = coef(model)["karnof_noisy"],
    cd40_coef = coef(model)["cd40_noisy"],
    r_squared = summary(model)$r.squared
  )
}) %>%
  bind_rows()

# Re-assign simple numeric row names
rownames(results_no_intercept ) = 1:nrow(results_no_intercept )
print(results_no_intercept )

# Noise Coefficient estimates vs. noise level plot
coef_plot = results_no_intercept %>%
  pivot_longer(cols = c(age_coef, wtkg_coef, karnof_coef, cd40_coef),
               names_to = "predictor",
               values_to = "coefficient") %>%
  ggplot(aes(x = noise_sd, y = coefficient, color = predictor)) +
  geom_line(linewidth = 1) +  # Replace size with linewidth
  geom_point(size = 2) +
  labs(
    title = "Impact of Noise on Coefficient Estimates",
    x = "Noise Standard Deviation",
    y = "Coefficient Estimate",
    color = "Predictor"
  ) +
  theme_minimal()

print(coef_plot)

# Plot of R-squared vs. noise level
r2_plot = ggplot(results_no_intercept, aes(x = noise_sd, y = r_squared)) +
  geom_line(size = 1, color = "blue") +
  geom_point(size = 2, color = "blue") +
  labs(
    title = "Impact of Noise on Model R-Squared",
    x = "Noise Standard Deviation",
    y = "R-Squared"
  ) +
  theme_minimal()

print(r2_plot)

# Add the true coefficients to the results
results_long = results_no_intercept %>%
  pivot_longer(cols = c(age_coef, wtkg_coef, karnof_coef, cd40_coef),
               names_to = "predictor",
               values_to = "coefficient") %>%
  mutate(
    true_coefficient = case_when(
      predictor == "age_coef" ~ true_coefs["age"],
      predictor == "wtkg_coef" ~ true_coefs["wtkg"],
      predictor == "karnof_coef" ~ true_coefs["karnof"],
      predictor == "cd40_coef" ~ true_coefs["cd40"]
    )
  )

print(results_long, n= 24)

# Plot the noisy estimates along with the true coefficients
coef_with_true_plot = ggplot(results_long, aes(x = noise_sd, y = coefficient, color = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true_coefficient, color = predictor), linetype = "dashed", size = 0.8) +
  labs(
    title = "Noise on Coefficient Estimates with True Values",
    x = "Noise Standard Deviation",
    y = "Coefficient Estimate",
    color = "Predictor"
  ) +
  theme_minimal()

print(coef_with_true_plot)

# Function to compute unbiased estimates
compute_unbiased_estimates = function(data, noise_sd, true_data) {
  n = nrow(data)
  
  # Matrix X (predictors) and y (response)
  X = as.matrix(data %>% select(age_noisy, wtkg_noisy, karnof_noisy, cd40_noisy))
  y = as.matrix(data$cd420)
  
  # Calculate X'X / n and X'y / n
  XtX_n = t(X) %*% X / n
  Xty_n = t(X) %*% y / n
  
  S2 = diag(noise_sd^2, ncol(X)) # Variance of the noise
  
  # Correct for noise: (X'X / n - S2 * I)
  XtX_n_corrected = XtX_n - S2 #* diag(ncol(X))
  
  # Calculate unbiased coefficients
  beta_unbiased = solve(XtX_n_corrected) %*% Xty_n
  colnames(beta_unbiased) = "corrected_estimates"
  
  # Add predictor names
  beta_df = as.data.frame(beta_unbiased) %>%
    mutate(
      predictor = c("age_coef", "wtkg_coef", "karnof_coef", "cd40_coef"),
      noise_sd = noise_sd
    )
  
  return(beta_df)
}

# Calculate unbiased estimates for each noise level
unbiased_results = lapply(noise_levels, function(noise_sd) {
  noisy_data = add_noise(actg_data, noise_sd)
  compute_unbiased_estimates(noisy_data, noise_sd) # Removed `true_data`
}) %>%
  bind_rows()

unbiased_results

# Combine with true coefficients
unbiased_results = unbiased_results %>%
  mutate(
    true_coefficient = case_when(
      predictor == "age_coef" ~ true_coefs["age"],
      predictor == "wtkg_coef" ~ true_coefs["wtkg"],
      predictor == "karnof_coef" ~ true_coefs["karnof"],
      predictor == "cd40_coef" ~ true_coefs["cd40"]
    )
  )

unbiased_results


# Combine with true coefficients
unbiased_results = unbiased_results %>%
  mutate(
    true_coefficient = case_when(
      predictor == "age_coef" ~ true_coefs["age"],
      predictor == "wtkg_coef" ~ true_coefs["wtkg"],
      predictor == "karnof_coef" ~ true_coefs["karnof"],
      predictor == "cd40_coef" ~ true_coefs["cd40"]
    )
  )

# Transform to wide format (desired structure)
wide_table = unbiased_results %>%
  select(noise_sd, predictor, corrected_estimates) %>%
  pivot_wider(
    names_from = predictor,  # Column names will be the predictor values (e.g., age_coef)
    values_from = corrected_estimates  # Values for the new columns will be the corrected_estimates
  )

print(wide_table)


# Plot corrected estimates vs. true coefficients
corrected_plot = ggplot(unbiased_results, aes(x = noise_sd, y = corrected_estimates, color = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true_coefficient, color = predictor), linetype = "dashed", size = 0.8) +
  labs(
    title = "Unbiased Estimates Using Correction Formula",
    x = "Noise Standard Deviation",
    y = "Coefficient Estimate",
    color = "Predictor"
  ) +
  theme_minimal()

print(corrected_plot)

# Combine unbiased estimates with noisy estimates and true values
comparison_results = unbiased_results %>%
  select(predictor, noise_sd, corrected_estimates, true_coefficient) %>%
  left_join(
    results_long %>%
      select(predictor, noise_sd, coefficient) %>%
      rename(noisy_estimates = coefficient),
    by = c("predictor", "noise_sd")
  )

# Pivot data for plotting
comparison_long = comparison_results %>%
  pivot_longer(
    cols = c(corrected_estimates, noisy_estimates, true_coefficient),
    names_to = "estimate_type",
    values_to = "value"
  )

# Plot the comparison
comparison_plot = ggplot(comparison_long, aes(x = noise_sd, y = value, color = predictor, linetype = estimate_type)) +
  geom_line(size = 1) +
  geom_point(size = 2, aes(shape = estimate_type)) +
  labs(
    title = " Unbiased, Noisy, and True Coefficient Estimates",
    x = "Noise Standard Deviation",
    y = "Coefficient Estimate",
    color = "Predictor",
    linetype = "Estimate Type",
    shape = "Estimate Type"
  ) +
  theme_minimal()

print(comparison_plot)

# Function to calculate R-squared using corrected coefficients
compute_corrected_r2 = function(data, corrected_coefs) {
  n = nrow(data)
  
  # Matrix X (predictors) and y (response)
  X = as.matrix(data %>% select(age_noisy, wtkg_noisy, karnof_noisy, cd40_noisy))
  y = as.matrix(data$cd420)
  
  # Predicted values using corrected coefficients
  y_pred = X %*% corrected_coefs
  
  # Calculate R-squared
  ssr = sum((y - y_pred)^2)  # Sum of squared residuals
  sst = sum((y - mean(y))^2) # Total sum of squares
  r2 = 1 - ssr / sst
  
  return(r2)
}

# Calculate R-squared for corrected estimates at each noise level
r2_corrected = lapply(noise_levels, function(noise_sd) {
  noisy_data = add_noise(actg_data, noise_sd)
  corrected_coefs = compute_unbiased_estimates(noisy_data, noise_sd, actg_data) %>%
    select(corrected_estimates) %>%
    as.matrix()
  r2 = compute_corrected_r2(noisy_data, corrected_coefs)
  
  data.frame(noise_sd = noise_sd, r_squared_corrected = r2)
}) %>%
  bind_rows()

# Plot R-squared for corrected estimates without rounding
r2_corrected_plot = ggplot(r2_corrected, aes(x = noise_sd, y = r_squared_corrected)) +
  geom_line(size = 1, color = "green") +
  geom_point(size = 2, color = "green") +
  labs(
    title = "R-Squared of Corrected Estimates",
    x = "Noise Standard Deviation",
    y = "R-Squared (Corrected)"
  ) +
  theme_minimal()

print(r2_corrected_plot)


# Combine corrected and noisy R-squared values
# Join the R-squared values from noisy and corrected results
r2_combined = results_no_intercept %>%
  select(noise_sd, r_squared) %>%
  rename(r_squared_noisy = r_squared) %>%
  left_join(r2_corrected, by = "noise_sd") %>%
  pivot_longer(
    cols = c(r_squared_noisy, r_squared_corrected),
    names_to = "type",
    values_to = "r_squared"
  )

# Print r2_combined without rounding
options(digits = 15)  # Ensuring full precision display
print(r2_combined)

# Plot corrected R-squared and noisy R-squared without rounding
r2_combined_plot = ggplot(r2_combined, aes(x = noise_sd, y = r_squared, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Comparison of Corrected and Noisy R-Squared",
    x = "Noise Standard Deviation",
    y = "R-Squared",
    color = "Estimate Type"
  ) +
  theme_minimal()

print(r2_combined_plot)


# Function to calculate standard errors using variance-covariance matrix
compute_corrected_se = function(data, corrected_coefs) {
  n = nrow(data)
  
  # Matrix X (predictors) and y (response)
  X = as.matrix(data %>% select(age_noisy, wtkg_noisy, karnof_noisy, cd40_noisy))
  y = as.matrix(data$cd420)
  
  # Predicted values and residuals
  y_pred = X %*% corrected_coefs
  residuals = y - y_pred
  
  # Variance of residuals (sigma^2)
  sigma_sq = sum(residuals^2) / (n - ncol(X)) # Adjusted for degrees of freedom
  
  # Variance-covariance matrix
  var_covar_matrix = sigma_sq * solve(t(X) %*% X)
  
  # Standard errors (square root of diagonal elements)
  se = sqrt(diag(var_covar_matrix))
  
  return(se)
}

# Calculate SEs for each noise level
se_corrected = lapply(noise_levels, function(noise_sd) {
  noisy_data = add_noise(actg_data, noise_sd)
  corrected_coefs = compute_unbiased_estimates(noisy_data, noise_sd, actg_data) %>%
    select(corrected_estimates) %>%
    as.matrix()
  se = compute_corrected_se(noisy_data, corrected_coefs)
  
  data.frame(
    noise_sd = noise_sd,
    predictor = c("age_coef", "wtkg_coef", "karnof_coef", "cd40_coef"), # Explicit predictor names
    se_corrected = se
  )
}) %>%
  bind_rows()

head(se_corrected)


# Calculate SEs for each noise level
se_corrected = lapply(noise_levels, function(noise_sd) {
  # Add noise to the data
  noisy_data = add_noise(actg_data, noise_sd)
  
  # Compute corrected coefficients
  corrected_coefs = compute_unbiased_estimates(noisy_data, noise_sd, actg_data) %>%
    select(corrected_estimates) %>%
    as.matrix()
  
  # Compute standard errors
  se = compute_corrected_se(noisy_data, corrected_coefs)
  
  # Validate the length of 'se' (expecting 4)
  if (length(se) != 4) {
    stop("Standard errors vector length is not 4. Check the compute_corrected_se function.")
  }
  
  # Create a data frame to store SEs for each noise level
  se_df = data.frame(
    noise_sd = noise_sd,
    age_se = se[1],
    wtkg_se = se[2],
    karnof_se = se[3],
    cd40_se = se[4]
  )
  
  return(se_df)
}) %>%
  bind_rows()

head(se_corrected)

# Calculate SEs for each noise level
se_corrected = lapply(noise_levels, function(noise_sd) {
  noisy_data = add_noise(actg_data, noise_sd)
  corrected_coefs = compute_unbiased_estimates(noisy_data, noise_sd, actg_data) %>%
    select(corrected_estimates) %>%
    as.matrix()
  se = compute_corrected_se(noisy_data, corrected_coefs)
  
  # Create a data frame to store SEs for each noise level
  se_df = data.frame(
    noise_sd = noise_sd,
    age_se = se[1],
    wtkg_se = se[2],
    karnof_se = se[3],
    cd40_se = se[4]
  )
  
  return(se_df)
}) %>%
  bind_rows()

# Plot corrected standard errors
se_plot = ggplot(se_corrected, aes(x = noise_sd)) +
  geom_line(aes(y = age_se, color = "age"), size = 1) +
  geom_line(aes(y = wtkg_se, color = "wtkg"), size = 1) +
  geom_line(aes(y = karnof_se, color = "karnof"), size = 1) +
  geom_line(aes(y = cd40_se, color = "cd40"), size = 1) +
  geom_point(aes(y = age_se, color = "age"), size = 2) +
  geom_point(aes(y = wtkg_se, color = "wtkg"), size = 2) +
  geom_point(aes(y = karnof_se, color = "karnof"), size = 2) +
  geom_point(aes(y = cd40_se, color = "cd40"), size = 2) +
  labs(
    title = "Corrected Standard Errors vs. Noise Level",
    x = "Noise Standard Deviation",
    y = "Standard Error",
    color = "Predictor"
  ) +
  theme_minimal()

print(se_plot)


# Function to calculate SEs for the noisy estimates
compute_noisy_se = function(model) {
  # Extract standard errors from the variance-covariance matrix of the model
  sqrt(diag(vcov(model)))
}

# Calculate SEs for each noise level
se_noisy_results = lapply(noise_levels, function(noise_sd) {
  noisy_data = add_noise(actg_data, noise_sd)
  model = fit_model(noisy_data)
  se = compute_noisy_se(model)
  data.frame(
    noise_sd = noise_sd,
    age_se = se["age_noisy"],
    wtkg_se = se["wtkg_noisy"],
    karnof_se = se["karnof_noisy"],
    cd40_se = se["cd40_noisy"]
  )
}) %>%
  bind_rows()

# Pivot data for visualization
se_noisy_long = se_noisy_results %>%
  pivot_longer(cols = c(age_se, wtkg_se, karnof_se, cd40_se),
               names_to = "predictor",
               values_to = "standard_error")

# Plot SEs for noisy estimates
se_noisy_plot = ggplot(se_noisy_long, aes(x = noise_sd, y = standard_error, color = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Standard Errors of Noisy Estimates vs. Noise Level",
    x = "Noise Standard Deviation",
    y = "Standard Error",
    color = "Predictor"
  ) +
  theme_minimal()

print(se_noisy_plot)

# Combine noisy and corrected results for comparison
se_comparison = bind_rows(
  se_noisy_results %>% mutate(type = "Noisy"),
  se_corrected %>% mutate(type = "Corrected")
)

# Pivot data for plotting
se_comparison_long = se_comparison %>%
  pivot_longer(cols = c(age_se, wtkg_se, karnof_se, cd40_se),
               names_to = "predictor",
               values_to = "standard_error")

# Plot comparison of SEs(Corrected and Noisy)
ggplot(se_comparison_long, aes(x = noise_sd, y = standard_error, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~predictor, scales = "free_y") +
  labs(
    title = "Comparison of SEs: Noisy vs Corrected Estimates",
    x = "Noise Standard Deviation",
    y = "Standard Error",
    color = "Estimate Type"
  ) +
  theme_minimal()


