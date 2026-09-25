# Optional: Load environment from previous step
# load("02_modeling_output.RData")

# Utworzenie folderu na wykresy, jeśli nie istnieje
if(!dir.exists("plots")) dir.create("plots")

### 1. General Deviance Results
res_freq_test <- data_set_test$ClaimNb / data_set_test$Exposure
exp_freq_test <- data_set_test$Exposure

dev_null <- poisson_deviance_function(cbind(res_freq_test, pred_test_null, exp_freq_test))
dev_glm  <- poisson_deviance_function(cbind(res_freq_test, pred_test_glm, exp_freq_test))
dev_rf   <- poisson_deviance_function(cbind(res_freq_test, pred_test_rf, exp_freq_test))
dev_rb   <- poisson_deviance_function(cbind(res_freq_test, pred_test_rb, exp_freq_test))
dev_gb   <- poisson_deviance_function(cbind(res_freq_test, pred_test_gb, exp_freq_test))
dev_ens  <- poisson_deviance_function(cbind(res_freq_test, pred_test_ensemble, exp_freq_test))
dev_xgb  <- poisson_deviance_function(cbind(res_freq_test, pred_test_xgb, exp_freq_test))

deviance_results <- data.frame(
  Model = c("Null Model (Baseline)", 
            "Poisson GLM", 
            "Random Forest (Parametric)", 
            "Response Boosting", 
            "Gradient Boosting", 
            "Ensemble (Average)", 
            "XGBoost"),
  Poisson_Deviance_10e2 = c(dev_null, dev_glm, dev_rf, dev_rb, dev_gb, dev_ens, dev_xgb) * 100 
)

deviance_results <- deviance_results[order(deviance_results$Poisson_Deviance_10e2), ]
print("Poisson Deviance on Test Set:")
print(deviance_results)

### 2. Forecast Dominance Analysis
p_seq <- seq(1, 3, 0.1)
dev_table <- matrix(0, nrow = length(p_seq), ncol = 7)
preds_list <- list(pred_test_null, pred_test_glm, pred_test_rf, pred_test_rb, pred_test_gb, pred_test_ensemble, pred_test_xgb)

for (i in seq_along(p_seq)) {
  p <- p_seq[i]
  devs <- numeric(7)
  
  for (m_idx in 1:7) {
    if (p == 1) {
      devs[m_idx] <- poisson_deviance_function(cbind(res_freq_test, preds_list[[m_idx]], exp_freq_test))
    } else if (p == 2) {
      devs[m_idx] <- gamma_deviance_function(cbind(res_freq_test, preds_list[[m_idx]], exp_freq_test))
    } else {
      devs[m_idx] <- tweedie_deviance_function(p, cbind(res_freq_test, preds_list[[m_idx]], exp_freq_test))
    }
  }
  dev_table[i, ] <- devs
}

# Ratio to GLM (Model 2 is GLM)
ratio_matrix <- sweep(dev_table[, c(1, 3, 4, 5, 6, 7)], 1, dev_table[, 2], "/")

plot_titles <- c("Null Model", "Random Forest", "Response Boosting", 
                 "Gradient Boosting", "Ensemble", "XGBoost")
plot_colors <- c("black", "darkgreen", "orange", "purple", "brown", "cyan")

# --- ZAPIS WYKRESU: Forecast Dominance Grid ---
png("plots/forecast_dominance_grid.png", width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3))

for (j in 1:6) {
  y_limit <- c(min(0.95, min(ratio_matrix[, j])), max(1.05, max(ratio_matrix[, j])))
  plot(p_seq, ratio_matrix[, j], type = "l", col = plot_colors[j], lwd = 2, 
       ylim = y_limit, main = plot_titles[j], ylab = "Ratio to GLM", xlab = "Parameter p")
  abline(h = 1, col = "red", lwd = 2, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()
# ----------------------------------------------


# Dominance among top algorithms (Ratio to XGBoost)
ratio_rb  <- dev_table[, 4] / dev_table[, 7]
ratio_gb  <- dev_table[, 5] / dev_table[, 7]
ratio_ens <- dev_table[, 6] / dev_table[, 7]

y_limit_zoom <- c(min(c(ratio_rb, ratio_gb, ratio_ens)) - 0.0005, max(c(ratio_rb, ratio_gb, ratio_ens)) + 0.0005)

# --- ZAPIS WYKRESU: Forecast Dominance Zoom ---
png("plots/forecast_dominance_leading.png", width = 1200, height = 800, res = 150)
plot(p_seq, ratio_gb, type = "l", col = "purple", lwd = 2, 
     ylim = y_limit_zoom, main = "Forecast Dominance - Leading Algorithms", 
     ylab = "Deviance Ratio to XGBoost", xlab = "Parameter p (Tweedie)")
lines(p_seq, ratio_rb, col = "orange", lwd = 2)
lines(p_seq, ratio_ens, col = "brown", lwd = 2)
abline(h = 1, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Gradient Boosting", "Response Boosting", "Ensemble"),
       col = c("purple", "orange", "brown"), lty = 1, lwd = 2, cex = 0.9)
dev.off()
# ----------------------------------------------


### 3. Actuarial Predictor Evaluation

model_predictions <- list(
  "Null_Model"        = pred_test_null,
  "GLM"               = pred_test_glm,
  "Random_Forest"     = pred_test_rf,
  "Response_Boosting" = pred_test_rb,
  "Gradient_Boosting" = pred_test_gb,
  "Ensemble"          = pred_test_ensemble,
  "XGBoost"           = pred_test_xgb
)

evaluation_metrics <- data.frame(
  Model = character(),
  Global_Balance_Expected = numeric(),
  Murphy_Score = numeric(),
  Murphy_Reliability = numeric(),
  Murphy_Resolution = numeric(),
  ABC = numeric(),
  Gini_ML_Original = numeric(),
  Gini_ML_Calibrated = numeric(),
  stringsAsFactors = FALSE
)

# Suppress plot generation during calculation
pdf(file = tempfile())

for (model_name in names(model_predictions)) {
  pred <- model_predictions[[model_name]]
  
  expected_freq <- sum(exp_freq_test * pred) / sum(exp_freq_test)
  murphy <- score_decomposition(res_freq_test, pred, exp_freq_test, "Poisson")
  
  z <- concentration_curve(res_freq_test, pred, exp_freq_test)
  
  fit_autocalibrate <- auto_calibrate(res_freq_test, pred, exp_freq_test)
  pred_calibrated <- fit_autocalibrate(pred) + 10^(-6)
  
  z_calibrated <- concentration_curve(res_freq_test, pred_calibrated, exp_freq_test)
  
  new_row <- data.frame(
    Model = model_name,
    Global_Balance_Expected = expected_freq * 100,
    Murphy_Score = murphy[1],
    Murphy_Reliability = murphy[2],
    Murphy_Resolution = murphy[3],
    ABC = z$abc,
    Gini_ML_Original = z$gini_ml,
    Gini_ML_Calibrated = z_calibrated$gini_ml
  )
  
  evaluation_metrics <- rbind(evaluation_metrics, new_row)
}

dev.off()

observed_frequency <- sum(data_set_test$ClaimNb) / sum(exp_freq_test)
cat(sprintf("Observed portfolio frequency (%%): %.4f\n", observed_frequency * 100))
print("Evaluation Metrics Summary:")
print(evaluation_metrics)


### 4. Graphical Evaluation

# --- ZAPIS WYKRESU: CORP / Lift Plots ---
png("plots/corp_lift_plots.png", width = 2400, height = 1200, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (model_name in names(model_predictions)) {
  pred <- model_predictions[[model_name]]
  model_name_clean <- gsub("_", " ", model_name)
  
  exp_vs_obs_plot(res_freq_test, pred, exp_freq_test, 10)
  title(main = model_name_clean, line = 0.15, cex.main = 1)
}
par(mfrow = c(1, 1))
dev.off()
# ----------------------------------------


# --- ZAPIS WYKRESU: Autocalibration Plots ---
png("plots/autocalibration_plots.png", width = 2400, height = 1200, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (model_name in names(model_predictions)) {
  pred <- model_predictions[[model_name]]
  model_name_clean <- gsub("_", " ", model_name)
  
  fit_autocalibrate <- auto_calibrate(res_freq_test, pred, exp_freq_test)
  pred_sort <- sort(pred)
  pred_calibrated <- fit_autocalibrate(pred_sort) + 10^(-6)
  
  plot(pred_calibrated ~ pred_sort, type = "s", col = "royalblue", lwd = 2,
       main = paste("Calibration:", model_name_clean), 
       xlab = "Original Predictor", ylab = "Calibrated Predictor")
  abline(0, 1, col = "firebrick", lwd = 2, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()
# --------------------------------------------


# Collect Data for CC and LC curves
results_original <- list()
results_calibrated <- list()

pdf(file = tempfile()) 
for (model_name in names(model_predictions)) {
  pred <- model_predictions[[model_name]]
  
  results_original[[model_name]] <- concentration_curve(res_freq_test, pred, exp_freq_test)
  fit_auto <- auto_calibrate(res_freq_test, pred, exp_freq_test)
  pred_cal <- fit_auto(pred) + 10^(-6)
  
  results_calibrated[[model_name]] <- concentration_curve(res_freq_test, pred_cal, exp_freq_test)
}
dev.off() 


# --- ZAPIS WYKRESU: Concentration Curves ---
png("plots/concentration_curves.png", width = 2400, height = 1200, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (model_name in names(model_predictions)) {
  model_name_clean <- gsub("_", " ", model_name)
  
  z_orig <- results_original[[model_name]]
  z_cal  <- results_calibrated[[model_name]]
  
  lorenz_before <- round(z_orig$gini_eco / 2, 2)
  conc_before   <- round(z_orig$abc + lorenz_before, 2)
  
  lorenz_after <- round(z_cal$gini_eco / 2, 2)
  conc_after   <- round(z_cal$abc + lorenz_after, 2)
  
  x_axis <- seq(1, 0, length.out = length(z_orig$concentration))
  
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Cumulative Portfolio", ylab = "Cumulative Claims",
       main = paste("Concentration:", model_name_clean), cex.main = 1.1)
  
  abline(0, 1, col = "black", lwd = 1.5) 
  lines(x_axis, z_orig$concentration, col = "royalblue", lwd = 2, lty = 2)
  lines(x_axis, z_cal$concentration, col = "darkorange", lwd = 2, lty = 1)
  
  legend("topleft", 
         legend = c(sprintf("Before Calibration (%.2f)", conc_before), 
                    sprintf("After Calibration (%.2f)", conc_after)),
         col = c("royalblue", "darkorange"), lty = c(2, 1), lwd = 2, cex = 0.85, bty = "n")
}
par(mfrow = c(1, 1))
dev.off()
# -------------------------------------------


# --- ZAPIS WYKRESU: Lorenz Curves ---
png("plots/lorenz_curves.png", width = 2400, height = 1200, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (model_name in names(model_predictions)) {
  model_name_clean <- gsub("_", " ", model_name)
  
  z_orig <- results_original[[model_name]]
  z_cal  <- results_calibrated[[model_name]]
  
  lorenz_before <- round(z_orig$gini_eco / 2, 2)
  lorenz_after  <- round(z_cal$gini_eco / 2, 2)
  
  x_axis <- seq(1, 0, length.out = length(z_orig$lorenz))
  
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Cumulative Portfolio", ylab = "Cumulative Predictions",
       main = paste("Lorenz:", model_name_clean), cex.main = 1.1)
  
  abline(0, 1, col = "black", lwd = 1.5) 
  lines(x_axis, z_orig$lorenz, col = "forestgreen", lwd = 2, lty = 2)
  lines(x_axis, z_cal$lorenz, col = "purple", lwd = 2, lty = 1)
  
  legend("topleft", 
         legend = c(sprintf("Before Calibration (%.2f)", lorenz_before), 
                    sprintf("After Calibration (%.2f)", lorenz_after)),
         col = c("forestgreen", "purple"), lty = c(2, 1), lwd = 2, cex = 0.85, bty = "n")
}
par(mfrow = c(1, 1))
dev.off()
# ------------------------------------


# --- ZAPIS WYKRESU: XGBoost vs GLM Predictions ---
png("plots/xgboost_vs_glm.png", width = 1200, height = 1000, res = 150)
plot(pred_test_xgb ~ pred_test_glm, col = "blue", cex = 0.5,
     xlab = "Prediction: XGBoost Model", ylab = "Prediction: GLM Model",
     main = "Predictor Comparison: XGBoost vs GLM")
abline(0, 1, col = "red", lwd = 2)
abline(v = quantile(pred_test_xgb, 0.95), col = "green", lty = 2, lwd = 2)

add_loc_reg(y = pred_test_xgb, x = pred_test_glm, sp = 0.1, d = 1)

legend("topleft", 
       legend = c("Ideal match (y=x)", "95th percentile XGBoost", "Mean trend"),
       col = c("red", "green", "green"), lty = c(1, 2, 1), lwd = c(2, 2, 2), cex = 0.8)
dev.off()
# -------------------------------------------------