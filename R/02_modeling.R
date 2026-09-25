# Optional: Load environment from previous step
# load("01_data_prep_output.RData")

library(rpart)

### 1. Null Model
model_poisson_null <- glm(ClaimNb ~ 1 + offset(log(Exposure)), 
                          family = poisson(link = "log"), 
                          data = data_set_train_full)

pred_test_null <- predict(model_poisson_null, newdata = data_set_test, type = "response") / data_set_test$Exposure

### 2. GLM Model
model_poisson_glm <- glm(ClaimNb ~ offset(log(Exposure)) + AreaGLM + VehPowerGLM + Region + 
                           VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
                         family = poisson(link = "log"), data = data_set_train_full)

pred_test_glm <- predict(model_poisson_glm, newdata = data_set_test, type = "response") / data_set_test$Exposure

### 3. Random Forest (Parametric Bootstrap)
all_features <- c("AreaGLM", "VehPowerGLM", "Region", "VehAgeGLM", 
                  "DrivAgeGLM", "BonusMalusGLM", "DensityGLM", "VehBrand", "VehGas")

rf_grid <- expand.grid(
  maxdepth = c(5, 10),
  minbucket = c(10, 50),
  m_features = c(5, 7)
)

n_cv_groups <- 5
n_trees <- 200

set.seed(123)
cv_group <- sample(rep(1:n_cv_groups, length.out = nrow(data_set_train)))
set.seed(NULL)

cv_results_rf <- numeric(nrow(rf_grid))

for (g in 1:nrow(rf_grid)) {
  current_maxdepth <- rf_grid$maxdepth[g]
  current_minbucket <- rf_grid$minbucket[g]
  current_m <- rf_grid$m_features[g]
  
  err_group <- numeric(n_cv_groups)
  
  for (k in 1:n_cv_groups) {
    idx_val <- which(cv_group == k)
    idx_train <- setdiff(1:nrow(data_set_train), idx_val)
    
    expected_claims_cv <- predict(model_poisson_glm, newdata = data_set_train[idx_train, ], type = "response")
    forest_predictions <- rep(0, length(idx_val))
    
    for (m in 1:n_trees) {
      simulated_ClaimNb <- rpois(n = length(idx_train), lambda = expected_claims_cv)
      
      parametric_train_data <- data_set_train[idx_train, ]
      parametric_train_data$ClaimNb <- simulated_ClaimNb
      
      sampled_features <- sample(all_features, current_m, replace = FALSE)
      tree_formula <- as.formula(paste("cbind(Exposure, ClaimNb) ~", paste(sampled_features, collapse = " + ")))
      
      model_tree <- rpart(
        tree_formula,
        data = parametric_train_data,
        method = "poisson",
        control = rpart.control(xval = 0, maxdepth = current_maxdepth, minbucket = current_minbucket, cp = 10^(-4))
      )
      
      forest_predictions <- forest_predictions + predict(model_tree, newdata = data_set_train[idx_val, ])
    }
    
    forest_predictions_avg <- forest_predictions / n_trees
    
    data_deviance <- cbind(
      data_set_train$ClaimNb[idx_val] / data_set_train$Exposure[idx_val],
      forest_predictions_avg,
      data_set_train$Exposure[idx_val]
    )
    err_group[k] <- poisson_deviance_function(data_deviance)
  }
  cv_results_rf[g] <- mean(err_group)
  cat(sprintf("Random Forest CV - Combination %d/%d - Loss: %.4f\n", g, nrow(rf_grid), cv_results_rf[g] * 100))
}

best_rf_params <- rf_grid[which.min(cv_results_rf), ]

# Final RF Estimation
final_forest_predictions <- 0
expected_claims_final <- predict(model_poisson_glm, newdata = data_set_train_full, type = "response")

for (m in 1:n_trees) {
  simulated_ClaimNb_final <- rpois(n = nrow(data_set_train_full), lambda = expected_claims_final)
  
  parametric_data <- data_set_train_full
  parametric_data$ClaimNb <- simulated_ClaimNb_final
  
  sampled_features <- sample(all_features, best_rf_params$m_features, replace = FALSE)
  tree_formula <- as.formula(paste("cbind(Exposure, ClaimNb) ~", paste(sampled_features, collapse = " + ")))
  
  model_tree <- rpart(
    tree_formula,
    data = parametric_data,
    method = "poisson",
    control = rpart.control(xval = 0, maxdepth = best_rf_params$maxdepth, minbucket = best_rf_params$minbucket, cp = 10^(-4))
  )
  final_forest_predictions <- final_forest_predictions + predict(model_tree, newdata = data_set_test)
}
pred_test_rf <- final_forest_predictions / n_trees

### 4. Response Boosting

rb_grid <- expand.grid(maxdepth = c(2, 4), minbucket = c(10, 50), alpha = c(0.05, 0.1))
M_max <- 100 
best_cv_loss <- Inf
best_rb_params <- NULL

for (g in 1:nrow(rb_grid)) {
  cv_loss_rb <- matrix(0, nrow = n_cv_groups, ncol = M_max)
  
  for (k in 1:n_cv_groups) {
    idx_val <- which(cv_group == k)  
    idx_train <- setdiff(1:nrow(data_set_train), idx_val)
    
    cv_pred <- sum(data_set_train$ClaimNb[idx_train]) / sum(data_set_train$Exposure[idx_train])
    exposure_boost <- data_set_train$Exposure[idx_train] * cv_pred
    
    for (m in 1:M_max) {
      model_boost <- rpart(
        cbind(exposure_boost, ClaimNb) ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
        data = data_set_train[idx_train, ], 
        method = "poisson",
        control = rpart.control(xval = 0, maxdepth = rb_grid$maxdepth[g], minbucket = rb_grid$minbucket[g], cp = 10^(-6))
      )
      
      cv_pred <- cv_pred * (predict(model_boost, newdata = data_set_train[idx_val, ]))^rb_grid$alpha[g]
      exposure_boost <- exposure_boost * (predict(model_boost))^rb_grid$alpha[g]
      
      data_deviance <- cbind(data_set_train$ClaimNb[idx_val] / data_set_train$Exposure[idx_val], cv_pred, data_set_train$Exposure[idx_val])
      cv_loss_rb[k, m] <- poisson_deviance_function(data_deviance)
    }
  }
  
  mean_cv_loss <- min(apply(cv_loss_rb, 2, mean))
  if (mean_cv_loss < best_cv_loss) {
    best_cv_loss <- mean_cv_loss
    best_rb_params <- rb_grid[g, ]
  }
}

# Early Stopping (Validation Set)
init_freq <- sum(data_set_train$ClaimNb) / sum(data_set_train$Exposure)
exp_train_es <- data_set_train$Exposure * init_freq
pred_valid_es <- rep(init_freq, nrow(data_set_valid))
valid_loss_rb <- numeric(M_max)

for (m in 1:M_max) {
  model_boost_es <- rpart(
    cbind(exp_train_es, ClaimNb) ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
    data = data_set_train, method = "poisson",
    control = rpart.control(xval = 0, maxdepth = best_rb_params$maxdepth, minbucket = best_rb_params$minbucket, cp = 10^(-6))
  )
  exp_train_es <- exp_train_es * (predict(model_boost_es))^best_rb_params$alpha
  pred_valid_es <- pred_valid_es * (predict(model_boost_es, newdata = data_set_valid))^best_rb_params$alpha
  
  data_deviance_val <- cbind(data_set_valid$ClaimNb / data_set_valid$Exposure, pred_valid_es, data_set_valid$Exposure)
  valid_loss_rb[m] <- poisson_deviance_function(data_deviance_val)
}
opt_M_rb <- which.min(valid_loss_rb)

# Final RB Estimation
final_freq <- sum(data_set_train_full$ClaimNb) / sum(data_set_train_full$Exposure)
final_exp_rb <- data_set_train_full$Exposure * final_freq
pred_test_rb <- rep(final_freq, nrow(data_set_test))

for (m in 1:opt_M_rb) {
  model_boost_final <- rpart(
    cbind(final_exp_rb, ClaimNb) ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
    data = data_set_train_full, method = "poisson",
    control = rpart.control(xval = 0, maxdepth = best_rb_params$maxdepth, minbucket = best_rb_params$minbucket, cp = 10^(-6))
  )
  final_exp_rb <- final_exp_rb * (predict(model_boost_final))^best_rb_params$alpha
  pred_test_rb <- pred_test_rb * (predict(model_boost_final, newdata = data_set_test))^best_rb_params$alpha
}

### 5. Gradient Boosting
gb_grid <- expand.grid(maxdepth = c(2, 4), minbucket = c(10, 50), alpha = c(0.05, 0.1))
best_gb_loss <- Inf
best_gb_params <- NULL

for (g in 1:nrow(gb_grid)) {
  cv_loss_gb <- matrix(0, nrow = n_cv_groups, ncol = M_max)
  
  for (k in 1:n_cv_groups) {
    idx_val <- which(cv_group == k)  
    idx_train <- setdiff(1:nrow(data_set_train), idx_val)
    
    mean_freq <- sum(data_set_train$ClaimNb[idx_train]) / sum(data_set_train$Exposure[idx_train])
    cv_pred_train <- rep(mean_freq, length(idx_train))
    cv_pred_val <- rep(mean_freq, length(idx_val))
    
    for (m in 1:M_max) {
      grad <- (data_set_train$ClaimNb[idx_train] / data_set_train$Exposure[idx_train]) - cv_pred_train
      
      model_gb <- rpart(
        grad ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
        data = data_set_train[idx_train, ], method = "anova", weights = Exposure,
        control = rpart.control(xval = 0, maxdepth = gb_grid$maxdepth[g], minbucket = gb_grid$minbucket[g], cp = 10^(-6))
      )
      
      leaves_train <- predict(model_gb)
      leaves_val <- predict(model_gb, newdata = data_set_train[idx_val, ])
      
      unique_leaves <- sort(unique(leaves_train))
      est_leaves <- numeric(length(unique_leaves))
      
      for (i in seq_along(unique_leaves)) {
        l <- unique_leaves[i]
        idx_l <- which(leaves_train == l)
        num <- sum(data_set_train$ClaimNb[idx_train[idx_l]])
        den <- sum(cv_pred_train[idx_l] * data_set_train$Exposure[idx_train[idx_l]])
        est_leaves[i] <- num / den
      }
      
      est_leaves_mat <- cbind(unique_leaves, est_leaves + 10^(-6))
      
      update_train <- est_leaves_mat[match(leaves_train, est_leaves_mat[,1]), 2]
      cv_pred_train <- cv_pred_train * (update_train)^gb_grid$alpha[g]
      
      update_val <- est_leaves_mat[match(leaves_val, est_leaves_mat[,1]), 2]
      cv_pred_val <- cv_pred_val * (update_val)^gb_grid$alpha[g]
      
      data_deviance <- cbind(data_set_train$ClaimNb[idx_val] / data_set_train$Exposure[idx_val], cv_pred_val, data_set_train$Exposure[idx_val])
      cv_loss_gb[k, m] <- poisson_deviance_function(data_deviance)
    }
  }
  
  mean_cv_loss <- min(apply(cv_loss_gb, 2, mean))
  if (mean_cv_loss < best_gb_loss) {
    best_gb_loss <- mean_cv_loss
    best_gb_params <- gb_grid[g, ]
  }
}

# Early Stopping
gb_train_es <- rep(init_freq, nrow(data_set_train))
gb_valid_es <- rep(init_freq, nrow(data_set_valid))
valid_loss_gb <- numeric(M_max)

for (m in 1:M_max) {
  grad_es <- (data_set_train$ClaimNb / data_set_train$Exposure) - gb_train_es
  
  model_gb_es <- rpart(
    grad_es ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
    data = data_set_train, method = "anova", weights = Exposure,
    control = rpart.control(xval = 0, maxdepth = best_gb_params$maxdepth, minbucket = best_gb_params$minbucket, cp = 10^(-6))
  )
  
  leaves_train_es <- predict(model_gb_es)
  leaves_val_es <- predict(model_gb_es, newdata = data_set_valid)
  
  unique_leaves_es <- sort(unique(leaves_train_es))
  est_leaves_es <- numeric(length(unique_leaves_es))
  
  for (i in seq_along(unique_leaves_es)) {
    l <- unique_leaves_es[i]
    idx_l <- which(leaves_train_es == l)
    est_leaves_es[i] <- sum(data_set_train$ClaimNb[idx_l]) / sum(gb_train_es[idx_l] * data_set_train$Exposure[idx_l])
  }
  
  est_mat_es <- cbind(unique_leaves_es, est_leaves_es + 10^(-6))
  
  update_train_es <- est_mat_es[match(leaves_train_es, est_mat_es[,1]), 2]
  gb_train_es <- gb_train_es * (update_train_es)^best_gb_params$alpha
  
  update_val_es <- est_mat_es[match(leaves_val_es, est_mat_es[,1]), 2]
  gb_valid_es <- gb_valid_es * (update_val_es)^best_gb_params$alpha
  
  data_deviance_val <- cbind(data_set_valid$ClaimNb / data_set_valid$Exposure, gb_valid_es, data_set_valid$Exposure)
  valid_loss_gb[m] <- poisson_deviance_function(data_deviance_val)
}
opt_M_gb <- which.min(valid_loss_gb)

# Final GB Estimation
gb_train_final <- rep(final_freq, nrow(data_set_train_full))
pred_test_gb <- rep(final_freq, nrow(data_set_test))

for (m in 1:opt_M_gb) {
  grad_final <- (data_set_train_full$ClaimNb / data_set_train_full$Exposure) - gb_train_final
  
  model_gb_final <- rpart(
    grad_final ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
    data = data_set_train_full, method = "anova", weights = Exposure,
    control = rpart.control(xval = 0, maxdepth = best_gb_params$maxdepth, minbucket = best_gb_params$minbucket, cp = 10^(-6))
  )
  
  leaves_train_fin <- predict(model_gb_final)
  leaves_test_fin <- predict(model_gb_final, newdata = data_set_test)
  
  unique_leaves_fin <- sort(unique(leaves_train_fin))
  est_leaves_fin <- numeric(length(unique_leaves_fin))
  
  for (i in seq_along(unique_leaves_fin)) {
    l <- unique_leaves_fin[i]
    idx_l <- which(leaves_train_fin == l)
    est_leaves_fin[i] <- sum(data_set_train_full$ClaimNb[idx_l]) / sum(gb_train_final[idx_l] * data_set_train_full$Exposure[idx_l])
  }
  
  est_mat_fin <- cbind(unique_leaves_fin, est_leaves_fin + 10^(-6))
  
  update_train_fin <- est_mat_fin[match(leaves_train_fin, est_mat_fin[,1]), 2]
  gb_train_final <- gb_train_final * (update_train_fin)^best_gb_params$alpha
  
  update_test_fin <- est_mat_fin[match(leaves_test_fin, est_mat_fin[,1]), 2]
  pred_test_gb <- pred_test_gb * (update_test_fin)^best_gb_params$alpha
}

### 6. XGBoost 
xgb_grid <- expand.grid(maxdepth = c(2, 4), lambda = c(1, 10), alpha = c(0.05, 0.1))
best_xgb_loss <- Inf
best_xgb_params <- NULL

for (g in 1:nrow(xgb_grid)) {
  cv_loss_xgb <- matrix(0, nrow = n_cv_groups, ncol = M_max)
  
  for (k in 1:n_cv_groups) {
    idx_val <- which(cv_group == k)  
    idx_train <- setdiff(1:nrow(data_set_train), idx_val)
    
    mean_freq <- sum(data_set_train$ClaimNb[idx_train]) / sum(data_set_train$Exposure[idx_train])
    eta_train <- rep(log(mean_freq), length(idx_train))
    eta_val <- rep(log(mean_freq), length(idx_val))
    
    for (m in 1:M_max) {
      mu_train <- exp(eta_train)
      R_train <- mu_train - (data_set_train$ClaimNb[idx_train] / data_set_train$Exposure[idx_train])
      H_train <- mu_train
      
      target_train <- -R_train / H_train
      weights_train <- data_set_train$Exposure[idx_train] * H_train
      
      model_xgb <- rpart(
        target_train ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
        data = data_set_train[idx_train, ], method = "anova", weights = weights_train,
        control = rpart.control(xval = 0, maxdepth = xgb_grid$maxdepth[g], minbucket = 10, cp = 10^(-6))
      )
      
      leaves_train <- predict(model_xgb)
      leaves_val <- predict(model_xgb, newdata = data_set_train[idx_val, ])
      
      unique_leaves <- sort(unique(leaves_train))
      est_leaves <- numeric(length(unique_leaves))
      
      for (i in seq_along(unique_leaves)) {
        l <- unique_leaves[i]
        idx_l <- which(leaves_train == l)
        sum_vR <- sum(data_set_train$Exposure[idx_train[idx_l]] * R_train[idx_l])
        sum_vH <- sum(data_set_train$Exposure[idx_train[idx_l]] * H_train[idx_l])
        est_leaves[i] <- -sum_vR / (sum_vH + 0.5 * xgb_grid$lambda[g])
      }
      
      est_mat <- cbind(unique_leaves, est_leaves)
      
      eta_train <- eta_train + xgb_grid$alpha[g] * est_mat[match(leaves_train, est_mat[,1]), 2]
      eta_val <- eta_val + xgb_grid$alpha[g] * est_mat[match(leaves_val, est_mat[,1]), 2]
      
      data_deviance <- cbind(data_set_train$ClaimNb[idx_val] / data_set_train$Exposure[idx_val], exp(eta_val), data_set_train$Exposure[idx_val])
      cv_loss_xgb[k, m] <- poisson_deviance_function(data_deviance)
    }
  }
  
  mean_cv_loss <- min(apply(cv_loss_xgb, 2, mean))
  if (mean_cv_loss < best_xgb_loss) {
    best_xgb_loss <- mean_cv_loss
    best_xgb_params <- xgb_grid[g, ]
  }
}

# Early Stopping
eta_train_es <- rep(log(init_freq), nrow(data_set_train))
eta_valid_es <- rep(log(init_freq), nrow(data_set_valid))
valid_loss_xgb <- numeric(M_max)

for (m in 1:M_max) {
  mu_train_es <- exp(eta_train_es)
  R_train_es <- mu_train_es - (data_set_train$ClaimNb / data_set_train$Exposure)
  H_train_es <- mu_train_es
  
  target_train_es <- -R_train_es / H_train_es
  weights_train_es <- data_set_train$Exposure * H_train_es
  
  model_xgb_es <- rpart(
    target_train_es ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
    data = data_set_train, method = "anova", weights = weights_train_es,
    control = rpart.control(xval = 0, maxdepth = best_xgb_params$maxdepth, minbucket = 10, cp = 10^(-6))
  )
  
  leaves_train_es <- predict(model_xgb_es)
  leaves_valid_es <- predict(model_xgb_es, newdata = data_set_valid)
  
  unique_leaves_es <- sort(unique(leaves_train_es))
  est_leaves_es <- numeric(length(unique_leaves_es))
  
  for (i in seq_along(unique_leaves_es)) {
    l <- unique_leaves_es[i]
    idx_l <- which(leaves_train_es == l)
    sum_vR <- sum(data_set_train$Exposure[idx_l] * R_train_es[idx_l])
    sum_vH <- sum(data_set_train$Exposure[idx_l] * H_train_es[idx_l])
    est_leaves_es[i] <- -sum_vR / (sum_vH + 0.5 * best_xgb_params$lambda)
  }
  
  est_mat_es <- cbind(unique_leaves_es, est_leaves_es)
  
  eta_train_es <- eta_train_es + best_xgb_params$alpha * est_mat_es[match(leaves_train_es, est_mat_es[,1]), 2]
  eta_valid_es <- eta_valid_es + best_xgb_params$alpha * est_mat_es[match(leaves_valid_es, est_mat_es[,1]), 2]
  
  data_deviance_val <- cbind(data_set_valid$ClaimNb / data_set_valid$Exposure, exp(eta_valid_es), data_set_valid$Exposure)
  valid_loss_xgb[m] <- poisson_deviance_function(data_deviance_val)
}
opt_M_xgb <- which.min(valid_loss_xgb)

# Final XGBoost Estimation
eta_train_final <- rep(log(final_freq), nrow(data_set_train_full))
eta_test_final <- rep(log(final_freq), nrow(data_set_test))

for (m in 1:opt_M_xgb) {
  mu_train_fin <- exp(eta_train_final)
  R_train_fin <- mu_train_fin - (data_set_train_full$ClaimNb / data_set_train_full$Exposure)
  H_train_fin <- mu_train_fin
  
  target_train_fin <- -R_train_fin / H_train_fin
  weights_train_fin <- data_set_train_full$Exposure * H_train_fin
  
  model_xgb_fin <- rpart(
    target_train_fin ~ AreaGLM + VehPowerGLM + Region + VehAgeGLM + DrivAgeGLM + BonusMalusGLM + DensityGLM + VehBrand + VehGas,
    data = data_set_train_full, method = "anova", weights = weights_train_fin,
    control = rpart.control(xval = 0, maxdepth = best_xgb_params$maxdepth, minbucket = 10, cp = 10^(-6))
  )
  
  leaves_train_fin <- predict(model_xgb_fin)
  leaves_test_fin <- predict(model_xgb_fin, newdata = data_set_test)
  
  unique_leaves_fin <- sort(unique(leaves_train_fin))
  est_leaves_fin <- numeric(length(unique_leaves_fin))
  
  for (i in seq_along(unique_leaves_fin)) {
    l <- unique_leaves_fin[i]
    idx_l <- which(leaves_train_fin == l)
    sum_vR <- sum(data_set_train_full$Exposure[idx_l] * R_train_fin[idx_l])
    sum_vH <- sum(data_set_train_full$Exposure[idx_l] * H_train_fin[idx_l])
    est_leaves_fin[i] <- -sum_vR / (sum_vH + 0.5 * best_xgb_params$lambda)
  }
  
  est_mat_fin <- cbind(unique_leaves_fin, est_leaves_fin)
  
  eta_train_final <- eta_train_final + best_xgb_params$alpha * est_mat_fin[match(leaves_train_fin, est_mat_fin[,1]), 2]
  eta_test_final <- eta_test_final + best_xgb_params$alpha * est_mat_fin[match(leaves_test_fin, est_mat_fin[,1]), 2]
}

pred_test_xgb <- exp(eta_test_final)

### 7. Ensemble Model
pred_test_ensemble <- (pred_test_rf + pred_test_rb + pred_test_gb) / 3

# Optional: Save workspace to carry over predictions
# save.image("02_modeling_output.RData")