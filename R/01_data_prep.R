options(scipen = 999)

# Source helper functions (formerly moje_funkcje_2.R)
# source("00_helpers.R")

library(CASdatasets)
library(monotone)
library(data.table)

### 1. Data Loading & Merging

data(freMTPL2freq)
data_set <- freMTPL2freq[, -2]
data_set$VehGas <- factor(data_set$VehGas)

data(freMTPL2sev)
sev <- freMTPL2sev
sev$ClaimNb <- 1

sev <- aggregate(cbind(ClaimAmount, ClaimNb) ~ IDpol, data = sev, FUN = sum)
names(sev)[2] <- "ClaimTotal"

data_set <- merge(x = data_set, y = sev, by = "IDpol", all.x = TRUE)
data_set[is.na(data_set)] <- 0

# Capping ClaimNb and Exposure
data_set <- data_set[which(data_set$ClaimNb <= 5), ]
data_set$Exposure <- pmin(data_set$Exposure, 1)

# Factor leveling
data_set$VehBrand <- factor(data_set$VehBrand,
                            levels = c("B1", "B2", "B3", "B4", "B5", "B6",
                                       "B10", "B11", "B12", "B13", "B14"))

levels(data_set$Region) <- c("R1", "R2", "R3", "R4", "R5", "R6",
                             "R7", "R8", "R9", "R10", "R11",
                             "R12", "R13", "R14", "R15", "R16", 
                             "R17", "R18", "R19", "R20", "R21", "R22")

### 2. Train, Validation, and Test Split

no_subsample <- 200000
frac_train   <- 0.70
frac_valid   <- 0.15
frac_test    <- 0.15

set.seed(100)
index <- sample(seq_len(nrow(data_set)), no_subsample, replace = FALSE)
data_set <- data_set[index, ]

shuffled_indices <- sample(seq_len(no_subsample))
cut_train <- floor(no_subsample * frac_train)
cut_valid <- cut_train + floor(no_subsample * frac_valid)

index_train <- shuffled_indices[1:cut_train]
index_valid <- shuffled_indices[(cut_train + 1):cut_valid]
index_test  <- shuffled_indices[(cut_valid + 1):no_subsample]
set.seed(NULL)

# Portfolio frequencies summary
cat("Train set frequency:", sum(data_set$ClaimNb[index_train]) / sum(data_set$Exposure[index_train]), "\n")
cat("Valid set frequency:", sum(data_set$ClaimNb[index_valid]) / sum(data_set$Exposure[index_valid]), "\n")
cat("Test set frequency:",  sum(data_set$ClaimNb[index_test]) / sum(data_set$Exposure[index_test]), "\n")

### 3. Feature Engineering

data_set$AreaGLM       <- as.integer(data_set$Area)
data_set$DensityGLM    <- log(data_set$Density)
data_set$VehPowerGLM   <- as.factor(pmin(data_set$VehPower, 9))
data_set$VehAgeGLM     <- as.factor(cut(data_set$VehAge, c(0, 5, 12, 101),
                                        labels = c("0-5", "6-12", "12+"),
                                        include.lowest = TRUE))
data_set$DrivAgeGLM    <- as.factor(cut(data_set$DrivAge, c(18, 20, 25, 30, 40, 50, 70, 101),
                                        labels = c("18-20", "21-25", "26-30", "31-40", "41-50", "51-70", "71+"),
                                        include.lowest = TRUE))
data_set$BonusMalusGLM <- pmin(data_set$BonusMalus, 150)

# Final Dataset Extraction
data_set_train      <- data_set[index_train, ]
data_set_valid      <- data_set[index_valid, ]
data_set_test       <- data_set[index_test, ]
data_set_train_full <- rbind(data_set_train, data_set_valid)

# Optional: Save environment to pass to modeling script
# save.image("01_data_prep_output.RData")