# Actuarial Modeling of Claim Frequency in Motor Insurance (MTPL)

## Project Overview
This project focuses on predictive modeling of insurance claim frequency using the `freMTPL2freq` dataset. The primary goal is to develop and thoroughly evaluate robust actuarial pricing models. The workflow utilizes custom loss functions (e.g., Poisson and Gamma deviance), manual hyperparameter optimization with early stopping, and advanced actuarial evaluation metrics.

## Documentation & Presentation
The repository includes a comprehensive presentation (`presentation/Actuarial_Modeling_Presentation.pdf`) that provides:
* **Theoretical Background:** Mathematical formulations of the evaluated algorithms.
* **Implementation Details:** A breakdown of the cross-validation, early stopping, and parametric bootstrap procedures.
* **Results Analysis:** A visual and quantitative comparison of the models.

## Key Features & Methodology
* **Algorithms Implemented:** Null Model, Poisson GLM, Random Forest (with Parametric Bootstrap), Response Boosting, Gradient Boosting, XGBoost, and an Ensemble model.
* **Custom Implementations:** 
  * Custom optimization routines using Poisson deviance.
  * Tree building directed by analytical negative gradients (Taylor expansion) tailored specifically for Poisson distributions.
* **Actuarial Evaluation:** Models were evaluated using Forecast Dominance (Tweedie distribution parameter space), Murphy's Score Decomposition, Lift Plots, Concentration Curves, and Gini indices (classic and ML-adjusted).

## Data
The dataset contains 200,000 observations of motor third-party liability (MTPL) policies, sourced from the `CASdatasets` package.
* **Train / Valid / Test split:** 70% / 15% / 15%.
* **Feature Engineering:** Log transformations (Density), discretization of continuous variables (VehAge, DrivAge), and exposure adjustments.

## Key Results
The XGBoost algorithm outperformed other models across key metrics (Poisson Deviance, Discrimination, Gini index). It demonstrated superior capability in discriminating extreme portfolio segments compared to the baseline GLM.
