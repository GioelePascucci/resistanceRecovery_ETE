# ==============================================================================
# Script title: main.R
# Description:  The script includes all the main analysis of the article,
#               included PGLS, plot generation, literature analysis, diagnostic 
#               test and sensitivity analyses.
# Manuscript title: Correlated evolution of thermal tolerance and egg size 
#                links resistance and recovery to extreme heat in amphibians.
# Manuscript authors: Gioele Pascucci, Miguel Tejedo, Lumír Gvoždík
# Script author: Gioele Pascucci
# Contact:  gioelepascucci@live.it
# Date: August 30, 2026
# R Version:  4.5.3
# ==============================================================================

rm(list = ls())
options(warn = -1)

library(readxl)     
library(dplyr)     
library(tidyr)      
library(tibble)     
library(ape)       
library(phytools)
library(phangorn)
library(geiger)     
library(nlme)       
library(MuMIn)      
library(mvMORPH)    
library(car)        
library(ggplot2)    
library(patchwork)  
library(gghalves)   
library(ggtree)     
library(sf)         
library(terra)     
library(rnaturalearth)     
library(rnaturalearthdata)

# ------------------------------------------------------------------------------
# Load and prepare the dataset and phylogeny
# ------------------------------------------------------------------------------

# Load the trait dataset, adjust to your path
main_data <- read_excel("../data/datasetTraits.xlsx")
main_data <- as.data.frame(main_data)
# Ensure missing values are stored as NA and not as string "NA"
main_data[main_data == "NA"] <- NA

# Ensure species names have underscores instead of spaces for phylogeny matching
main_data$Species <- gsub(" ", "_", main_data$Species)
rownames(main_data) <- main_data$Species

# Security check for missing data (NA)
df_trait <- main_data %>%
  drop_na(CT_max, Maturity, Abs_Latitude, Adult_SVL, Clutch_size, Egg_size, Altitude)

# Set categorical variables as factors
df_trait$Family <- as.factor(df_trait$Family)
df_trait$Order <- as.factor(df_trait$Order)
df_trait$Species <- as.factor(df_trait$Species)

# Load the Jetz & Pyron (2018) phylogeny, adjust to your path (Download link:)
trees <- read.nexus("../data/Jetz&Pyron2018.nex")

# Create a Maximum Clade Credibility (MCC) tree
consensus_tree <- maxCladeCred(trees)

# Replace spaces with underscores in tip labels
consensus_tree$tip.label <- gsub(" ", "_", consensus_tree$tip.label)

# Align the tree's taxonomy with Frost (2026)
consensus_tree$tip.label <- gsub("Litoria_chloris", "Chlorohyla_chloris", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Litoria_aurea", "Ranoidea_aurea", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Litoria_lesueurii", "Rhyaconastes_lesueurii", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Hyla_andersonii", "Dryophytes_andersonii", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Hyla_chrysoscelis", "Dryophytes_chrysoscelis", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Hyla_cinerea", "Dryophytes_cinereus", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Hyla_versicolor", "Dryophytes_versicolor", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Rana_catesbeiana", "Lithobates_catesbeianus", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Rana_sylvatica", "Lithobates_sylvaticus", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Rana_virgatipes", "Lithobates_virgatipes", consensus_tree$tip.label)
consensus_tree$tip.label <- gsub("Rhinella_schneideri", "Rhinella_diptycha", consensus_tree$tip.label)

# Check for mismatches between the tree and the dataset
check <- name.check(consensus_tree, df_trait)

# Prune the tree
if (length(check$tree_not_data) > 0) {
  pruned_tree <- drop.tip(consensus_tree, check$tree_not_data)
} else {
  pruned_tree <- consensus_tree
}

# Force the tree to be ultrametric
pruned_tree <- force.ultrametric(pruned_tree, method = "extend")

cat("Tree successfully pruned, Remaining species:", Ntip(pruned_tree), "\n")




# ------------------------------------------------------------------------------
# Phylogenetic signal and evolutionary model selection
# ------------------------------------------------------------------------------

#' Run Evolutionary Models and phylogenetic signal on a trait
#' @param trait_name Trait under analysis, matching the column name in the dataframe
#' @param df Dataframe containing the trait and species names
#' @param phylo Pruned phylogenetic tree
analyze_trait <- function(trait_name, df, phylo) {
  
  cat(sprintf("\n=========================================================\n"))
  cat(sprintf("   ANALYZING TRAIT: %s", toupper(trait_name)))
  cat(sprintf("\n=========================================================\n"))
  
  # Extract named vector for the trait
  trait_vector <- df[[trait_name]]
  names(trait_vector) <- df$Species
  
  # Fit evolutionary models
  cat("Fitting evolutionary models (BM, EB, OU, Lambda)...\n")
  
  # Suppress geiger's warnings during model fitting
  suppressWarnings({
    fitBM <- fitContinuous(phylo, trait_vector, model = "BM")
    fitEB <- fitContinuous(phylo, trait_vector, model = "EB")
    # Setting alpha bounds to accommodate both physiology and life history traits
    fitOU <- fitContinuous(phylo, trait_vector, model = "OU", bounds = list(alpha = c(1e-5, 100))) 
    fitL  <- fitContinuous(phylo, trait_vector, model = "lambda")
  })
  
  # Calculate AICc weights
  aicc_scores <- setNames(
    c(fitBM$opt$aicc, fitEB$opt$aicc, fitOU$opt$aicc, fitL$opt$aicc),
    c("BM", "EB", "OU", "Lambda")
  )
  weights <- aic.w(aicc_scores)
  
  # Compile model results into a dataframe
  model_df <- data.frame(
    Model  = c("Brownian Motion (BM)", "Early Burst (EB)", "Ornstein-Uhlenbeck (OU)", "Lambda"),
    logLik = c(fitBM$opt$lnL, fitEB$opt$lnL, fitOU$opt$lnL, fitL$opt$lnL),
    AICc   = c(fitBM$opt$aicc, fitEB$opt$aicc, fitOU$opt$aicc, fitL$opt$aicc),
    Weight = as.numeric(weights)
  ) %>% arrange(AICc) # Sort best model to top
  
  # Calculate phylogenetic signal
  cat("Calculating Phylogenetic Signal...\n")
  k_sig <- phylosig(phylo, trait_vector, test = TRUE, nsim = 10000)
  lambda_sig <- phylosig(phylo, trait_vector, method = "lambda", test = TRUE)
  
  # Compile signal results into a dataframe
  sig_df <- data.frame(
    Metric  = c("Blomberg's K", "Pagel's Lambda"),
    Value   = c(k_sig$K, lambda_sig$lambda),
    P_Value = c(k_sig$P, lambda_sig$P)
  )
  
  # Print result
  cat("\n--- Best Fitting Evolutionary Model ---\n")
  print(model_df, digits = 5)
  
  cat("\n--- Phylogenetic Signal Testing ---\n")
  print(sig_df, digits = 5)
  cat("\n")
}

# Pass each trait through the custom function just defined. 
analyze_trait("Adult_SVL", df_trait, pruned_tree)

analyze_trait("CT_max", df_trait, pruned_tree)

analyze_trait("Clutch_size", df_trait, pruned_tree)

analyze_trait("Egg_size", df_trait, pruned_tree)

analyze_trait("Maturity", df_trait, pruned_tree)


# ------------------------------------------------------------------------------
# PGLS model selection and estimates
# ------------------------------------------------------------------------------
cat("\n--- Fitting PGLS Models (Maximum Likelihood) ---\n")

# Define Correlation Structures
cor_lambda <- corPagel(value = 0.5, phy = pruned_tree, form = ~Species, fixed = FALSE)
cor_BM     <- corBrownian(phy = pruned_tree, form = ~Species)
cor_OU     <- corMartins(value = 1e-08, phy = pruned_tree, form = ~Species, fixed = FALSE)
cor_EB     <- corBlomberg(value = 0.5, phy = pruned_tree, form = ~Species, fixed = TRUE)

# Fit models using Maximum Likelihood (ML) for model comparison
base_formula <- CT_max ~ Adult_SVL + Clutch_size + Egg_size + Maturity + Abs_Latitude + Altitude

pgls_lambda <- gls(base_formula, data = df_trait, correlation = cor_lambda, method = "ML")
pgls_OU     <- gls(base_formula, data = df_trait, correlation = cor_OU, method = "ML")
pgls_EB     <- gls(base_formula, data = df_trait, correlation = cor_EB, method = "ML")
pgls_BM     <- gls(base_formula, data = df_trait, correlation = cor_BM, method = "ML")

# Compare Models
aicc_evolution_average <- setNames(
  c(MuMIn::AICc(pgls_BM), MuMIn::AICc(pgls_lambda), MuMIn::AICc(pgls_EB), MuMIn::AICc(pgls_OU)), 
  c("BM", "Lambda", "EB", "OU")
)

# Calculate weights and print output using MuMIn
weights_evolution_average <- MuMIn::Weights(aicc_evolution_average)
model_selection_df <- data.frame(
  Model = names(aicc_evolution_average),
  AICc = aicc_evolution_average,
  Weight = as.numeric(weights_evolution_average)
) %>% arrange(AICc)

# Lambda resulted best model
cat("\n--- Best phylogenetic correlation structure ---\n")
print(model_selection_df, digits = 4)

# Refit the best model (Lambda) using REML 
pgls_final<- gls(base_formula, 
                       data = df_trait, 
                       correlation = cor_lambda, 
                       method = "REML")

# Function to proper formatting CI, handling small coefficient values
format_ci <- function(est, lower, upper, threshold = 0.1, digits = 2) {
  # Return NA if any input value is missing
  if (is.na(est) || is.na(lower) || is.na(upper)) {
    return(NA_character_)
  }
  
  # Check the maximum absolute value among the inputs
  max_abs <- max(abs(c(est, lower, upper)), na.rm = TRUE)
  
  # Format with x10^-2 multiplier if maximum value is below the threshold
  if (max_abs < threshold && max_abs > 0) {
    fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f] \u00D710\u207B\u00B2")
    sprintf(fmt, est * 100, lower * 100, upper * 100)
  } else {
    fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
    sprintf(fmt, est, lower, upper)
  }
}

# Extract coefficients table
reml_summary <- as.data.frame(summary(pgls_final)$tTable)

# Calculate 95% CI (Estimate +/- 1.96 * SE)
reml_summary$Lower <- reml_summary$Value - (1.96 * reml_summary$Std.Error)
reml_summary$Upper <- reml_summary$Value + (1.96 * reml_summary$Std.Error)

# Extract R2-lik 
pgls_full_ML <- gls(CT_max ~ Adult_SVL + Clutch_size + Egg_size + Maturity + Abs_Latitude + Altitude, 
                    data = df_trait, 
                    correlation = cor_lambda, 
                    method = "ML")

# Null model 
pgls_null_ML <- gls(CT_max ~ 1, 
                    data = df_trait, 
                    correlation = cor_lambda, 
                    method = "ML")
r2_fixed_only <- r.squaredLR(pgls_full_ML, null = pgls_null_ML)

# Extract Pagel's lambda
lambda_val <- coef(pgls_final$modelStruct$corStruct, unconstrained = FALSE)
# Get Egg_size coefficient
egg_coef <- coef(pgls_final)["Egg_size"]
# Calculate Egg_size range across the dataset
egg_range <- max(df_trait$Egg_size, na.rm = TRUE) - min(df_trait$Egg_size, na.rm = TRUE)
# Calculate CTmax variation across Egg_size range
total_decline <- egg_coef * egg_range

# Format output
reml_summary$formatted_effect <- mapply(format_ci, 
                                        reml_summary$Value, 
                                        reml_summary$Lower, 
                                        reml_summary$Upper)

# Print formatted table
pgls_output <- data.frame(
  Term = rownames(reml_summary),
  Effect_95CI = reml_summary$formatted_effect,
  t_value = round(reml_summary$`t-value`, 2),
  P_Value = format.pval(reml_summary$`p-value`, eps = 0.001, digits = 3)
)

cat("\n--- PGLS results ---\n")
cat("\n-------------------------------\n")
cat(sprintf("\nFixed effects (excluding phylogeny) explained %.1f%% of the variation (R2lik).\n", r2_fixed_only * 100))
cat(sprintf("Phylogenetic signal (λ = %.2f).\n", lambda_val))
cat(sprintf("Predicted CTmax change across Egg Size range: %.2f °C\n", total_decline))
print(pgls_output)



# ------------------------------------------------------------------------------
# Multivariate evolutionary variance and pairwise evolutionary correlation (Jackknife 1000, 95% CIs)
# ------------------------------------------------------------------------------
cat("\n--- Multivariate Evolutionary Variance (1000 Jackknife Iterations) ---\n")

traits_orig <- c("Adult_SVL", "Maturity", "Clutch_size", "Egg_size", "CT_max")
traits_new  <- c("SVL", "Age at maturity", "Clutch size", "Egg size", "CTmax") 

# Z-scale trait values to ensure comparability
multi_matrix <- as.matrix(scale(df_trait[, traits_orig]))
colnames(multi_matrix) <- traits_new
# Keep species as row names, ensure correspondence between matrix row names and tree tips
rownames(multi_matrix) <- rownames(df_trait)
# Set number of iterations
n_boot <- 1000
# Set number of species to keep at each iteration
n_tips_keep <- floor(0.90 * Ntip(pruned_tree)) 

# Define a matrix to store the evolutionary variance, of each trait, for every iteration
boot_multi_sigma <- matrix(NA, nrow = n_boot, ncol = length(traits_new))
colnames(boot_multi_sigma) <- traits_new

# Define a three-dimensional array to store the correlation between all traits, for every iteration
boot_multi_cor <- array(NA, dim = c(length(traits_new), length(traits_new), n_boot))
dimnames(boot_multi_cor) <- list(traits_new, traits_new, NULL)


for (i in 1:n_boot) {
  if (i %% 100 == 0) cat("Iteration", i, "of", n_boot, "...\n")
  
  # Define species sub-sample our tree, keep only 90% of original species
  boot_tips <- sample(pruned_tree$tip.label, size = n_tips_keep, replace = FALSE)
  # Prune the tree to our sub-sample
  boot_tree <- keep.tip(pruned_tree, boot_tips)
  # Filter our trait dataset keeping only species in our sub-sample
  boot_data_matrix <- multi_matrix[boot_tips, ]
  
  tryCatch({
    
    # Fit the multivariate generalized least squares
    boot_fit <- mvgls(boot_data_matrix ~ 1, tree = boot_tree, model = "lambda", method = "LL")
    
    # Extracts the estimated evolutionary covariance matrix
    b_sigma <- boot_fit$sigma$Pinv
    
    # Extracts the diagonal elements of the covariance matrix (individual trait variances)
    boot_multi_sigma[i, ] <- diag(b_sigma)
    
    # Converts the evolutionary covariance matrix into a correlation matrix
    boot_multi_cor[, , i] <- cov2cor(b_sigma)
    
  }, error = function(e) {
    # Handle non-invertible covariance matrices caused by resampling artifacts
  })
}

# Calculate Means and 95% CI for variances
multi_df <- as.data.frame(boot_multi_sigma) %>%
  pivot_longer(cols = everything(), names_to = "Trait", values_to = "SigmaSq") %>%
  drop_na() 

multi_summary <- multi_df %>%
  group_by(Trait) %>%
  summarize(
    Mean = mean(SigmaSq, na.rm = TRUE),
    Lower = quantile(SigmaSq, probs = 0.025, na.rm = TRUE),
    Upper = quantile(SigmaSq, probs = 0.975, na.rm = TRUE)
  ) %>%
  arrange(Mean)

# Format result for variances
multi_summary$formatted_CI <- mapply(format_ci, 
                                     multi_summary$Mean, 
                                     multi_summary$Lower, 
                                     multi_summary$Upper)



# Extract pairwise evolutionary correlation, from correlation matrix
cor_pairs_list <- list()

for (r in 2:length(traits_new)) {
  for (c in 1:(r - 1)) {
    trait1 <- traits_new[r]
    trait2 <- traits_new[c]
    
    # Extract the 1000 jackknife estimates for a specific trait pair
    pair_cor_vals <- boot_multi_cor[r, c, ]
    pair_cor_vals <- pair_cor_vals[!is.na(pair_cor_vals)] # Drop failed iterations
    
    if (length(pair_cor_vals) > 0) {
      cor_pairs_list[[length(cor_pairs_list) + 1]] <- data.frame(
        Trait_1 = trait1,
        Trait_2 = trait2,
        Mean = mean(pair_cor_vals),
        Lower = as.numeric(quantile(pair_cor_vals, probs = 0.025)),
        Upper = as.numeric(quantile(pair_cor_vals, probs = 0.975))
      )
    }
  }
}

cor_summary <- bind_rows(cor_pairs_list)

# Format result
cor_summary$formatted_CI <- mapply(format_ci, 
                                   cor_summary$Mean, 
                                   cor_summary$Lower, 
                                   cor_summary$Upper)

# Print final output
cat("\n--- Evolutionary Variance Estimates ---\n")
print(multi_summary %>% select(Trait, formatted_CI))
cat("\n----------------------------------------------------------------------\n")
cat("\n--- Pairwise Evolutionary Correlations (95% CI) ---\n")
print(cor_summary %>% select(Trait_1, Trait_2, formatted_CI))



# ------------------------------------------------------------------------------
# Diasgnostic test and sensitivity analysis
# ------------------------------------------------------------------------------

# Check multicollinearity (VIF)
cat("\n--- Checking Phylogenetic Multicollinearity (Phylo-VIF) ---\n")

phylo_vif_values <- car::vif(pgls_final)

phylo_vif_data <- as.data.frame(phylo_vif_values) %>%
  rownames_to_column(var = "term") %>%
  rename(Phylo_VIF = phylo_vif_values) %>%
  arrange(Phylo_VIF) %>%
  mutate(
    term = factor(term, levels = term),
    is_collinear = Phylo_VIF > 2 # Flag values > 2
  )

print(phylo_vif_data)


# PGLS diagnostics
cat("\n--- Running Diagnostic Tests on PGLS Model ---\n")

# Extract normalized residuals and fitted values
model_res <- as.numeric(residuals(pgls_final, type = "normalized"))
model_fit <- as.numeric(fitted(pgls_final))

# Create a dataframe for ggplot
diag_data <- data.frame(
  Fitted = model_fit,
  Residuals = model_res
)

# Plot A: Residuals vs Fitted
p_res_fit <- ggplot(diag_data, aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6, size = 2, color = "black") +
  geom_hline(yintercept = 0, color = "#b30000", linetype = "dashed", linewidth = 1) +
  theme_bw(base_family = "sans") +
  labs(title = "Residuals vs Fitted", 
       x = "Fitted Values", 
       y = "Normalized Residuals") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Plot B: Q-Q Plot
p_qq <- ggplot(diag_data, aes(sample = Residuals)) +
  stat_qq(alpha = 0.6, size = 2) +
  stat_qq_line(color = "#b30000", linewidth = 1) +
  theme_bw(base_family = "sans") +
  labs(title = "Normal Q-Q Plot", 
       x = "Theoretical Quantiles", 
       y = "Sample Quantiles") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Combine plot
diagnostic_patchwork <- p_res_fit + p_qq + 
  plot_annotation(tag_levels = 'A')

# Print the patchwork plot
print(diagnostic_patchwork)


# Phylogenetic signal within residual
# Attach species names to residuals vector 
names(model_res) <- rownames(df_trait)
# Test for phylogenetic signal (using Pagel's Lambda)
resid_lambda <- phylosig(tree = pruned_tree, x = model_res, method = "lambda", test = TRUE)
# Print the results
cat("\n--- Checking Phylogenetic Signal in Residuals ---\n")
print(resid_lambda)


# Sensitivity analysis (Consider only species standardized at the same experimental conditions)

# Keep only species with same experimental conditions
# (Acclimation time=7 days; Acclimation temperature=20 °C: endpoint= LRR)
df_stdExp <- main_data %>%
  drop_na(CT_max, Maturity, Abs_Latitude, Adult_SVL, Clutch_size, Egg_size, Altitude, Std_CTmax) %>% 
  filter(Std_CTmax==TRUE)
# Prune the tree
check <- name.check(consensus_tree, df_stdExp)
if (length(check$tree_not_data) > 0) {
  pruned_sensitivity <- drop.tip(consensus_tree, check$tree_not_data)
} else {
  pruned_sensitivity <- consensus_tree
}
# Force the tree to be ultrametric
pruned_sensitivity <- force.ultrametric(pruned_sensitivity, method = "extend")

# Calculate correlation structures and PGLS
cor_lambda_s <- corPagel(value = 0.5, phy = pruned_sensitivity, form = ~Species, fixed = FALSE)
pgls_stdExp <- gls(base_formula, data = df_stdExp, correlation = cor_lambda_s, method = "REML")

# Extract coefficients table
reml_summary_s <- as.data.frame(summary(pgls_stdExp)$tTable)

# Calculate 95% CI (Estimate +/- 1.96 * SE)
reml_summary_s$Lower <- reml_summary_s$Value - (1.96 * reml_summary$Std.Error)
reml_summary_s$Upper <- reml_summary_s$Value + (1.96 * reml_summary$Std.Error)

# Format output
reml_summary_s$formatted_effect <- mapply(format_ci, 
                                        reml_summary_s$Value, 
                                        reml_summary_s$Lower, 
                                        reml_summary_s$Upper)

# Print formatted table
pgls_output_s<- data.frame(
  Term = rownames(reml_summary_s),
  Effect_95CI = reml_summary_s$formatted_effect,
  t_value = round(reml_summary_s$`t-value`, 2),
  P_Value = format.pval(reml_summary_s$`p-value`, eps = 0.001, digits = 3)
)

cat("\n--- Sensitivity analysis: only species with CTmax standardized at the same experimental conditions ---\n")
cat("Tree successfully pruned, Remaining species:", Ntip(pruned_sensitivity), "\n")
cat("\n--- PGLS Coefficients & 95% CIs ---\n")
print(pgls_output_s)




# ------------------------------------------------------------------------------
# Plots (Figure 1, Figure 2, Figure 3, Figure S2)
# ------------------------------------------------------------------------------

# Figure 1 (Literature analysis)
# Load literature analysis result, adjust to your path
lit_data <- read_excel("../data/literatureAnalysis.xlsx")

# Set explicit factor levels and calculate row counts per category
df_counts <- lit_data %>%
  mutate(Category = factor(Category, levels = c("Recovery", "Resistance", "Both"))) %>%
  group_by(Category, .drop = FALSE) %>%
  summarize(Count = n(), .groups = "drop")

# Define plot features
plot_literature <- ggplot(df_counts, aes(x = Category, y = Count)) +
  geom_col(fill = "#31688E", width = 0.6) + 
  labs(
    x = "Trait category", 
    y = "Number of papers"
  ) +
  theme_classic(base_size = 20, base_family = "sans") +
  theme(
    axis.title.y = element_text(face = "plain"),
    axis.text.x  = element_text(face = "plain", color = "black"),
    axis.text.y  = element_text(face = "plain", color = "black"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Display result
print(df_counts)
# Display plot
print(plot_literature)
# Export plot
ggsave(
  "../figure/figure1.pdf", 
  plot = plot_literature, 
  width = 6, 
  height = 5, 
  dpi = 600
)
ggsave(
  "../figure/figure1.png", 
  plot = plot_literature, 
  width = 6, 
  height = 5, 
  dpi = 600
)

#Figure 2 (Partial residual)
# Define standardized theme for all partial effect plots
partial_theme <- theme_classic(base_family = "sans", base_size = 14) +
  theme(
    axis.title.x = element_text(size = 16, color = "black", face = "plain"),
    axis.text    = element_text(color = "grey20", size = 12),
    plot.margin  = margin(10, 15, 10, 10)
  )

# Define color palette for plot elements
col_ribbon <- "#31688E" # Color for the 95% Confidence Interval ribbon
col_dots <- "#153347" # Color for the partial residual points

# Define a function to calculates partial residuals and fitted values for a specific predictor, 
# holding all other predictors in the model at their mean.
# pred_var: predictor variable, model: pgls model for coefficient extraction, data: dataset with predictor values
partial_effect <- function(pred_var, model, data) {
  
  # Extract model matrix, coefficients, variance-covariance matrix, and residuals
  X   <- model.matrix(formula(model), data = data)
  b   <- coef(model)
  V   <- vcov(model) 
  res <- as.numeric(residuals(model, type = "response"))
  
  # Generate a prediction grid, selected predictor varies, others are held at their mean
  # Generate a sequence of x-coordinates within our predictor range
  grid <- seq(min(X[, pred_var]), max(X[, pred_var]), length.out = 5000)
  # Define a new matrix to store the future fitted line coordinates
  # All predictors are set at their mean values except our selected predictor
  newX <- matrix(colMeans(X), nrow = length(grid), ncol = ncol(X), byrow = TRUE,
                 dimnames = list(NULL, colnames(X)))
  newX[, "(Intercept)"] <- 1
  # Remove mean values for the selected predictor and substitute with x-coordinates
  newX[, pred_var] <- grid
  
  # Calculate predicted fit (y-coordinates) using model slopes (b)
  fit <- as.numeric(newX %*% b)
  # Add standard error
  se <- sqrt(rowSums((newX %*% V) * newX))
  
  # Calculate the baseline effect of all other predictors
  others <- setdiff(names(b), pred_var)
  mu_others <- sum(b[others] * colMeans(X)[others])
  
  # Return dataframes for plotting the fitted line and partial residual
  list(
    # Fitted line and 95% confidence interval
    line = data.frame(
      X = grid, 
      predicted = fit,
      ymin = fit - 1.96 * se, 
      ymax = fit + 1.96 * se
    ),
    # Point coordinates (partial residual)
    pts  = data.frame(
      X = X[, pred_var],
      # Residual + baseline of other predictor+ effect of our selected predictor
      presid = res + mu_others + b[pred_var] * X[, pred_var]
    )
  )
}

# Define a function to apply same style and plotting
# pred_var: predictor variable, xlab: label for x-axis, panel_label: custom letter to identify the panel
make_panel <- function(pred_var, xlab, panel_label, model = pgls_final, data = df_trait) {
  pe <- partial_effect(pred_var, model, data)
  
  ggplot() +
    # Plot with jitter to reduce overlapping
    geom_jitter(data = pe$pts, aes(x = X, y = presid),
                width = 0.0, height = 0.2, size = 0.8, alpha = 0.5, color = col_dots) +
    # Plot 95% Confidence Interval ribbon
    geom_ribbon(data = pe$line, aes(x = X, ymin = ymin, ymax = ymax),
                fill = col_ribbon, alpha = 0.2) +
    # Plot fitted partial effect line
    geom_line(data = pe$line, aes(x = X, y = predicted),
              color = col_ribbon, linewidth = 1.5) +
    # Add plot annotation (A, B, C, D) in the top left corner
    annotate("text", x = -Inf, y = Inf, label = panel_label,
             hjust = -2, vjust = 1.5, size = 6) +
    labs(x = xlab, y = NULL) +
    partial_theme
}

# Call the previous define functions to create a plot(panel) for every predictor
plot_egg    <- make_panel("Egg_size", "Egg size (mm)", "A")
plot_clutch <- make_panel("Clutch_size", "Clutch size (No. of eggs)", "B")
plot_mat    <- make_panel("Maturity", "Age at Maturity (years)", "D")
plot_svl <- make_panel("Adult_SVL", "Snout-vent length (mm)", "C") 

# Create a unique y-axis label for the plot
y_axis_plot <- ggplot() +
  annotate("text", x = 1, y = 1, label = "Critical thermal maximum (°C)",
           angle = 90, size = 6, fontface = "plain") +
  theme_void()

# Combine the 4 plots into a 2x2 grid
# Apply y-axis limits and number formatting
data_panels <- ((plot_egg | plot_clutch) / (plot_svl | plot_mat)) & 
  coord_cartesian(ylim = c(28, 42)) & 
  scale_y_continuous(labels = function(x) round(x))

# Merge the shared y-axis and panel grid
plot_combined <- y_axis_plot + data_panels + plot_layout(widths = c(1, 25))

print(plot_combined)
ggsave(
  filename = "../figure/figure2.png", 
  plot = plot_combined, 
  width = 8, height = 6.5, dpi = 600, bg = "white"
)
ggsave(
  filename = "../figure/figure2.pdf", 
  plot = plot_combined, 
  width = 8, height = 6.5, dpi = 600, bg = "white"
)


# Figure 3 (Evolutionary variances)
# Order factor levels by increasing evolutionary variance
trait_order <- multi_summary$Trait
multi_df$Trait      <- factor(multi_df$Trait, levels = trait_order)
multi_summary$Trait <- factor(multi_summary$Trait, levels = trait_order)

# Color palette mapped to traits
distinct_colors <- c(
  "Egg size"        = "#E69F00",
  "CTmax"           = "#56B4E9",
  "SVL"             = "#009E73",
  "Clutch size"     = "#CC79A7",
  "Age at maturity" = "#D55E00"
)

# Custom x-axis labels
trait_labels <- expression(
  "Egg size"        = "Egg\nsize",
  "Clutch size"     = "Clutch\nsize",
  "SVL"             = "Snout-vent\nlength",
  "Age at maturity" = "Age at\nmaturity",
  "CTmax"           = CT[max]
)

# Set plot features
pub_theme <- theme_classic(base_family = "sans", base_size = 20) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 20, margin = margin(t = 15)),
    axis.text.y      = element_text(color = "black", size = 20),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(face = "plain", size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position  = "none", 
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank() 
  )

# Display plot
plot_multi <- ggplot() +
  # Half-violin plot showing posterior/jackknife distributions
  geom_half_violin(
    data = multi_df, 
    aes(x = Trait, y = SigmaSq, fill = Trait),
    side = "r", 
    position = position_nudge(x = 0.1), # Nudging set to 0.1
    color = NA, 
    alpha = 0.7, 
    trim = FALSE
  ) +
  # Point range for Mean and 95% CIs
  geom_pointrange(
    data = multi_summary, 
    aes(x = Trait, y = Mean, ymin = Lower, ymax = Upper),
    position = position_nudge(x = -0.1), 
    size = 0.8, 
    linewidth = 1, 
    color = "black"
  ) +
  scale_fill_manual(values = distinct_colors) +
  scale_x_discrete(labels = trait_labels) +
  labs(y = "Evolutionary variance (SD)") +
  pub_theme

print(plot_multi)

# Export plot
ggsave(
  filename = "../figure/figure3.png", 
  plot = plot_multi, 
  width = 8, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)
ggsave(
  filename = "../figure/figure3.pdf", 
  plot = plot_multi, 
  width = 8, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)

# Figure S2-A
# Group the tree tips by ORDER
order_list <- split(main_data$Species, main_data$Order)
tree_grouped <- groupOTU(pruned_tree, order_list)
order_colors <- c(
  "Anura"   = "#E69F00", # Warm Orange
  "Caudata" = "#0072B2", # Dark Blue
  "0"       = "grey75"   # Neutral grey for deep backbone
)

# Plot the basic tree
plot_phylo <- ggtree(
  tree_grouped,
  layout = "circular",
  aes(color = group), 
  linewidth = 0.6
) +
  scale_color_manual(values = order_colors) +
  theme_void() +
  theme(
    legend.position = "none",        
    text = element_text(family = "sans"),
    plot.margin = margin(40, 55, 80, 55, "pt") # Expanded margins to prevent clipping
  )

# Add external family ring and radial labelling
# Group species by family
family_groups <- split(main_data$Species, main_data$Family)

plot_df <- plot_phylo$data
for (fam_name in names(family_groups)) {
  sp_in_fam <- family_groups[[fam_name]]
  
  target_node <- if (length(sp_in_fam) > 1) {
    MRCA(pruned_tree, sp_in_fam)
  } else {
    which(pruned_tree$tip.label == sp_in_fam[1])
  }
  
  if (!is.null(target_node) && target_node %in% plot_df$node) {
    
    # Extract arc-midpoint from family tips 
    fam_tip_angles <- plot_df$angle[plot_df$label %in% sp_in_fam]
    mid_angle      <- mean(fam_tip_angles)
    
    # Set label inclination and orientation
    if (mid_angle > 90 && mid_angle < 270) {
      label_angle <- mid_angle + 180
      label_hjust <- 1        # Right-align toward the ring
    } else {
      label_angle <- mid_angle
      label_hjust <- 0        # Left-align away from the ring
    }
    
    plot_phylo <- plot_phylo +
      geom_cladelabel(
        node = target_node,
        label = fam_name,
        align = TRUE,
        offset = 6,       # Distance from tip ends to the black arc bars
        offset.text = 3,  # Distance between the arc bar and the text label
        barsize = 1,      # Thickness of the outer family arc
        fontsize = 5,   # Font size for family names 
        family = "sans",   
        angle = label_angle,
        hjust = label_hjust,
        horizontal = FALSE,
        color = "black"
      )
  }
}


# Figure S2-B
# Load IUCN elaborated polygons, adjust to your path (Be sure to run 'speciesRange.r' before)
IUCN <- st_read("../data/IUCN_final.gpkg")
# Import world boundaries map
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Set template for raster plotting
r_template <- rast(
  xmin = -180, xmax = 180,
  ymin = -90,  ymax = 90,
  resolution = 1
)

#Converts IUCN polygons into a terra vector object
IUCN_vect <- vect(IUCN)
# Assigns a value of 1 to every species range polygon
IUCN_vect$presence <- 1
# Calulcate species richness: number of overlapping polygons
richness_raster <- rasterize(
  IUCN_vect,
  r_template,
  field = "presence",
  fun = "sum",
  background = NA
)

richness_df <- as.data.frame(richness_raster, xy = TRUE, na.rm = TRUE)
colnames(richness_df)[3] <- "richness"

# Set plot features
plot_range <- ggplot() +
  
  geom_sf(
    data = world,
    fill = "grey80",
    colour = "white",
    linewidth = 0.15
  ) +
  
  geom_tile(
    data = richness_df,
    aes(x = x, y = y, fill = richness),
    alpha = 0.90
  ) +
  
  scale_fill_viridis_c(
    option = "cividis",
    trans = "log1p",
    breaks = c(3, 10, 30),
    name = "Richness of\nthe study species"
  ) +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal(base_size = 20) +
  
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 18, colour = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 18, vjust = 0.5),
    legend.text  = element_text(size = 16),
    
    legend.key.height = unit(5, "mm"),
    legend.key.width  = unit(20, "mm"),
    
    legend.box.margin = margin(t = 10, r = 0, b = 0, l = 0),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    
    plot.margin = margin(t = 10, r = 10, b = -5, l = 10)
  ) +
  
  labs(x = "Longitude", y = "Latitude")

# Display plots
print(plot_phylo)
print(plot_range)
# Combine and export the plot
plot_phylo <- plot_phylo + 
  theme(plot.margin = margin(t = 50, r = 60, b = 80, l = 10, unit = "pt"))

plot_range <- plot_range + 
  theme(plot.margin = margin(t = 10, r = 60, b = 10, l = 10, unit = "pt"))

plot_phyloRange <- (plot_phylo / plot_range) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 28,family = "sans"))
# Export the plot
ggsave("../figure/figureS2.png", plot = plot_phyloRange, width = 10, height = 12, dpi = 600, bg = "white")
ggsave("../figure/figureS2.pdf", plot = plot_phyloRange, width = 10, height = 12, dpi = 600, bg = "white")
