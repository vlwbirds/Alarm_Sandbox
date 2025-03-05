################################################################################################################################################
################################################################################################################################################
################################################ Species Level GLM with Weighted AAIC ##########################################################
################################################################################################################################################
################################################################################################################################################

# Load necessary packages
library(dplyr)
library(readr)
library(here)
library(broom)
library(tidyverse)

# Load dataset
HwkGLM <- read_csv(here("data/HwkGLM.csv"))

df <- HwkGLM %>%
  mutate(across(c("Habitat_Primary",
                  "Habitat_Under_Dom1", 
                  "Habitat_Under_Dom2", 
                  "SocialGroup", 
                  "TarSpp",
                  "Alarm_Spp",
                  "Alarm_Sex",
                  "Tar_Sex",
                  "Audio_React",
                  "Pre_Vocal",
                  "Mid_Vocal",
                  "Pre_Behav"
  ), as.factor)) %>%
  mutate(across(c("Obs1_bird_dist",
                  "Shoot_bird_dist",
                  "Alarm_Qty",
                  "Alarm_Conspecific",
                  "abundance",
                  "richness",
                  "Proximity",
                  "Flight_dist",
                  "Flight_apex",
                  "Flight_proximity",
                  "Pre_Bird_Height",
                  "Bird_Dense_Avg",
                  "Post_Dist",
                  "Post_Bird_Height",
                  "Post_Dense_Avg",
                  "Veg_under_height",
                  "Veg_mid_low",
                  "Veg_mid_high",
                  "Veg_can_low",
                  "Veg_can_high",
                  "Rad_10_Avg",
                  "Veg_Under_Avg",
                  "Forage_Tree_height"
  ), as.numeric)) %>%
  mutate(Height_Difference = Flight_apex - Pre_Bird_Height) %>%
  mutate(across(where(is.numeric) & !c(Alarm_Presence), scale))  # Exclude Alarm_Presence

# Load dataset into inputdata
inputdata <- df

# Define predictor variables
predictors <- c("Tar_Sex", 
                "Pre_Bird_Height",
                "Flight_apex",
                "Flight_proximity",
                "Veg_under_height",
                "Bird_Dense_Avg",
                "Height_Difference")

# Define species of interest (only haru)
species <- "haru"

# Filter data for haru
species_data <- inputdata %>%
  filter(TarSpp == species) %>%
  filter(rowSums(is.na(select(., all_of(predictors)))) == 0) %>%
  mutate(Alarm_Presence = as.numeric(Alarm_Presence)) %>%  # Ensure numeric
  filter(Alarm_Presence %in% c(0, 1))  # Ensure binary response

# Initialize list for model summaries
model_summaries <- list()

# Debugging: Output sample size for the species
cat("Processing species:", species, "- Sample size:", nrow(species_data), "\n")

# Iterate over predictor combinations
for (num_predictors in 1:2) {
  predictor_combinations <- combn(predictors, num_predictors, simplify = FALSE)
  
  for (combination in predictor_combinations) {
    valid_combination <- all(sapply(combination, function(var) length(unique(species_data[[var]])) > 1))
    if (!valid_combination) next
    
    formula_str <- paste("Alarm_Presence ~", paste(combination, collapse = " + "))
    model <- glm(as.formula(formula_str), data = species_data, family = binomial)
    
    k <- length(combination) + 1  # Number of parameters (including intercept)
    n <- nrow(species_data)
    aic <- AIC(model)
    aaic <- aic + ((2 * k * (k + 1)) / (n - k - 1))
    
    model_anova <- anova(model, test = "Chisq")
    f_stat <- ifelse(ncol(model_anova) >= 4, model_anova$Deviance[2] / model_anova$`Resid. Df`[2], NA)
    global_p <- model_anova$`Pr(>Chi)`[2]
    
    model_summaries[[paste(species, paste(combination, collapse = "_"), sep = "_")]] <- tidy(model) %>%
      mutate(species = species,
             predictors = paste(combination, collapse = " + "),
             AAIC = aaic,
             sample_size = nrow(species_data),
             F_statistic = f_stat,
             Global_p_value = global_p,
             R_squared = 1 - (model$deviance / model$null.deviance))
  }
}

# Combine all summaries into a single data frame
all_summaries <- bind_rows(model_summaries)

# Compute Delta AAIC and Akaike Weights (only for haru)
all_summaries <- all_summaries %>%
  mutate(Delta_AAIC = AAIC - min(AAIC))

# Export results
write_csv(all_summaries, here("output/GLM_summaries_AAIC.csv"))
cat("Model summaries saved to 'output/GLM_summaries_AAIC.csv'\n")

# Create summary table with weights
coef_table <- all_summaries %>%
  select(species, predictors, term, p.value, AAIC, Delta_AAIC) %>%
  pivot_wider(names_from = term, values_from = c(p.value), names_prefix = "p_") %>%
  mutate(model_index = row_number(),
         Akaike_weight = exp(-0.5 * Delta_AAIC) / sum(exp(-0.5 * Delta_AAIC))) %>%
  arrange(AAIC)

# Save summary table
write_csv(coef_table, here("output/GLM_Summary_Table.csv"))

# Select top 10 best-fitted models
top_10_models <- coef_table %>%
  slice_min(order_by = AAIC, n = 10)

write_csv(top_10_models, here("output/Top_10_GLM_Models.csv"))

# Merge weighted Akaike weights with the summaries
all_summaries <- all_summaries %>%
  mutate(model_index = match(predictors, coef_table$predictors)) %>%
  left_join(select(coef_table, model_index, Akaike_weight), by = "model_index")

# Fill missing predictors with 0 estimates
all_summaries_filled <- all_summaries %>%
  complete(model_index, term, fill = list(estimate = 0, std.error = 0))

# Compute weighted means
weighted_SE <- all_summaries_filled %>% 
  filter(Delta_AAIC < 4) %>% 
  group_by(term) %>% 
  summarise(
    weighted_estimate = weighted.mean(estimate, Akaike_weight, na.rm = TRUE),
    weighted_SE = weighted.mean(std.error, Akaike_weight, na.rm = TRUE)
  )

weighted_SE$z.score <- weighted_SE$weighted_estimate / weighted_SE$weighted_SE
weighted_SE$p.value <- pnorm(abs(weighted_SE$z.score), lower.tail = F) # calculate p value given a z.score 
weighted_SE$t.dist <- pt(abs(weighted_SE$z.score), nrow(all_summaries_filled %>% filter(Delta_AAIC < 4))) # calculate p value from a t score, Student t Distribution

# Save final weighted results
write_csv(weighted_SE, here("output/GLM_Weighted_Estimates_HARU.csv"))

cat("Weighted estimates saved to 'output/GLM_Weighted_Estimates_HARU.csv'\n")
