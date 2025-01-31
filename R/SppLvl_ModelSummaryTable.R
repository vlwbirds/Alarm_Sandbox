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
                  "Alarm_Qty",  # Exclude Alarm_Presence from scaling
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
predictors <- c("Tar_Sex", "Pre_Bird_Height", "Flight_apex", "Proximity", "Veg_under_height", "Bird_Dense_Avg", "Height_Difference")

# Set minimum sample size threshold
min_sample_size <- 30  

# Filter data for the current species
species_data <- inputdata %>%
  filter(TarSpp == species) %>%
  filter(rowSums(is.na(select(., all_of(predictors)))) == 0) %>%
  mutate(Alarm_Presence = as.numeric(Alarm_Presence)) %>%  # Ensure numeric
  filter(Alarm_Presence %in% c(0, 1))  # Ensure binary response

# Recheck sample size after filtering
if (nrow(species_data) > min_sample_size) {
  cat("Skipping species:", species, "- Sample size after filtering:", nrow(species_data), "\n")
  next
}


# Initialize a list to store model summaries
model_summaries <- list()

# Iterate over each species in the filtered list
for (species in species_list) {
  # Filter data for the current species
  species_data <- inputdata %>%
    filter(TarSpp == species) %>%
    filter(rowSums(is.na(select(., all_of(predictors)))) == 0)
  
  # Recheck sample size after filtering
  if (nrow(species_data) < min_sample_size) {
    cat("Skipping species:", species, "- Sample size after filtering:", nrow(species_data), "\n")
    next
  }
  
  # Debugging: Output sample size for the species
  cat("Processing species:", species, "- Sample size:", nrow(species_data), "\n")
  
  # Iterate over the number of predictors in the model
  for (num_predictors in 1:3) {
    # Generate all combinations of predictors for the current number
    predictor_combinations <- combn(predictors, num_predictors, simplify = FALSE)
    
    for (combination in predictor_combinations) {
      # Check if all predictors in the combination have sufficient levels
      valid_combination <- all(sapply(combination, function(var) length(unique(species_data[[var]])) > 1))
      if (!valid_combination) next
      
      # Construct the model formula
      predictor_terms <- paste(combination, collapse = " + ")
      formula_str <- paste("Alarm_Presence ~", predictor_terms)
      model <- glm(as.formula(formula_str), data = species_data, family = binomial)
      
      # Compute AAIC
      k <- length(combination) + 1  # Number of parameters (including intercept)
      n <- nrow(species_data)
      aic <- AIC(model)
      aaic <- aic + ((2 * k * (k + 1)) / (n - k - 1))
      
      # Compute F-statistic and p-value for global model
      model_anova <- anova(model, test = "Chisq")
      f_stat <- ifelse(ncol(model_anova) >= 4, model_anova$Deviance[2] / model_anova$`Resid. Df`[2], NA)
      global_p <- model_anova$`Pr(>Chi)`[2]
      
      # Store the summary
      model_summaries[[paste(species, paste(combination, collapse = "_"), sep = "_")]] <- tidy(model) %>%
        mutate(species = species,
               predictors = predictor_terms,
               AAIC = aaic,
               sample_size = nrow(species_data),
               F_statistic = f_stat,
               Global_p_value = global_p,
               R_squared = 1 - (model$deviance / model$null.deviance))
    }
  }
}

# Combine all summaries into one data frame
all_summaries <- bind_rows(model_summaries)

# Compute Delta AAIC and Akaike Weights
all_summaries <- all_summaries %>%
  group_by(species) %>%
  mutate(Delta_AAIC = AAIC - min(AAIC),
         Akaike_weight = exp(-0.5 * Delta_AAIC) / sum(exp(-0.5 * Delta_AAIC))) %>%
  ungroup()

str(all_summaries)

# Export the combined summaries to a CSV file
write_csv(all_summaries, here("output/GLM_summaries_AAIC.csv"))

cat("Model summaries saved to 'output/GLM_summaries_AAIC.csv'\n")

view(read_csv(here("output/GLM_summaries_AAIC.csv")))

# Reshape the model summaries into a table where each row represents a single model
coef_table <- all_summaries %>%
  select(species, predictors, term, estimate, p.value, AAIC, Delta_AAIC, Akaike_weight, std.error) %>%
  pivot_wider(names_from = term, values_from = c(estimate, p.value), names_prefix = "p_") %>%
  mutate(model_index = row_number()) %>%
  arrange(AAIC)  # Sorting by best model (lowest AAIC)
view(coef_table)
# Select the top 10 best-fitted models
top_10_models <- coef_table %>%
  slice_min(order_by = AAIC, n = 20)

# Save to CSV
write_csv(top_10_models, here("output/Top_10_GLM_Models.csv"))

view(read_csv(here("output/Top_10_GLM_Models.csv")))
# Print the table
print(top_10_models)

###### TO DO
# long output, one predictor per row indexed by specific model
# questions: how to identify a good model, which predictors are important?
# how do we summarise effects across multiple models?
# model weight and index needed to link summary table to long table, left join by model index
# recalculate aaic weight from model summary table "coef_table", right now it's from the long table
# 

inner_join(all_summaries, coef_table %>% select(c(Akaike_weight, Delta_AAIC, predictors, model_index)), predictors)

