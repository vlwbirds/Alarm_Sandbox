########################################################################################################
########################## Spp Level Dredge ############################################################
########################################################################################################

# Load necessary packages
library(dplyr)
library(readr)
library(here)
library(broom)
library(tidyverse)

# Load dataset
HwkGLM <- read_csv(here("data/HwkGLM.csv"))

# Assigning appropriate data types
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
                  "Alarm_Presence",
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
  mutate(Height_Difference = Flight_apex - Pre_Bird_Height)

# Load dataset into inputdata
inputdata <- df

# Define predictor variables
predictors <- c("Tar_Sex", "Pre_Bird_Height", "Flight_apex", "Proximity", "Veg_under_height", "Bird_Dense_Avg")

# Set minimum sample size threshold
min_sample_size <- 15  # Adjust this threshold as needed

# Get a list of unique species that meet the sample size threshold
species_list <- inputdata %>%
  group_by(TarSpp) %>%
  summarize(sample_size = n(), .groups = "drop") %>%
  filter(sample_size >= min_sample_size) %>%
  pull(TarSpp)

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
      
      # Always include the predictors in the model
      predictor_terms <- paste(combination, collapse = " + ")
      
      # Test the model with predictors only
      formula_str <- paste("Alarm_Presence ~", predictor_terms)
      model <- glm(as.formula(formula_str), data = species_data, family = binomial)
      
      # Store the summary
      model_summaries[[paste(species, paste(combination, collapse = "_"), sep = "_")]] <- tidy(model) %>%
        mutate(species = species,
               predictors = predictor_terms,
               interactions = "None",
               AIC = AIC(model),
               BIC = BIC(model),
               sample_size = nrow(species_data))
      
      # For 2 predictors, include interaction terms
      if (num_predictors == 2) {
        interaction_terms <- combn(combination, 2, function(x) paste(x, collapse = ":"), simplify = TRUE)
        interaction_terms_str <- paste(interaction_terms, collapse = " + ")
        
        # Test the model with predictors + interactions
        formula_str_with_interactions <- paste("Alarm_Presence ~", paste(predictor_terms, interaction_terms_str, sep = " + "))
        model_with_interactions <- glm(as.formula(formula_str_with_interactions), data = species_data, family = binomial)
        
        # Store the summary
        model_summaries[[paste(species, paste(combination, collapse = "_"), "interactions", sep = "_")]] <- tidy(model_with_interactions) %>%
          mutate(species = species,
                 predictors = predictor_terms,
                 interactions = interaction_terms_str,
                 AIC = AIC(model_with_interactions),
                 BIC = BIC(model_with_interactions),
                 sample_size = nrow(species_data))
      }
    }
  }
}

# Combine all summaries into one data frame
all_summaries <- bind_rows(model_summaries)

# Export the combined summaries to a CSV file
write_csv(all_summaries, here("output/GLM_summaries_3Predictors_SppLvl.csv"))

cat("Model summaries saved to 'output/GLM_summaries_3Predictors_SppLvl.csv'\n")

view(read_csv(here("output/GLM_summaries_3Predictors_SppLvl.csv")))

# add columns to df: delta AICC, AICC, Index models for groupby(), AICC weight - number from 0 to 1 to describe influence of model
# take a weighted average of coefficient values
# explore dredge packages
# groupby(weighted average) %>% 
# summarise()
# calculate a new table where each row is a model singularly to calculate the weight of each model
## identify plausible models, do averaging on best fits
## look up examples of model averaging
## explore MuMIn package
# coef_table %>% filter(delta<4) %>% group_by(param) %>% summarise(mod.avg_coef=weighted.mean(estimate, weight))
## scale parameters based on standard deviation