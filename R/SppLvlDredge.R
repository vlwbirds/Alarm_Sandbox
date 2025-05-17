########################################################################################################
########################## Spp Level Dredge ############################################################
########################################################################################################

# Load necessary packages
library(AICcmodavg)
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


# Store summaries and aictabs
model_summaries <- list()
aictab_outputs <- list()

for (species in species_list) {
  species_data <- inputdata %>%
    filter(TarSpp == species) %>%
    filter(rowSums(is.na(select(., all_of(predictors)))) == 0)
  
  if (nrow(species_data) < min_sample_size) {
    cat("Skipping species:", species, "- Sample size after filtering:", nrow(species_data), "\n")
    next
  }
  
  cat("Processing species:", species, "- Sample size:", nrow(species_data), "\n")
  
  # Lists to store models and names
  species_models <- list()
  model_names <- c()
  
  for (num_predictors in 1:3) {
    predictor_combinations <- combn(predictors, num_predictors, simplify = FALSE)
    
    for (combination in predictor_combinations) {
      if (!all(sapply(combination, function(var) length(unique(species_data[[var]])) > 1))) next
      
      predictor_terms <- paste(combination, collapse = " + ")
      formula_str <- paste("Alarm_Presence ~", predictor_terms)
      model <- glm(as.formula(formula_str), data = species_data, family = binomial)
      
      model_id <- paste(species, paste(combination, collapse = "_"), sep = "_")
      species_models[[model_id]] <- model
      model_names <- c(model_names, model_id)
      
      model_summaries[[model_id]] <- tidy(model) %>%
        mutate(species = species,
               predictors = predictor_terms,
               interactions = "None",
               AIC = AIC(model),
               BIC = BIC(model),
               sample_size = nrow(species_data))
      
      # Include interactions if 2 predictors
      if (num_predictors == 2) {
        interaction_terms <- combn(combination, 2, function(x) paste(x, collapse = ":"), simplify = TRUE)
        interaction_terms_str <- paste(interaction_terms, collapse = " + ")
        formula_str_with_interactions <- paste("Alarm_Presence ~", paste(predictor_terms, interaction_terms_str, sep = " + "))
        model_interact <- glm(as.formula(formula_str_with_interactions), data = species_data, family = binomial)
        
        interact_id <- paste(species, paste(combination, collapse = "_"), "interactions", sep = "_")
        species_models[[interact_id]] <- model_interact
        model_names <- c(model_names, interact_id)
        
        model_summaries[[interact_id]] <- tidy(model_interact) %>%
          mutate(species = species,
                 predictors = predictor_terms,
                 interactions = interaction_terms_str,
                 AIC = AIC(model_interact),
                 BIC = BIC(model_interact),
                 sample_size = nrow(species_data))
      }
    }
  }
  
  # Generate species-safe variable name
  species_varname <- gsub("[^A-Za-z0-9]", "_", species)
  assign(paste0("models_", species_varname),
         aictab(cand.set = species_models, modnames = model_names))
  aictab_outputs[[species]] <- aictab(cand.set = species_models, modnames = model_names)
}

# Combine all summaries into one data frame
all_summaries <- bind_rows(model_summaries)

# Export the combined summaries to a CSV file
write_csv(all_summaries, here("output/GLM_summaries_3Predictors_SppLvl_DeltaAIC.csv"))

cat("Model summaries saved to 'output/GLM_summaries_3Predictors_SppLvl.csv'\n")

view(read_csv(here("output/GLM_summaries_3Predictors_SppLvl.csv")))
GLM_sum <- read_csv(here("output/GLM_summaries_3Predictors_SppLvl.csv"))

# calculate AIC and model averages
aictab(cand.set = GLM_sum, modnames = mod.names)

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