########################################################################################################
########################## Spp Level Dredge ############################################################
########################################################################################################

# Load necessary packages
library(dplyr)
library(readr)
library(here)
library(broom)  # For tidy model summaries
library(purrr)  # For working with lists

# Loading Reduced GLM Data
HwkGLM <- read_csv(here("data/HwkGLM.csv"))

# Assigning Appropriat Data Types
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
  ), as.factor)) %>% # Convert multiple columns to factors
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
  ), as.numeric)) # Convert multiple columns to numeric

# Create the new column by calculating the difference
df <- df %>%
  mutate(Height_Difference = Flight_apex - Pre_Bird_Height)

# Load dataset
inputdata <- df

# Define predictor variables and interactions as character vectors
predictors <- c("Pre_Bird_Height",
                "Bird_Dense_Avg",
                "Flight_proximity", 
                "Height_Difference",
                "Veg_under_height",
                "Veg_mid_low", 
)  # Specify predictor column names
interactions <- c("Flight_proximity:Pre_Bird_Height", "Height_Difference:Flight_proximity")  # Define interaction terms if needed

# Set minimum sample size threshold
min_sample_size <- 10  # Adjust this threshold as needed

# Get a list of unique species that meet the sample size threshold
species_list <- inputdata %>%
  group_by(TarSpp) %>%
  filter(n() >= min_sample_size) %>%
  pull(TarSpp) %>%
  unique()

# Initialize a list to store model summaries for each species and combination
model_summaries <- list()

# Iterate over each species and fit the GLM individually with all combinations of predictors
for (species in species_list) {
  # Filter dataframe for the current species
  df <- inputdata %>% filter(TarSpp == species)
  
  # Remove rows with NA in any of the specified predictor columns to retain maximum values
  df <- df %>% filter(rowSums(is.na(select(., all_of(predictors)))) == 0)
  
  # Check if there are enough rows after filtering for NAs
  if (nrow(df) < min_sample_size) {
    message(paste("Skipping", species, "- insufficient data after NA filtering."))
    next
  }
  
  # Generate all possible non-empty combinations of predictors
  for (k in 1:length(predictors)) {
    predictor_combinations <- combn(predictors, k, simplify = FALSE)
    
    # Iterate over each combination of predictors
    for (combination in predictor_combinations) {
      # Check if each predictor in the combination has more than one level
      valid_combination <- all(sapply(combination, function(var) length(unique(df[[var]])) > 1))
      
      if (!valid_combination) {
        message(paste("Skipping combination", paste(combination, collapse = ", "), "for species", species, "- insufficient levels."))
        next
      }
      
      # Create the formula string for the current predictor combination
      predictor_terms <- paste(combination, collapse = " + ")
      interaction_terms <- if (length(interactions) > 0) {
        paste(interactions, collapse = " + ")
      } else {
        NULL
      }
      
      # Combine predictors and interactions into the final formula
      formula_str <- paste("Alarm_Presence ~", paste(c(predictor_terms, interaction_terms), collapse = " + "))
      formula_obj <- as.formula(formula_str)
      
      # Fit the binomial GLM if there are predictor terms or interactions
      model <- glm(formula_obj, data = df, family = binomial)
      
      # Get a tidy summary of the model using broom::tidy
      model_summary <- tidy(model) %>%
        mutate(species = species,
               predictors = paste(combination, collapse = ", "),
               AIC = AIC(model),
               BIC = BIC(model))  # Add species, predictors, AIC, and BIC to each row
      
      # Append the model summary to the list
      model_summaries[[paste(species, paste(combination, collapse = "_"), sep = "_")]] <- model_summary
    }
  }
}

# Combine all summaries into one data frame
all_summaries <- bind_rows(model_summaries)

# Export the combined summaries to a CSV file
write_csv(all_summaries, here("output/GLM_summaries_by_species_combinations.csv"))

cat("Model summaries saved to 'output/GLM_summaries_by_species_combinations.csv'\n")
