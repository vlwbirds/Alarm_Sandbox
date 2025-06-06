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
view(df)

# Load dataset into inputdata
inputdata <- df

# Define predictor variables
predictors <- c("Tar_Sex", 
                "Pre_Bird_Height",
                "Flight_apex",
                "Flight_proximity",
                "Veg_under_height",
                "Bird_Dense_Avg",
                "Height_Difference"
                )

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
  mutate(Delta_AAIC = AAIC - min(AAIC)) %>%
  ungroup()

str(all_summaries)

# Export the combined summaries to a CSV file
write_csv(all_summaries, here("output/GLM_summaries_AAIC.csv"))

cat("Model summaries saved to 'output/GLM_summaries_AAIC.csv'\n")

view(read_csv(here("output/GLM_summaries_AAIC.csv")))

# Reshape the model summaries into a table where each row represents a single model
coef_table <- all_summaries %>%
  select(species, predictors, term, p.value, AAIC, Delta_AAIC) %>%
  pivot_wider(names_from = term, values_from = c( p.value), names_prefix = "p_") %>%
  mutate(model_index = row_number(),
         Akaike_weight = exp(-0.5 * Delta_AAIC) / sum(exp(-0.5 * Delta_AAIC))) %>%
  arrange(AAIC)  # Sorting by best model (lowest AAIC)
view(coef_table)

# save to csv
write_csv(coef_table, here("output/GLM_Summary_Table.csv"))

# Select the top 10 best-fitted models
top_10_models <- coef_table %>%
  slice_min(order_by = AAIC, n = 20)

# Save to CSV
write_csv(top_10_models, here("output/Top_10_GLM_Models.csv"))

view(read_csv(here("output/Top_10_GLM_Models.csv")))
# Print the table
print(top_10_models)

######
# questions: how to identify a good model, which predictors are important?
# how do we summarise effects across multiple models?

# Add model index to all_summaries based on predictors
all_summaries <- all_summaries %>%
  mutate(model_index = match(predictors, coef_table$predictors))

# Merge the long-format all_summaries with coef_table to include AAIC weight
all_summaries_weighted <- all_summaries %>%
  left_join(select(coef_table, model_index, Akaike_weight), by = "model_index")

# View the updated table
view(all_summaries_weighted)

# Save to CSV
write_csv(all_summaries_weighted, here("output/GLM_summaries_with_weights.csv"))

# Take the table number 2 %>%
#   Filter to delta aicc below 4%>%
#   Group by term %>%
#   Summarise to take weighted means of estimates and SE, weighting Parameter estimates by model weight
# 
# That should return one estimate and se per term.
#
# Might need to backfill 0s to have the weighting work correctly but that's for after getting this to work 
# I.e. a model with only predictor A also technically can be thought of having predictor B and C but their estimates and SE are 0 and only A has a nonzero estimate in this case.


# Backfilling missing predictors with 0 estimates and SE
all_summaries_filled <- all_summaries_weighted %>%
  complete(model_index, term, fill = list(estimate = 0, std.error = 0))

weighted.var.se <- function(x, w, na.rm=FALSE)
  #  Computes the variance of a weighted mean following Cochran 1977 definition
{
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  n = length(w)
  xWbar = weighted.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}

# Calculate weighted means
weighted_SE <- all_summaries_filled %>% 
  filter(Delta_AAIC < 4) %>% 
  group_by(term) %>% 
  summarise(
    weighted_estimate = weighted.mean(estimate, Akaike_weight, na.rm = T),
    weighted_SE = weighted.mean(std.error, Akaike_weight, na.rm = T)
  ) 

# View results
view(weighted_SE)

# Save to CSV
write_csv(weighted_SE, here("output/weighted_parameter_estimates.csv"))


# So you are missing one step in the call to summarise()
# weighted_estimate = sum(estimate * Akaike_weight, na.rm = TRUE)
# ...needs to be standardized by the sum of weights, such that ....
# weighted_estimate = sum(estimate * Akaike_weight, na.rm = TRUE)/sum(Akaike_weight, na.rm=T)
# 
# Alternatively, you should be able to use >summarise(we=weighted.mean(x=estimate, w=Akaike_weight)) to get the same result. 
# Does it? If not, we might be touching on exactly that issue you identified earlier 
# about the weights not summing to one *after* the filtering to models with delta aicc below 4. 
# Could you explore that and explain that to me when we meet and have code to show me how it works?
