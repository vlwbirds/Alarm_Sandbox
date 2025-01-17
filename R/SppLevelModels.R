#2024-12-10 
#vince is about to be awol for a month
#i'm gonna try to get the mcmcglmm or the brms running
#goal is the be able to start the model simplification
#2024-01-08 restart this

library(MCMCglmm)
library(coda)
library(ape)
library(castor)
library(phytools)
library(here)
library(tidyverse)
library(magrittr)
library(performance)
library(GGally)

#load tree
tree<-read.nexus(here("data/trees/AllBirdsHackett1_summary.tre"))
# Loading Reduced GLM Data
df <- read_csv(here("data/HwkAll.csv"))
# Replace '#N/A' with NA
df[df == "#N/A"] <- NA
#remove NA from BirdTree 
df <- df[!is.na(df$BirdTree), ]
#some data organization is missing
#to do: grou pby timestamp , slice 1. this filters out data to only include target individuals
target_individual_responses<-df %>% group_by(Timestamp) %>% slice(1)
target_individual_responses$BirdTree<-gsub(" ", "_", target_individual_responses$BirdTree)

cleaned_data<-target_individual_responses

pruned_tree <- keep.tip(tree, cleaned_data$BirdTree)
all(pruned_tree$tip.label %in% cleaned_data$BirdTree) #T
all(cleaned_data$BirdTree %in% pruned_tree$tip.label) #T



###########################
#STEP1 Models: ssp-levels 
###########################



#Audio_React [alarm/call/quiet] : alarm=[1], call|quiet|w==[0]
#Y: alarm [1/0]
cleaned_data$Audio_React %>% table #there's a W - ask Vince about this
cleaned_data$alarm#is NAs, manually reconstruct as above
cleaned_data$alarm<-ifelse(cleaned_data$Audio_React=="a", 1, 0)
cleaned_data$alarm %>% table #ok

#############vars for global models: species-specific models
#Spp, Sex, 
# (also use in interspecific model), 
#Flight_apex (maximum height of model hawk)
#Pre_Bird_Height
#[ignore this height difference for now but it might be important]
#height_difference_between_focalbird_and model = Pre_Bird_Height - Flight_apex
#Proximity: of hawk-focal bird
#Veg_under_height; avg height of understory in 5m radius
#Pre_Bird_Dense: relative density of vegetation of focal bird at moment of 

#only for spp-specific models
#THAR & HARU
##########################ecological model with only THAR######################
df_sp_thar<-cleaned_data %>% filter(BirdTree=="Thamnomanes_ardesiacus")
#choose vars
names(df_sp_thar)
df_thar<-df_sp_thar %>% ungroup() %>% select(c("Timestamp", "alarm", "Audio_React",
                                               "TarSpp", "BirdTree", "Tar_Sex", "Pre_Bird_Height", "Flight_apex", "Proximity", "Veg_under_height", "Pre_Bird_Dense")) 
str(df_thar)
#tibbles are atrocious to work with
df_thar<-as.data.frame(df_thar)
#5 variables that should not be characters are - some NAs, verify in not un full dataset above
for(i in 7:10){df_thar[,i]<-as.numeric(df_thar[,i])}
#double check if there is row-level info in the full data set that could help salvage the NAs in the reduced
df_thar %>% filter(!complete.cases(.)) #NAs in flight apex and Vegunder
bad_thar_experiments<-df_thar %>% filter(!complete.cases(.)) %$% Timestamp
View(df %>% filter(Timestamp %in% bad_thar_experiments))
#to me it looks like the NAs are legit.... Vince should confirm this
#work with the NAs for now
df_thar<-df_thar%>%mutate(height_difference=Pre_Bird_Height-Flight_apex)
#make density tractable - sum?
df_thar$Pre_Bird_Dense_numeric<-sapply(strsplit(df_thar$Pre_Bird_Dense, ""), FUN=function(x){mean(as.numeric(x))})
#we now have a viable analytic dataframe
#first, we check if all continuous variables require transformation
#to make a pairplot, drop charactervalues with too many or uninformative levels 
ggpairs(df_thar %>% select(-which(names(df_thar) %in% c("Timestamp", "TarSpp", "BirdTree", "Pre_Bird_Dense"))))
#flight apex, Proximity and height difference seem like they have long tails and could benefit from log-transformation
#why? because variables must be Z-transformed before regression, this means that their SDs need to mean something
#a one unit increase in SD is hard to interpret when the distribution is wildly non-normal
df_thar_df_clean<-df_thar %>% transmute(alarm=alarm, sex=Tar_Sex, prebirdheight=Pre_Bird_Height, flight_apex_ln=log(Flight_apex), proximity_ln=log(Proximity),
                                        veg_under_height=Veg_under_height, height_difference=height_difference, pre_bird_density=Pre_Bird_Dense_numeric)
#now we scale and center all numeric predictor variables
for(i in 3:8){
  df_thar_df_clean[,i]<-as.numeric(scale(df_thar_df_clean[,i]))
  names(df_thar_df_clean)[i]<-paste0(names(df_thar_df_clean)[i], "_std")
}
#give a last look at the true, scaled, predictor-level correlation - any visible big ones here are going to make us look out for VIFs
ggpairs(df_thar_df_clean)
#no surprise here, proximity and height difference are highly correlated, but that's right as we defined one as a function of the other
#so for the regression analyses we should just pick one
#ill arbitrarily ignore the height difference for now as that is perhaps harder to interpret
full_thar<-glm(alarm~sex+prebirdheight_std+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial")
#the model fits
check_collinearity(full_thar)
#we are ok, only now do we know that this is a model we can actually work with
#now let's assume we care about potential pairwise interactions (forgot to add sex but point is moot)
full_thar<-glm(alarm~sex+
                 prebirdheight_std+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std+
                 
                 prebirdheight_std:flight_apex_ln_std+prebirdheight_std:proximity_ln_std+prebirdheight_std:veg_under_height_std+prebirdheight_std:pre_bird_density_std+ 
                 flight_apex_ln_std:proximity_ln_std+flight_apex_ln_std:veg_under_height_std+flight_apex_ln_std:pre_bird_density_std+
                 proximity_ln_std:veg_under_height_std+proximity_ln_std:pre_bird_density_std+
                 veg_under_height_std:pre_bird_density_std,
               data = df_thar_df_clean, family = "binomial")
#the model wont converge - why would that be? well in this model there are 34 datapoints and 16 parameters
#you want at least 10N per model parameter
#we are underpowered to consider this kind of model
#hence, we need a simpler model
#ideally, we would have a priori, strict set of interactions that we want to consider
#barring that, in order to restrict our model space in an way that is less arbitrary, we can first search for any significant, additive effects (via backwards model selection)
#and then probe for any pairwise interactions involving the remaining, significant variables.
#so we go back to the additive model 
full_thar<-glm(alarm~sex+prebirdheight_std+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial")
summary(full_thar)
#remove the factor with highest p-value and refit
#that's this
#                   Estimate Std. Error z value Pr(>|z|)
#prebirdheight_std     0.04041    0.39994   0.101    0.920

#that info goes in the table (excel and word) - this info is relevant for the paper
summary(glm(alarm~sex+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
#                   Estimate Std. Error z value Pr(>|z|)
#sexm                  -0.1851     0.7967  -0.232    0.816
summary(glm(alarm~flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
#                   Estimate Std. Error z value Pr(>|z|)
#flight_apex_ln_std     0.1733     0.5541   0.313    0.754
summary(glm(alarm~proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
#                   Estimate Std. Error z value Pr(>|z|)
#proximity_ln_std      -0.1818     0.4095  -0.444    0.657
summary(glm(alarm~veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
#                   Estimate Std. Error z value Pr(>|z|)
#veg_under_height_std  -0.4762     0.4024  -1.183    0.237
summary(glm(alarm~pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
#                   Estimate Std. Error z value Pr(>|z|)
#pre_bird_density_std   0.4129     0.3655   1.130    0.259

#ok - i didn't expect this, but this is interesting in the sense that it suggests that for this species, these ecological values don't matter
#interesting because we could now predict that instead, phylogeny matters, and we have a lot of data to support this idea

############point of discussion: do we want to do anything to explore any pairwise interactions? controlled dredging would be an option, possibly quite justifiable given the 
#low power, but we'd have to have clear guidelines about how much model space to explore
#full dredging is not an option because the model space exhausts our sampling effort quickly
#one idea could be to explore all possible, plausible models with up to 3 model predictors: 3 because 3 model parameters (4 with intercept) with N34 is pushing it
#could be 3 additive effects, could be 2 additive effects and their interaction

#discussion: 2025-01-08
#carfeully construct model space: only those models with max 3predictors (all models with 1 predictor, all models with 2 predictors, all models with 3 predictors, all models with 2 predictors and their interaction)
#maake model averaging table (aicc, delta aicc, waicc, one column per model coefficient and one more for its se)
#note: as AIC is sensitive to sample size, all models that you make for here have to be fitted to the same dataset, so those 2 rows with NAs have to be excluded here
#note that this doesnt apply above as we were doing model selection not by AIC, but by significance of parameters, which is not sensitive in the same way as AIC





##########################ecological model with only HARU######################


###########################
#STEP2: multispecies models with phylo:
#predictors; include the important stuff from step1 AND

#BodyMass
#SOCIAL GRADIENT
#foraging strategy: Remsen classification
#ADD: species-level traits
#(diet maybe)
#eltontraits: foraging stratum. #weighted average like will sweet's project?
#residual eye size?
###########################



















# # Check for matching species between cleaned data and tree
# species_to_keep <- unique(cleaned_data$BirdTree)
# matching_species <- intersect(tree$tip.label, species_to_keep)
# 
# if (length(matching_species) == 0) {
#   stop("No matching species found between cleaned data and tree. Check species names for consistency.")
# }
# 
# # Prune the tree to include only matching species
# pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, matching_species))
# 
# # Add branch lengths to the pruned tree if they are missing
# if (is.null(pruned_tree$edge.length)) {
#   pruned_tree$edge.length <- rep(1, nrow(pruned_tree$edge))  # Assign a default branch length of 1
# }
# 
# # Root the tree if necessary
# if (!is.rooted(pruned_tree)) {
#   # Choose a suitable outgroup if you have one, or root it by midpoint if none
#   rooted_tree <- midpoint.root(pruned_tree)  # Midpoint rooting as a general approach
# } else {
#   rooted_tree <- pruned_tree
# }
# 
# # Set node labels (optional, if you need labeled nodes for further analyses)
# rooted_tree$node.label <- paste0("node", 1:length(rooted_tree$node.label))
# 
# # Assign arbitrary branch lengths if they are missing
# if (is.null(rooted_tree$edge.length)) {
#   rooted_tree$edge.length <- rep(1, nrow(rooted_tree$edge))
# }
# 
# # Make the tree ultrametric using `chronos` with a simple constant rate model
# ultrametric_tree <- chronos(rooted_tree, lambda = 1)  # lambda = 1 applies a constant rate

# Verify the tree is now ultrametric
is.ultrametric(ultrametric_tree)  # Should return TRUE

# Ensure BirdTree contains only species present in the tree
tree_species <- ultrametric_tree$tip.label
cleaned_data <- cleaned_data %>% filter(BirdTree %in% tree_species)

# Convert BirdTree to a factor with levels matching tree species, then drop unused levels
cleaned_data$BirdTree <- factor(cleaned_data$BirdTree, levels = tree_species)
cleaned_data$BirdTree <- droplevels(cleaned_data$BirdTree)

# Confirm that levels of BirdTree in cleaned_data match the row names of inv.phylo$Ainv
inv.phylo <- inverseA(ultrametric_tree, nodes = "ALL", scale = TRUE)

# Get species names in ginverse and check alignment with cleaned_data's BirdTree
ginverse_species <- rownames(inv.phylo$Ainv)
data_species <- levels(cleaned_data$BirdTree)

# Check for any mismatch
missing_in_ginverse <- setdiff(data_species, ginverse_species)
missing_in_data <- setdiff(ginverse_species, data_species)

# Display messages if there are mismatches
if (length(missing_in_ginverse) > 0) {
  warning("Species in data but not in ginverse: ", paste(missing_in_ginverse, collapse = ", "))
}
if (length(missing_in_data) > 0) {
  warning("Species in ginverse but not in data: ", paste(missing_in_data, collapse = ", "))
}

# Get the species labels in the Ainv matrix
Ainv_species <- rownames(inv.phylo$Ainv)

# Identify species in tree_species that are not in Ainv_species
missing_species <- setdiff(tree_species, Ainv_species)
if (length(missing_species) > 0) {
  cat("Warning: The following species are in tree_species but not in Ainv:\n")
  print(missing_species)
}

# Filter Ainv to include only species in both tree_species and Ainv_species
common_species <- intersect(tree_species, Ainv_species)

# Check the dimensions of Ainv and the row/column names
cat("Dimensions of Ainv:", dim(inv.phylo$Ainv), "\n")
cat("Species in Ainv:", rownames(inv.phylo$Ainv), "\n")

# Check the species in cleaned_data and the pruned tree
cat("Species in cleaned_data:", unique(cleaned_data$BirdTree), "\n")
cat("Species in pruned tree:", pruned_tree$tip.label, "\n")

species_only_Ainv <- inv.phylo$Ainv[common_species, common_species, drop = FALSE]

# Now continue with model fitting using the filtered species_only_Ainv
prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
              R = list(V = 1, nu = 0.002))

mcmcglmm_mod <- list()
for (chain in 1:3) {
  mcmcglmm_mod[[chain]] <- MCMCglmm(
    Alarm_Presence ~ Pre_Bird_Height + SocialGroup,
    random = ~BirdTree,
    data = cleaned_data,
    family = "categorical",
    ginverse = list(BirdTree = species_only_Ainv),
    prior = prior,
    nitt = 25000,
    burnin = 5000,
    thin = 1000,
    verbose = TRUE,
    pr = TRUE
  )
}

# Save the model
saveRDS(mcmcglmm_mod, here("output/models/mcmcglmm_mod.rds"))




# Diagnostic checks
par(mfrow = c(3, 2), mar = c(2, 2, 1, 2))
plot(do.call(mcmc.list, lapply(mcmcglmm_mod, function(m) m$Sol)), ask = FALSE)
gelman_diag <- gelman.diag(do.call(mcmc.list, lapply(mcmcglmm_mod, function(m) m$Sol)))
print(gelman_diag$psrf[, 1] %>% range)  # Check Rhat values
summary(mcmcglmm_mod[[1]])  # Effective sample sizes

# Combine chains for inference if convergence is good
all_chains <- mcmcglmm_mod[[1]]
all_chains$Sol <- as.mcmc(rbind(mcmcglmm_mod[[1]]$Sol, mcmcglmm_mod[[2]]$Sol, mcmcglmm_mod[[3]]$Sol))
all_chains$VCV <- as.mcmc(rbind(mcmcglmm_mod[[1]]$VCV, mcmcglmm_mod[[2]]$VCV, mcmcglmm_mod[[3]]$VCV))
summary(all_chains)  # Final summary across all chains

inv.phylo <- inverseA(rooted_tree, nodes = "ALL", scale = TRUE)

df <- df[df$BirdTree %in% rownames(inv.phylo$Ainv), ]

# Define prior
prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
              R = list(V = 1, nu = 0.002))

# Create list to store model chains
mcmcglmm_mod <- list()

# Run MCMCglmm for 3 chains
for (chain in 1:3) {
  mcmcglmm_mod[[chain]] <- MCMCglmm(
    Alarm_Presence~ Pre_Bird_Height + SocialGroup, 
    random = ~ BirdTree, 
    family = "categorical",
    ginverse = list(BirdTree = inv.phylo$Ainv),  # `ginverse` is outside `data`
    prior = prior,  # Prior is also outside `data`
    data = df, 
    nitt = 25000, 
    burnin = 5000, 
    thin = 1000, 
    verbose = TRUE, 
    pr = TRUE
  )
}

# Save the model
saveRDS(mcmcglmm_mod, here("output/models/mcmcglmm_mod.rds"))

# Diagnose model, assess satisfactory convergence
par(mfrow = c(3, 2), mar = c(2, 2, 1, 2))  # Adjust the plotting window

# Plot traces for each model chain
plot(do.call(mcmc.list, lapply(mcmcglmm_mod, function(m) m$Sol)), ask = FALSE)

# Check convergence with Gelman-Rubin diagnostic
gelman_diag_sol <- gelman.diag(do.call(mcmc.list, lapply(mcmcglmm_mod, function(m) m$Sol)))
gelman_diag_sol

# Ensure Rhat values are below 1.1
range(gelman_diag_sol$psrf[, 1])

# Inspect model summary (Neff values, etc.)
summary(mcmcglmm_mod[[1]])

# Perform inference on a single model's chain
summary(mcmcglmm_mod[[1]])

matched_species <- levels(cleaned_data$BirdTree) %in% rownames(inv.phylo$Ainv)
cleaned_data <- cleaned_data[cleaned_data$BirdTree %in% rownames(inv.phylo$Ainv), ]
cleaned_data$BirdTree <- droplevels(cleaned_data$BirdTree)


mcmcglmm_mod[[chain]] <- MCMCglmm(Alarm_Presence ~ Pre_Bird_Height + SocialGroup, 
                                  random = ~BirdTree, 
                                  data = cleaned_data, 
                                  family = "categorical",
                                  ginverse = list(BirdTree = inv.phylo$Ainv), 
                                  prior = prior, 
                                  nitt = 25000, 
                                  burnin = 5000, 
                                  thin = 1000, 
                                  verbose = TRUE, 
                                  pr = TRUE)
