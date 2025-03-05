#2024-12-10 
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
library(phylolm)
library(gridExtra)
lines<-theme(axis.line.x.bottom = element_line(),axis.line.y.left = element_line())
#theme_tufte2<-theme_tufte(base_family = "Helvetica")

#load tree
tree<-read.nexus(here("data/trees/AllBirdsHackett1_summary.tre"))
# Loading Reduced GLM Data
df <- read_csv(here("data/HwkGLM.csv"))

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

# #only for spp-specific models
# #THAR & HARU
# ##########################ecological model with only THAR######################
# df_sp_thar<-cleaned_data %>% filter(BirdTree=="Thamnomanes_ardesiacus")
# #choose vars
# names(df_sp_thar)
# df_thar<-df_sp_thar %>% ungroup() %>% select(c("Timestamp", "alarm", "Audio_React",
#                         "TarSpp", "BirdTree", "Tar_Sex", "Pre_Bird_Height", "Flight_apex", "Proximity", "Veg_under_height", "Pre_Bird_Dense")) 
# str(df_thar)
# #tibbles are atrocious to work with
# df_thar<-as.data.frame(df_thar)
# #5 variables that should not be characters are - some NAs, verify in not un full dataset above
# for(i in 7:10){df_thar[,i]<-as.numeric(df_thar[,i])}
# #double check if there is row-level info in the full data set that could help salvage the NAs in the reduced
# df_thar %>% filter(!complete.cases(.)) #NAs in flight apex and Vegunder
# bad_thar_experiments<-df_thar %>% filter(!complete.cases(.)) %$% Timestamp
# View(df %>% filter(Timestamp %in% bad_thar_experiments))
# #to me it looks like the NAs are legit.... Vince should confirm this
# #work with the NAs for now
# df_thar<-df_thar%>%mutate(height_difference=Pre_Bird_Height-Flight_apex)
# #make density tractable - sum?
# df_thar$Pre_Bird_Dense_numeric<-sapply(strsplit(df_thar$Pre_Bird_Dense, ""), FUN=function(x){mean(as.numeric(x))})
# #we now have a viable analytic dataframe
# #first, we check if all continuous variables require transformation
# #to make a pairplot, drop charactervalues with too many or uninformative levels 
# ggpairs(df_thar %>% select(-which(names(df_thar) %in% c("Timestamp", "TarSpp", "BirdTree", "Pre_Bird_Dense"))))
# #flight apex, Proximity and height difference seem like they have long tails and could benefit from log-transformation
# #why? because variables must be Z-transformed before regression, this means that their SDs need to mean something
# #a one unit increase in SD is hard to interpret when the distribution is wildly non-normal
# df_thar_df_clean<-df_thar %>% transmute(alarm=alarm, sex=Tar_Sex, prebirdheight=Pre_Bird_Height, flight_apex_ln=log(Flight_apex), proximity_ln=log(Proximity),
#                                         veg_under_height=Veg_under_height, height_difference=height_difference, pre_bird_density=Pre_Bird_Dense_numeric)
# #now we scale and center all numeric predictor variables
# for(i in 3:8){
#   df_thar_df_clean[,i]<-as.numeric(scale(df_thar_df_clean[,i]))
#   names(df_thar_df_clean)[i]<-paste0(names(df_thar_df_clean)[i], "_std")
# }
# #give a last look at the true, scaled, predictor-level correlation - any visible big ones here are going to make us look out for VIFs
# ggpairs(df_thar_df_clean)
# #no surprise here, proximity and height difference are highly correlated, but that's right as we defined one as a function of the other
# #so for the regression analyses we should just pick one
# #ill arbitrarily ignore the height difference for now as that is perhaps harder to interpret
# full_thar<-glm(alarm~sex+prebirdheight_std+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial")
# #the model fits
# check_collinearity(full_thar)
# #we are ok, only now do we know that this is a model we can actually work with
# #now let's assume we care about potential pairwise interactions (forgot to add sex but point is moot)
# full_thar<-glm(alarm~sex+
#                  prebirdheight_std+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std+
#                  
#                  prebirdheight_std:flight_apex_ln_std+prebirdheight_std:proximity_ln_std+prebirdheight_std:veg_under_height_std+prebirdheight_std:pre_bird_density_std+ 
#                flight_apex_ln_std:proximity_ln_std+flight_apex_ln_std:veg_under_height_std+flight_apex_ln_std:pre_bird_density_std+
#                  proximity_ln_std:veg_under_height_std+proximity_ln_std:pre_bird_density_std+
#                  veg_under_height_std:pre_bird_density_std,
#                data = df_thar_df_clean, family = "binomial")
# #the model wont converge - why would that be? well in this model there are 34 datapoints and 16 parameters
# #you want at least 10N per model parameter
# #we are underpowered to consider this kind of model
# #hence, we need a simpler model
# #ideally, we would have a priori, strict set of interactions that we want to consider
# #barring that, in order to restrict our model space in an way that is less arbitrary, we can first search for any significant, additive effects (via backwards model selection)
# #and then probe for any pairwise interactions involving the remaining, significant variables.
# #so we go back to the additive model 
# full_thar<-glm(alarm~sex+prebirdheight_std+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial")
# summary(full_thar)
# #remove the factor with highest p-value and refit
# #that's this
# #                   Estimate Std. Error z value Pr(>|z|)
# #prebirdheight_std     0.04041    0.39994   0.101    0.920
# 
# #that info goes in the table (excel and word) - this info is relevant for the paper
# summary(glm(alarm~sex+flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
# #                   Estimate Std. Error z value Pr(>|z|)
# #sexm                  -0.1851     0.7967  -0.232    0.816
# summary(glm(alarm~flight_apex_ln_std+proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
# #                   Estimate Std. Error z value Pr(>|z|)
# #flight_apex_ln_std     0.1733     0.5541   0.313    0.754
# summary(glm(alarm~proximity_ln_std+veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
# #                   Estimate Std. Error z value Pr(>|z|)
# #proximity_ln_std      -0.1818     0.4095  -0.444    0.657
# summary(glm(alarm~veg_under_height_std+pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
# #                   Estimate Std. Error z value Pr(>|z|)
# #veg_under_height_std  -0.4762     0.4024  -1.183    0.237
# summary(glm(alarm~pre_bird_density_std, data = df_thar_df_clean, family = "binomial"))
# #                   Estimate Std. Error z value Pr(>|z|)
# #pre_bird_density_std   0.4129     0.3655   1.130    0.259
# 
# #ok - i didn't expect this, but this is interesting in the sense that it suggests that for this species, these ecological values don't matter
# #interesting because we could now predict that instead, phylogeny matters, and we have a lot of data to support this idea
# 
# ############point of discussion: do we want to do anything to explore any pairwise interactions? controlled dredging would be an option, possibly quite justifiable given the 
# #low power, but we'd have to have clear guidelines about how much model space to explore
# #full dredging is not an option because the model space exhausts our sampling effort quickly
# #one idea could be to explore all possible, plausible models with up to 3 model predictors: 3 because 3 model parameters (4 with intercept) with N34 is pushing it
# #could be 3 additive effects, could be 2 additive effects and their interaction
# 
# #discussion: 2025-01-08
# #carfeully construct model space: only those models with max 3predictors (all models with 1 predictor, all models with 2 predictors, all models with 3 predictors, all models with 2 predictors and their interaction)
# #maake model averaging table (aicc, delta aicc, waicc, one column per model coefficient and one more for its se)
# #note: as AIC is sensitive to sample size, all models that you make for here have to be fitted to the same dataset, so those 2 rows with NAs have to be excluded here
# #note that this doesnt apply above as we were doing model selection not by AIC, but by significance of parameters, which is not sensitive in the same way as AIC
# 
# 









##########################ecological model with only HARU######################


###########################
#STEP2: multispecies models with phylo:
#predictors; include the important stuff from step1 AND


#as of 2025-01-09 well assume just for example that the only variable we care about is the Pre_Bird_Height
#as more variables come in, add these to the selection df below

#prepare ecological data
df_msp<- cleaned_data %>% ungroup() %>% select(c("Timestamp", "alarm", "Audio_React",
                                                 "TarSpp", "BirdTree", "Tar_Sex", "Pre_Bird_Height",
                                                 "SocialGroup")) 
str(df_msp)
#tibbles are atrocious to work with
df_msp<-as.data.frame(df_msp)
#variables that should not be characters are - some NAs, verify in not un full dataset above
for(i in 7){df_msp[,i]<-as.numeric(df_msp[,i])}
#double check if there is row-level info in the full data set that could help salvage the NAs in the reduced
df_msp %>% filter(!complete.cases(.)) #NA in prebird height in ts 240709-095534
bad_msp_experiments<-df_msp %>% filter(!complete.cases(.)) %$% Timestamp
View(df %>% filter(Timestamp %in% bad_msp_experiments))
#to me it looks like the NA is legit.... Vince should confirm this
#exclude NAs as the mcmcglmm models cannot tolerate NAs (unlike the GLMs used above)
df_msp<-df_msp %>% filter(complete.cases(.)) #NA in prebird height in ts 240709-095534

###########################ADD: species-level traits##########################
#now we assemble the 'interspecific' variables from other datasources 
#SOCIAL GRADIENT is already in the data

#BodyMass & diet
#get this from AVONET
bm<-read.csv(here("data/ELEData/TraitData/AVONET1_BirdLife.csv"), stringsAsFactors = F)
crosswalk<-read.csv(here("data/ELEData/PhylogeneticData/BirdLife-BirdTree crosswalk.csv"), stringsAsFactors = F)
bm<-left_join(bm, crosswalk)
bm$Species1<-gsub(" ", "_", bm$Species1)
bm$Species3<-gsub(" ", "_", bm$Species3)
#check species coverage
my_species<-unique(df_msp$BirdTree)
all(my_species %in% bm$Species1)#no
all(my_species %in% bm$Species3)#yes
bm$BirdTree<-bm$Species3
df_msp<-left_join(df_msp, bm %>% group_by(Species3) %>% slice(1), by="BirdTree")
#eltontraits: foraging stratum. #weighted average like will sweet's project?
el<-read.csv(here("data/EltonTraits.csv"), stringsAsFactors = F)
el$Scientific<-gsub(" ", "_", el$Scientific)
all(my_species %in% el$Scientific)#T
el$BirdTree<-el$Scientific
df_msp<-left_join(df_msp, el %>% select(c(41, 26:29)))
#assuming level 1/4 is ground, level 2/4 is understory, level 3/4 is midcanopy and level 4/4 is canopy, 
#we can calculate the weighted average of time spent at each height
df_msp$foraging_stratum<- (df_msp$ForStrat.ground/100*0.25+df_msp$ForStrat.understory/100*0.5+df_msp$ForStrat.midhigh/100*0.75+df_msp$ForStrat.canopy/100*1)/4
#residual eye size - ausprey 2021
eyes<-read.csv(here("data/Ausprey_ProcB/Data/Ritland eyes final 20210330.csv"), stringsAsFactors = F)
eyes$Species1<-eyes$species_jetz
eyes<-left_join(eyes, crosswalk)
eyes$Species1<-gsub(" ", "_", eyes$Species1)
eyes$Species3<-gsub(" ", "_", eyes$Species3)
table(my_species %in% eyes$Species1)#F
table(my_species %in% eyes$Species3)#F
#ok this can be cleaned up but lets assume that's all we get - 17 missing data
eyes<-eyes %>% filter(Species1 %in% my_species) %>% group_by(Species1) %>% slice(1) %>% as.data.frame()
#construct residual eye size using phylogenetic regression and allometry from only these species
row.names(eyes)<-eyes$Species1
resid_ad<-resid(phylolm(ad~mass_final, data = eyes, model="lambda", phy= keep.tip(tree, eyes$Species1)))
resid_td<-resid(phylolm(td~mass_final, data = eyes, model="lambda", phy= keep.tip(tree, eyes$Species1)))
resid_eyes_df=data.frame(res_ad=resid_ad, res_td=resid_td, Species1=names(resid_ad))
eyes<-left_join(eyes, resid_eyes_df)
df_msp<-left_join(df_msp , eyes %>% transmute(Species1=Species1, res_ad=res_ad, res_td=res_td, stratum_ausprey=stratum))
#foraging strategy: Remsen classification - perhaps get this from Ari
###########################

#first to make sure we have good predictors that are not correlated we can fit an exploratory glm
df_msp_clean<-df_msp %>% select(c("Species1", "Species3", "alarm", "Tar_Sex", "Pre_Bird_Height", "SocialGroup", "Mass", "Trophic.Niche", "foraging_stratum", "res_ad","res_td"))
#do we need transformations before we do standardization?
ggpairs(df_msp_clean %>% select(-c(which(names(df_msp_clean) %in% c("Species1", "Species3"))))  )
#yep, mass and pre bird heigt
df_msp_clean$Pre_Bird_Height_ln<-log(df_msp_clean$Pre_Bird_Height+1) #plus one to avoid Inf bc you cant log 0
df_msp_clean$Mass_ln<-log(df_msp_clean$Mass)
#now standardize the continuous variables
df_msp_clean$Pre_Bird_Height_ln_std<-as.numeric(scale(df_msp_clean$Pre_Bird_Height_ln))
df_msp_clean$Mass_ln_std<-as.numeric(scale(df_msp_clean$Mass_ln))
df_msp_clean$foraging_stratum_std<-as.numeric(scale(df_msp_clean$foraging_stratum))
df_msp_clean$res_ad_std<-as.numeric(scale(df_msp_clean$res_ad))
df_msp_clean$res_td_std<-as.numeric(scale(df_msp_clean$res_td))
#ggpairs(df_msp_clean %>% select(-c(which(names(df_msp_clean) %in% c("Species1", "Species3")))))
ggpairs(df_msp_clean %>% select(c("alarm", "Tar_Sex", "Pre_Bird_Height_ln_std","SocialGroup","Trophic.Niche","foraging_stratum_std","Mass_ln_std","res_ad_std")))


full_model<-glm(alarm~Tar_Sex+Pre_Bird_Height_ln_std+
    SocialGroup+Trophic.Niche+foraging_stratum_std+Mass_ln_std+res_ad_std, data =df_msp_clean )
#doesnt fit, maybe too many (99) NAs from eye variable
full_model_noeyes<-glm(alarm~Tar_Sex+Pre_Bird_Height_ln_std+
                  SocialGroup+Trophic.Niche+foraging_stratum_std+Mass_ln_std, data =df_msp_clean )
#models fit, now evaluate collinearity assumption
check_collinearity(full_model_noeyes)
# 
# Low Correlation
# 
# Term  VIF   VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# Tar_Sex 1.33 [1.17, 1.62]         1.15      0.75     [0.62, 0.85]
# Pre_Bird_Height_ln_std 1.73 [1.48, 2.11]         1.31      0.58     [0.47, 0.68]
# foraging_stratum_std 2.11 [1.77, 2.59]         1.45      0.47     [0.39, 0.57]
# Mass_ln_std 1.38 [1.21, 1.67]         1.17      0.73     [0.60, 0.83]
# 
# Moderate Correlation
# 
# Term  VIF   VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# SocialGroup 5.88 [4.70, 7.44]         2.43      0.17     [0.13, 0.21]
# Trophic.Niche 5.65 [4.52, 7.14]         2.38      0.18     [0.14, 0.22]

#basically either choose trophic niche or social group
full_model_noeyes<-glm(alarm~Tar_Sex+Pre_Bird_Height_ln_std+
                         SocialGroup+foraging_stratum_std+Mass_ln_std, data =df_msp_clean )
check_collinearity(full_model_noeyes)
#this is what we want. good to go.



















#######################now from that plausible glm, we now adjust the correlation structure to account for phylogenetic relatedness and repeated sampling within species
#so we go from fitting glms to pglmms
#we can first try this in mcmcglmm which people like because it gives you bayesian pvalues in the output
#if that doesnt work we can do it in brms
inv.phylo <- inverseA(keep.tip(tree, df_msp_clean$Species3),nodes="TIPS",scale=TRUE)

mcmcglmm_mod <- list()
for (chain in 1:3) {
  mcmcglmm_mod[[chain]] <- MCMCglmm(
    alarm ~ Tar_Sex+Pre_Bird_Height_ln_std + SocialGroup,
    random=~Species3,
    data = df_msp_clean,
    family = "categorical",
    ginverse = list(Species3 = inv.phylo$Ainv),
    prior = prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                          R = list(V = 1, nu = 0.002)),
    nitt = 60000,
    burnin = 10000,
    thin = 50,
    verbose = TRUE
  )
}

saveRDS(mcmcglmm_mod, here("output/models/mcmcglmm_mod.rds"))
par(mfrow=c(13,2), mar=c(2,2,1,2))
plot(do.call(mcmc.list, lapply(mcmcglmm_mod, function(m) m$Sol)), ask=F) #yep
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod, function(m) m$Sol)))
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod, function(m) m$Sol)))$psrf[,1] %>% range # all below 1.0033758
summary(mcmcglmm_mod[[1]])#Nef ~1000

# Iterations = 10001:59951
# Thinning interval  = 50
# Sample size  = 1000 
# 
# DIC: 2.275227 
# 
# G-structure:  ~Species3
# 
#          post.mean l-95% CI u-95% CI eff.samp
# Species3     61672 0.001211   204831    236.7
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     90870    28258   177610    167.2
# 
# Location effects: alarm ~ Tar_Sex + Pre_Bird_Height_ln_std + SocialGroup 
# 
#                        post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)            -16380.76 -29149.21  -1353.75    2.111  0.004 ** 
# Tar_Sexm                  -85.31   -265.20     65.51 1000.000  0.312    
# Tar_Sexu                 -183.10   -450.11     34.02  845.180  0.128    
# Pre_Bird_Height_ln_std     25.50    -49.83    134.05  648.055  0.568    
# SocialGrouplek           7127.85 -14589.33  24067.43    1.363  0.700    
# SocialGroupmsf          16280.24   1365.04  29042.65    2.131  0.002 ** 
# SocialGrouppair         16033.16   1135.23  28852.19    2.128  0.004 ** 
# SocialGroupsolo         15895.46   1310.01  28880.88    2.157  0.006 ** 
# SocialGroupssf          16354.49   2135.41  29902.14    2.122 <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

df_msp_clean %>% count(SocialGroup, alarm)
# 
# SocialGroup alarm  n
# 1          ant     0  6
# 2          lek     0 15
# 3          msf     0 47
# 4          msf     1 36
# 5         pair     0 36
# 6         pair     1  6
# 7         solo     0 44
# 8         solo     1  3
# 9          ssf     0  7
# 10         ssf     1  7
#full separation visible in ant and lek, so we should lump

df_msp_clean$SocialGroup_lumped<-df_msp_clean$SocialGroup
df_msp_clean$SocialGroup_lumped[df_msp_clean$SocialGroup_lumped=="ant"]<-"msf"
df_msp_clean$SocialGroup_lumped[df_msp_clean$SocialGroup_lumped=="lek"]<-"ssf"
df_msp_clean %>% count(SocialGroup_lumped, alarm)
# SocialGroup_lumped alarm  n
# 1                msf     0 53
# 2                msf     1 36
# 3               pair     0 36
# 4               pair     1  6
# 5               solo     0 44
# 6               solo     1  3
# 7                ssf     0 22
# 8                ssf     1  7
#we have solved the issue of full separation

mcmcglmm_mod2 <- list()
for (chain in 1:4) {
  mcmcglmm_mod2[[chain]] <- MCMCglmm(
    alarm ~ Tar_Sex+Pre_Bird_Height_ln_std + SocialGroup_lumped,
    random=~Species3,
    data = df_msp_clean,
    family = "categorical",
    ginverse = list(Species3 = inv.phylo$Ainv),
    prior = prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                          R = list(V = 1, nu = 0.002)),
    nitt = 60000,
    burnin = 10000,
    thin = 50,
    verbose = TRUE,
  )
}
saveRDS(mcmcglmm_mod2, here("output/models/mcmcglmm_mod2.rds"))
par(mfrow=c(5,2), mar=c(2,2,1,2))
plot(do.call(mcmc.list, lapply(mcmcglmm_mod2, function(m) m$Sol)), ask=F) #yep
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod2, function(m) m$Sol)))
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod2, function(m) m$Sol)))$psrf[,1] %>% range # all below 1.0033758
summary(mcmcglmm_mod2[[1]])#Nef ~1000 or at least >>100

# Iterations = 10001:59951
# Thinning interval  = 50
# Sample size  = 1000 
# 
# DIC: 2.272355 
# 
# G-structure:  ~Species3
# 
#           post.mean l-95% CI u-95% CI eff.samp
# Species3    208885     2720   518835    206.2
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units     82332    26400   151497      228
# 
# Location effects: alarm ~ Tar_Sex + Pre_Bird_Height_ln_std + SocialGroup_lumped 
# 
#                        post.mean  l-95% CI  u-95% CI eff.samp pMCMC   
# (Intercept)            -322.6771 -952.4835  215.1067    499.2 0.212   
# Tar_Sexm               -116.5188 -295.1277   43.3023   1000.0 0.166   
# Tar_Sexu               -134.6686 -395.8136  190.5823   1000.0 0.342   
# Pre_Bird_Height_ln_std   40.2275  -59.5990  145.8167    879.8 0.412   
# SocialGroup_lumpedpair -198.5215 -437.0681    0.3518    872.1 0.054 . 
# SocialGroup_lumpedsolo -316.6627 -592.4230  -49.2901    483.0 0.004 **
# SocialGroup_lumpedssf   -18.4429 -325.3699  263.6359   1117.8 0.878   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#technically, this is interpretable and makes biological sense -- there's a strong effect of social strategy on response probability, 
#even though in this model, we don't have any meaningful ecological variables...

#however, the SEs of the parameters are terrible, usually indicating poor mixing, correlated variables but that isn't the case, hmm....

df_msp_clean %>% count(Tar_Sex, alarm)
#good separation on sex
#my gut feeling is that the categorical terms are heavily phylogenetically structured, so those two pieces could be competing for variance.
#first we try a ranked social gradient
#then we could try a two-step modelling approach, more on that later
#solo-pair-ssf-msf

#can we first justify it with the -albeit controversial - groupsize/abundance data?
cleaned_data$abundance %>% hist#log
cleaned_data$abundance_log<-log(cleaned_data$abundance)
cleaned_data$SocialGroup_factor<-factor(cleaned_data$SocialGroup, levels = c("solo", "ant" , "lek","pair", "ssf", "msf"))
ggplot(cleaned_data, aes(y=abundance_log, x=SocialGroup_factor, color=SocialGroup_factor))+
  geom_point()+
  geom_boxplot()
#yeah this tracks
cleaned_data$SocialGroup_factor_lumped<-cleaned_data$SocialGroup_factor
cleaned_data$SocialGroup_factor_lumped[cleaned_data$SocialGroup_factor_lumped=="ant"]<-"msf"
cleaned_data$SocialGroup_factor_lumped[cleaned_data$SocialGroup_factor_lumped=="lek"]<-"ssf"
cleaned_data$SocialGroup_factor_lumped<-factor(cleaned_data$SocialGroup_factor_lumped, levels = c("solo", "pair", "ssf", "msf"))
#are we justified in using the rank of these lumped groups as a social gradient?
grid.arrange(
  ggplot(cleaned_data, aes(y=abundance_log, x=SocialGroup_factor_lumped, color=SocialGroup_factor))+
  geom_point()+
  geom_violin()+
  geom_boxplot()+
  #theme_tufte2+
    lines,
  ggplot(cleaned_data, aes(y=abundance_log, x=SocialGroup_factor_lumped, color=SocialGroup_factor_lumped))+
    geom_point()+
    geom_violin()+
    geom_boxplot()+
    #theme_tufte2+
    lines
)
#yes we are justified
cleaned_data$SocialGroup_factor_lumped_numeric<-as.numeric(cleaned_data$SocialGroup_factor_lumped)
cor.test(cleaned_data$SocialGroup_factor_lumped_numeric, cleaned_data$abundance_log)
#Pearson's product-moment correlation
# 
# data:  cleaned_data$SocialGroup_factor_lumped_numeric and cleaned_data$abundance_log
# t = 22.313, df = 206, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.7960738 0.8767589
# sample estimates:
#   cor 
# 0.8410322 


#significant correlation so yes, we can now say that species in our experiments that were part of these social groups that are high in our rank 
#indeed were found in groups that seemed to have more individuals, even though the abundance estimates could potentially be flawed
#nevertheless, we are now justified in using this rank of social group in the model

#ok but the model above ran and that worked, so i think we might not need to use the rank as a predictor variable.
#the figure above would be a good supplemental figure
df_msp_clean$SocialGroup_lumped<-factor(df_msp_clean$SocialGroup_lumped, levels = c("solo", "pair", "ssf", "msf"))
df_msp_clean$SocialGroup_lumped_numeric<-as.numeric(df_msp_clean$SocialGroup_lumped)

mcmcglmm_mod3 <- list()
for (chain in 1:4) {
  mcmcglmm_mod3[[chain]] <- MCMCglmm(
    alarm ~ Tar_Sex+Pre_Bird_Height_ln_std + SocialGroup_lumped_numeric,
    random=~Species3,
    data = df_msp_clean,
    family = "categorical",
    ginverse = list(Species3 = inv.phylo$Ainv),
    prior = prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                          R = list(V = 1, nu = 0.002)),
    nitt = 60000,
    burnin = 10000,
    thin = 50,
    verbose = TRUE,
  )
}
saveRDS(mcmcglmm_mod3, here("output/models/mcmcglmm_mod3.rds"))
par(mfrow=c(13,2), mar=c(2,2,1,2))
plot(do.call(mcmc.list, lapply(mcmcglmm_mod3, function(m) m$Sol)), ask=F) #yep
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod3, function(m) m$Sol)))
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod3, function(m) m$Sol)))$psrf[,1] %>% range # all below 1.0033758
summary(mcmcglmm_mod3[[1]])#Nef ~1000 or at least >>100

# Iterations = 10001:59951
# Thinning interval  = 50
# Sample size  = 1000 
# 
# DIC: 2.529092 
# 
# G-structure:  ~Species3
# 
# post.mean l-95% CI u-95% CI eff.samp
# Species3    137739     6594   352979    279.6
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units     69369    17137   135058    143.8
# 
# Location effects: alarm ~ Tar_Sex + Pre_Bird_Height_ln_std + SocialGroup_lumped_numeric 
# 
#                            post.mean l-95% CI u-95% CI eff.samp pMCMC   
# (Intercept)                  -601.76 -1142.11  -173.80    373.7 0.004 **
# Tar_Sexm                     -105.58  -256.43    37.86    999.9 0.122   
# Tar_Sexu                     -124.12  -402.76   114.49   1000.0 0.308   
# Pre_Bird_Height_ln_std         41.13   -42.50   138.78    658.3 0.348   
# SocialGroup_lumped_numeric     93.30    32.32   172.69    394.1 0.004 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#could run for longer but this result seems decent and plausible
#increases in social group complexity ranking make it more likely for species to respond to alarm calls






############################################################################
######################The two step approach######################################
############################################################################



#other alternative, recommendable only if the ecological variables reveal little effect (like Thar seems to show..)..
#..is to fit an intercept only model (with phylo and rep measure), then extract the random intercepts, and 
#...do a phylolm of those random slopes -- that reduces the complexity and lets us perhaps avoid some of the ungainly large effect sizes
#and have a clean phylogeny output where we map the random intercepts onto the phylogeny and color them

mcmcglmm_mod4 <- list()
for (chain in 1:4) {
  mcmcglmm_mod4[[chain]] <- MCMCglmm(
    alarm ~ 1,
    random=~Species3,
    data = df_msp_clean,
    family = "categorical",
    ginverse = list(Species3 = inv.phylo$Ainv),
    prior = prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                          R = list(V = 1, nu = 0.002)),
    nitt = 110000,
    burnin = 10000,
    thin = 100,
    verbose = TRUE,
    pr=T
  )
}
saveRDS(mcmcglmm_mod4, here("output/models/mcmcglmm_mod4.rds"))
par(mfrow=c(10,2), mar=c(2,2,1,2))
plot(do.call(mcmc.list, lapply(mcmcglmm_mod4, function(m) m$Sol)), ask=F) #nope
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod4, function(m) m$Sol)))
# Potential scale reduction factors:
#   
#                                       Point est. Upper C.I.
# (Intercept)                                1.57       2.62
# Species3.Crotophaga_ani                    1.09       1.14
# Species3.Brotogeris_cyanoptera             1.12       1.21
# Species3.Todirostrum_chrysocrotaphum       1.10       1.21
# Species3.Platyrinchus_coronatus            1.25       1.77
# Species3.Platyrinchus_platyrhynchos        1.37       2.10
# Species3.Camptostoma_obsoletum             1.13       1.35
# Species3.Attila_bolivianus                 1.42       2.24
# Species3.Myiozetetes_cayanensis            1.44       2.27
# Species3.Myiodynastes_maculatus            1.25       1.79
# Species3.Tyrannus_melancholicus            1.40       2.18
# Species3.Pyrocephalus_rubinus              1.14       1.42
# Species3.Laniocera_hypopyrra               1.11       1.27
# Species3.Pachyramphus_minor                1.42       2.24
# Species3.Terenotriccus_erythrurus          1.31       1.93
# Species3.Lipaugus_vociferans               1.09       1.11
# Species3.Pipra_fasciicauda                 1.32       2.18
# Species3.Pipra_mentalis                    1.20       1.57
# Species3.Machaeropterus_pyrocephalus       1.22       1.64
# Species3.Lepidothrix_coronata              1.25       1.82
# Species3.Tyranneutes_stolzmanni            1.17       1.42
# Species3.Myrmotherula_longipennis          1.35       2.03
# Species3.Myrmotherula_axillaris            1.17       1.53
# Species3.Myrmotherula_menetriesii          1.11       1.30
# Species3.Dichrozona_cincta                 1.12       1.31
# Species3.Cymbilaimus_lineatus              1.34       2.03
# Species3.Thamnomanes_ardesiacus            1.46       2.33
# Species3.Thamnomanes_schistogynus          1.48       2.39
# Species3.Thamnophilus_schistaceus          1.11       1.29
# Species3.Thamnophilus_aethiops             1.11       1.28
# Species3.Myrmeciza_hemimelaena             1.18       1.54
# Species3.Phlegopsis_nigromaculata          1.10       1.15
# Species3.Rhegmatorhina_melanosticta        1.09       1.13
# Species3.Gymnopithys_salvini               1.08       1.12
# Species3.Willisornis_poecilinotus          1.09       1.16
# Species3.Myrmeciza_hyperythra              1.13       1.28
# Species3.Hypocnemoides_maculicauda         1.10       1.21
# Species3.Myrmeciza_fortis                  1.11       1.29
# Species3.Myrmoborus_leucophrys             1.52       2.49
# Species3.Myrmoborus_myotherinus            1.25       1.78
# Species3.Synallaxis_gujanensis             1.12       1.27
# Species3.Philydor_pyrrhodes                1.10       1.22
# Species3.Nasica_longirostris               1.10       1.22
# Species3.Xiphorhynchus_guttatus            1.13       1.33
# Species3.Sittasomus_griseicapillus         1.11       1.24
# Species3.Dendrocincla_merula               1.12       1.29
# Species3.Formicarius_analis                1.08       1.11
# Species3.Cyanocorax_violaceus              1.45       2.30
# Species3.Thryothorus_leucotis              1.08       1.10
# Species3.Cyphorhinus_arada                 1.09       1.11
# Species3.Ramphocaenus_melanurus            1.09       1.11
# Species3.Ammodramus_aurifrons              1.16       1.47
# Species3.Cacicus_cela                      1.48       2.39
# Species3.Habia_rubica                      1.65       2.82
# Species3.Volatinia_jacarina                1.14       1.40
# Species3.Ramphocelus_carbo                 1.52       2.48
# Species3.Sporophila_caerulescens           1.17       1.52
# Species3.Tangara_schrankii                 1.16       1.50
# Species3.Tangara_chilensis                 1.17       1.51
# Species3.Trogon_curucui                    1.21       1.64
# Species3.Trogon_melanurus                  1.19       1.56
# Species3.Trogon_collaris                   1.19       1.58
# Species3.Baryphthengus_martii              1.11       1.21
# Species3.Monasa_nigrifrons                 1.09       1.15
# Species3.Monasa_morphoeus                  1.34       2.01
# Species3.Chelidoptera_tenebrosa            1.10       1.11
# Species3.Galbula_cyanescens                1.08       1.16
# Species3.Eubucco_richardsoni               1.10       1.22
# Species3.Veniliornis_affinis               1.13       1.32
# Species3.Campephilus_rubricollis           1.13       1.33
# 
# Multivariate psrf
# 
# 1.87
gelman.diag(do.call(mcmc.list,lapply(mcmcglmm_mod4, function(m) m$Sol)))$psrf[,1] %>% range #not converging
summary(mcmcglmm_mod4[[1]])#Nef very low

# Iterations = 10001:109901
# Thinning interval  = 100
# Sample size  = 1000 
# 
# DIC: 193.5605 
# 
# G-structure:  ~Species3
# 
# post.mean l-95% CI u-95% CI eff.samp
# Species3     28.34    1.825    84.67     18.1
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units     3.988 0.000309    21.43    18.18
# 
# Location effects: alarm ~ 1 
# 
# post.mean l-95% CI u-95% CI eff.samp pMCMC   
# (Intercept)   -5.6306 -12.3236  -0.7499    40.84 0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########################can we get brms on this?###############
source(here("R/2025-02-13 files for phylo models/my_brms_funs.R"))
source(here("R/2025-02-13 files for phylo models/my_model_simp_figure.R"))
source(here("R/2025-02-13 files for phylo models/my_logihist.R"))


library(brms)
alarm_use_m1 <- brm(bf(alarm ~ 1 + (1|gr(Species3, cov = A)) + (1|species)) 
                   + bernoulli(), 
                   data = df_msp_clean %>% mutate(species=Species3),
                   data2 = list(A=ape::vcv.phylo(keep.tip(tree, df_msp_clean$Species3))),
                   control = list(adapt_delta=adapt, max_treedepth=treedepth),
                   sample_prior = TRUE, save_pars(group = FALSE),
                   chains = chains, cores=cores, iter = iter, warmup = warmup, thin = thin,
                   file =  here("output/models/brms_intercept_mod.rds"))
plot(alarm_use_m1)
summary(alarm_use_m1)
#this is converging a whole lot better than the mcmcglmm
#probably due to the recommended variance partitioning into phylo and non-phylo

# Family: bernoulli 
# Links: mu = logit 
# Formula: alarm ~ 1 + (1 | gr(Species3, cov = A)) + (1 | species) 
# Data: df_msp_clean %>% mutate(species = Species3) (Number of observations: 207) 
# Draws: 4 chains, each with iter = 25000; warmup = 15000; thin = 10;
# total post-warmup draws = 4000
# 
# Multilevel Hyperparameters:
#   ~species (Number of levels: 69) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.64      0.89     0.12     3.52 1.00     3202     3775
# 
# ~Species3 (Number of levels: 69) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.27      0.18     0.02     0.69 1.00     2913     3690
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -2.88      1.41    -6.37    -0.54 1.00     3918     3818
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

#calc phylo signal
#https://discourse.mc-stan.org/t/phylogenetic-signal-for-bernoulli-families-by-brms/19608/12
#bernoulli pglmms don't have a sigma term (as Gaussian do), but the principal seems the same
#calc the proportion of total variance due to the phylogenetically pooled variable
hyp <- "sd_Species3__Intercept^2 / (sd_Species3__Intercept^2 + sd_species__Intercept^2) = 0"
(hyp <- hypothesis(alarm_use_m1, hyp, class = NULL))
#0.12 +- 0.23
plot(hyp)
#12% of variance is attributable to phylogeny, based on Estimate

#ok we could easily extract those species level phylogenetic random intercepts and then plot them in a two-step approach
posterior_summary(alarm_use_m1) 

ri<-posterior_summary(alarm_use_m1) %>% as.data.frame() %>% mutate(param=row.names(.))
alpha=ri$Estimate[1]
output_p_alarm=ri %>% filter(grepl(ri$param, pattern="_Species3[", fixed=T)) %>% 
  mutate(RI=Estimate+alpha) %>% mutate(p_alarm=boot::inv.logit(RI)) %>%
  mutate(Species3=gsub(param, pattern="r_Species3[", replacement="", fixed=T)) %>%
  transmute(Species3=gsub(Species3, pattern=",Intercept]", replacement="", fixed=T), 
            p_alarm=p_alarm)
row.names(output_p_alarm)<-NULL

write_csv(output_p_alarm, here("output/analysis/PhyloAlarmProbability.csv"))
#                       Species3    p_alarm
# 1         Ammodramus_aurifrons 0.10470449
# 2            Attila_bolivianus 0.29755353
# 3         Baryphthengus_martii 0.01858282
# 4        Brotogeris_cyanoptera 0.02690657
# 5                 Cacicus_cela 0.24105210
# 6      Campephilus_rubricollis 0.01434371
# 7        Camptostoma_obsoletum 0.09169042
# 8       Chelidoptera_tenebrosa 0.02850721
# 9               Crotophaga_ani 0.02328796
# 10        Cyanocorax_violaceus 0.35812837
# 11        Cymbilaimus_lineatus 0.18803780
# 12           Cyphorhinus_arada 0.04057341
# 13         Dendrocincla_merula 0.01450167
# 14           Dichrozona_cincta 0.07157980
# 15         Eubucco_richardsoni 0.01579567
# 16          Formicarius_analis 0.01990713
# 17          Galbula_cyanescens 0.01564006
# 18         Gymnopithys_salvini 0.04765233
# 19                Habia_rubica 0.30184260
# 20   Hypocnemoides_maculicauda 0.04865621
# 21         Laniocera_hypopyrra 0.07377959
# 22        Lepidothrix_coronata 0.02003363
# 23         Lipaugus_vociferans 0.04501443
# 24 Machaeropterus_pyrocephalus 0.01983389
# 25            Monasa_morphoeus 0.10691829
# 26           Monasa_nigrifrons 0.04127033
# 27      Myiodynastes_maculatus 0.16176516
# 28      Myiozetetes_cayanensis 0.28185586
# 29            Myrmeciza_fortis 0.07100400
# 30       Myrmeciza_hemimelaena 0.08859437
# 31        Myrmeciza_hyperythra 0.04682334
# 32       Myrmoborus_leucophrys 0.18081317
# 33      Myrmoborus_myotherinus 0.10674859
# 34      Myrmotherula_axillaris 0.09122697
# 35    Myrmotherula_longipennis 0.13302254
# 36    Myrmotherula_menetriesii 0.07098313
# 37         Nasica_longirostris 0.01436339
# 38          Pachyramphus_minor 0.25952711
# 39          Philydor_pyrrhodes 0.01584583
# 40    Phlegopsis_nigromaculata 0.05013844
# 41           Pipra_fasciicauda 0.01799260
# 42              Pipra_mentalis 0.01933028
# 43      Platyrinchus_coronatus 0.13892305
# 44  Platyrinchus_platyrhynchos 0.21509485
# 45        Pyrocephalus_rubinus 0.09930925
# 46      Ramphocaenus_melanurus 0.04427064
# 47           Ramphocelus_carbo 0.22947886
# 48  Rhegmatorhina_melanosticta 0.04908809
# 49   Sittasomus_griseicapillus 0.01486703
# 50     Sporophila_caerulescens 0.11563006
# 51       Synallaxis_gujanensis 0.01527343
# 52           Tangara_chilensis 0.10139366
# 53           Tangara_schrankii 0.10385541
# 54    Terenotriccus_erythrurus 0.16408731
# 55      Thamnomanes_ardesiacus 0.17661756
# 56    Thamnomanes_schistogynus 0.25353416
# 57       Thamnophilus_aethiops 0.07176510
# 58    Thamnophilus_schistaceus 0.07015118
# 59        Thryothorus_leucotis 0.04099675
# 60 Todirostrum_chrysocrotaphum 0.07026660
# 61             Trogon_collaris 0.01086816
# 62              Trogon_curucui 0.01133943
# 63            Trogon_melanurus 0.01067179
# 64      Tyranneutes_stolzmanni 0.02477242
# 65      Tyrannus_melancholicus 0.22008822
# 66         Veniliornis_affinis 0.01384590
# 67          Volatinia_jacarina 0.09834256
# 68    Willisornis_poecilinotus 0.05309964
# 69      Xiphorhynchus_guttatus 0.01384750

####################
#nice now we plot the probability of alarming against intersp variables and put it on a tree

species_level_data<-df_msp_clean %>% group_by(Species3) %>% slice(1) %>% select(c(1,2,6:11))
output_p_alarm<-left_join(output_p_alarm, species_level_data)
output_p_alarm$SocialGroup_factor<-factor(output_p_alarm$SocialGroup, levels = c("solo", "ant" , "lek","pair", "ssf", "msf"))
output_p_alarm$SocialGroup_factor_lumped<-output_p_alarm$SocialGroup_factor
output_p_alarm$SocialGroup_factor_lumped[output_p_alarm$SocialGroup_factor_lumped=="ant"]<-"msf"
output_p_alarm$SocialGroup_factor_lumped[output_p_alarm$SocialGroup_factor_lumped=="lek"]<-"ssf"
output_p_alarm$SocialGroup_factor_lumped<-factor(output_p_alarm$SocialGroup_factor_lumped, levels = c("solo", "pair", "ssf", "msf"))

ggplot(output_p_alarm, aes(y=p_alarm, x=SocialGroup_factor, color=SocialGroup_factor))+
  geom_violin()+geom_boxplot()+geom_point()
ggplot(output_p_alarm, aes(y=p_alarm, x=SocialGroup_factor_lumped, color=SocialGroup_factor))+geom_violin()+geom_boxplot()+geom_point()
row.names(output_p_alarm)<-output_p_alarm$Species3
output_p_alarm$p_alarm_ln<-log(output_p_alarm$p_alarm)
output_p_alarm$SocialGroup_factor_lumped_numeric<-as.numeric(output_p_alarm$SocialGroup_factor_lumped)
#what intersp factors matter?
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3)), model="lambda")
# Call:
#   phylolm(formula = p_alarm_ln ~ SocialGroup_factor_lumped_numeric, 
#           data = output_p_alarm, phy = keep.tip(tree, output_p_alarm$Species3))
# 
# AIC logLik 
# 119.96 -56.98 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -1.2742 -0.4623  0.6207  1.3423  2.2391 
# 
# Mean tip height: 83.00898
# Parameter estimate(s) using ML:
#   sigma2: 0.01089198 
# 
# Coefficients:
#                                    Estimate    StdErr t.value   p.value    
# (Intercept)                       -3.780452  0.384227 -9.8391 1.232e-14 ***
# SocialGroup_factor_lumped_numeric  0.128620  0.044605  2.8835  0.005283 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.1104	Adjusted R-squared: 0.09712 

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3)), model="lambda")
# Call:
#   phylolm(formula = p_alarm_ln ~ SocialGroup_factor_lumped, data = output_p_alarm, 
#           phy = keep.tip(tree, output_p_alarm$Species3))
# 
# AIC logLik 
# 123.74 -56.87 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -1.2863 -0.4737  0.6085  1.3508  2.2270 
# 
# Mean tip height: 83.00898
# Parameter estimate(s) using ML:
#   sigma2: 0.01085707 
# 
# Coefficients:
#                               Estimate   StdErr t.value  p.value    
# (Intercept)                   -3.66029  0.38749 -9.4462 8.17e-14 *** lumpedsolo
# SocialGroup_factor_lumpedpair  0.18796  0.15606  1.2044 0.232808    
# SocialGroup_factor_lumpedssf   0.22367  0.24521  0.9122 0.365048    
# SocialGroup_factor_lumpedmsf   0.40646  0.14311  2.8401 0.006015 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.1132	Adjusted R-squared: 0.07232 
# still missing are msf different from ssf, and pairs etc.

## Output we want: one single table that summarize the effects of each parameter, phylogeny plot, summary table for each sig variable
## fig panel: ranked social gradient = x, y = prob alarm logged
## scale and center global parameters
## paper ready output: summary table, one row per predictor, one fig panel per important predictor

output_p_alarm$mass_ln_std<-as.numeric(scale(log(output_p_alarm$Mass)))
output_p_alarm$foraging_stratum_std<-as.numeric(scale(output_p_alarm$foraging_stratum))

m1_intersp<-phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+mass_ln_std, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda")
check_collinearity(m1_intersp)
car::vif(lm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+mass_ln_std, 
          data = output_p_alarm))
#yeah we are ok
summary(m1_intersp)
# Call:
#   phylolm(formula = p_alarm_ln ~ SocialGroup_factor_lumped_numeric + 
#             Trophic.Niche + foraging_stratum_std + mass_ln_std, data = output_p_alarm, 
#           phy = keep.tip(tree, output_p_alarm$Species3), model = "lambda")
# 
# AIC logLik 
# 115.74 -48.87 
# 
# Raw residuals:
#   Min       1Q   Median       3Q      Max 
# -1.06535  0.04828  0.71939  1.43856  2.08354 
# 
# Mean tip height: 83.00898
# Parameter estimate(s) using ML:
#   lambda : 1
# sigma2: 0.008610318 
# 
# Coefficients:
#                                     Estimate     StdErr  t.value   p.value    
# (Intercept)                       -4.4441867  0.4105230 -10.8257 6.364e-16 ***
# SocialGroup_factor_lumped_numeric  0.1247163  0.0421189   2.9610  0.004340 ** 
# Trophic.NicheGranivore             0.4834120  0.4074847   1.1863  0.240018    
# Trophic.NicheInvertivore           0.6770238  0.2641659   2.5629  0.012822 *  
# Trophic.NicheOmnivore              1.1187162  0.3254500   3.4374  0.001054 ** 
# foraging_stratum_std               0.0848540  0.0707410   1.1995  0.234896    
# mass_ln_std                        0.0040068  0.0838214   0.0478  0.962028    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.2968	Adjusted R-squared: 0.2287 
#remove mass
#                                    Estimate     StdErr  t.value   p.value    
#mass_ln_std                        0.0013299  0.0838617   0.0159 0.9873984    
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std, 
        data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#remove foraging stratum
#                                   Estimate     StdErr  t.value   p.value    
#foraging_stratum_std               0.083336  0.069400   1.2008 0.2343185    
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#the final model has arrived
#same model parameters from 2024-12-10 as 2025-02-20 but different output. Why?
#
# 2025-02-20
# Call:
#   phylolm(formula = p_alarm_ln ~ SocialGroup_factor_lumped_numeric + 
#             Trophic.Niche, data = output_p_alarm, phy = keep.tip(tree, 
#                                                                  output_p_alarm$Species3), model = "lambda")
# 
# AIC logLik 
# 113.34 -49.67 
# 
# Raw residuals:
#   Min       1Q   Median       3Q      Max 
# -1.02304 -0.01439  0.85717  1.47885  2.05876 
# 
# Mean tip height: 83.00898
# Parameter estimate(s) using ML:
#   lambda : 1
# sigma2: 0.008812688 
# 
# Coefficients:
#                                    Estimate    StdErr  t.value   p.value    
# (Intercept)                       -4.387427  0.403500 -10.8734 3.511e-16 ***
# SocialGroup_factor_lumped_numeric  0.128009  0.041444   3.0887  0.002972 ** 
# Trophic.NicheGranivore             0.239180  0.342433   0.6985  0.487413    
# Trophic.NicheInvertivore           0.618782  0.257851   2.3998  0.019326 *  
# Trophic.NicheOmnivore              1.075928  0.315226   3.4132  0.001118 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.2802	Adjusted R-squared: 0.2352 
# 
# Note: p-values and R-squared are conditional on lambda=1.

# 2024-12-10
# Call:
#   phylolm(formula = p_alarm_ln ~ SocialGroup_factor_lumped_numeric + 
#             Trophic.Niche, data = output_p_alarm, phy = keep.tip(tree, 
#                                                                  output_p_alarm$Species3), model = "lambda")
# 
# AIC logLik 
# 113.37 -49.68 
# 
# Raw residuals:
#   Min       1Q   Median       3Q      Max 
# -1.07797 -0.00685  0.82396  1.50207  2.09584 
# 
# Mean tip height: 83.00898
# Parameter estimate(s) using ML:
#   lambda : 1
# sigma2: 0.008815904 
# 
# Coefficients:
#                                       Estimate    StdErr  t.value   p.value    
# (Intercept)                       -4.381155  0.403573 -10.8559 3.756e-16 ***
# SocialGroup_factor_lumped_numeric  0.126525  0.041451   3.0524 0.0033044 ** 
# Trophic.NicheGranivore             0.197042  0.342495   0.5753 0.5670958    
# Trophic.NicheInvertivore           0.607644  0.257898   2.3561 0.0215398 *  
# Trophic.NicheOmnivore              1.096349  0.315284   3.4773 0.0009156 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.2869	Adjusted R-squared: 0.2423 
# 
# Note: p-values and R-squared are conditional on lambda=1.










#oops forgot eyes
output_p_alarm$res_ad_std<-as.numeric(scale(output_p_alarm$res_ad))
output_p_alarm$res_td_std<-as.numeric(scale(output_p_alarm$res_td))
#gotta take out nas before model will fit
output_p_alarm_nona<-output_p_alarm %>% filter(complete.cases(.))
m1_intersp_eyes<-phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+mass_ln_std+res_td_std+res_ad_std, 
                    data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda")
check_collinearity(m1_intersp_eyes)
#pick one
m1_intersp_eyes<-phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+mass_ln_std+res_td_std, 
                         data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda")
check_collinearity(m1_intersp_eyes)
#yep
car::vif(lm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+mass_ln_std+res_td_std, 
            data = output_p_alarm_nona))
#yep
summary(m1_intersp_eyes)
#remove mass
#                                       Estimate    StdErr  t.value   p.value    
#mass_ln_std                       -0.057783  0.110712 -0.5219  0.604754    
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+res_td_std, 
        data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda"))
# Coefficients:
#                                    Estimate    StdErr t.value   p.value    
# (Intercept)                       -4.722115  0.486355 -9.7092 5.863e-12 ***
# SocialGroup_factor_lumped_numeric  0.194381  0.058112  3.3449  0.001829 ** 
# Trophic.NicheGranivore             0.670541  0.501034  1.3383  0.188545    
# Trophic.NicheInvertivore           0.848537  0.306557  2.7680  0.008585 ** 
# Trophic.NicheOmnivore              1.017104  0.377449  2.6947  0.010340 *  
# foraging_stratum_std               0.161710  0.094191  1.7168  0.093945 .  
# res_td_std                         0.162279  0.095066  1.7070  0.095771 . 

#marginal effect of eyesize, remove, and use big data
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#remove foraging stratum
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#final model indicates effects of trophic level and social(numeric)



################what if we treat social as discrete from the get go?#############################
mx_model=phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche+foraging_stratum_std+res_td_std+mass_ln_std, 
        data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda")
check_collinearity(mx_model)
#good
summary(mx_model)
#remove mass
#                              Estimate    StdErr  t.value   p.value    
#mass_ln_std                   -0.057316  0.113954 -0.5030  0.618046    
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche+foraging_stratum_std+res_td_std, 
                data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda"))
#remove foraging stratum
#                              Estimate    StdErr  t.value   p.value
#foraging_stratum_std           0.157994  0.097337  1.6232  0.113043    
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche+res_td_std, 
                data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda"))
#remove eyesize, remove, and use big data
#                              Estimate    StdErr  t.value   p.value
#res_td_std                     0.122919  0.096537  1.2733  0.210651    
summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#final model

# Coefficients:
#                                 Estimate   StdErr  t.value   p.value    
# (Intercept)                   -4.26427  0.40773 -10.4585 2.565e-15 ***
# SocialGroup_factor_lumpedpair  0.18810  0.14462   1.3007  0.198183    
# SocialGroup_factor_lumpedssf   0.21693  0.23338   0.9295  0.356228    
# SocialGroup_factor_lumpedmsf   0.40082  0.13289   3.0162  0.003708 ** 
# Trophic.NicheGranivore         0.20073  0.35065   0.5724  0.569099    
# Trophic.NicheInvertivore       0.60658  0.26330   2.3038  0.024601 *  
# Trophic.NicheOmnivore          1.10213  0.32133   3.4299  0.001079 **


output_p_alarm %>% count(SocialGroup_factor_lumped, Trophic.Niche) %>% 
  pivot_wider(values_from = n, names_from  =SocialGroup_factor_lumped )
# A tibble: 4 × 5
#Trophic.Niche    solo  pair   ssf   msf
#   <chr>         <int> <int> <int> <int>
#1 Frugivore         2     1     5     3
#2 Invertivore      18    13     1    17
#3 Granivore        NA     1     1     1
#4 Omnivore         NA     2     2     2