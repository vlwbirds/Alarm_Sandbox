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

#match df spp to tree spp format
df$BirdTree<-gsub(" ", "_", df$BirdTree)

cleaned_data<-df

pruned_tree <- keep.tip(tree, cleaned_data$BirdTree)
all(pruned_tree$tip.label %in% cleaned_data$BirdTree) #T
all(cleaned_data$BirdTree %in% pruned_tree$tip.label) #T
plot(pruned_tree)

###########################
#STEP2: multispecies models with phylo:
#predictors; include the important stuff from step1 AND


#prepare ecological data
df_msp<- cleaned_data
str(df_msp)
#tibbles are atrocious to work with
df_msp<-as.data.frame(df_msp)
#variables that should not be characters are - some NAs, verify in not un full dataset above
for(i in 7){df_msp[,i]<-as.numeric(df_msp[,i])}
#double check if there is row-level info in the full data set that could help salvage the NAs in the reduced
# df_msp %>% filter(!complete.cases(.)) #NA in prebird height in ts 240709-095534
# bad_msp_experiments<-df_msp %>% filter(!complete.cases(.)) %$% Timestamp
#View(df %>% filter(Timestamp %in% bad_msp_experiments))
#to me it looks like the NA is legit.... Vince should confirm this
#exclude NAs as the mcmcglmm models cannot tolerate NAs (unlike the GLMs used above)
# df_msp<-df_msp %>% filter(complete.cases(.)) #NA in prebird height in ts 240709-095534

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
df_msp<-left_join(df_msp, bm %>% group_by(Species3), by="BirdTree")
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

inv.phylo <- inverseA(keep.tip(tree, df_msp$BirdTree),nodes="TIPS",scale=TRUE)

##########################can we get brms on this?###############
source(here("R/2025-02-13 files for phylo models/my_brms_funs.R"))
source(here("R/2025-02-13 files for phylo models/my_model_simp_figure.R"))
source(here("R/2025-02-13 files for phylo models/my_logihist.R"))

#troubleshooting small tree to big tree
# colSums(is.na(df_msp)) # SocialGroup, Trophic.Niche, and Species 3 = 0 NAs
# all(df_msp$Species3 %in% tree$tip.label) # Species3 and tree$tip.label = TRUE
# A <- ape::vcv.phylo(keep.tip(tree, df_msp$Species3))
# print(dim(A)) # same 114 and 114

df_msp <- df_msp %>%
  mutate(Alarm1 = ifelse(Audio_React == "a", 1, 0))

library(brms)
alarm_use_m1 <- brm(bf(Alarm1 ~ 1 + (1|gr(Species3, cov = A)) + (1|species)) 
                    + bernoulli(), 
                    data = df_msp %>% mutate(species=Species3),
                    data2 = list(A=ape::vcv.phylo(keep.tip(tree, df_msp$Species3))),
                    control = list(adapt_delta=adapt, max_treedepth=treedepth),
                    sample_prior = TRUE, save_pars(group = FALSE),
                    chains = chains, cores=cores, iter = iter, warmup = warmup, thin = thin,
                    file =  here("output/models/brms_intercept_mod_all.rds"))
plot(alarm_use_m1)
summary(alarm_use_m1)
#this is converging a whole lot better than the mcmcglmm
#probably due to the recommended variance partitioning into phylo and non-phylo


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

write_csv(output_p_alarm, here("output/analysis/PhyloAlarmProbability_All.csv"))

p_alarm_tar <- read_csv(here("output/analysis/PhyloAlarmProbability.csv"))
output_p_alarm <- read_csv(here("output/analysis/PhyloAlarmProbability_All.csv"))
foraging_strategy <- read_csv(here("data/ForagingStrategy.csv"))
foraging_strategy$BirdTree<-gsub(" ", "_", foraging_strategy$BirdTree)

foraging_strategy <- foraging_strategy %>% 
  rename(Species3 = BirdTree)

# Filter out species in p_alarm_all that are not in p_alarm_tar
filtered_p_alarm_all <- p_alarm_all %>%
  filter(Species3 %in% p_alarm_tar$Species3)

# Merge the two datasets by species to compare p_alarm values
merged_data <- left_join(p_alarm_tar, filtered_p_alarm_all, by = "Species3", suffix = c("_tar", "_all"))

# Run a linear model to compare p_alarm values between the two datasets
lm_model <- lm(p_alarm_tar ~ p_alarm_all, data = merged_data)

# Print the summary of the linear model to get significance
summary(lm_model)

#nice now we plot the probability of alarming against intersp variables and put it on a tree
view(df_msp)
species_level_data<-df_msp %>% group_by(Species3) %>% select(c(1,2,6:11,84,94,108))
output_p_alarm<-left_join(output_p_alarm, species_level_data)
output_p_alarm<-left_join(output_p_alarm, foraging_strategy)
output_p_alarm$SocialGroup_factor<-factor(output_p_alarm$SocialGroup, levels = c("solo", "ant" , "lek","pair", "ssf", "msf"))
output_p_alarm$SocialGroup_factor_lumped<-output_p_alarm$SocialGroup_factor
output_p_alarm$SocialGroup_factor_lumped[output_p_alarm$SocialGroup_factor_lumped=="ant"]<-"msf"
output_p_alarm$SocialGroup_factor_lumped[output_p_alarm$SocialGroup_factor_lumped=="lek"]<-"ssf"
output_p_alarm$SocialGroup_factor_lumped<-factor(output_p_alarm$SocialGroup_factor_lumped, levels = c("solo", "pair", "ssf", "msf"))

ggplot(output_p_alarm, aes(y=p_alarm, x=SocialGroup_factor, color=SocialGroup_factor))+
  geom_violin()+geom_boxplot()+geom_point()
ggplot(output_p_alarm, aes(y=p_alarm, x=SocialGroup_factor_lumped, color=SocialGroup_factor))+geom_violin()+geom_boxplot()+geom_point()
output_p_alarm <- output_p_alarm %>%
  distinct(Species3, .keep_all = TRUE) # keep only distict Species to remove duplicate rows
row.names(output_p_alarm) <- output_p_alarm$Species3
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
# 105.52 -49.76 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -1.2315 -0.2256  0.6154  1.0203  1.9202 
# 
# Mean tip height: 108.8036
# Parameter estimate(s) using ML:
#   sigma2: 0.005636675 
# 
# Coefficients:
#                                      Estimate    StdErr  t.value p.value    
#   (Intercept)                       -4.330636  0.381176 -11.3612  <2e-16 ***
#   SocialGroup_factor_lumped_numeric  0.008601  0.025954   0.3314   0.741    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.0009796	Adjusted R-squared: -0.00794 

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3)), model="lambda")

# Call:
#   phylolm(formula = p_alarm_ln ~ SocialGroup_factor_lumped, data = output_p_alarm, 
#           phy = keep.tip(tree, output_p_alarm$Species3))
# 
# AIC logLik 
# 108    -49 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -1.2315 -0.2733  0.6154  1.0350  1.9462 
# 
# Mean tip height: 108.8036
# Parameter estimate(s) using ML:
#   sigma2: 0.005562502 
# 
# Coefficients:
#                                Estimate    StdErr  t.value p.value    
# (Intercept)                   -4.369540  0.380637 -11.4796  <2e-16 ***
# SocialGroup_factor_lumpedpair  0.146291  0.121739   1.2017  0.2321    
# SocialGroup_factor_lumpedssf   0.038727  0.140039   0.2765  0.7826    
# SocialGroup_factor_lumpedmsf   0.073289  0.088602   0.8272  0.4099    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.01413	Adjusted R-squared: -0.01276 

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

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#remove foraging stratum

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))
#the final model has arrived
#same model parameters from 2024-12-10 as 2025-02-20 but different output. Why?











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

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped_numeric+Trophic.Niche+foraging_stratum_std+res_td_std, 
                data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda"))


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

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche+foraging_stratum_std+res_td_std, 
                data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda"))

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche+res_td_std, 
                data = output_p_alarm_nona, phy=keep.tip(tree, output_p_alarm_nona$Species3), model="lambda"))

summary(phylolm(p_alarm_ln~SocialGroup_factor_lumped+Trophic.Niche, 
                data = output_p_alarm, phy=keep.tip(tree, output_p_alarm$Species3), model="lambda"))



output_p_alarm %>% count(SocialGroup_factor_lumped, Trophic.Niche) %>% 
  pivot_wider(values_from = n, names_from  =SocialGroup_factor_lumped )
