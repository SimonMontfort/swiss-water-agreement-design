### Swiss Water Coooperation
## Authors: Simon Montfort, Manuel Fischer, James Hollway, Nicolas Jager
# R.version
# platform       x86_64-apple-darwin17.0     
# arch           x86_64                      
# os             darwin17.0                  
# system         x86_64, darwin17.0          
# status                                     
# major          4                           
# minor          0.4                         
# year           2021                        
# month          02                          
# day            15                          
# svn rev        80002                       
# language       R                           
# version.string R version 4.0.4 (2021-02-15)
# nickname       Lost Library Book   
######################################
### Load packages
######################################

rm(list = ls())

library("dplyr")
library("ggplot2")
library("summarytools")
library("survival")
library("survminer")
library("igraph")
library("ggrepel")
library("Hmisc")
library("sandwich")

######################################
# description of code
######################################

# The code proceeds in steps:
# (1) loading treaty data and subsetting to bilaterals
# (2) definition of risk set to contiguity over land
# (3) variables and covariates 
# (4) descriptives
# (5) analysis
# (6) robustness checks

######################################
# (1) treaty design data
######################################

# change directory here to replicate 
setwd("/Users/simon/Documents/repo/swiss-water-agreement-design/ModelInput")

load("dyadicdat_2.RData")

# subset to bilateral agreements
dyadicdat <- dyadicdat %>% filter(bilateral == 1)

# # subset to those issues included
dyadicdat <- dyadicdat %>%
  filter(pollution == 1 | fish == 1 | shipping == 1 | construction == 1 )


##########
# (2) risk set
##########

# risk set rivers
load("risk_set.RData")

# join datasets
dyadicdat <- left_join(dyads, dyadicdat %>% dplyr::select(-canton1, -canton2), by = c("cantons" = "cantons"))

# for those cantons which did not experience an event at the end of the study time, i.e. empty risk sets should have the end of the study period in the relevant time variable
dyadicdat$time[is.na(dyadicdat$time)] <- as.numeric(as.Date("2020-01-01")) # in numeric format, 18262 is the end of the study period, this is lower than 

# replace NAs with zero because these are the observations that did not experience the event. For all those observations, the data is NA although the mechanisms (DVs) should be zero as
# to event did not occur
dyadicdat[is.na(dyadicdat)] <- 0

############
# (3) variables
############

# I. sequencing
dyadicdat <- dyadicdat %>%
  group_by(cantons) %>%
  arrange(cantons, time) %>%
  # cumulative number of prior agreements
  mutate(agr_cum = row_number() -1,
         conflict_cum = cumsum(conflict) -1,
         monitoring_cum = cumsum(monitoring) -1,
         commission_cum = cumsum(commission) -1) %>%
  mutate(agr_cum = ifelse(agr_cum == -1, agr_cum +1, agr_cum), 
         conflict_cum = ifelse(conflict_cum == -1, conflict_cum +1, conflict_cum),
         monitoring_cum = ifelse( monitoring_cum == -1, monitoring_cum +1, monitoring_cum),
         commission_cum = ifelse(commission_cum == -1, commission_cum +1, commission_cum)) %>% 
  # at least one prior agreement
  mutate(agr_one = ifelse(agr_cum >= 1, 1, 0),
         conflict_one = ifelse(conflict_cum >= 1, 1, 0),
         monitoring_one = ifelse(monitoring_cum >= 1, 1, 0),
         commission_one = ifelse(commission_cum >= 1, 1, 0)) %>% 
  # square root of the cumulative number of prior agreements
  mutate(agreement_sqrt = sqrt(agr_cum),
         conflict_sqrt = sqrt(conflict_cum),
         monitoring_sqrt = sqrt(monitoring_cum),
         commission_sqrt =sqrt(commission_cum)) 

# inspect data setup, check sequencing variables
dyadicdat %>%
  arrange(cantons, time) %>%
  dplyr::select(cantons, time, treaty_yes, conflict, monitoring, commission, agr_cum, conflict_cum, monitoring_cum, commission_cum) %>%
  as.data.frame()

# II. external mechanisms
dyadicdat <- dyadicdat %>%
  mutate(any_ext_yes = ifelse(mon_comm == 1 | conf_body == 1 | conf_fed == 1 | commission == 1, 1, 0),
         any_mech_yes = ifelse(commission == 1 | monitoring == 1 | conflict == 1, 1, 0))

## look at structure of the data to check that cummulative sums have been calculated correctly
dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(id = cur_group_id()) %>% 
  dplyr::select(id, cantons, time, enforcement, monitoring, conflict, commission_cum, monitoring_cum, conflict_cum, agr_cum) %>% as.data.frame()

# add dyad id to cluster s.e. later for modelling
dyadicdat <- dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(dyad_id = cur_group_id())

## III. upstream/downstream operationalisation 
load("US_DS_matrix_polygons_within_2.RData")
for (i in 1:nrow(dyadicdat)){
  dyadicdat$drain_1_to_2[i] <- mat[dyadicdat$canton1[i] == rownames(mat), dyadicdat$canton2[i] == colnames(mat)]
  dyadicdat$drain_2_to_1[i] <- mat[dyadicdat$canton2[i] == rownames(mat), dyadicdat$canton1[i] == colnames(mat)]
}

# ## Water Symmetry based on trade interdependence by Crescenzi but adapted to fit water interdependence, p.30
dyadicdat$total_wat <- (dyadicdat$drain_1_to_2 + dyadicdat$drain_2_to_1)
dyadicdat$share_drain_1 <- dyadicdat$drain_1_to_2/dyadicdat$total_wat
dyadicdat$share_drain_2 <- dyadicdat$drain_2_to_1/dyadicdat$total_wat
dyadicdat$symmetry <- 1 - abs(dyadicdat$share_drain_1 - dyadicdat$share_drain_2)

# adapt risk set becuse cantons cannot cooperate when no water flows from one to the other
dyadicdat <- dyadicdat[!is.nan(dyadicdat$symmetry),]

# IV. language
load("bi_lingue.RData")
dyadicdat <- left_join(dyadicdat, lang, by = c("dyad_id"))
rm(bi_lingue)

#################
# (4) descriptives
#################

## for descriptives, subset to those dyads that experienced an agreement after 1980 or did not experience an agreement at all
desc_dat <- dyadicdat[dyadicdat$year >= 1980 & dyadicdat$year != 0,]

# I. correlation between design principles
chi <- as.data.frame(matrix(nrow = 6, ncol = 3))
chi[1,2] <- round(chisq.test(table(desc_dat$conflict, desc_dat$monitoring))$statistic, 3)
chi[2,2] <- paste0( "(", round(chisq.test(table(desc_dat$conflict, desc_dat$monitoring))$p.value, 3), ")")
chi[1,3] <- round(chisq.test(table(desc_dat$conflict, desc_dat$commission))$statistic, 3)
chi[2,3] <- paste0( "(", round(chisq.test(table(desc_dat$conflict, desc_dat$commission))$p.value, 3), ")")
chi[3,3] <- round(chisq.test(table(desc_dat$monitoring, desc_dat$commission))$statistic, 3)
chi[4,3] <- paste0( "(", round(chisq.test(table(desc_dat$monitoring, desc_dat$commission))$p.value, 3), ")")
chi[is.na(chi)] <- "-"
colnames(chi) <- c("Conflict", "Monitoring",  "Commission")
rownames(chi) <- c("Conflict", "", "Monitoring", " ",  "Commission", "  ")
sjPlot::tab_df(chi, show.rownames = T, footnote = "P-value in parentheses", show.footnote = T)

# II. simple descriptives of co-appearance
co_app <- as.data.frame(matrix(nrow = 4, ncol = 3))
co_app[1,1] <- table(desc_dat$pollution, desc_dat$commission)[2,2]
co_app[1,2] <- table(desc_dat$pollution, desc_dat$monitoring)[2,2]
co_app[1,3] <- table(desc_dat$pollution, desc_dat$conflict)[2,2]
co_app[2,1] <- table(desc_dat$shipping, desc_dat$commission)[2,2]
co_app[2,2] <- table(desc_dat$shipping, desc_dat$monitoring)[2,2]
co_app[2,3] <- table(desc_dat$shipping, desc_dat$conflict)[2,2]
co_app[3,1] <- table(desc_dat$fish, desc_dat$commission)[2,2]
co_app[3,2] <- table(desc_dat$fish, desc_dat$monitoring)[2,2]
co_app[3,3] <- table(desc_dat$fish, desc_dat$conflict)[2,2]
co_app[4,1] <- table(desc_dat$construction, desc_dat$commission)[2,2]
co_app[4,2] <- table(desc_dat$construction, desc_dat$monitoring)[2,2]
co_app[4,3] <- table(desc_dat$construction, desc_dat$conflict)[2,2]
colnames(co_app) <- c("Commission", "Monitoring",  "Conflict")
rownames(co_app) <- c("Pollution", "Shipping", "Fish", "Construction")
sjPlot::tab_df(co_app, show.rownames = T)

## III. standard summary table
library(summarytools)
sum_dat <- dyadicdat[dyadicdat$year >= 1980,] %>% ungroup() %>% 
  dplyr::select(treaty_yes, commission, monitoring, conflict, agr_cum, commission_cum, monitoring_cum, conflict_cum, pollution, shipping, fish, construction, bi_lingue, symmetry) 
sum_tab <- summarytools::descr(sum_dat, transpose = T, stats = c("n.valid", "mean", "sd", "min", "q1", "med", "q3", "max"), order = "preserve")
sum_tab <- round(sum_tab, 2)
rownames(sum_tab) <- c("Agreement (yes)", "Commission (yes)", "Monitoring (yes)", "Conflict Resolution (yes)", "Prior Agreements", "Prior Commission", "Prior Monitoring", 
                       "Prior Conflict Resolution", "Pollution (yes)", "Shipping (yes)", "Fish (yes)", "Construction (yes)", "Bi-lingue (yes)", "Symmetry")
sjPlot::tab_df(sum_tab, show.rownames = T)

## III. plot design principles and issue area over time
# create data first
p <- dyadicdat %>%
  group_by(ID_SM) %>%
  slice(1) %>%
  as_tibble() %>%
  mutate(date = as.Date(time, origin = "1970-01-01"),
         year = as.numeric(format(date,'%Y')),
         agreement_yes = ifelse(!is.na(ID_SM), 1, 0),
         Period = NA,
         Period = replace(Period,  year< 1930, "<1930"), 
         Period = replace(Period, year>= 1930 & year< 1940, "1930-1939"), 
         Period = replace(Period, year>= 1940 & year< 1950, "1940-1949"), 
         Period = replace(Period, year>= 1950 & year< 1960, "1950-1959"), 
         Period = replace(Period, year>= 1960 & year< 1970, "1960-1969"), 
         Period = replace(Period, year>= 1970 & year< 1980, "1970-1979"), 
         Period = replace(Period, year>= 1980 & year< 1990, "1980-1989"), 
         Period = replace(Period, year>= 1990 & year< 2000, "1990-1999"),
         Period = replace(Period, year>= 2000 & year< 2010, "2000-2009"),
         Period = replace(Period, year>= 2010 & year<= 2020, "2010-2020")) 

# create a custom theme for the layout
theme_SM <- function () { 
  theme_classic() +
    theme(
      text = element_text(size=12),
      axis.text=element_text(size=12, angle = 90),
      axis.ticks.length.x = unit(.2, "cm"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
      panel.grid.major.y = element_line(colour = "grey40"),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(3, "lines"),
      strip.background = element_rect(color = "white"),
      strip.text.x = element_text(size = 15),
      strip.placement = "inside",
      legend.position = "none"
    )
}

# Figure 1: Bilateral Agreements and Design principles over time
p %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c("agreement_yes", "conflict", "monitoring", "commission")) %>%
  group_by(Period, variable) %>% 
  summarise(Frequency = sum(value)) %>% 
  ggplot() + 
  geom_bar(mapping = aes(x=Period, y=Frequency, fill = Period), width=.3, stat="identity", position = "dodge") +
  facet_wrap(~variable, ncol=4, scale = "free", 
             labeller = labeller(variable = c("agreement_yes" = "Agreement (yes)", "commission" = "Commission (yes)", "monitoring" = "Monitoring (yes)", "conflict" = "Conflict (yes)"))) +
  theme_SM() +
  ylab("New Agreements") +   xlab("") +
  scale_fill_manual(values =  c(rep("grey", 6), rep("black", 4))) +
  scale_alpha_manual( guide = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, NA), color = "black") ) ) +
  scale_y_continuous(limits=c(0,20), expand = c(0, 0))
dev.off()

# Figure Appendix: Bilateral Agreements and Issue Area over Time
p %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c("pollution", "shipping", "fish", "construction")) %>%
  group_by(Period, variable) %>% 
  summarise(Frequency = sum(value)) %>% 
  ggplot() + 
  geom_bar(  mapping = aes(x=Period, y=Frequency, fill = Period) ,width=.3,
             stat="identity", position = "dodge") +
  facet_wrap(~variable, ncol=4, scale = "free",
             labeller = labeller(variable = c("pollution" = "Pollution", "shipping" = "Shipping", "fish" = "Fish", "construction" = "Construction"))) +
  theme_SM() +
  ylab("New Agreements") +   xlab("") +
  scale_fill_manual(values =  c(rep("grey", 6), rep("black", 4))) +
  scale_alpha_manual( guide = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, NA), color = "black") ) ) +
  scale_y_continuous(limits=c(0,20), expand = c(0, 0))
dev.off()

############
# (5) analysis of repeated events: main model + robustness checks
############

attach(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

# Main Model: >=1
cox_1 <- coxph(Surv(time, treaty_yes) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
cox_2 <- coxph(Surv(time, treaty_yes) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
cox_3 <- coxph(Surv(time, any_mech_yes) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
cox_4 <- coxph(Surv(time, any_mech_yes) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_5 <- coxph(Surv(time, conflict) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
cox_6 <- coxph(Surv(time, conflict) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_7 <- coxph(Surv(time, monitoring) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
cox_8 <- coxph(Surv(time, monitoring) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_9 <- coxph(Surv(time, commission) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
cox_10 <- coxph(Surv(time, commission) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

# Robustness Checks 1: Cumulative number
rob_1.1 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_2.1 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.1 <-  coxph(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_4.1 <-  coxph(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.1 <-  coxph(Surv(time, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_6.1 <-  coxph(Surv(time, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.1 <-  coxph(Surv(time, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_8.1 <-  coxph(Surv(time, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_9.1 <-  coxph(Surv(time, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum,  cluster = dyad_id)
rob_10.1 <-  coxph(Surv(time, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

# Robustness Checks 2: Sqrt
rob_1.2 <- coxph(Surv(time, treaty_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_2.2 <- coxph(Surv(time, treaty_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.2 <- coxph(Surv(time, any_mech_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_4.2 <- coxph(Surv(time, any_mech_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.2 <- coxph(Surv(time, conflict) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_6.2 <- coxph(Surv(time, conflict) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.2 <- coxph(Surv(time, monitoring) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_8.2 <- coxph(Surv(time, monitoring) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_9.2 <- coxph(Surv(time, commission) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt,  cluster = dyad_id)
rob_10.2 <- coxph(Surv(time, commission) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

# Robustness Checks 3: Weibull parametric model
weib_1 <- survreg(Surv(time, treaty_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "weibull")
weib_2 <- survreg(Surv(time, treaty_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_3 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "weibull")
weib_4 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_5 <- survreg(Surv(time, conflict) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "weibull")
weib_6 <- survreg(Surv(time, conflict) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping  + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_7 <- survreg(Surv(time, monitoring) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "weibull")
weib_8 <- survreg(Surv(time, monitoring) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_9 <- survreg(Surv(time, commission) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "weibull")
weib_10 <- survreg(Surv(time, commission) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")

# Robustness Checks 4: Exponential parametric model
exp_1 <- survreg(Surv(time, treaty_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "exponential")
exp_2 <- survreg(Surv(time, treaty_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
exp_3 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "exponential")
exp_4 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
exp_5 <- survreg(Surv(time, conflict) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "exponential")
exp_6 <- survreg(Surv(time, conflict) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping  + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
exp_7 <- survreg(Surv(time, monitoring) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "exponential")
exp_8 <- survreg(Surv(time, monitoring) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
exp_9 <- survreg(Surv(time, commission) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + cluster(dyad_id), dist = "exponential")
exp_10 <- survreg(Surv(time, commission) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")

model_list_main <- list(cox_1, cox_2, cox_3, cox_4, cox_5, cox_6, cox_7, cox_8, cox_9, cox_10)
model_list_rob_1 <- list(rob_1.1, rob_2.1, rob_3.1, rob_4.1, rob_5.1, rob_6.1, rob_7.1, rob_8.1, rob_9.1, rob_10.1)
model_list_rob_2 <- list(rob_1.2, rob_2.2, rob_3.2, rob_4.2, rob_5.2, rob_6.2,rob_7.2, rob_8.2, rob_9.2, rob_10.2)
model_list_rob_3 <- list(weib_1, weib_2, weib_3, weib_4, weib_5, weib_6, weib_7, weib_8, weib_9, weib_10)
model_list_rob_4 <- list(exp_1, exp_2, exp_3, exp_4, exp_5, exp_6, exp_7, exp_8, exp_9, exp_10)

texreg::texreg(model_list_main,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Cox, Historical Effects: >= 1 after 1980")
texreg::texreg(model_list_rob_1,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Cox, Historical Effects: cumulative after 1980")
texreg::texreg(model_list_rob_2,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Cox, Historical Effects: sqrt(cumulative) after 1980")
texreg::texreg(model_list_rob_3,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Weibull, Historical Effects: cumulative after 1980")
texreg::texreg(model_list_rob_4,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Exponential, Historical Effects: cumulative after 1980")

# does not work well for me with the robust standard errors, not shown correctly
sjPlot::tab_model(cox_1, cox_2, cox_3, cox_4, cox_5, cox_6, cox_7, cox_8, cox_9, cox_10,
                  show.aic = T, show.loglik = T, show.r2 = F, show.ci = F, show.p = F, collapse.se = T, p.style = "stars",
                  vcov.fun = "HC",  vcov.type = "HC3"
                  # se = lapply(list(rob_2_tt, rob_4_tt, rob_6_tt, rob_10_tt), function(x) summary(x)$coefficients[, 4]) # does not work
                  # file = "model1.doc
)

detach(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

attach(dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])

# Robustness Checks 1: Cumulative number
rob_1.5 <-  coxph(Surv(time, treaty_yes) ~  agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_2.5 <-  coxph(Surv(time, treaty_yes) ~  agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.5 <-  coxph(Surv(time, any_mech_yes) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_4.5 <-  coxph(Surv(time, any_mech_yes) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.5 <-  coxph(Surv(time, conflict) ~  agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_6.5 <-  coxph(Surv(time, conflict) ~  agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.5 <-  coxph(Surv(time, monitoring) ~  agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_8.5 <-  coxph(Surv(time, monitoring) ~  agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_9.5 <-  coxph(Surv(time, commission) ~  agr_one + conflict_one + monitoring_one + commission_one,  cluster = dyad_id)
rob_10.5 <-  coxph(Surv(time, commission) ~  agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

model_list_rob_5 <- list(rob_1.5, rob_2.5, rob_3.5, rob_4.5, rob_5.5, rob_6.5, rob_7.5, rob_8.5, rob_9.5, rob_10.5)
texreg::texreg(model_list_rob_5,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Cox, Historical Effects: >=1 after 1945")

# Robustness Checks 1: Cumulative number
rob_1.5 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_2.5 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.5 <-  coxph(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_4.5 <-  coxph(Surv(time, any_mech_yes) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.5 <-  coxph(Surv(time, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_6.5 <-  coxph(Surv(time, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.5 <-  coxph(Surv(time, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_8.5 <-  coxph(Surv(time, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_9.5 <-  coxph(Surv(time, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum,  cluster = dyad_id)
rob_10.5 <-  coxph(Surv(time, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

model_list_rob_5 <- list(rob_1.5, rob_2.5, rob_3.5, rob_4.5, rob_5.5, rob_6.5, rob_7.5, rob_8.5, rob_9.5, rob_10.5)
texreg::texreg(model_list_rob_5,
               custom.model.names=c("Agreement", "Agreement", "Any Mech", "Any Mech", "Conflict", "Conflict", "Monitoring", "Monitoring", "Commission", "Commission"),
               caption="Model: Cox, Historical Effects: cumulative after 1945")

detach(dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])

############
# (6) output tables and result plot
############

library(splines)
library(Greg)

b <- ggpredict(cox_4, c("fish"), type = "survival")
plot(b)

cox_1 <- coxph(Surv(time, treaty_yes) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)

rob_1.1 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

library(ggplot2)


cox.zph(rob_1.1)
termplot(cox_2, term=1, se=TRUE, col.term=1, col.se=1)

org_par <- par(xaxs = "i", ask = TRUE)
plotHR(rob_1.1, term = "agr_cum", plot.bty = "o", xlim = c(0, 7), xlab = "Age")


agr_df_1 <- with(dyadicdat,
                 data.frame(agr_one = rep(mean(agr_one, na.rm = TRUE), 2),
                            commission_one = rep(mean(commission_one, na.rm = TRUE), 2),
                            monitoring_one = rep(mean(monitoring_one, na.rm = TRUE), 2),
                            conflict_one = c(0, 1),  
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            shipping = rep(mean(shipping, na.rm = TRUE), 2),
                            # fish = rep(mean(fish, na.rm = TRUE), 2),
                            construction = rep(mean(construction, na.rm = TRUE), 2),
                            symmetry = rep(mean(symmetry, na.rm = TRUE), 2),
                            bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2)
                 )
)

fit_1 <- survfit(cox_6, newdata = agr_df_1)
ggsurvplot(fit_1,
           conf.int = TRUE,
           legend.labs=c("No Prior Commission", "One Prior Commission"),
           xlab = "Time in days",
           ylab = "Probability without Additional Mechanism",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           data = agr_df_1)

agr_df_2 <- with(dyadicdat,
                 data.frame(agr_one = rep(mean(agr_one, na.rm = TRUE), 2),
                            commission_one = rep(mean(commission_one, na.rm = TRUE), 2),
                            monitoring_one = rep(mean(monitoring_one, na.rm = TRUE), 2),
                            conflict_one = c(0, 1),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            shipping = rep(mean(shipping, na.rm = TRUE), 2),
                            fish = rep(mean(fish, na.rm = TRUE), 2),
                            construction = rep(mean(construction, na.rm = TRUE), 2),
                            symmetry = rep(mean(symmetry, na.rm = TRUE), 2),
                            bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2)
                 )
)

fit_1 <- survfit(cox_4, newdata = agr_df_2)
ggsurvplot(fit_1,
           conf.int = TRUE,
           legend.labs=c("No Prior Commission", "One Prior Commission"),
           xlab = "Time in days",
           ylab = "Probability without Additional Mechanism",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           data = agr_df_2)
fit_2 <- survfit(cox_6, newdata = agr_df_2)
ggsurvplot(fit_2,
           conf.int = TRUE,
           legend.labs=c("No Prior Commission", "One Prior Commission"),
           xlab = "Time in days",
           ylab = "Probability without Additional Conflict Resolution",
           # break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           data = agr_df_2)


############
# (6) Assumptions and Fixes
############

# Models 
cox.zph(cox_2)
cox.zph(cox_4)
cox.zph(cox_6)
cox.zph(cox_8)
cox.zph(cox_10)

attach(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
# https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf p. 20
cox_2_tt <- coxph(Surv(time, treaty_yes) ~ agr_one + tt(agr_one) + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + pollution	+ shipping +	fish + tt(fish) + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_4_tt <- coxph(Surv(time, any_mech_yes) ~ agr_one + tt(agr_one) + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + pollution	+ shipping + tt(shipping) +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_6_tt <- coxph(Surv(time, commission) ~ agr_one + tt(agr_one) + commission_one + tt(commission_one) + monitoring_one + tt(monitoring_one) + conflict_one + tt(conflict_one) + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_10_tt <- coxph(Surv(time, conflict) ~ agr_one  + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + conflict_one + pollution	+ shipping + construction + symmetry + bi_lingue, cluster = dyad_id)

cox_2_tt
cox_4_tt
cox_6_tt
cox_10_tt
detach(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

# does not work well for me with the robust standard errors, not shown correctly
sjPlot::tab_model(rob_2_tt, rob_4_tt, rob_6_tt, rob_10_tt,
                  show.aic = T, show.loglik = T, show.r2 = F, show.ci = F, show.p = F, collapse.se = T, p.style = "stars"
                  # se = lapply(list(rob_2_tt, rob_4_tt, rob_6_tt, rob_10_tt), function(x) summary(x)$coefficients[, 4]) # does not work
                  # file = "model1.doc
                  )

############
# (6) Plot Results
############

# plot results:
agr_df_1 <- with(dyadicdat,
                 data.frame(agr_one = rep(mean(agr_one, na.rm = TRUE), 2),
                            commission_one = c(0, 1),
                            monitoring_one = rep(mean(monitoring_one, na.rm = TRUE), 2),
                            conflict_one = rep(mean(conflict_one, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            shipping = rep(mean(shipping, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            fish = rep(mean(fish, na.rm = TRUE), 2),
                            construction = rep(mean(construction, na.rm = TRUE), 2),
                            symmetry = rep(mean(symmetry, na.rm = TRUE), 2),
                            bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2)
                 )
)

fit_agr_1 <- survfit(cox_6, newdata = agr_df_1)
fit_agr_1 <- survfit(cox_10, newdata = agr_df_1)

ggsurvplot(fit_agr_1,
           conf.int = TRUE,
           legend.labs=c("0 Commissions", "1 Commissions"),
           xlab = "Time in days",
           ylab = "Probability without Additional Agreement",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           data = agr_df_1)
ggsurvplot(fit_agr_1,
           conf.int = TRUE,
           legend.labs=c("0 Commissions", "1 Commissions"),
           xlab = "Time in days",
           ylab = "Probability without Additional Agreement",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           data = agr_df_1)

# Plot Cox Models
results <- function(cox_1){
  names <- names(summary(cox_1)$coefficients[,2])
  yhat  <- summary(cox_1)$coefficients[, 2]
  se    <- summary(cox_1)$coefficients[, 4]
  z <- qnorm(.975)
  lower <- yhat - z*se
  upper <- yhat + z*se
  index <- 1:length(names)
  
  dv <- if (grepl("treaty_yes", cox_1$formula[2]) == T) {
        "Agreement"
      } else if (grepl("any_mech_yes", cox_1$formula[2]) == T) {
        "Any Mechanism"
      } else if (grepl("conflict", cox_1$formula[2]) == T) {
        "Conflict"
      } else if (grepl("monitoring", cox_1$formula[2]) == T) {
        "Montoring"
      } else if (grepl("commission", cox_1$formula[2]) == T) {
        "Commission"}
  
  dvs <- rep(dv, length(names))
  
  df <- data.frame(dvs, names, yhat, se, lower, upper, index)[1:4,]
  
  return(df)
}

model_list <- list(cox_1, cox_2, cox_3, cox_4, cox_5, cox_6, cox_7, cox_8, cox_9, cox_10,
     rob_1.1, rob_2.1, rob_3.1, rob_4.1, rob_5.1, rob_6.1, rob_7.1, rob_8.1, rob_9.1, rob_10.1)

results_list <- lapply(model_list, function(x) results(x))
results_list <- lapply(seq_along(results_list), function(x) cbind(results_list[[x]], unique.id=x))
results <- do.call("rbind", results_list)

results$shape[results$unique.id %in% c(1:10)] <-  "round"
results$shape[results$unique.id %in% c(11:20)] <-  "triangle"
results$shape[results$unique.id %in% c(21:30)] <-  "square"
results$shape[results$unique.id %in% c(31:40)] <-  "diamond"
results$shape <- factor(results$shape, levels= unique(results$shape))

is.even <- function(x) x %% 2 == 0
results$model <- ifelse(is.even(results$unique.id), "Main Effects", "Full Model")
results$model <- factor(results$model, levels= unique(results$model))

results$names <- rep(c("Prior Agreements","Prior Conflict Resolution", "Prior Monitoring", "Prior Commissions"), length(results$names)/4)

# to reverse order of the labels coherent with the table
results$names <- factor(results$names,levels=rev(unique(results$names)))

ggplot(results) + 
  geom_pointrange(position = position_dodge(width = 0.5),  aes(x = names, y = yhat, ymin = lower, ymax = upper, group = unique.id, col = shape, shape = model)) + 
  scale_x_discrete(limits = rev(levels("names"))) +
  coord_flip() +
  geom_hline(yintercept = 1, lty = 3) +
  ylab("95% Confidence Intervals around the Hazards Ratio") +
  xlab("Explanatory Variables") + #switch because of the coord_flip() above
  ggtitle("") +
  theme_minimal() +
  # scale_y_discrete(limits="hazard ratio") + 
  # scale_shape_manual(name = "model",
  #                    labels = c("Main Effects", "Full Model"),
  #                    values = c(1,2)) +
  scale_colour_manual(name = "shape",
                      labels = c(">= 1", "Cumulative"),
                      values = c("grey","black")) +
  theme(text=element_text(size=12, color="black"),
        strip.text.x = element_text(size = 12),
        panel.spacing = unit(3, "lines"),
        axis.text.x = element_text(size=12, vjust=0.5, color = 'black'),  # x-axis labels
        axis.text.y = element_text(size=12, vjust=0.5, color = 'black'),  # y-axis labels
        axis.title.x = element_text(size=15, vjust=0.1),                # x-title justification
        axis.title.y = element_text(size=15, vjust=1.5),                 # y-title justification
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14.5)
  ) +
  facet_wrap(~dvs, ncol=5,
             labeller = labeller(variable = c("treaty_yes" = "Agreement (yes)",
                                              "commission" = "Commission (yes)",
                                              "monitoring" = "Monitoring (yes)",
                                              "conflict" = "Conflict (yes)")))

