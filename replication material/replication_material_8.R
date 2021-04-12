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

load("dyadicdat.RData")

# create DV for agreements
dyadicdat <- dyadicdat %>%
  mutate(treaty_yes = ifelse(!is.na(ID_SM), 1, 0))

dyadicdat <- dyadicdat %>%
  mutate(year = as.numeric(format(as.Date(time, origin = "1970-01-01"),'%Y')),
         date = as.numeric(format(as.Date(time, origin = "1970-01-01"),'%Y%m%d')))

# subset to bilateral agreements
dyadicdat <- dyadicdat %>% filter(bilateral == 1)

# # subset to those issues included
dyadicdat <- dyadicdat %>%
  filter(pollution == 1 | fish == 1 | shipping == 1 | construction == 1 )

# dyadic data: also make sure that it is orderered by name
dyadicdat <- dyadicdat %>%
  mutate(canton1 = as.character(canton1),
         canton2 = as.character(canton2),
         cantons = ifelse(canton1 < canton2, paste(canton1, canton2), paste(canton2, canton1)))

##########
# (2) risk set
##########

# risk set rivers
load("risk_set.RData")

# join datasets
dyadicdat <- left_join(dyads, dyadicdat %>% dplyr::select(-canton1, -canton2), by = c("cantons" = "cantons"))

# for those cantons which did not experience an event at the end of the study time, i.e. empty risk sets should have the end of the study period in the relevant time variable
dyadicdat$time[is.na(dyadicdat$time)] <- as.numeric(as.Date("2020-01-01")) # in numeric format, 18262 is the end of the study period, this is lower than 

# bilateral agreements are only contiguous, this is just to double check that they. should be empty
dyadicdat$cantons[!dyadicdat$cantons %in% dyads$cantons] # perfect, it is empty

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
  mutate(agr_cum = row_number() -1,
         conflict_cum = cumsum(conflict) -1,
         monitoring_cum = cumsum(monitoring) -1,
         commission_cum = cumsum(commission) -1) %>%
  mutate(agr_cum = ifelse(agr_cum == -1, agr_cum +1, agr_cum),
         conflict_cum = ifelse(conflict_cum == -1, conflict_cum +1, conflict_cum),
         monitoring_cum = ifelse( monitoring_cum == -1, monitoring_cum +1, monitoring_cum),
         commission_cum = ifelse(commission_cum == -1, commission_cum +1, commission_cum))

dyadicdat <- dyadicdat %>%
  mutate(agr_one = ifelse(agr_cum >= 1, 1, 0),
         conflict_one = ifelse(conflict_cum >= 1, 1, 0),
         monitoring_one = ifelse(monitoring_cum >= 1, 1, 0),
         commission_one = ifelse(commission_cum >= 1, 1, 0)) 

dyadicdat <- dyadicdat %>%
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

# ## Water Salience, Symmetry and Interdependence based on trade interdependence by Crescenzi but adapted to fit water interdependence, p.30
dyadicdat$total_wat <- (dyadicdat$drain_1_to_2 + dyadicdat$drain_2_to_1)
dyadicdat$share_drain_1 <- dyadicdat$drain_1_to_2/dyadicdat$total_wat
dyadicdat$share_drain_2 <- dyadicdat$drain_2_to_1/dyadicdat$total_wat
dyadicdat$symmetry <- 1 - abs(dyadicdat$share_drain_1 - dyadicdat$share_drain_2)

# adapt risk set becuse cantons cannot cooperate when no water flows from one to the other
dyadicdat <- dyadicdat[!is.nan(dyadicdat$symmetry),]

# IV. language
load("bi_lingue.RData")
dyadicdat <- left_join(dyadicdat, bi_lingue, by = c("cantons" = "cantons", "ID_SM" = "ID_SM"))
rm(bi_lingue)

#################
# (4) descriptives
#################

# for descriptives, subset to those dyads that experienced an agreement after 1980 or did not experience an agreement at all
plot_dat <- dyadicdat[dyadicdat$year >= 1980 &
                        dyadicdat$year != 0,]

# I. correlation between mechanisms
chi <- matrix(nrow = 6, ncol = 3)
colnames(chi) <- c("Commission", "Monitoring",  "Conflict")
rownames(chi) <- c("Commission", "", "Monitoring", "",  "Conflict", "")

chi[1,2] <- round(chisq.test(table(plot_dat$commission, plot_dat$monitoring))$statistic, 3)
chi[2,2] <- paste0( "(", round(chisq.test(table(plot_dat$commission, plot_dat$monitoring))$p.value, 3), ")")
chi[1,3] <- round(chisq.test(table(plot_dat$commission, plot_dat$conflict))$statistic, 3)
chi[2,3] <- paste0( "(", round(chisq.test(table(plot_dat$commission, plot_dat$conflict))$p.value, 3), ")")
chi[3,3] <- round(chisq.test(table(plot_dat$conflict, plot_dat$monitoring))$statistic, 3)
chi[4,3] <- paste0( "(", round(chisq.test(table(plot_dat$commission, plot_dat$monitoring))$p.value, 3), ")")
chi[is.na(chi)] <- "-"

latex(chi, file="")

# II. simple descriptives of co-appearance
co_app <- matrix(nrow = 4, ncol = 3)

co_app[1,1] <- table(plot_dat$pollution, plot_dat$commission)[2,2]
co_app[1,2] <- table(plot_dat$pollution, plot_dat$monitoring)[2,2]
co_app[1,3] <- table(plot_dat$pollution, plot_dat$conflict)[2,2]

co_app[2,1] <- table(plot_dat$shipping, plot_dat$commission)[2,2]
co_app[2,2] <- table(plot_dat$shipping, plot_dat$monitoring)[2,2]
co_app[2,3] <- table(plot_dat$shipping, plot_dat$conflict)[2,2]

co_app[3,1] <- table(plot_dat$fish, plot_dat$commission)[2,2]
co_app[3,2] <- table(plot_dat$fish, plot_dat$monitoring)[2,2]
co_app[3,3] <- table(plot_dat$fish, plot_dat$conflict)[2,2]

co_app[4,1] <- table(plot_dat$construction, plot_dat$commission)[2,2]
co_app[4,2] <- table(plot_dat$construction, plot_dat$monitoring)[2,2]
co_app[4,3] <- table(plot_dat$construction, plot_dat$conflict)[2,2]

colnames(co_app) <- c("Commission", "Monitoring",  "Conflict")
rownames(co_app) <- c("Pollution", "Shipping", "Fish", "Construction")

latex(co_app, file="")

## III. plot design principles and issue area over time
# create data first
p <- dyadicdat %>%
  group_by(ID_SM) %>%
  slice(1) %>%
  as_tibble() %>%
  # filter(year >= 1980) %>%
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

# Figure 1: Bilateral Agreements and Design principles over time
p %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c("agreement_yes", "commission", "monitoring", "conflict")) %>%
  group_by(Period, variable) %>% 
  summarise(Frequency = sum(value)) %>% 
  ggplot() + 
  geom_bar(  mapping = aes(x=Period, y=Frequency, fill = Period), width=.3,
             stat="identity", position = "dodge") +
  facet_wrap(~variable, ncol=4, scale = "free",
             labeller = labeller(variable = c("agreement_yes" = "Agreement (yes)", 
                                              "commission" = "Commission (yes)",  
                                              "monitoring" = "Monitoring (yes)", 
                                              "conflict" = "Conflict (yes)"))) +
  theme_classic() +
  # coord_flip() +
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
  ) +
  ylab("New Agreements") +   xlab("") +
  scale_fill_manual(values =  c(rep("grey", 6), rep("black", 4))) +
  scale_alpha_manual( guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                               shape = c(16, NA),
                                                               color = "black") ) ) +
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
             labeller = labeller(variable = c("pollution" = "Pollution", 
                                              "shipping" = "Shipping",  
                                              "fish" = "Fish", 
                                              "construction" = "Construction"))) +
  theme_classic() +
  # coord_flip() +
  theme(
    text = element_text(size=12),
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
  ) +
  ylab("New Agreements") +   xlab("") +
  scale_fill_manual(values =  c(rep("grey", 6), rep("black", 4))) +
  scale_alpha_manual( guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                               shape = c(16, NA),
                                                               color = "black") ) ) +
  scale_y_continuous(limits=c(0,20), expand = c(0, 0))
dev.off()

## IV. standard summary table
library(summarytools)
sum_dat <- dyadicdat[dyadicdat$year >= 1980,] %>%  ungroup() %>%
  dplyr::select(treaty_yes, commission, monitoring, conflict, agr_cum, commission_cum, monitoring_cum, conflict_cum, pollution, shipping, fish, construction, 
                bi_lingue, symmetry) 

sum_tab <- summarytools::descr(sum_dat, transpose = T, stats = c("n.valid", "mean", "sd", "min", "q1", "med", "q3", "max"), order = "preserve")
sum_tab <- round(sum_tab, 2)
rownames(sum_tab) <- c("Agreement (yes)",
                       "Commission (yes)",
                       "Monitoring (yes)",
                       "Conflict Resolution (yes)",
                       "Prior Agreements",
                       "Prior Commission",
                       "Prior Monitoring",
                       "Prior Conflict Resolution",
                       "Pollution (yes)",
                       "Shipping (yes)",
                       "Fish (yes)",
                       "Construction (yes)",
                       "Bi-lingue (yes)",
                       "Symmetry")
latex(sum_tab, file="")

# are the data set up correctly? Triple check
dyadicdat %>%
  arrange(cantons, time) %>%
  dplyr::select(ID_SM, cantons, time, year, date, commission, monitoring, conflict, commission_cum, monitoring_cum, conflict_cum, agr_cum, treaty_yes) %>%
  as.data.frame()

############
# (5) analysis of repeated events: main model + robustness checks
############

attach(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

# # Main Model: >=1
cox_1 <- coxph(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one, cluster = dyad_id)
cox_2 <- coxph(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_3 <- coxph(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one, cluster = dyad_id)
cox_4 <- coxph(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_5 <- coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one, cluster = dyad_id)
cox_6 <- coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_7 <- coxph(Surv(time, mon_cant) ~ agr_one + commission_one + monitoring_one + conflict_one, cluster = dyad_id)
cox_8 <- coxph(Surv(time, mon_cant) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
cox_9 <- coxph(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one, cluster = dyad_id)
cox_10 <- coxph(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)

# Main Model
rob_1.1 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum, cluster = dyad_id)
rob_2.1 <-  coxph(Surv(time, treaty_yes) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.1 <-  coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum, cluster = dyad_id)
rob_4.1 <-  coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.1 <-  coxph(Surv(time, commission) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum,  cluster = dyad_id)
rob_6.1 <-  coxph(Surv(time, commission) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.1 <-  coxph(Surv(time, mon_cant) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum, cluster = dyad_id)
rob_8.1 <-  coxph(Surv(time, mon_cant) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_9.1 <-  coxph(Surv(time, conflict) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum, cluster = dyad_id)
rob_10.1 <- coxph(Surv(time, conflict) ~  agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)

# Robustness Checks 2: Sqrt
rob_1.2 <- coxph(Surv(time, treaty_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt, cluster = dyad_id)
rob_2.2 <- coxph(Surv(time, treaty_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.2 <- coxph(Surv(time, any_mech_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt, cluster = dyad_id)
rob_4.2 <- coxph(Surv(time, any_mech_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.2 <- coxph(Surv(time, commission) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt,  cluster = dyad_id)
rob_6.2 <- coxph(Surv(time, commission) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.2 <- coxph(Surv(time, mon_cant) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt, cluster = dyad_id)
rob_8.2 <- coxph(Surv(time, mon_cant) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_9.2 <- coxph(Surv(time, conflict) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt, cluster = dyad_id)
rob_10.2 <- coxph(Surv(time, conflict) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)

# Robustness Checks 3: Weibull parametric model
weib_1 <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + cluster(dyad_id), dist = "weibull")
weib_2 <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_3 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + cluster(dyad_id), dist = "weibull")
weib_4 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_5 <- survreg(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + cluster(dyad_id), dist = "weibull")
weib_6 <- survreg(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_7 <- survreg(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + cluster(dyad_id), dist = "weibull")
weib_8 <- survreg(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")
weib_9 <- survreg(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + cluster(dyad_id), dist = "weibull")
weib_10 <- survreg(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum + pollution	+ shipping  + construction + symmetry + bi_lingue + cluster(dyad_id), dist = "weibull")

reg <- function(model_list, model){
  coef_names_cox <- c("Prior Agreements", "Prior Commission", "Prior Monitoring", "Prior Conflict Resolution", "Pollution", "Shipping", "Fish", "Construction","Symmetry", "Bi-lingue")
  coef_names_weib <- c("Intercept", "Prior Agreements", "Prior Commission", "Prior Monitoring", "Prior Conflict Resolution", "Log(scale)", "Pollution", "Shipping",
                       "Fish", "Construction", "Symmetry", "Bi-lingue")
  texreg::texreg(model_list,
                  # file = file_name,
                  custom.header = list("Agreement (yes)" = 1:2, "Any Mechanism (yes)" = 3:4, "Commission (yes)" = 5:6, "Monitoring (yes)" = 7:8, "Conflict (yes)" = 9:10),
                  booktabs = T,
                  custom.coef.names = if(model == "Weibull") coef_names_weib else coef_names_cox,
                  caption = "Survival Model for Bilateral Water Agreements after 1980 based on a Weibull Distribution",
  )
}

model_list_main <- list(cox_1, cox_2, cox_3, cox_4, cox_5, cox_6, cox_7, cox_8, cox_9, cox_10)
model_list_rob_1 <- list(rob_1.1, rob_2.1, rob_3.1, rob_4.1, rob_5.1, rob_6.1, rob_7.1, rob_8.1, rob_9.1, rob_10.1)
model_list_rob_2 <- list(rob_1.2, rob_2.2, rob_3.2, rob_4.2, rob_5.2, rob_6.2,rob_7.2, rob_8.2, rob_9.2, rob_10.2)
model_list_rob_3 <- list(weib_1, weib_2, weib_3, weib_4, weib_5, weib_6, weib_7, weib_8, weib_9, weib_10)

sjPlot::tab_model(cox_1, cox_2, cox_3, cox_4, cox_5, cox_6, cox_7, cox_8, cox_9, cox_10,
                  show.aic = TRUE, show.loglik = TRUE, show.r2 = FALSE, show.ci = FALSE,
                  vcov.fun = "CL", 
                  vcov.type = "HC4",
                  show.se = T
                  # file = "model1.doc"
                  )
# sjPlot::plot_model(cox_1)

reg(model_list_main, "Cox")

# Models 
cox.zph(rob_2.1)
cox.zph(rob_4.1)
cox.zph(rob_6.1)
cox.zph(rob_8.1)
cox.zph(rob_10.1)

rob_2_tt <- coxph(Surv(time, treaty_yes) ~ agr_one + tt(agr_one) + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_4_tt <- coxph(Surv(time, any_mech_yes) ~ agr_one + tt(agr_one) + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_6_tt <- coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_10_tt <- coxph(Surv(time, conflict) ~ agr_one + tt(agr_one) + commission_one + monitoring_one + tt(monitoring_one) + conflict_one + tt(conflict_one) + pollution	+ shipping + construction + symmetry + bi_lingue, cluster = dyad_id)

sjPlot::tab_model(rob_2_tt, rob_4_tt, rob_6_tt, rob_10_tt,
                  show.aic = T, show.loglik = T, show.r2 = F, show.ci = F, show.p = F, collapse.se = T, p.style = "stars"
                  # se = lapply(list(rob_2_tt, rob_4_tt, rob_6_tt, rob_10_tt), function(x) summary(x)$coefficients[, 4])
                  # file = "model1.doc
                  )



texreg::texreg(list(rob_2_tt, rob_4_tt, rob_6_tt, rob_10_tt),
               # file = file_name,
               # custom.header = list("Agreement (yes)" = 1:2, "Any Mechanism (yes)" = 3:4, "Commission (yes)" = 5:6, "Monitoring (yes)" = 7:8, "Conflict (yes)" = 9:10),
               booktabs = T,
               # custom.coef.names = if(model == "Weibull") coef_names_weib else coef_names_cox,
               # caption = "Survival Model for Bilateral Water Agreements after 1980 based on a Weibull Distribution",
)


vfit3 <-  coxph(Surv(time, status) ~ trt + prior + karno + tt(karno),
                data=veteran,
                tt = function(x, t, ...) x * log(t+20))

# change directory here to replicate 
setwd("/Users/simon/Documents/repo/swiss-water-agreement-design/ModelOutput")
reg(model_list_main, "Cox")
reg(model_list_rob_1, "Cox")
reg(model_list_rob_2, "Cox")
reg(model_list_rob_3, "Weibull")

detach(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

# plot results:
agr_df_1 <- with(dyadicdat,
                 data.frame(agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
                            commission_cum = c(0, 2),
                            monitoring_cum = rep(mean(monitoring_cum, na.rm = TRUE), 2),
                            conflict_cum = rep(mean(conflict_cum, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            shipping = rep(mean(shipping, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            fish = rep(mean(fish, na.rm = TRUE), 2),
                            construction = rep(mean(construction, na.rm = TRUE), 2),
                            symmetry = rep(mean(symmetry, na.rm = TRUE), 2),
                            bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2)
                 )
)

fit_agr <- survfit(cox_10, newdata = agr_df_1)
ggsurvplot(fit_agr,
           conf.int = TRUE,
           legend.labs=c("0 Commissions", "2 Commissions"),
           xlab = "Time in days",
           ylab = "Probability without Additional Concflict Resolution",
           break.time.by = 4000,
           ggtheme = theme_light(), 
           palette = "grey",
           # risk.table = "abs_pct",  # absolute number and percentage at risk.
           # risk.table.y.text.col = T,# colour risk table text annotations.
           data = agr_df_1)

fit_agr <- survfit(cox_6, newdata = agr_df_1)
ggsurvplot(fit_agr,
           conf.int = TRUE,
           legend.labs=c("0 Commissions", "2 Commissions"),
           xlab = "Time in days",
           ylab = "Probability without Additional Commission",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           palette = "grey",
           # risk.table = "abs_pct",  # absolute number and percentage at risk.
           # risk.table.y.text.col = T,# colour risk table text annotations.
           # risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           # surv.median.line = "hv",  # add the median survival pointer.
           # xlim = c(0, 2000),
           # ggtheme = theme_minimal(),
           data = agr_df_1)

cox.zph(cox_2)
cox.zph(cox_4)
cox.zph(cox_6)
cox.zph(cox_8)
cox.zph(cox_10)

############
# (5) robustness checks
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

agr_df_2 <- with(dyadicdat,
                 data.frame(agr_one = rep(mean(agr_one, na.rm = TRUE), 2),
                            commission_one = rep(mean(commission_one, na.rm = TRUE), 2),
                            monitoring_one = rep(mean(monitoring_one, na.rm = TRUE), 2),
                            conflict_one = c(0, 1),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            shipping = rep(mean(shipping, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            fish = rep(mean(fish, na.rm = TRUE), 2),
                            construction = rep(mean(construction, na.rm = TRUE), 2),
                            symmetry = rep(mean(symmetry, na.rm = TRUE), 2),
                            bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2)
                 )
)

fit_agr_1 <- survfit(rob_4.4, newdata = agr_df_1)
ggsurvplot(fit_agr_1,
           conf.int = TRUE,
           legend.labs=c("0 Commissions", "2 Commissions"),
           xlab = "Time in days",
           ylab = "Probability without Additional Agreement",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           # risk.table = "abs_pct",  # absolute number and percentage at risk.
           # risk.table.y.text.col = T,# colour risk table text annotations.
           # risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           # surv.median.line = "hv",  # add the median survival pointer.
           # xlim = c(0, 2000),
           # ggtheme = theme_minimal(),
           
           data = agr_df_1)

fit_agr_2 <- survfit(rob_10.4, newdata = agr_df_2)
ggsurvplot(fit_agr_2,
           conf.int = TRUE,
           legend.labs=c("0 Commissions", "2 Commissions"),
           xlab = "Time in days",
           ylab = "Probability without Additional Agreement with Conflict Resolution Provision",
           break.time.by = 4000,
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           # risk.table = "abs_pct",  # absolute number and percentage at risk.
           # risk.table.y.text.col = T,# colour risk table text annotations.
           # risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           # surv.median.line = "hv",  # add the median survival pointer.
           # xlim = c(0, 2000),
           # ggtheme = theme_minimal(),
           data = agr_df_2)


# Cox Model
rob_1.4 <-    coxph(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_2.4 <-    coxph(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_3.4 <-    coxph(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_4.4 <-    coxph(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_5.4 <-   coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_6.4 <-   coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_7.4 <-   coxph(Surv(time, mon_cant) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_8.4 <-   coxph(Surv(time, mon_cant) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_9.4 <-     coxph(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                     cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_10.4 <-     coxph(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping  + construction + symmetry + bi_lingue,
                      cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

model_list <- list(rob_1.4, rob_2.4, 
                   rob_3.4, rob_4.4,
                   rob_5.4, rob_6.4, 
                   rob_7.4, rob_8.4, 
                   rob_9.4, rob_10.4)

cox.zph(rob_2.4)
cox.zph(rob_4.4)
cox.zph(rob_6.4)
cox.zph(rob_8.4)
cox.zph(rob_10.4)

dyadicdat <- dyadicdat %>%
  mutate(agreement_sqrt = sqrt(agr_cum),
         conflict_sqrt = sqrt(conflict_cum),
         monitoring_sqrt = sqrt(monitoring_cum),
         commission_sqrt =sqrt(commission_cum)) 



# Cox Model
rob_1.4 <-    coxph(Surv(time, treaty_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_2.4 <-    coxph(Surv(time, treaty_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_3.4 <-    coxph(Surv(time, any_mech_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_4.4 <-    coxph(Surv(time, any_mech_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_5.4 <-   coxph(Surv(time, commission) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_6.4 <-   coxph(Surv(time, commission) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_7.4 <-   coxph(Surv(time, mon_cant) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_8.4 <-   coxph(Surv(time, mon_cant) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                   cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_9.4 <-     coxph(Surv(time, conflict) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt,
                     cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_10.4 <-     coxph(Surv(time, conflict) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt + pollution	+ shipping  + construction + symmetry + bi_lingue,
                      cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

model_list <- list(rob_1.4, rob_2.4, 
                   rob_3.4, rob_4.4,
                   rob_5.4, rob_6.4, 
                   rob_7.4, rob_8.4, 
                   rob_9.4, rob_10.4)

cox.zph(rob_2.4)
cox.zph(rob_4.4)
cox.zph(rob_6.4)
cox.zph(rob_8.4)
cox.zph(rob_10.4)


# Weibull distribution because the event rate is increasing over time
weib_1 <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_2 <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_3 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_4 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_5 <- survreg(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_6 <- survreg(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_7 <- survreg(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_8 <- survreg(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_9 <- survreg(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
weib_10 <- survreg(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping  + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")

texreg::texreg(list(weib_1, weib_2, weib_3, weib_4, weib_5, weib_6, weib_7, weib_8, weib_9, weib_10),
               custom.header = list("Agreement (yes)" = 1:2, "Any Mechanism (yes)" = 3:4, "Commission (yes)" = 5:6, "Monitoring (yes)" = 7:8, "Conflict (yes)" = 9:10),
               booktabs = T,
               sideways = T,
               custom.coef.names = c("Intercept",
                                     "Prior Agreements",
                                     "Prior Commission",
                                     "Prior Monitoring",
                                     "Prior Conflict Resolution",
                                     "Log (scale)",
                                     "Pollution",
                                     "Shipping",
                                     "Fish",
                                     "Construction",
                                     "Symmetry",
                                     "Bi-lingue"),
               caption = "Survival Model for Bilateral Water Agreements after 1980 based on a Weibull Distribution",
               )

## Parametric estimation with Weibull distribution


library(survminer)
library(tidyr)

s <- with(lung,Surv(time,status))
weibull.null <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")

fKM <- survfit(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
               + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

pred.sex1 = predict(weibull.null, newdata=list(agr_cum=0,
                                               commission_cum = mean(dyadicdat$commission_cum, na.rm = TRUE),
                                               monitoring_cum = mean(dyadicdat$monitoring_cum, na.rm = TRUE),
                                               conflict_cum = mean(dyadicdat$conflict_cum, na.rm = TRUE),
                                               pollution = mean(dyadicdat$pollution, na.rm = TRUE),                                               
                                               shipping = mean(dyadicdat$shipping, na.rm = TRUE),
                                               fish = mean(dyadicdat$fish, na.rm = TRUE),
                                               construction = mean(dyadicdat$construction, na.rm = TRUE),
                                               symmetry = mean(dyadicdat$symmetry, na.rm = TRUE),
                                               bi_lingue = mean(dyadicdat$bi_lingue, na.rm = TRUE),
                                               type="quantile"))
pred.sex2 = predict(weibull.null, newdata=list(agr_cum=2,
                                               commission_cum = mean(dyadicdat$commission_cum),
                                               monitoring_cum = mean(dyadicdat$monitoring_cum),
                                               conflict_cum = mean(dyadicdat$conflict_cum, na.rm = TRUE),
                                               pollution = mean(dyadicdat$pollution, na.rm = TRUE),                                               
                                               shipping = mean(dyadicdat$shipping, na.rm = TRUE),
                                               fish = mean(dyadicdat$fish, na.rm = TRUE),
                                               construction = mean(dyadicdat$construction, na.rm = TRUE),
                                               symmetry = mean(dyadicdat$symmetry, na.rm = TRUE),
                                               bi_lingue = mean(dyadicdat$bi_lingue, na.rm = TRUE),
                                               type="quantile"))

df = data.frame(y=seq(1,7,by=1), sex1=pred.sex1, sex2=pred.sex2)
df_long = gather(df, key= "sex", value="time", -y)

ggsurvplot(fKM[c(1,3)], data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], risk.table = T,
           legend = "none")





############
# (6) robustness checks
############

dyadicdat <- dyadicdat %>%
  mutate(agr_one = ifelse(agr_cum >= 1, 1, 0),
         conflict_one = ifelse(conflict_cum >= 1, 1, 0),
         monitoring_one = ifelse(monitoring_cum >= 1, 1, 0),
         commission_one = ifelse(commission_cum >= 1, 1, 0)) 
  
# Weibull distribution because the event rate is increasing over time
rob_1.1 <- survreg(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one 
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_2.1 <- survreg(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_cum 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_3.1 <- survreg(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one 
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_4.1 <- survreg(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_5.1 <- survreg(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one 
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_6.1 <- survreg(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_7.1 <- survreg(Surv(time, monitoring) ~ agr_one + commission_one + monitoring_one + conflict_one
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_8.1 <- survreg(Surv(time, monitoring) ~ agr_one + commission_one + monitoring_one + conflict_one 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_9.1 <- survreg(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_10.1 <- survreg(Surv(time, conflict) ~ agr_cum + commission_one + monitoring_one + conflict_one 
                   + pollution	+ shipping  + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")


texreg::texreg(list(weib_1, weib_2, weib_3, weib_4, weib_5, weib_6, weib_7, weib_8, weib_9, weib_10),
               custom.header = list("Agreement (yes)" = 1:2, "Any Mechanism (yes)" = 3:4, "Commission (yes)" = 5:6, "Monitoring (yes)" = 7:8, "Conflict (yes)" = 9:10),
               booktabs = T,
               sideways = T,
               custom.coef.names = c("Intercept",
                                     "Prior Agreements",
                                     "Prior Commission",
                                     "Prior Monitoring",
                                     "Prior Conflict Resolution",
                                     "Log (scale)",
                                     "Pollution",
                                     "Shipping",
                                     "Fish",
                                     "Construction",
                                     "Symmetry",
                                     "Bi-lingue"),
               caption = "Survival Model for Bilateral Water Agreements after 1980 based on a Weibull Distribution",
)

dyadicdat <- dyadicdat %>%
  mutate(agreement_sqrt = sqrt(agr_cum),
         conflict_sqrt = sqrt(conflict_cum),
         monitoring_sqrt = sqrt(monitoring_cum),
         commission_sqrt =sqrt(commission_cum)) 

# Weibull distribution because the event rate is increasing over time
rob_1.2 <- survreg(Surv(time, treaty_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_2.2 <- survreg(Surv(time, treaty_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_3.2 <- survreg(Surv(time, any_mech_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_4.2 <- survreg(Surv(time, any_mech_yes) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_5.2 <- survreg(Surv(time, commission) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_6.2 <- survreg(Surv(time, commission) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_7.2 <- survreg(Surv(time, monitoring) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_8.2 <- survreg(Surv(time, monitoring) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                  + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_9.2 <- survreg(Surv(time, conflict) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt
                  + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")
rob_10.2 <- survreg(Surv(time, conflict) ~ agreement_sqrt + commission_sqrt + monitoring_sqrt + conflict_sqrt 
                   + pollution	+ shipping  + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "weibull")

texreg::texreg(list(weib_1, weib_2, weib_3, weib_4, weib_5, weib_6, weib_7, weib_8, weib_9, weib_10),
               custom.header = list("Agreement (yes)" = 1:2, "Any Mechanism (yes)" = 3:4, "Commission (yes)" = 5:6, "Monitoring (yes)" = 7:8, "Conflict (yes)" = 9:10),
               booktabs = T, 
               sideways = T,
               custom.coef.names = c("Intercept",
                                     "Prior Agreements",
                                     "Prior Commission",
                                     "Prior Monitoring",
                                     "Prior Conflict Resolution",
                                     "Log (scale)",
                                     "Pollution",
                                     "Shipping",
                                     "Fish",
                                     "Construction",
                                     "Symmetry",
                                     "Bi-lingue"),
               caption = "Survival Model for Bilateral Water Agreements after 1980 based on a Weibull Distribution",
)

# Exponential assumes that the rate at which events occur is constant --> fit is lower
rob_1.3 <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_2.3 <- survreg(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_3.3 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_4.3 <- survreg(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_5.3 <- survreg(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_6.3 <- survreg(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_7.3 <- survreg(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_8.3 <- survreg(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_9.3 <- survreg(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum
                     + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")
rob_10.3 <- survreg(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                     + pollution	+ shipping     + construction + symmetry + bi_lingue + cluster(dyad_id), data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], dist = "exponential")

texreg::texreg(model_list2 <- list(exp_1, exp_2, exp_3, exp_4, exp_5, exp_6, exp_7, exp_8, exp_9, exp_10),
               custom.header = list("Agreement (yes)" = 1:2, "Any Mechanism (yes)" = 3:4, "Commission (yes)" = 5:6, "Monitoring (yes)" = 7:8, "Conflict (yes)" = 9:10),
               booktabs = T,
               sideways = T,
               # custom.coef.names = c("Intercept",     
               #                    "Prior Agreements",
               #                    "Prior Commission",
               #                    "Prior Monitoring",
               #                    "Prior Conflict Resolution",
               #                    "Pollution",
               #                    "Shipping",
               #                    "Fish",
               #                    "Construction",
               #                    "Symmetry",
               #                    "Bi-lingue"),
               caption = "Survival Model for Bilateral Water Agreements after 1980 based on an Exponential Distribution"
)


# Plot effect sizes of the Cox Models
results <- function(res_agreement_1){
  names <- rownames(summary(res_agreement_1)$table)
  yhat  <- summary(res_agreement_1)$table[, 1]
  se    <- summary(res_agreement_1)$table[, 2]
  z <- qnorm(.975)
  lower <- yhat - z*se
  upper <- yhat + z*se
  index <- 1:length(names)
    
    df <- data.frame(names, yhat, se, lower, upper, index)
  
  return(df)
}

model_list3 <- list(output1.1, output1.2, output2.1, output2.2, output3.1, output3.2, output4.1, output4.2, output5.1, output5.2,
                    output1.3, output1.4, output2.3, output2.4, output3.3, output3.4, output4.3, output4.4, output5.3, output5.4)

results_list <- lapply(model_list3, function(x) results(x))
results_list <- lapply(seq_along(results_list), function(x) cbind(results_list[[x]], unique.id=x))
results <- do.call("rbind", results_list)
results$dv[results$unique.id %in% c(1:2, 11:12)] <- "Agreement"
results$dv[results$unique.id %in% c(3:4, 13:14)] <- "Any Mechanism"
results$dv[results$unique.id %in% c(5:6, 15:16)] <- "Commission"
results$dv[results$unique.id %in% c(7:8, 17:18)] <- "Montoring"
results$dv[results$unique.id %in% c(9:10, 19:20)] <- "Conflict"

results$shape <- ifelse(results$unique.id %in% c(3,4,7,8,11,12,15,16,19,20), "round", "triangle")
results$shape <- factor(results$shape, levels= unique(results$shape))
results$model <- ifelse(results$unique.id %in% c(1,3,5,7,9,11,13,15,17,19), "Main Effects", "Full Model")
results$model <- factor(results$model, levels= unique(results$model))

results <- results[results$names %in% c("agr_cum", "commission_cum", "monitoring_cum", "conflict_cum"),]
results$names[results$names == "agr_cum"] <- "Prior Agreements"
results$names[results$names == "commission_cum"] <- "Prior Commissions"
results$names[results$names == "monitoring_cum"] <- "Prior Monitoring"
results$names[results$names == "conflict_cum"] <- "Prior Conflict Resolution"

# to reverse order of the labels coherent with the table
results$names <- factor(results$names, levels = rev(unique(results$names)))

ggplot(results) + 
  geom_pointrange(position = position_dodge(width = 0.5),  aes(x = names, y = yhat, ymin = lower, ymax = upper, group = unique.id, col = shape, shape = model)) + 
  scale_x_discrete(limits = rev(levels("names"))) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 3) +
  ylab("95% Confidence Intervals around the Hazards Ratio") +
  xlab("Explanatory Variables") + #switch because of the coord_flip() above
  ggtitle("") +
  theme_minimal() +
  # scale_y_discrete(limits="hazard ratio") + 
  scale_shape_manual(name = "model",
                     labels = c("Main Effects", "Full Model"),
                     values = c(1,2)) +
  scale_colour_manual(name = "shape",
                      labels = c("Weibull", "Exponential"),
                      values = c("grey","black")) +
  # scale_shape_manual(values = c(16,17,16,17, 16,17)) +
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
  facet_wrap(~dv, ncol=5,
             labeller = labeller(variable = c("treaty_yes" = "Agreement (yes)", 
                                              "commission" = "Commission (yes)",  
                                              "monitoring" = "Monitoring (yes)", 
                                              "conflict" = "Conflict (yes)"))) 

dyadicdat <- dyadicdat %>%
  mutate(agr_one = ifelse(agr_cum >= 1, 1, 0),
         conflict_one = ifelse(conflict_cum >= 1, 1, 0),
         monitoring_one = ifelse(monitoring_cum >= 1, 1, 0),
         commission_one = ifelse(commission_cum >= 1, 1, 0)) 

summary(dyadicdat$commission_one)

# Cox Model
rob_1.4 <-    coxph(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_2.4 <-    coxph(Surv(time, treaty_yes) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_3.4 <-    coxph(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_4.4 <-    coxph(Surv(time, any_mech_yes) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_5.4 <-   coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_6.4 <-   coxph(Surv(time, commission) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_7.4 <-   coxph(Surv(time, mon_cant) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_8.4 <-   coxph(Surv(time, mon_cant) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_9.4 <-     coxph(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
rob_10.4 <-     coxph(Surv(time, conflict) ~ agr_one + commission_one + monitoring_one + conflict_one + pollution	+ shipping  + construction + symmetry + bi_lingue,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

model_list <- list(rob_1.4, rob_2.4, 
                   rob_3.4, rob_4.4,
                   rob_5.4, rob_6.4, 
                   rob_7.4, rob_8.4, 
                   rob_9.4, rob_10.4)

cox.zph(rob_2.4)
cox.zph(rob_4.4)
cox.zph(rob_6.4)
cox.zph(rob_8.4)
cox.zph(rob_10.4)

library(stargazer)
stargazer(model_list,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements after 1980-01-01",
          column.separate = c(2, 1),
          dep.var.labels=c("Agreement (yes)", "Any Mechanism (yes)", "Commission (yes)", "Monitoring (yes)", "Conflict (yes)"),
          covariate.labels=c("Prior Agreements",
                             "Prior Commission",
                             "Prior Monitoring",
                             "Prior Conflict Resolution",
                             "Pollution",
                             "Shipping",
                             "Fish",
                             "Construction",
                             "Symmetry",
                             "Bi-lingue"),
          column.sep.width = "-20pt",
          notes = "standard errors clustered by dyad, $^{\\cdot}$ P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001",
          notes.append = FALSE,
          align=TRUE,
          no.space = TRUE,
          df = FALSE,
          float.env = "sidewaystable",
          stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list, function(x) summary(x)$coefficients[, 4]),
          star.char = c("\\cdot", "*", "**", "***"),
          star.cutoffs = c(.1, .05, .01, .001),
          style = "ajs",
          omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("AIC   ", as.character(round(unlist(lapply(model_list, function(x) AIC(x))), 3))),
                           c("BIC   ", as.character(round(unlist(lapply(model_list, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list, function(x) (x)$nevent))))
          )
)


dyadicdat$agr_cum <- sqrt(dyadicdat$agr_cum)
dyadicdat$commission_cum <- sqrt(dyadicdat$commission_cum)
dyadicdat$monitoring_cum <- sqrt(dyadicdat$monitoring_cum)
dyadicdat$conflict_cum <- sqrt(dyadicdat$conflict_cum)


# Cox Model
res_agreement_main_1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_agreement_main_2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 +	pollution	+ shipping + fish + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_any_mech_main_1 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_any_mech_main_2 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                +	pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                                cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_main_1 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_main_2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 +	pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_main_1 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_main_2 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 +	pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_main_1 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_main_2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 +	pollution	+ shipping  + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

model_list <- list(res_agreement_main_1, res_agreement_main_2, 
                   res_any_mech_main_1, res_any_mech_main_2,
                   res_commission_main_1, res_commission_main_2, 
                   res_monitoring_main_1, res_monitoring_main_2, 
                   res_conflict_main_1, res_conflict_main_2)


cox.zph(res_agreement_main_2)
cox.zph(res_commission_main_2)
cox.zph(res_monitoring_main_2)
cox.zph(res_conflict_main_2)

res_agreement_1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat)
res_agreement_2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat)
res_any_mech_1 <-     coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                cluster = dyad_id, data = dyadicdat)
res_any_mech_2 <-     coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                +	pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                                cluster = dyad_id, data = dyadicdat)
res_commission_1 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat)
res_commission_2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat)
res_monitoring_1 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat)
res_monitoring_2 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 + pollution	+ shipping +	fish + construction + symmetry + bi_lingue ,
                                 cluster = dyad_id, data = dyadicdat)
res_conflict_1 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                                 cluster = dyad_id, data = dyadicdat)
res_conflict_2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                                 + pollution	+ shipping  + construction + symmetry + bi_lingue,
                                 cluster = dyad_id, data = dyadicdat)

model_list <- list(res_agreement_main_1, res_agreement_main_2, 
                   res_any_mech_main_1, res_any_mech_main_2,
                   res_commission_main_1, res_commission_main_2, 
                   res_monitoring_main_1, res_monitoring_main_2, 
                   res_conflict_main_1, res_conflict_main_2)

library(stargazer)
stargazer(model_list,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements after 1980-01-01",
          column.separate = c(2, 1),
          dep.var.labels=c("Agreement (yes)", "Any Mechanism (yes)", "Commission (yes)", "Monitoring (yes)", "Conflict (yes)"),
          covariate.labels=c("Prior Agreements",
                             "Prior Commission",
                             "Prior Monitoring",
                             "Prior Conflict Resolution",
                             "Pollution",
                             "Shipping",
                             "Fish",
                             "Construction",
                             "Symmetry",
                             "Bi-lingue"),
          column.sep.width = "-20pt",
          notes = "standard errors clustered by dyad, $^{\\cdot}$ P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001",
          notes.append = FALSE,
          align=TRUE,
          no.space = TRUE,
          df = FALSE,
          float.env = "sidewaystable",
          stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list, function(x) summary(x)$coefficients[, 4]),
          star.char = c("\\cdot", "*", "**", "***"),
          star.cutoffs = c(.1, .05, .01, .001),
          style = "ajs",
          omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("AIC   ", as.character(round(unlist(lapply(model_list, function(x) AIC(x))), 3))),
                           c("BIC   ", as.character(round(unlist(lapply(model_list, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list, function(x) (x)$nevent))))
          )
)


res_agreement_1.1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum ,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_agreement_1.2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                              +	pollution	+ shipping +	fish + construction
                              + symmetry + bi_lingue,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_any_mech_1.1 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                             cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_any_mech_1.2 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                             +	pollution	+ shipping +	fish + construction
                             + symmetry + bi_lingue,
                             cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_commission_1.1 <-   coxph(Surv(time, commission) ~  
                                + agr_cum + commission_cum + monitoring_cum + conflict_cum,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_commission_1.2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                              +	pollution	+ shipping +	fish + construction
                              + symmetry + bi_lingue  
                              ,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_monitoring_1.1 <-   coxph(Surv(time, mon_cant) ~ 
                                + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                              ,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_monitoring_1.2 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                              +	pollution	+ shipping +	fish + construction
                              + symmetry + bi_lingue  
                              ,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_conflict_1.1 <-     coxph(Surv(time, conflict) ~    
                                + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                              ,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_conflict_1.2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                              +	pollution	+ shipping  + construction
                              + symmetry + bi_lingue  
                              ,
                              cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])

model_list <- list(res_agreement_1.1, res_agreement_1.2, 
                   res_any_mech_1.1, res_any_mech_1.2,
                   res_commission_1.1, res_commission_1.2, 
                   res_monitoring_1.1, res_monitoring_1.2, 
                   res_conflict_1.1, res_conflict_1.2)


stargazer(model_list,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, all bilateral agreements after 1945-01-01",
          column.separate = c(3, 1),
          dep.var.labels=c("Agreement (yes)", "Commission (yes)", "Monitoring (yes)", "Conflict (yes)"),
          covariate.labels=c("Prior Agreements",
                             "Prior Commission",
                             "Prior Monitoring",
                             "Prior Conflict Resolution",
                             "Pollution",
                             "Shipping",
                             "Fish",
                             "Construction",
                             "Symmetry",
                             "Bi-lingue"),
          column.sep.width = "-20pt",
          notes = "standard errors clustered by dyad, $^{\\cdot}$ P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001",
          notes.append = FALSE,
          align=TRUE,
          no.space = TRUE,
          df = FALSE,
          float.env = "sidewaystable",
          stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list, function(x) summary(x)$coefficients[, 4]),
          star.char = c("\\cdot", "*", "**", "***"),
          star.cutoffs = c(.1, .05, .01, .001),
          style = "ajs",
          omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("AIC", as.character(round(unlist(lapply(model_list, function(x) AIC(x))), 3))),
                           c("BIC", as.character(round(unlist(lapply(model_list, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list, function(x) (x)$nevent))))
          )
)

model_list <- list(res_agreement_1, res_agreement_2, res_agreement_1.1, res_agreement_1.2, 
                   res_any_mech_main_1, res_any_mech_main_2, res_any_mech_1.1, res_any_mech_1.2,
                   res_commission_1, res_commission_2, res_commission_1.1, res_commission_1.2, 
                   res_monitoring_1, res_monitoring_2, res_monitoring_1.1, res_monitoring_1.2, 
                   res_conflict_1, res_conflict_2, res_conflict_1.1, res_conflict_1.2)

cox.zph(res_agreement_1.2)
cox.zph(res_commission_1.2)
cox.zph(res_monitoring_1.2)
cox.zph(res_conflict_1.2)
