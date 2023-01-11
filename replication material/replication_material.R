#### Trajectories of Institutional Design: Intergovernmental Cooperation in Federal Systems
### Authors: Simon Montfort, Manuel Fischer, James Hollway, Nicolas Jager
## Script compiled by Simon Montfort
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
library("tidyr")
library("igraph")
library("ggrepel")
library("Hmisc")
library("sandwich")
library("broom")
library("texreg")
library("scales")
library("stargazer")
library("wordreg")

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
setwd("/Users/simon/Documents/repo/swiss-water-agreement-design")

load("ModelInput/dyadicdat_2.RData")

# subset to bilateral agreements
dyadicdat <- dyadicdat %>% filter(bilateral == 1)

# subset to those issues included
dyadicdat <- dyadicdat %>% filter(pollution == 1 | fish == 1 | shipping == 1 | construction == 1 )

##########
# (2) risk set
##########

# risk set rivers
load("ModelInput/risk_set.RData")

# join datasets
dyadicdat <- left_join(dyads, dyadicdat %>% dplyr::select(-canton1, -canton2), by = c("cantons" = "cantons"))

# for those cantons which did not experience an event at the end of the study time, i.e. empty risk sets should have the end of the study period in the relevant time variable
dyadicdat$time[is.na(dyadicdat$time)] <- as.numeric(as.Date("2020-01-01")) # in numeric format, 18262 is the end of the study period, this is lower than 

# replace NAs of variables such as agreement, conflict, monitoring, commission with zero because the variables of these booleans are not present
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

# II. any mechanisms
dyadicdat <- dyadicdat %>%
  mutate(any_mech_yes = ifelse(commission == 1 | monitoring == 1 | conflict == 1, 1, 0))

## look at structure of the data to check that cummulative sums have been calculated correctly
dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(id = cur_group_id()) %>% 
  dplyr::select(id, cantons, time, enforcement, monitoring, conflict, commission_cum, monitoring_cum, conflict_cum, agr_cum) %>% as.data.frame()

# add dyad id to cluster s.e. later for modelling
dyadicdat <- dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(dyad_id = cur_group_id())

## III. upstream/downstream operationalisation 
load("ModelInput/US_DS_matrix_polygons_within_2.RData")
for (i in 1:nrow(dyadicdat)){
  dyadicdat$drain_1_to_2[i] <- mat[dyadicdat$canton1[i] == rownames(mat), dyadicdat$canton2[i] == colnames(mat)]
  dyadicdat$drain_2_to_1[i] <- mat[dyadicdat$canton2[i] == rownames(mat), dyadicdat$canton1[i] == colnames(mat)]
}

# IV. Water Symmetry based on trade interdependence by Crescenzi but adapted to fit water interdependence, p.30
dyadicdat$total_wat <- (dyadicdat$drain_1_to_2 + dyadicdat$drain_2_to_1)
dyadicdat$share_drain_1 <- dyadicdat$drain_1_to_2/dyadicdat$total_wat
dyadicdat$share_drain_2 <- dyadicdat$drain_2_to_1/dyadicdat$total_wat
dyadicdat$symmetry <- 1 - abs(dyadicdat$share_drain_1 - dyadicdat$share_drain_2)

# V. language
load("ModelInput/lang.RData")
dyadicdat <- left_join(dyadicdat, lang, by = c("dyad_id"))

# adapt risk set because cantons do not cooperate when no water flows from one to the other
dyadicdat <- dyadicdat[!is.nan(dyadicdat$symmetry),]

# number of cantons in the risk set:
length(unique(dyadicdat$cantons))
# number of repeated events:
sum(duplicated(dyadicdat$cantons[dyadicdat$treaty_yes == 1 & dyadicdat$year >= 1980]))

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
# sjPlot::tab_df(chi, show.rownames = T, footnote = "P-value in parentheses", show.footnote = T)
latex(chi, file="", caption = "Chi-Squared Correlation between Design Principles")

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
sjPlot::tab_df(co_app, show.rownames = T, title = "Design Principles by Issue Area", file = "ModelOutput/tables/co-appearance_table.doc")
# latex(co_app, file="", caption = "Design Principles by Issue Area", booktabs = T)

t.test(desc_dat$symmetry[desc_dat$commission == 1], desc_dat$symmetry[desc_dat$monitoring == 1], var.equal = F)
t.test(desc_dat$symmetry[desc_dat$commission == 1], desc_dat$symmetry[desc_dat$conflict == 1], var.equal = F)
t.test(desc_dat$symmetry[desc_dat$conflict == 1], desc_dat$symmetry[desc_dat$monitoring == 0], var.equal = F)


## III. standard summary table
library(summarytools)
sum_dat <- dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,] %>% 
  ungroup() %>% 
  dplyr::select(treaty_yes, conflict, monitoring, commission, agr_cum, conflict_cum, monitoring_cum, commission_cum, pollution, shipping, fish, construction, bi_lingue, symmetry) 
sum_tab <- summarytools::descr(sum_dat, transpose = T, stats = c("n.valid", "mean", "sd", "min", "med", "max"), order = "p")
rownames(sum_tab) <- c("Agreement (yes)","Conflict Resolution (yes)",  "Monitoring (yes)", "Commission (yes)",  "Prior Agreements", "Prior Conflict Resolution", "Prior Monitoring", 
                       "Prior Commission", "Pollution (yes)", "Shipping (yes)", "Fish (yes)", "Construction (yes)", "Bi-lingue (yes)", "Symmetry")
sum_tab <- round(sum_tab, 2)
sum_tab <- cbind(sum_tab, c(rep("no", 4), rep("yes", 4), rep("no", 6)))
colnames(sum_tab)[7] <- "Time-varying"
sjPlot::tab_df(sum_tab, show.rownames = T)

# numbers in the data:
length(unique(dyadicdat$ID_SM[dyadicdat$year >= 1980])) # number of agreements
nrow(dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == "0",]) # number of observations
unique(dyadicdat$cantons[duplicated(dyadicdat$cantons)]) # dyads that sign two or more agreements
length(unique(dyadicdat$cantons[duplicated(dyadicdat$cantons)])) # how many dyads sign two or more agreements
length(unique(dyadicdat$cantons[duplicated(dyadicdat$cantons)])) # how many dyads sign two or more agreements
length(unique(dyadicdat$cantons)) # risk set
# length(unique(dyadicdat$cantons[duplicated(dyadicdat$cantons)])) # how many dyads sign two or more agreements
length(unique(dyadicdat$ID_SM[dyadicdat$commission == 1 & dyadicdat$year >= 1980])) # number of commissions
length(unique(dyadicdat$ID_SM[dyadicdat$conflict == 1 & dyadicdat$year >= 1980])) # number of conflict
length(unique(dyadicdat$ID_SM[dyadicdat$monitoring == 1 & dyadicdat$year >= 1980])) # number of conflict

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
p1 <- p %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c("agreement_yes", "conflict", "monitoring", "commission")) %>%
  group_by(Period, variable) %>% 
  summarise(Frequency = sum(value)) %>% 
  ggplot() + 
  geom_bar(mapping = aes(x=Period, y=Frequency, fill = Period), width=.5, stat="identity", position = "dodge") +
  facet_wrap(~variable, ncol=4, scale = "free", 
             labeller = labeller(variable = c("agreement_yes" = "Agreement (yes)", "commission" = "Commission (yes)", "monitoring" = "Monitoring (yes)", "conflict" = "Conflict (yes)"))) +
  theme_SM() +
  ylab("New Agreements") +   xlab("") +
  scale_fill_manual(values =  c(rep("grey", 6), rep("black", 4))) +
  scale_alpha_manual( guide = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, NA), color = "black") ) ) +
  scale_y_continuous(limits=c(0,18), expand = c(0, 0))
ggsave(p1, filename = "ModelOutput/Fig_1.pdf",height = 5, width = 10)

# Figure Appendix: Bilateral Agreements and Issue Area over Time
p2 <- p %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c("pollution", "shipping", "fish", "construction")) %>%
  group_by(Period, variable) %>% 
  summarise(Frequency = sum(value)) %>% 
  ggplot() + 
  geom_bar(mapping = aes(x=Period, y=Frequency, fill = Period) ,width=.5, stat="identity", position = "dodge") +
  facet_wrap(~variable, ncol=4, scale = "free",
             labeller = labeller(variable = c("pollution" = "Pollution", "shipping" = "Shipping", "fish" = "Fish", "construction" = "Construction"))) +
  theme_SM() +
  ylab("New Agreements") +   xlab("") +
  scale_fill_manual(values =  c(rep("grey", 6), rep("black", 4))) +
  scale_alpha_manual( guide = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, NA), color = "black") ) ) +
  scale_y_continuous(limits=c(0,18), expand = c(0, 0))
ggsave(p2, filename = "ModelOutput/Fig_2.pdf",height = 5, width = 10)

############
# (5) setup data to model time-varying co-variates
############

# prepare the data for the analysis of time dependent covariates as suggested by Therneau et al. 2021 https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
dyadicdat2 <- dyadicdat %>% 
  group_by(cantons) %>% 
  mutate(tstart = sort(unique(dyadicdat$time))[1],
         tstop = time,
         tstart = ifelse(is.na(lag(time)), tstart, lag(time))) %>% 
  as.data.frame()

# disadvantage fo the time-varying co-variates: cannot model events that occur at the same time
# move tiny bit back as suggested by Therneau https://rdrr.io/cran/survival/f/inst/doc/survival.pdf p. 14
dyadicdat2$tstart[dyadicdat2$tstart == dyadicdat2$tstop] <- (dyadicdat2$tstart[dyadicdat2$tstart == dyadicdat2$tstop] -.01)

example_set_up <- dyadicdat2 %>%
  dplyr::select(ID_SM, year, cantons, tstart, tstop, commission, commission_cum, monitoring, monitoring_cum, conflict, conflict_cum) %>% 
  slice(1:30) %>% 
  dplyr::rename(id = ID_SM, 
                Year = year, 
                Dyad = cantons,
                Commission = commission,
                `Prior Commissions` = commission_cum,
                Monitoring = monitoring,
                `Prior Monitoring` = monitoring_cum,
                `Conflict Resolution` = conflict,
                `Prior Conflict Resolution` = conflict_cum) 
rownames(example_set_up) <- rep("", nrow(example_set_up))
  # rownames_to_column(var = "row") %>% 
  # dplyr::select(-example)
latex(example_set_up, file="", caption = "Data Setup", booktabs = T, landscape = F)

# number of cantons in the risk set
length(unique(dyadicdat2$cantons))
# number of cantons in the risk set that do not sign an agreement
length(unique(dyadicdat2$cantons[dyadicdat2$ID_SM == 0]))
# number of cantons in the risk set that sign an agreement
length(unique(dyadicdat2$cantons[dyadicdat2$ID_SM != 0]))
# number of agreements by these cantons
length(unique(dyadicdat2$ID_SM[dyadicdat2$ID_SM != 0 & dyadicdat2$year >= 1980]))
# number of agreements by these cantons
length(unique(dyadicdat2$ID_SM[dyadicdat2$ID_SM != 0 & dyadicdat2$year >= 1980 & dyadicdat2$commission == 1]))
length(unique(dyadicdat2$ID_SM[dyadicdat2$ID_SM != 0 & dyadicdat2$year >= 1980 & dyadicdat2$monitoring == 1]))
length(unique(dyadicdat2$ID_SM[dyadicdat2$ID_SM != 0 & dyadicdat2$year >= 1980 & dyadicdat2$conflict == 1]))

############
# (6) analysis of repeated events: main model + robustness checks
############

attach(dyadicdat2[dyadicdat2$year >= 1980 | dyadicdat2$year == 0,])

# Sqrt transformation of the legacy effects (all those with sqrt transformations)
rob_1.1 <- coxph(Surv(tstart, tstop, treaty_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_2.1 <- coxph(Surv(tstart, tstop, treaty_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.1 <- coxph(Surv(tstart, tstop, conflict) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_4.1 <- coxph(Surv(tstart, tstop, conflict) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.1 <- coxph(Surv(tstart, tstop, monitoring) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt, cluster = dyad_id)
rob_6.1 <- coxph(Surv(tstart, tstop, monitoring) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.1 <- coxph(Surv(tstart, tstop, commission) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt,  cluster = dyad_id)
rob_8.1 <- coxph(Surv(tstart, tstop, commission) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

# Robustness Checks 1: Cumulative number
cox_1 <-  coxph(Surv(tstart, tstop, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id, ties = "breslow")
cox_2 <-  coxph(Surv(tstart, tstop, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, ties = "breslow")
cox_3 <-  coxph(Surv(tstart, tstop, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id, ties = "breslow")
cox_4 <-  coxph(Surv(tstart, tstop, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id, ties = "breslow")
cox_5 <-  coxph(Surv(tstart, tstop, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id, ties = "breslow")
cox_6 <-  coxph(Surv(tstart, tstop, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, ties = "breslow")
cox_7 <-  coxph(Surv(tstart, tstop, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum,  cluster = dyad_id, ties = "breslow")
cox_8 <-  coxph(Surv(tstart, tstop, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, ties = "breslow")

# Main Model: >=1
rob_1.2 <- coxph(Surv(tstart, tstop, treaty_yes) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_2.2 <- coxph(Surv(tstart, tstop, treaty_yes) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.2 <- coxph(Surv(tstart, tstop, conflict) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_4.2 <- coxph(Surv(tstart, tstop, conflict) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.2 <- coxph(Surv(tstart, tstop, monitoring) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_6.2 <- coxph(Surv(tstart, tstop, monitoring) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.2 <- coxph(Surv(tstart, tstop, commission) ~ agr_one + conflict_one + monitoring_one + commission_one, cluster = dyad_id)
rob_8.2 <- coxph(Surv(tstart, tstop, commission) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

model_list_rob_1 <- list(rob_1.1, rob_2.1, rob_3.1, rob_4.1, rob_5.1, rob_6.1, rob_7.1, rob_8.1)
model_list_main <- list(cox_1, cox_2, cox_3, cox_4, cox_5, cox_6, cox_7, cox_8)
model_list_rob_2 <- list(rob_1.2, rob_2.2, rob_3.2, rob_4.2, rob_5.2, rob_6.2,rob_7.2, rob_8.2)

texreg::wordreg(model_list_rob_1, custom.model.names=c("1", "2", "3", "4", "5", "6", "7", "8"), custom.coef.names=c("Sqrt(Prior Agreements)", "Sqrt(Prior Conflict Resolution)", "Sqrt(Prior Monitoring)", "Sqrt(Prior Commission)", "Pollution", "Shipping", "Fish", "Construction", "Symmetry", "Bi-lingue"),
                title="Cox model with the square root of cumulative legacy effects for data after 1980", file = "ModelOutput/main.doc")
texreg::wordreg(model_list_rob_1, booktabs = T, scalebox = 0.9, custom.model.names=c("1", "2", "3", "4", "5", "6", "7", "8"), float.pos = "!htb", label = "table:main",
                custom.coef.names=c("Sqrt(Prior Agreements)", "Sqrt(Prior Conflict Resolution)", "Sqrt(Prior Monitoring)", "Sqrt(Prior Commission)", "Pollution", "Shipping", "Fish", "Construction", "Symmetry", "Bi-lingue"),
                caption="Cox model with the square root of cumulative legacy effects for data after 1980", file = "ModelOutput/square_root.doc")
texreg::wordreg(model_list_main, booktabs = T, scalebox = 0.9, custom.model.names=c("1", "2", "3", "4", "5", "6", "7", "8"), float.pos = "!htb", label = "table:rob1",
                custom.coef.names=c("Prior Agreements", "Prior Conflict Resolution", "Prior Monitoring", "Prior Commission", "Pollution", "Shipping", "Fish", "Construction", "Symmetry", "Bi-lingue"),
                caption="Cox model with cumulative legacy effects for data after 1980", file = "ModelOutput/cumulative.doc")
texreg::wordreg(model_list_rob_2, booktabs = T, scalebox = 0.9, custom.model.names=c("1", "2", "3", "4", "5", "6", "7", "8"), float.pos = "!htb", label = "table:rob2",
                custom.coef.names=c("Prior Agreements", "Prior Conflict Resolution", "Prior Monitoring", "Prior Commission", "Pollution", "Shipping", "Fish", "Construction", "Symmetry", "Bi-lingue"),
                caption="Cox model witht legacy effects $>=1$ for data after 1980", file = "ModelOutput/larger_or_equal_than_one.doc")


### check for the cox proportional hazards assumption, outliers, residuals and non-linearity for the main model
## Agreement
# test proportional-hazards (PH) assumption
(zph_agr <- cox.zph(rob_2.1))
# plot to investigate relationship with time
pdf("ModelOutput/prop_haz_agr.pdf")
ggcoxzph(zph_agr, ggtheme = theme_minimal(base_size = 6)) 
dev.off()
# outliers
dev.off()
pdf("ModelOutput/pol_bil_outliers_agr.pdf")
ggcoxdiagnostics(rob_2.1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# residuals
pdf("ModelOutput/pol_bil_residuals_agr.pdf")
ggcoxdiagnostics(rob_2.1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# non-linearity
pdf("ModelOutput/pol_bil_non_lin_agr.pdf")
ggcoxfunctional(Surv(tstart, tstop, treaty_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                cluster = dyad_id, data = dyadicdat2[dyadicdat2$year >= 1980 | dyadicdat2$year == 0,])
dev.off()

## i) Conflict
# test proportional-hazards (PH) assumption
(zph_con <- cox.zph(rob_4.1))
# plot to investigate relationship with time
pdf("ModelOutput/prop_haz_con.pdf")
ggcoxzph(zph_con, ggtheme = theme_minimal(base_size = 6))
dev.off()
# outliers
pdf("ModelOutput/pol_bil_outliers_con.pdf")
ggcoxdiagnostics(rob_4.1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# residuals
pdf("ModelOutput/pol_bil_residuals_con.pdf")
ggcoxdiagnostics(rob_4.1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# non-linearity
pdf("ModelOutput/pol_bil_non_lin_con.pdf")
ggcoxfunctional(Surv(tstart, tstop, conflict) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping  + construction + symmetry + bi_lingue,
                cluster = dyad_id, data = dyadicdat2[dyadicdat2$year >= 1980 | dyadicdat2$year == 0,])
dev.off()

## ii) Monitoring
# test proportional-hazards (PH) assumption
(zph_mon <- cox.zph(rob_6.1))
# plot to investigate relationship with time
pdf("ModelOutput/prop_haz_mon.pdf")
ggcoxzph(zph_mon, ggtheme = theme_minimal(base_size = 6))
dev.off()
# outliers
pdf("ModelOutput/pol_bil_outliers_mon.pdf")
ggcoxdiagnostics(rob_6.1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# residuals
pdf("ModelOutput/pol_bil_residuals_mon.pdf")
ggcoxdiagnostics(rob_6.1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# non-linearity
pdf("ModelOutput/pol_bil_non_lin_mon.pdf")
ggcoxfunctional(Surv(tstart, tstop, monitoring) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                cluster = dyad_id, data = dyadicdat2[dyadicdat2$year >= 1980 | dyadicdat2$year == 0,])
dev.off()

## iii) Commission
# test proportional-hazards (PH) assumption
(zph_com <- cox.zph(rob_8.1))
# plot to investigate relationship with time
pdf("ModelOutput/prop_haz_com.pdf")
ggcoxzph(zph_com, ggtheme = theme_minimal(base_size = 6))
dev.off()
# outliers
pdf("ModelOutput/pol_bil_outliers_com.pdf")
ggcoxdiagnostics(rob_8.1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# residuals
pdf("ModelOutput/pol_bil_residuals_com.pdf")
ggcoxdiagnostics(rob_8.1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
dev.off()
# non-linearity
pdf("ModelOutput/pol_bil_non_lin_com.pdf")
ggcoxfunctional(Surv(tstart, tstop, commission) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution	+ shipping +	fish + construction + symmetry + bi_lingue,
                cluster = dyad_id, data = dyadicdat2[dyadicdat2$year >= 1980 | dyadicdat2$year == 0,])
dev.off()

rob_2.1_tt <- coxph(Surv(tstart, tstop, treaty_yes) ~ agreement_sqrt + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution + shipping +	fish + tt(fish) + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
rob_4.1_tt <- coxph(Surv(tstart, tstop, conflict) ~ agreement_sqrt + tt(agreement_sqrt) + conflict_sqrt + monitoring_sqrt + commission_sqrt + pollution + shipping + fish + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
model_list_rob_1 <- list(rob_1.1, rob_2.1_tt, rob_3.1, rob_4.1_tt, rob_5.1, rob_6.1, rob_7.1, rob_8.1)

# Table 4 in the main text
stargazer(model_list_rob_1, type = "text", title = "Cox model with the square root of cumulative legacy effects for data after 1980. PH-Tests are conducted before time transformation", 
          dep.var.labels=c("Agreement", "Conflict", "Monitoring", "Commission"), column.sep.width = "-5pt",
          notes = "standard errors are clustered by dyad, P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001",
          notes.append = FALSE, align=TRUE, no.space = TRUE, digits = 2, stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list_rob_1, function(x) summary(x)$coefficients[, 4]), 
          star.char = c("*", "**", "***"), star.cutoffs = c(.05, .01, .001), style = "ajs", omit.stat = c("rsq", "max.rsq", "lr", "wald", "logrank", "ll"),
          add.lines = list(c("AIC   ", as.character(round(unlist(lapply(model_list_rob_1, function(x) AIC(x))), 3))),
                           c("BIC   ", as.character(round(unlist(lapply(model_list_rob_1, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list_rob_1, function(x) (x)$nevent)))),
                           c("Missings", as.character(c(sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$conflict)), sum(is.na(dyadicdat2$conflict)),
                                                        sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$commission)), sum(is.na(dyadicdat2$commission))))),
                           c("PH-Test", sapply(c(tail(cox.zph(rob_1.1)$table, 1)[3], tail(cox.zph(rob_2.1)$table, 1)[3], tail(cox.zph(rob_3.1)$table, 1)[3], tail(cox.zph(rob_4.1)$table, 1)[3],
                                                       tail(cox.zph(rob_5.1)$table, 1)[3], tail(cox.zph(rob_6.1)$table, 1)[3], tail(cox.zph(rob_7.1)$table, 1)[3], tail(cox.zph(rob_8.1)$table, 1)[3]
                                                       ), round, 3))
                           ),
          out = "ModelOutput/tables/main_table_square_root.html"
)

cox_2_tt <- coxph(Surv(tstart, tstop, treaty_yes) ~ agr_cum + tt(agr_cum) + conflict_cum + tt(conflict_cum) + monitoring_cum + tt(monitoring_cum) + commission_cum + tt(commission_cum) + pollution + shipping +	fish + tt(fish) + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
cox_3_tt <- coxph(Surv(tstart, tstop, conflict) ~ agr_cum + tt(agr_cum) + conflict_cum + tt(conflict_cum) + monitoring_cum + tt(monitoring_cum) + commission_cum + tt(commission_cum), cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
cox_7_tt <- coxph(Surv(tstart, tstop, commission) ~ agr_cum + conflict_cum + monitoring_cum + commission_cum + tt(conflict_cum), cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
cox_8_tt <- coxph(Surv(tstart, tstop, commission) ~ agr_cum + tt(agr_cum) + conflict_cum + tt(conflict_cum) + monitoring_cum + commission_cum + tt(commission_cum) + pollution + shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))

model_list_main <- list(cox_1, cox_2_tt, cox_3_tt, cox_4, cox_5, cox_6, cox_7_tt, cox_8_tt)

# Table 5 in the Appendix
stargazer(model_list_main, type = "latex", title = "Cox Proportional Hazards Model with cummulative legacy effects for prior agreements in the time-period from 1980 to 2020. PH-Tests are conducted before time transformation", 
          dep.var.labels=c("Agreement", "Conflict", "Monitoring", "Commission"), column.sep.width = "-5pt",
          notes = "standard errors are clustered by dyad, P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001, PH-Tests are conducted before time transformation",
          notes.append = FALSE, align=TRUE, no.space = TRUE, digits = 2, stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list_main, function(x) summary(x)$coefficients[, 4]),
          star.char = c("*", "**", "***"), star.cutoffs = c(.05, .01, .001), style = "ajs", omit.stat = c("rsq", "max.rsq", "lr", "wald", "logrank", "ll"),
          add.lines = list(c("AIC   ", as.character(round(unlist(lapply(model_list_main, function(x) AIC(x))), 3))),
                           c("BIC   ", as.character(round(unlist(lapply(model_list_main, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list_main, function(x) (x)$nevent)))),
                           c("Missings", as.character(c(sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$conflict)), sum(is.na(dyadicdat2$conflict)),
                                                        sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$commission)), sum(is.na(dyadicdat2$commission))))),
                           c("PH-Test", sapply(c(tail(cox.zph(cox_1)$table, 1)[3], tail(cox.zph(cox_2)$table, 1)[3], tail(cox.zph(cox_3)$table, 1)[3], tail(cox.zph(cox_4)$table, 1)[3],
                                                 tail(cox.zph(cox_5)$table, 1)[3], tail(cox.zph(cox_6)$table, 1)[3], tail(cox.zph(cox_7)$table, 1)[3], tail(cox.zph(cox_8)$table, 1)[3]
                           ), round, 3))
          ),
          out = "ModelOutput/tables/cumulative.html"
)

lapply(model_list_rob_2, cox.zph)

rob_2.2_tt <- coxph(Surv(tstart, tstop, treaty_yes) ~ agr_one + conflict_one + monitoring_one + commission_one + pollution + shipping +	fish + tt(fish) + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
rob_3.2_tt <- coxph(Surv(tstart, tstop, conflict) ~ agr_one + tt(agr_one) + conflict_one + monitoring_one + commission_one, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
rob_4.2_tt <- coxph(Surv(tstart, tstop, conflict) ~ agr_one + tt(agr_one) + conflict_one + monitoring_one + tt(monitoring_one) + commission_one + pollution + shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
rob_7.2_tt <- coxph(Surv(tstart, tstop, commission) ~ agr_one + tt(agr_one) + conflict_one + tt(conflict_one) + monitoring_one + tt(monitoring_one) + commission_one + tt(commission_one), cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))
rob_8.2_tt <- coxph(Surv(tstart, tstop, commission) ~ agr_one + conflict_one + monitoring_one + tt(monitoring_one) + commission_one + pollution + shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id, tt = function(x, t, ...) x * log(t+20))

model_list_rob_2 <- list(rob_1.2, rob_2.2_tt, rob_3.2_tt, rob_4.2_tt, rob_5.2, rob_6.2, rob_7.2_tt, rob_8.2_tt)

# Table 6 in the Appendix
stargazer(model_list_rob_2, type = "latex", title = "Cox model with legacy effects $>=1$ legacy effects for prior agreements in the time-period from 1980 to 2020. PH-Tests are conducted before time transformation",
          dep.var.labels=c("Agreement", "Conflict", "Monitoring", "Commission"), column.sep.width = "-5pt",
          notes = "standard errors are clustered by dyad, P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001",
          notes.append = FALSE, align=TRUE, no.space = TRUE, digits = 2, stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list_rob_2, function(x) summary(x)$coefficients[, 4]),
          star.char = c("*", "**", "***"), star.cutoffs = c(.05, .01, .001), style = "ajs", omit.stat = c("rsq", "max.rsq", "lr", "wald", "logrank", "ll"),
          add.lines = list(c("AIC   ", as.character(round(unlist(lapply(model_list_rob_2, function(x) AIC(x))), 3))),
                           c("BIC   ", as.character(round(unlist(lapply(model_list_rob_2, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list_rob_2, function(x) (x)$nevent)))),
                           c("Missings", as.character(c(sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$conflict)), sum(is.na(dyadicdat2$conflict)),
                                                        sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$commission)), sum(is.na(dyadicdat2$commission))))),
                           c("PH-Test", sapply(c(tail(cox.zph(rob_1.2)$table, 1)[3], tail(cox.zph(rob_2.2)$table, 1)[3], tail(cox.zph(rob_3.2)$table, 1)[3], tail(cox.zph(rob_4.2)$table, 1)[3],
                                                 tail(cox.zph(rob_5.2)$table, 1)[3], tail(cox.zph(rob_6.2)$table, 1)[3], tail(cox.zph(rob_7.2)$table, 1)[3], tail(cox.zph(rob_8.2)$table, 1)[3]
                           ), round, 3))
          ),
          out = "ModelOutput/tables/larger_or_equal_than_one.html"
)


detach(dyadicdat2[dyadicdat2$year >= 1980 | dyadicdat2$year == 0,])

attach(dyadicdat2[dyadicdat2$year >= 1945 | dyadicdat2$year == 0,])

# Robustness Checks 1: Cumulative number
rob_1.3 <-  coxph(Surv(tstart, tstop, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_2.3 <-  coxph(Surv(tstart, tstop, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_3.3 <-  coxph(Surv(tstart, tstop, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_4.3 <-  coxph(Surv(tstart, tstop, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.3 <-  coxph(Surv(tstart, tstop, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_6.3 <-  coxph(Surv(tstart, tstop, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_7.3 <-  coxph(Surv(tstart, tstop, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum,  cluster = dyad_id)
rob_8.3 <-  coxph(Surv(tstart, tstop, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue, cluster = dyad_id)

model_list_rob_3 <- list(rob_1.3, rob_2.3, rob_3.3, rob_4.3, rob_5.3, rob_6.3, rob_7.3, rob_8.3)
lapply(model_list_rob_3, cox.zph)

# Robustness Checks 1: Cumulative number
rob_1.3_tt <-  coxph(Surv(tstart, tstop, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + tt(commission_cum), cluster = dyad_id)
rob_2.3_tt <-  coxph(Surv(tstart, tstop, treaty_yes) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + tt(bi_lingue), cluster = dyad_id)
rob_3.3_tt <-  coxph(Surv(tstart, tstop, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_4.3_tt <-  coxph(Surv(tstart, tstop, conflict) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution + tt(pollution) + shipping  + construction + symmetry + bi_lingue, cluster = dyad_id)
rob_5.3_tt <-  coxph(Surv(tstart, tstop, monitoring) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum, cluster = dyad_id)
rob_6.3_tt <-  coxph(Surv(tstart, tstop, monitoring) ~  agr_cum + conflict_cum + tt(conflict_cum) + monitoring_cum + commission_cum + tt(commission_cum) + pollution	+ shipping +	fish + construction + symmetry + bi_lingue + tt(bi_lingue), cluster = dyad_id)
rob_7.3_tt <-  coxph(Surv(tstart, tstop, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum,  cluster = dyad_id)
rob_8.3_tt <-  coxph(Surv(tstart, tstop, commission) ~  agr_cum + conflict_cum + monitoring_cum + commission_cum + pollution + tt(pollution)	+ shipping + tt(shipping) +	fish + construction + symmetry + bi_lingue + tt(bi_lingue), cluster = dyad_id)

model_list_rob_3_tt <- list(rob_1.3_tt, rob_2.3_tt, rob_3.3_tt, rob_4.3_tt, rob_5.3_tt, rob_6.3_tt, rob_7.3_tt, rob_8.3_tt)

# Table 7 in the Appendix
stargazer(model_list_rob_3_tt, type = "latex", title = "Cox Proportional Hazards Model with cumulative legacy effects for prior agreements in the time-period from 1945 to 2020. PH-Tests are conducted before time transformation",
          dep.var.labels=c("Agreement", "Conflict", "Monitoring", "Commission"), column.sep.width = "-5pt",
          notes = "standard errors are clustered by dyad, P $<$ 0.1, $^{*}$ P $<$ 0.05, $^{**}$ P $<$ 0.01, $^{*}$ P $<$ 0.001",
          notes.append = FALSE, align=TRUE, no.space = TRUE, digits = 2, stargazer_stat_code_list = list(c("aic", "bic")),
          se = lapply(model_list_rob_3_tt, function(x) summary(x)$coefficients[, 4]),
          star.char = c("*", "**", "***"), star.cutoffs = c(.05, .01, .001), style = "ajs", omit.stat = c("rsq", "max.rsq", "lr", "wald", "logrank", "ll"),
          add.lines = list(c("AIC   ", as.character(round(unlist(lapply(model_list_rob_3_tt, function(x) AIC(x))), 3))),
                           c("BIC   ", as.character(round(unlist(lapply(model_list_rob_3_tt, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list_rob_3_tt, function(x) (x)$nevent)))),
                           c("Missings", as.character(c(sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$treaty_yes)), sum(is.na(dyadicdat2$conflict)), sum(is.na(dyadicdat2$conflict)),
                                                        sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$monitoring)), sum(is.na(dyadicdat2$commission)), sum(is.na(dyadicdat2$commission))))),
                           c("PH-Test", sapply(c(tail(cox.zph(rob_1.3)$table, 1)[3], tail(cox.zph(rob_2.3)$table, 1)[3], tail(cox.zph(rob_3.3)$table, 1)[3], tail(cox.zph(rob_4.3)$table, 1)[3],
                                                 tail(cox.zph(rob_5.3)$table, 1)[3], tail(cox.zph(rob_6.3)$table, 1)[3], tail(cox.zph(rob_7.3)$table, 1)[3], tail(cox.zph(rob_8.3)$table, 1)[3]
                           ), round, 3))
          )
)

attach(dyadicdat2[dyadicdat2$year >= 1945 | dyadicdat2$year == 0,])


############
# (6) Plot Results
############
## Plot Kaplan-Meier survival curves:
agr_df_1 <- with(dyadicdat,
                 data.frame(agreement_sqrt = rep(mean(agreement_sqrt, na.rm = TRUE), 2),
                            commission_sqrt = c(0, 1),
                            monitoring_sqrt = rep(mean(monitoring_sqrt, na.rm = TRUE), 2),
                            conflict_sqrt = rep(mean(conflict_sqrt, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            shipping = rep(mean(shipping, na.rm = TRUE), 2),
                            pollution = rep(mean(pollution, na.rm = TRUE), 2),
                            fish = rep(mean(fish, na.rm = TRUE), 2),
                            construction = rep(mean(construction, na.rm = TRUE), 2),
                            symmetry = rep(mean(symmetry, na.rm = TRUE), 2),
                            bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2)
                 )
)

fit_agr_1 <- survfit(rob_4.1, newdata = agr_df_1)
p3 <- ggsurvplot(fit_agr_1,
                 conf.int = TRUE,
                 legend.labs=c("0 Commissions", "1 Commission"),
                 xlab = "Time in days",
                 ylab = "Probability without Additional Conflict Resolution",
                 break.time.by = 4000,
                 ggtheme = theme_bw(base_size = 14), 
                 data = agr_df_1)
p3 # stored manually
