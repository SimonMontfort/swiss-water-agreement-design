### Swiss Water Coooperation
## Authors: Simon Montfort, Manuel Fischer, James Hollway, Nicolas Jager
# > R.version
# _                           
# platform       x86_64-apple-darwin17.0     
# arch           x86_64                      
# os             darwin17.0                  
# system         x86_64, darwin17.0          
# status                                     
# major          4                           
# minor          0.3                         
# year           2020                        
# month          10                          
# day            10                          
# svn rev        79318                       
# language       R                           
# version.string R version 4.0.3 (2020-10-10)
# nickname       Bunny-Wunnies Freak Out     
######################################
### Load packages
######################################

rm(list = ls())

library("foreign")
library("texreg")
library("lme4")
library("xtable")
library("dplyr")
library("gridExtra")
library("effects")
library("ggplot2")
library("rio")
library("summarytools")
library("survival")
library("survminer")
library("dplyr")
library("gdata")
library("igraph")
library("ggrepel")

######################################
# description of code
######################################

# The code proceeds in steps:
# (1) loading treaty data and subsetting to bilaterals
# (2) definition of risk set to contiguity over land
# (3) variables and covariates 
# (4) descriptives
# (5) generation of data sets for different DVs
# (6) analysis
# (7) assumptions and rectifications

######################################
# (1) treaty design data
######################################

# change directory here to replicate 
setwd("/Users/simon/Documents/repo/swiss-water-agreement-design/ModelInput")

load("dyadicdat.RData")

# subset to those issues included
dyadicdat <- dyadicdat %>%
  filter(pollution == 1 | shipping == 1 | fish ==  1 | power ==  1)

# subset to bilateral river agreements
dyadicdat <- dyadicdat %>% filter(bilateral == 1
                                  & river == 1
                                  )

# create DV for agreements
dyadicdat <- dyadicdat %>%
  mutate(treaty_yes = ifelse(!is.na(ID_SM), 1, 0))

##########
# (2) risk set
##########

# risk set rivers
load("cantons_border.RData")
cont_adj <- cont_adj[order(rownames(cont_adj)), order(colnames(cont_adj))]
cont_adj[cont_adj >0] <- 1
dyads <- as.data.frame(get.edgelist(graph_from_adjacency_matrix(cont_adj, mode = "undirected")))

# # risk set lakes
# load("risk_set_lakes.RData")
# risk_set_lakes[risk_set_lakes >0] <- 1
# dyads_2 <- as.data.frame(get.edgelist(graph_from_adjacency_matrix(risk_set_lakes, mode = "undirected")))
# dyads <- rbind.data.frame(dyads, dyads_2)

# make sure the dyads are ordered by name, these are the cantons in the risk set. Only contiguous cantons can sign agreements
colnames(dyads) <- c("canton1", "canton2")
dyads <- dyads %>%
  mutate(canton1 = as.character(canton1),
         canton2 = as.character(canton2),
         cantons = ifelse(canton1 < canton2, paste(canton1, canton2), paste(canton2, canton1))) %>%
  distinct(cantons, .keep_all = T)

############
# (3) variables and covariates 
############

## I. upstream/downstream operationalisation 

load("cantons_border.RData")
# make sure matrix is ordered
cont_adj <- cont_adj[order(rownames(cont_adj)), order(colnames(cont_adj))]
# merge to dyadic data
for (i in 1:nrow(dyads)){
  dyads$border_length[i] <- cont_adj[dyads$canton1[i] == rownames(cont_adj), dyads$canton2[i] == colnames(cont_adj)]
}

# load border length data and prepare df for merging to main dyadic data
load("bord_to_size.RData")
border_to_size <- as.data.frame(matrix(NA, 26, 2))
colnames(border_to_size) <- c("canton", "border_to_size")
border_to_size$canton <- rownames(cont_adj)
border_to_size$border_to_size <- bord_to_size
rm(bord_to_size)

left_join(dyadicdat, border_to_size, by = c("canton1" = "canton"))
left_join(dyadicdat, border_to_size, by = c("canton2" = "canton"))

# dyadic data: also make sure that it is orderered by name
dyadicdat <- dyadicdat %>%
  mutate(canton1 = as.character(canton1),
         canton2 = as.character(canton2),
         cantons = ifelse(canton1 < canton2, paste(canton1, canton2), paste(canton2, canton1)))

# join datasets
dyadicdat <- left_join(dyads, dyadicdat %>% dplyr::select(-canton1, -canton2), by = c("cantons" = "cantons"))

dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(id = cur_group_id()) %>% as.data.frame() %>% dplyr::select(id, cantons, time, enforcement, monitoring, conflict)

dyadicdat$time[is.na(dyadicdat$time)] <- as.numeric(as.Date("2020-01-01")) # in numeric format, 18262 is the end of the study period, this is lower than 

# bilateral agreements are only contiguous, this is just to double check that they. should be empty
dyadicdat$cantons[!dyadicdat$cantons %in% dyads$cantons] # perfect, it is empty

# replace NAs with zero because these are the observations that did not experience the event
dyadicdat[is.na(dyadicdat)] <- 0

sum(dyadicdat$enforcement)
sum(dyadicdat$monitoring)
sum(dyadicdat$conflict)
sum(dyadicdat$commission)

length(unique(dyadicdat$cantons))

## II. upstream/downstream operationalisation 
load("US_DS_matrix_polygons_within_1.RData")
drainage <- cbind(dyadicdat$canton1, dyadicdat$canton1)
colnames(drainage) <- c("canton1", "canton2")
colnames(mat)[colnames(mat) == "NI"] <- "NW"
rownames(mat)[rownames(mat) == "NI"] <- "NW"
for (i in 1:nrow(dyadicdat)){
  dyadicdat$drain_1_to_2[i] <- mat[dyadicdat$canton1[i] == rownames(mat), dyadicdat$canton2[i] == colnames(mat)]
  dyadicdat$drain_2_to_1[i] <- mat[dyadicdat$canton2[i] == rownames(mat), dyadicdat$canton1[i] == colnames(mat)]
}

## Water Salience, Symmetry and Interdependence based on trade interdependence by Crescenzi, p.30
# total water
total_wat <- as.data.frame(rowSums(mat) + colSums(mat)) 
total_wat$canton <- rownames(total_wat)
colnames(total_wat) <- c("total_water", "canton")
dyadicdat <- left_join(dyadicdat, total_wat, by = c("canton1" = "canton"))
dyadicdat <- left_join(dyadicdat, total_wat, by = c("canton2" = "canton"))

# water_share i: drainage area from i to j divided by total drainage of i
dyadicdat$water_share.x <- (dyadicdat$drain_1_to_2+ dyadicdat$drain_2_to_1)/dyadicdat$total_water.x
# water_share j: drainage area from j to i divided by total drainage of j
dyadicdat$water_share.y <- (dyadicdat$drain_1_to_2+ dyadicdat$drain_2_to_1)/dyadicdat$total_water.y
# salience  ij: square root of water_share i time water_share j
dyadicdat$salience <- sqrt(dyadicdat$water_share.x * dyadicdat$water_share.y)
# symmetry  ij: one minus the absolute difference in water_share i and water_share j
dyadicdat$symmetry_1 <- 1 - abs(dyadicdat$water_share.x - dyadicdat$water_share.y)

# our symmetry operationalisation (Manuel's suggestion) 
dyadicdat$total_wat <- (dyadicdat$drain_1_to_2 + dyadicdat$drain_2_to_1)
dyadicdat$share_drain_1 <- dyadicdat$drain_1_to_2/dyadicdat$total_wat
dyadicdat$share_drain_2 <- dyadicdat$drain_2_to_1/dyadicdat$total_wat
dyadicdat$symmetry <- 1 - abs(dyadicdat$share_drain_1 - dyadicdat$share_drain_2)

# inspect the data
tab <- cbind(round(dyadicdat$drain_1_to_2, 3), round(dyadicdat$drain_2_to_1, 3), round(dyadicdat$water_share.x, 3), round(dyadicdat$water_share.y, 3), round(dyadicdat$salience , 3), round(dyadicdat$symmetry, 3), round(dyadicdat$symmetry_1, 3), dyadicdat$cantons, dyadicdat$commission, dyadicdat$monitoring, dyadicdat$conflict)
colnames(tab) <- c("drain_1_to_2", "drain_2_to_1", "water_share.x","water_share.y", "symmetry", "symmetry_1", "interdpendence", "cantons", "commission", "monitoring", "conflict")
tab # This shows that the symmetry score by Manuel produces more sensible numbers

# adapt risk set becuse cantons cannot cooperate when no water flows from one to the other
dyadicdat <- dyadicdat[!is.nan(dyadicdat$symmetry),]

# III. shared lakes
load("shared_lakes_adj.RData")
shared_lakes_adj <- shared_lakes_adj[order(rownames(shared_lakes_adj)), order(colnames(shared_lakes_adj))]
shared_lakes_adj
for (i in 1:nrow(dyadicdat)){
  dyadicdat$lakes[i] <- shared_lakes_adj[dyadicdat$canton1[i] == rownames(shared_lakes_adj), dyadicdat$canton2[i] == colnames(shared_lakes_adj)]
}

# IV. border length el 
load("cantons_border_el.RData")
load("cantons_size.RData")
cantons_border_el <- as.data.frame(cantons_border_el)
cantons_size <- as.data.frame(cantons_size)
cantons_size <- cantons_size[,1:2] # drop spatial column, not needed
cantons_size$cantons_size <- as.numeric(cantons_size$cantons_size)/1000 # in meters

# order cantons names in the dyad to be able to correctly merge them below
cantons_border_el <- cantons_border_el %>%
  mutate(canton1 = as.character(canton1),
         canton2 = as.character(canton2),
         cantons = ifelse(canton1 < canton2, paste(canton1, canton2), paste(canton2, canton1))) %>%
  arrange(cantons)

# merge dyadic boundary length
dyadicdat <- left_join(dyadicdat, cantons_border_el %>% dplyr::select(total_boundary_length, cantons), by = c("cantons"))
# merge the size to each canton in the dyad
dyadicdat <- left_join(dyadicdat, cantons_size %>% as.data.frame() %>% dplyr::select(iso, cantons_size), by = c("canton1" = "iso"))
dyadicdat <- left_join(dyadicdat, cantons_size %>% as.data.frame() %>% dplyr::select(iso, cantons_size), by = c("canton2" = "iso"))

## Variable construction
# border length relative to size
dyadicdat$bord_length_to_size <- as.numeric(dyadicdat$total_boundary_length / (dyadicdat$cantons_size.x + dyadicdat$cantons_size.y))
# calculate the size size of both cantons in the dyad
dyadicdat$total_size <- as.numeric(dyadicdat$cantons_size.x + dyadicdat$cantons_size.y)/1000
# calculate the size size of both cantons in the dyad
dyadicdat$size <- as.numeric(abs(dyadicdat$cantons_size.x - dyadicdat$cantons_size.y )/(dyadicdat$cantons_size.x + dyadicdat$cantons_size.y))

dyadicdat$total_boundary_length <- as.numeric(dyadicdat$total_boundary_length)
dyadicdat$cantons_size.x <- as.numeric(dyadicdat$cantons_size.x)
dyadicdat$cantons_size.y <- as.numeric(dyadicdat$cantons_size.y)
dyadicdat$bord_length_to_size <- as.numeric(dyadicdat$bord_length_to_size)

# V. sequencing
dyadicdat <- dyadicdat %>%
  group_by(cantons) %>%
  arrange(cantons, time) %>%
  mutate(agr_cum = row_number() -1,
         enforcement_cum = cumsum(enforcement) -1,
         monitoring_cum = cumsum(monitoring) -1,
         conflict_cum = cumsum(conflict) -1,
         commission_cum = cumsum(commission) -1) %>%
  mutate(enforcement_cum = ifelse(enforcement_cum == -1, enforcement_cum +1, enforcement_cum),
         monitoring_cum = ifelse( monitoring_cum == -1, monitoring_cum +1, monitoring_cum),
         conflict_cum = ifelse(conflict_cum == -1, conflict_cum +1, conflict_cum),
         commission_cum = ifelse(commission_cum == -1, commission_cum +1, commission_cum))

dyadicdat %>%
  arrange(cantons, time) %>%
  dplyr::select(cantons, time, enforcement, monitoring, conflict, enforcement_cum, monitoring_cum, conflict_cum, agr_cum) %>%
  as.data.frame()

# VI. language
covariates <- read.table("Kantonsdaten von Manuel_2.csv", sep = ",", header = T)

dyadicdat <- left_join(dyadicdat, covariates %>% dplyr::select(Kantone, German), by = c("canton1" = "Kantone"))
dyadicdat <- left_join(dyadicdat, covariates %>% dplyr::select(Kantone, German), by = c("canton2" = "Kantone"))
dyadicdat <- left_join(dyadicdat, covariates %>% dplyr::select(Kantone, French), by = c("canton1" = "Kantone"))
dyadicdat <- left_join(dyadicdat, covariates %>% dplyr::select(Kantone, French), by = c("canton2" = "Kantone"))
dyadicdat <- left_join(dyadicdat, covariates %>% dplyr::select(Kantone, Italian), by = c("canton1" = "Kantone"))
dyadicdat <- left_join(dyadicdat, covariates %>% dplyr::select(Kantone, Italian), by = c("canton2" = "Kantone"))

dyadicdat$german <- ifelse(dyadicdat$German.x + dyadicdat$German.y == 2 & dyadicdat$French.x + dyadicdat$French.y < 1 & dyadicdat$Italian.x + dyadicdat$Italian.y < 1, 1, 0)
dyadicdat$french <- ifelse(dyadicdat$French.x + dyadicdat$French.y == 2 & dyadicdat$German.y + dyadicdat$German.y < 1 & dyadicdat$Italian.x + dyadicdat$Italian.y < 1, 1, 0)
dyadicdat$italian <- ifelse(dyadicdat$Italian.x + dyadicdat$Italian.y == 2  & dyadicdat$French.x + dyadicdat$French.y < 1 & dyadicdat$German.x + dyadicdat$German.y < 1, 1, 0)

dyadicdat$bi_lingue.x <- ifelse(dyadicdat$Italian.x + dyadicdat$French.x + dyadicdat$German.x > 1, 1, 0)
dyadicdat$bi_lingue.y <- ifelse(dyadicdat$Italian.y + dyadicdat$French.y + dyadicdat$German.y > 1, 1, 0)

dyadicdat$bi_lingue <- ifelse(dyadicdat$bi_lingue.x + dyadicdat$bi_lingue.y >= 2, 1, 0)
dyadicdat$bi_lingue <- ifelse(dyadicdat$bi_lingue.x + dyadicdat$bi_lingue.y >= 1, 1, 0)

# add dyad id
dyadicdat <- dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(dyad_id = cur_group_id())

## look at structure of the data to check that cummulative sums have been calculated correctly
dyadicdat %>% group_by(cantons) %>%  arrange(cantons, time) %>% mutate(id = cur_group_id()) %>% 
  dplyr::select(id, cantons, time, enforcement, monitoring, conflict, enforcement_cum, monitoring_cum, conflict_cum, agr_cum) %>% as.data.frame()

#################
# (4) descriptives
#################
par(mfrow = c(1,1), mar = c(5, 6, 2, 2 ))

## sequentiality
library(TraMineR)
dyadicdat %>%
  dplyr::select(cantons, time, monitoring, conflict, commission) %>%
  filter(sum( monitoring, conflict, commission) != 0) %>%
  reshape2::melt(., id.vars = c("cantons", "time"), measure.vars = c("commission", "monitoring", "conflict")) %>%
  filter(value == 1) %>%
  dplyr::select(-value) %>%
  summarise(id = cantons,
            time = as.Date(time, origin = "1970-01-01"),
            event = variable ) %>%
  arrange(id) %>%
  mutate(id = as.factor(id),
         time =as.numeric(time),
         event = as.character(event)) %>%
  seqecreate(., id = id, timestamp = time, event = event) %>%
  seqpcplot(.,
            # alphabet = alph,
            filter = list(type = "function",
                          value = "cumfreq",
                          level = .6),
            order.align = "last",
            # ltype = "non-embeddable",
            cex = 1.5, lwd = .4,
            lcourse = "downwards")

## correlation between mechanisms
chisq.test(table(dyadicdat$commission, dyadicdat$monitoring), correct=FALSE)
chisq.test(table(dyadicdat$commission, dyadicdat$conflict), correct=FALSE)
chisq.test(table(dyadicdat$conflict, dyadicdat$monitoring), correct=FALSE)

## correlation between DVs and main EVs
cor(dyadicdat$monitoring, dyadicdat$symmetry)
cor(dyadicdat$conflict, dyadicdat$symmetry)
cor(dyadicdat$commission, dyadicdat$symmetry)
cor(dyadicdat$treaty_yes, dyadicdat$symmetry)

## simple descriptives of co-appearance
table(dyadicdat$pollution,	dyadicdat$shipping, dyadicdat$fish, dyadicdat$power, dyadicdat$commission, dyadicdat$monitoring, dyadicdat$conflict)
table(dyadicdat$pollution, dyadicdat$commission)
table(dyadicdat$pollution, dyadicdat$monitoring)
table(dyadicdat$pollution, dyadicdat$enforcement)
table(dyadicdat$pollution, dyadicdat$conflict)

table(dyadicdat$shipping, dyadicdat$commission)
table(dyadicdat$shipping, dyadicdat$monitoring)
table(dyadicdat$shipping, dyadicdat$enforcement)
table(dyadicdat$shipping, dyadicdat$conflict)

table(dyadicdat$fish, dyadicdat$commission)
table(dyadicdat$fish, dyadicdat$monitoring)
table(dyadicdat$fish, dyadicdat$enforcement)
table(dyadicdat$fish, dyadicdat$conflict)

table(dyadicdat$power, dyadicdat$commission)
table(dyadicdat$power, dyadicdat$monitoring)
table(dyadicdat$power, dyadicdat$enforcement)
table(dyadicdat$power, dyadicdat$conflict)

# how many agreeement per issue?
sum(dyadicdat$pollution)
sum(dyadicdat$shipping)
sum(dyadicdat$fish)
sum(dyadicdat$power)


### plot 1
plot1 <- dyadicdat %>%
  group_by(cantons) %>%
  summarise(no_of_agr = sum(!is.na(title)),
            no_of_mechs = sum(commission, monitoring, conflict)/no_of_agr,
            symmetry = symmetry,
            bi_lingue = bi_lingue) %>%
  distinct() %>%
  ungroup() 

ggplot(plot1, aes(symmetry, no_of_agr, label = cantons)) +
  geom_point(aes(shape = factor(bi_lingue))) + 
  geom_text_repel(data          = subset(plot1, no_of_agr > 1),
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.5,
                  size          = 2,
                  segment.color = "grey10",
                  segment.alpha = .5,
                  # direction     = "y"
                  ) +
  stat_smooth(method = "lm",
              col = "black",
              se = FALSE,
              size = .5) +
  theme_classic(base_size = 16) +
  labs(y = "Number of Agreements", x = "Symmetry", shape = "Bi-Lingue") 

## plot 2
dyadicdat %>%
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
         Period = replace(Period, year>= 2010 & year<= 2020, "2010-2020")) %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c("agreement_yes", "commission", "monitoring", "conflict")) %>%
  group_by(Period, variable) %>%
  summarise(Frequency = sum(value)) %>%
  ggplot() + 
  geom_bar(aes(y = Frequency, x = Period, width=.5, color = "black"), color = NA,
           stat="identity", position = "dodge") +
  facet_wrap(~variable) +
  theme_classic() +
  # coord_flip() +
  theme(
    axis.text=element_text(size=10, angle = 90),
    panel.grid.major.y = element_line(colour = "grey80"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1, "lines")) 
dev.off()

# we can probably depete this
dyadicdat %>%
  group_by(ID_SM) %>%
  slice(1) %>%
  as_tibble() %>%
  # filter(year >= 1980) %>%
  mutate(date = as.Date(time, origin = "1970-01-01"),
         year = as.numeric(format(date,'%Y')),
         agreement_yes = ifelse(!is.na(ID_SM), 1, 0),
         Period = NA,
         Period = replace(Period, year < 1940, "<1940"), 
         Period = replace(Period, year>= 1930 & year< 1940, "1930-1939"), 
         Period = replace(Period, year>= 1940 & year< 1950, "1940-1949"), 
         Period = replace(Period, year>= 1950 & year< 1960, "1950-1959"), 
         Period = replace(Period, year>= 1960 & year< 1970, "1960-1969"), 
         Period = replace(Period, year>= 1970 & year< 1980, "1970-1979"), 
         Period = replace(Period, year>= 1980 & year< 1990, "1980-1989"), 
         Period = replace(Period, year>= 1990 & year< 2000, "1990-1999"),
         Period = replace(Period, year>= 2000 & year< 2010, "2000-2009"),
         Period = replace(Period, year>= 2010 & year<= 2020, "2010-2020")) %>%
  reshape2::melt(., id.vars = "Period", measure.vars = c( "commission", "monitoring", "conflict")) %>%
  group_by(Period, variable) %>%
  summarise(Frequency = sum(value)) %>%
  arrange(variable) %>%
  filter(Frequency != 0) %>%
  ggplot(aes(y = factor(Frequency), x = Period, width=.5, fill = factor(variable))) + 
  geom_bar(aes(),
           stat="identity", position = "dodge") +
  # facet_wrap(~variable) +
  scale_fill_grey(start=0.8, end=0.2) +
  theme_bw() +
  # coord_flip() +
  theme(
    axis.text=element_text(size=10, angle = 90),
    panel.grid.major.y = element_line(colour = "grey80") ,
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1, "lines"))
dev.off()

## standard summary table
library(summarytools)
library(Hmisc) 
sum_dat <- dyadicdat %>%
  dplyr::select(treaty_yes, commission, monitoring, conflict, commission_cum, monitoring_cum, conflict_cum, salience, symmetry, border_length, total_size, bi_lingue) %>%
  ungroup()
sum_tab <- summarytools::descr(sum_dat, transpose = T, stats = c("mean", "sd", "min", "q1", "med", "q3", "max"))
sum_tab <- round(sum_tab, 2)
# rownames(sum_tab) <- c("Bi-lingue",
#                        "Border Length",
#                        "Commission",
#                        "Conflict",
#                        "Monitoring",
#                        "Salience",
#                        "Symmetry",
#                        "Dyad Size" )
latex(sum_tab, file="")


##############
# (5) generation of data sets for different DVs
##############
# first transform time variable so that it start with 0 to avoid problems with logarithmic transformation later when analysing assumptions and remedy violations
# dyadicdat$time <- dyadicdat$time + (min(dyadicdat$time)*(-1))

# # total time past in the study:
# start <- as.numeric(as.Date("1814-01-01")) # start
# end <- as.numeric(as.Date("2020-01-01")) # end

# it is only possible to run repeated event cox models whith data set up in a counting process (as opposed to multiple events as we currently have it) when all tstart are smaller
# than tstop. This is not given in this data because some dyads sign two agreements at the same day.
dyadicdat$tstart < dyadicdat$tstop

# dyadicdat %>%
#   mutate(tstart = start,
#          tstop = end) %>%
#   arrange(cantons, time) %>%
#   ungroup() %>%
#   mutate(tstop = ifelse(lag(cantons) != cantons & !is.na(lag(cantons)), time, tstop),
#          tstop = ifelse(lag(cantons) == cantons  & !is.na(lag(cantons)), time, tstop),
#          tstart = ifelse(lag(cantons) == cantons  & !is.na(lag(cantons)), lag(tstop), tstart),
#          tstart = ifelse(lag(cantons) == cantons  & !is.na(lag(cantons)), lag(tstop), tstart),
#          # tstart = ifelse(time != tstart, time, tstop)
#          ) %>%
#   dplyr::select(cantons, time, tstart, tstop, ID_SM) %>%
#   as.data.frame()

dyadicdat %>%
  select(cantons, time, ID_SM )

## I. Monitoring:
# (1) to subset to monitoring, we need to replace the value of the event when the agreement was signed as right censored, i.e. the end of the study time
mon_dat <- dyadicdat %>%
  arrange(dyad_id, time) %>%
  group_by(dyad_id) %>%
  mutate(sum = sum(monitoring)) %>%
  ungroup() %>%
  mutate(time = ifelse(sum == 0, as.numeric(as.Date("2020-01-01")), time))

# (2) keep in the dataset only those for which monitoring was not present at the end of the study period and those for which monitoring was an event that happened before the
# end of the study period:
mon_dat <- mon_dat[mon_dat$monitoring == 0 & mon_dat$time == 18262 | # monitoring was not present at the end of the study period -> right censored
                   mon_dat$monitoring == 1 & mon_dat$time != 18262,] # event observed
# now remove all the observations which did not have monitoring present and but another mechanism and are therefore still present in the data. These need be excluded as only
# right censored data have the end of the study time which is as.numeric(as.Date("2020-01-01")) = 18262 and those which experienced an event (which is here the signature of
# an agreement with the given mechanism)
mon_dat <- mon_dat %>% 
  group_by(dyad_id) %>%
  mutate(dup = ifelse(time == 18262 & duplicated(time, fromLast = TRUE), 1, 0)) %>%
  filter(dup != 1) %>%   dplyr::select(-dup)

mon_dat$time <- mon_dat$time + (min(dyadicdat$time)*(-1))

## II. Conflict Resolution:
# (1) to subset to conflict, we need to replace the value of the event when the agreement was signed as right censored, i.e. the end of the study time
con_dat <- dyadicdat %>%
  arrange(dyad_id, time) %>%
  group_by(dyad_id) %>%
  mutate(sum = sum(conflict)) %>%
  ungroup() %>%
  mutate(time = ifelse(sum== 0, as.numeric(as.Date("2020-01-01")), time)) 
# (2) keep in the dataset only those for which conflict was not present at the end of the study period and those for which conflict was an event that happened before the
# end of the study period:
con_dat <- con_dat[con_dat$conflict == 0 & con_dat$time == 18262 | # conflict was not present at the end of the study period -> right censored
                         con_dat$conflict == 1 & con_dat$time != 18262,] # event observed
# now remove all the observations which did not have conflict present and but another mechanism and are therefore present twice per dyad. These need be excluded as only
# right censored data have the end of the study time which is as.numeric(as.Date("2020-01-01")) = 18262 and those which experienced an event (which is here the signature of
# an agreement with the given mechanism)
con_dat <- con_dat %>% 
  group_by(dyad_id) %>%
  mutate(dup = ifelse(time == 18262 & duplicated(time, fromLast = TRUE), 1, 0)) %>%
  filter(dup != 1) %>%   dplyr::select(-dup)

con_dat$time <- con_dat$time + (min(dyadicdat$time)*(-1))

## III. Commission:
# (1) to subset to commission, we need to replace the value of the event when the agreement was signed as right censored, i.e. the end of the study time
com_dat <- dyadicdat %>%
  arrange(dyad_id, time) %>%
  group_by(dyad_id) %>%
  mutate(sum = sum(commission)) %>%
  ungroup() %>%
  mutate(time = ifelse(sum == 0, as.numeric(as.Date("2020-01-01")), time)) 
# (2) keep in the dataset only those for which commission was not present at the end of the study period and those for which commission was an event that happened before the
# end of the study period:
com_dat <- com_dat[com_dat$commission == 0 & com_dat$time == 18262 | # commission was not present at the end of the study period -> right censored
                     com_dat$commission == 1 & com_dat$time != 18262,] # event observed
# now remove all the observations which did not have commission present and but another mechanism and are therefore present twice per dyad. These need be excluded as only
# right censored data have the end of the study time which is as.numeric(as.Date("2020-01-01")) = 18262 and those which experienced an event (which is here the signature of
# an agreement with the given mechanism)
com_dat <- com_dat %>% 
  group_by(dyad_id) %>%
  mutate(dup = ifelse(time == 18262 & duplicated(time, fromLast = TRUE), 1, 0)) %>%
  filter(dup != 1) %>%   dplyr::select(-dup)

com_dat$time <- com_dat$time + (min(dyadicdat$time)*(-1))

# non-negative numeric time: shift to only positive values to allow time transformation tt() to rectify the violation of assumptions
dyadicdat$time <- dyadicdat$time + (min(dyadicdat$time)*(-1))

## look at objects:
mon_dat
con_dat
com_dat

# how many events per mechanism?
sum(mon_dat$monitoring)
sum(con_dat$conflict)
sum(com_dat$commission)

# are the data set up correctly?
dyadicdat %>%
  arrange(cantons, time) %>%
  dplyr::select(cantons, time, commission, monitoring, conflict, commission_cum, monitoring_cum, conflict_cum, agr_cum, treaty_yes) %>%
  as.data.frame()

# are the data set up correctly?
mon_dat %>%
  arrange(cantons, time) %>%
  dplyr::select(cantons, time, commission, monitoring, conflict, monitoring_cum, conflict_cum, agr_cum, treaty_yes) %>%
  as.data.frame()

# are the data set up correctly?
con_dat %>%
  arrange(cantons, time) %>%
  dplyr::select(cantons, time, commission, monitoring, conflict, monitoring_cum, conflict_cum, agr_cum, treaty_yes) %>%
  as.data.frame()

# are the data set up correctly?
com_dat %>%
  arrange(cantons, time) %>%
  dplyr::select(cantons, time, commission, monitoring, conflict, monitoring_cum, conflict_cum, agr_cum, treaty_yes) %>%
  as.data.frame()

############
# (6) analysis of repeated events
############

# which functional form should the coefficient of symmetry have? Outliers do not seem to be such an issue for this index between zero and one.
# theoretically, I would also not necessarily expect an effect which is decreasing at the margin, which is what sqrt() would imply. I simply tested for this 
# because we talked about this last time during the meeting with the others.
hist(dyadicdat$symmetry) 
hist(sqrt(dyadicdat$symmetry))

res_agreement_1 <-    coxph(Surv(time, treaty_yes) ~  salience + symmetry + border_length + total_size + bi_lingue 
                            + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish
                            + power
                           ,
                            cluster = dyad_id, data = dyadicdat)
res_agreement_2 <-    coxph(Surv(time, treaty_yes) ~  (salience) * symmetry + border_length + total_size + bi_lingue 
                            + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish
                            + power
                           ,
                            cluster = dyad_id, data = dyadicdat)
res_agreement_3 <-    coxph(Surv(time, treaty_yes) ~  sqrt(salience) + sqrt(symmetry) + border_length + total_size + bi_lingue 
                            + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish
                            + power
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_commission_1 <- coxph(Surv(time, commission) ~   (salience) + symmetry  + border_length + total_size + bi_lingue 
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish
                          + power
                          ,
                          cluster = dyad_id, data = com_dat)
res_commission_2 <- coxph(Surv(time, commission) ~   (salience) * symmetry   + border_length + total_size + bi_lingue + commission_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish
                          + power
                          ,
                          cluster = dyad_id, data = com_dat)
res_commission_3 <- coxph(Surv(time, commission) ~   sqrt(salience) + sqrt(symmetry)   + border_length + total_size + bi_lingue + commission_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          + pollution	+ shipping +	fish
                          + power
                          ,
                          cluster = dyad_id, data = com_dat)
res_monitoring_1 <- coxph(Surv(time, monitoring) ~  (salience) + symmetry    + border_length + total_size + bi_lingue + monitoring_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish
                          + power
                          ,
                          cluster = dyad_id, data = mon_dat)
res_monitoring_2 <- coxph(Surv(time, monitoring) ~   (salience) * symmetry   + border_length + total_size + bi_lingue + monitoring_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish
                          + power
                          ,
                          cluster = dyad_id, data = mon_dat)
res_monitoring_3 <- coxph(Surv(time, monitoring) ~  sqrt(salience) + sqrt(symmetry)   + border_length + total_size + bi_lingue + monitoring_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish
                          + power
                          ,
                          cluster = dyad_id, data = mon_dat)
res_conflict_1 <-   coxph(Surv(time, conflict) ~   (salience) + symmetry   + border_length + total_size + bi_lingue + conflict_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping  +	fish
                          + power
                          ,
                          cluster = dyad_id, data = con_dat)
res_conflict_2 <-   coxph(Surv(time, conflict) ~    (salience) * symmetry    + border_length + total_size + bi_lingue + conflict_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping + fish
                          + power
                          ,
                          cluster = dyad_id, data = con_dat)
res_conflict_3 <-   coxph(Surv(time, conflict) ~   sqrt(salience) + sqrt(symmetry)   + border_length + total_size + bi_lingue + conflict_cum  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping + fish
                          + power
                          ,
                          cluster = dyad_id, data = con_dat)


library(stargazer)
stargazer(res_agreement_1, res_agreement_2, res_agreement_3, res_commission_1, res_commission_2, res_commission_3, res_monitoring_1, res_monitoring_2, res_monitoring_3, res_conflict_1, res_conflict_2, res_conflict_3,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements",
          # column.labels   = c("Mechanisms", "Agreement"),
          column.separate = c(3, 1),
          dep.var.labels=c("Agreement", "Commission", "Monitoring", "Conflict"),
          # covariate.labels=c("Salience",
          #                    "Symmetry",
          #                    "Border Length",
          #                    "Dyad Size",
          #                    "Bi-lingue",
          #                    "Prior Commission",
          #                    "Prior Monitoring",
          #                    "Prior Conflict Resolution",
          #                    "Prior Agreements",
          #                    "Pollution",
          #                    "Shipping",
          #                    "Fish",
          #                    "Power"),
          column.sep.width = "-20pt",
          notes = "standard errors clustered by dyad",
          align=TRUE,
          no.space = TRUE,
          df = FALSE,
          float.env = "sidewaystable",
          stargazer_stat_code_list = list(c("aic", "bic")),
          star.char = c("°", "*", "**", "***"),
          star.cutoffs = c(.1, .05, .01, .001)
          # se=list(c(summary(res_enf_cant)$coefficients[, 2], summary(prob2)$coefficients[, 2], summary(prob3)$coefficients[, 2], summary(prob4)$coefficients[, 2] ))
)

res_agreement_1.1 <-    coxph(Surv(time, treaty_yes) ~ 
                            + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish 
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_agreement_2.1 <-    coxph(Surv(time, treaty_yes) ~  
                            + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_commission_1.1<- coxph(Surv(time, commission) ~  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish 
                          ,
                          cluster = dyad_id, data = com_dat)
res_commission_2.1 <- coxph(Surv(time, commission) ~  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          ,
                          cluster = dyad_id, data = com_dat)
res_monitoring_1.1 <- coxph(Surv(time, monitoring) ~  
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping  +	fish 
                          ,
                          cluster = dyad_id, data = mon_dat)
res_monitoring_2.1 <- coxph(Surv(time, monitoring) ~ 
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          ,
                          cluster = dyad_id, data = mon_dat)
res_conflict_1.1 <-   coxph(Surv(time, conflict) ~   
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping  +	fish 
                          ,
                          cluster = dyad_id, data = con_dat)
res_conflict_2.1 <-   coxph(Surv(time, conflict) ~    
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          ,
                          cluster = dyad_id, data = con_dat)


library(stargazer)
stargazer(res_agreement_1.1, res_agreement_2.1, res_commission_1.1, res_commission_2.1, res_monitoring_1.1, res_monitoring_2.1, res_conflict_1.1, res_conflict_2.1,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements",
          # column.labels   = c("Mechanisms", "Agreement"),
          column.separate = c(3, 1),
          dep.var.labels=c("Agreement", "Commission", "Monitoring", "Conflict"),
          # covariate.labels=c("Salience",
          #                    "Symmetry",
          #                    "Border Length",
          #                    "Dyad Size",
          #                    "Bi-lingue",
          #                    "Prior Commission",
          #                    "Prior Monitoring",
          #                    "Prior Conflict Resolution",
          #                    "Prior Agreements",
          #                    "Pollution",
          #                    "Shipping",
          #                    "Fish",
          #                    "Power"),
          column.sep.width = "-20pt",
          notes = "standard errors clustered by dyad",
          align=TRUE,
          no.space = TRUE,
          df = FALSE,
          float.env = "sidewaystable",
          stargazer_stat_code_list = list(c("aic", "bic")),
          star.char = c("°", "*", "**", "***"),
          star.cutoffs = c(.1, .05, .01, .001)
          # se=list(c(summary(res_enf_cant)$coefficients[, 2], summary(prob2)$coefficients[, 2], summary(prob3)$coefficients[, 2], summary(prob4)$coefficients[, 2] ))
)

############
# (7) assumptions 
############
## What to do if Assumptions are violated? https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf instructions by Therneau et al (2020)
# Cox Prop Hazard: Tbe effect of a variable is not constant over time
# (1) step function, split the data 
# (2) Continuous time-dependent coefficients

# test proportional-hazards (PH) assumption
(zph_agr <- cox.zph(res_agreement_1)) # p.vals below 0.05 indicate that there is an interaction with time
# plot to investigate relationship with time
# pdf("prop_haz_agr.pdf")
ggcoxzph(zph_agr)
dev.off()
# outliers
dev.off()
# pdf("pol_bil_outliers_agr.pdf")
ggcoxdiagnostics(res_agreement_1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
dev.off()

res_agreement_1.1 <-    coxph(Surv(time, treaty_yes) ~  salience + symmetry + border_length + tt(border_length) + total_size + tt(total_size) + bi_lingue 
                              + agr_cum + tt(agr_cum) + commission_cum + tt(commission_cum) + monitoring_cum + tt(monitoring_cum) + conflict_cum + tt(conflict_cum)
                              +	pollution + tt(pollution)	+ shipping +	fish 
                              # + river
                              ,
                              cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))

# test proportional-hazards (PH) assumption
(zph_agr <- cox.zph(res_commission_1)) # p.vals below 0.05 indicate that there is an interaction with time
# plot to investigate relationship with time
# pdf("prop_haz_agr.pdf")
ggcoxzph(zph_agr)
dev.off()
# outliers
dev.off()
# pdf("pol_bil_outliers_agr.pdf")
ggcoxdiagnostics(res_agreement_1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
dev.off()

res_commission_1.1 <- coxph(Surv(time, commission) ~   (salience) + symmetry  + border_length + total_size + bi_lingue 
                          + agr_cum + tt(agr_cum) + commission_cum + tt(commission_cum) + monitoring_cum + tt(monitoring_cum) + conflict_cum + tt(conflict_cum) 
                          +	pollution + tt(pollution)	+ shipping +	fish 
                          # + river
                          ,
                          cluster = dyad_id, data = com_dat, tt = function(x, t, ...) x * log(t+20))

# test proportional-hazards (PH) assumption
(zph_agr <- cox.zph(res_monitoring_1)) # p.vals below 0.05 indicate that there is an interaction with time
# plot to investigate relationship with time
# pdf("prop_haz_agr.pdf")
ggcoxzph(zph_agr)
dev.off()
# outliers
dev.off()
# pdf("pol_bil_outliers_agr.pdf")
ggcoxdiagnostics(res_monitoring_1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
dev.off()

res_monitoring_1.1 <- coxph(Surv(time, monitoring) ~  (salience) + symmetry    + border_length + tt(border_length) + total_size + bi_lingue + monitoring_cum + tt(monitoring_cum)  
                          + agr_cum + tt(agr_cum) + commission_cum + tt(commission_cum) + monitoring_cum + tt(monitoring_cum) + conflict_cum + tt(conflict_cum) 
                          +	pollution	+ shipping +	fish 
                          # + river
                          ,
                          cluster = dyad_id, data = mon_dat, tt = function(x, t, ...) x * log(t+20))

# test proportional-hazards (PH) assumption
(zph_agr <- cox.zph(res_conflict_1)) # p.vals below 0.05 indicate that there is an interaction with time
# plot to investigate relationship with time
# pdf("prop_haz_agr.pdf")
ggcoxzph(zph_agr)
dev.off()
# outliers
dev.off()
# pdf("pol_bil_outliers_agr.pdf")
ggcoxdiagnostics(res_conflict_1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

dev.off()

res_conflict_1.1 <-   coxph(Surv(time, conflict) ~   (salience) + symmetry   + border_length + total_size + bi_lingue
                          + agr_cum + tt(agr_cum) + commission_cum + tt(commission_cum) + monitoring_cum  + conflict_cum + tt(conflict_cum)  
                          +	pollution + tt(pollution)	+ shipping + tt(shipping)  + fish +	tt(fish) 
                          # + power
                          # + river
                          ,
                          cluster = dyad_id, data = con_dat, tt = function(x, t, ...) x * log(t+20))

library(stargazer)
stargazer(res_agreement_1.1, res_commission_1.1, res_monitoring_1.1, res_conflict_1.1, 
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements",
          # column.labels   = c("Mechanisms", "Agreement"),
          column.separate = c(3, 1),
          dep.var.labels=c("Agreement", "Commission", "Monitoring", "Conflict"),
          # covariate.labels=c("Salience",
          #                    "Symmetry",
          #                    "Border Length",
          #                    "Dyad Size",
          #                    "Bi-lingue",
          #                    "Prior Commission",
          #                    "Prior Monitoring",
          #                    "Prior Conflict Resolution",
          #                    "Prior Agreements",
          #                    "Pollution",
          #                    "Shipping",
          #                    "Fish",
          #                    "Power"),
          column.sep.width = "-20pt",
          notes = "standard errors clustered by dyad",
          align=TRUE,
          no.space = TRUE,
          df = FALSE,
          float.env = "sidewaystable",
          stargazer_stat_code_list = list(c("aic", "bic")),
          star.char = c(".", "*", "**", "***"),
          star.cutoffs = c(.1, .05, .01, .001)
          # se=list(c(summary(res_enf_cant)$coefficients[, 2], summary(prob2)$coefficients[, 2], summary(prob3)$coefficients[, 2], summary(prob4)$coefficients[, 2] ))
)


# 
# # plot results:
# agr_df <- with(dyadicdat,
#                data.frame(salience = c(summary(salience)[2], summary(salience)[5]),
#                           symmetry = c(summary(symmetry)[1], summary(symmetry)[6]),
#                           border_length = rep(mean(border_length, na.rm = TRUE), 2),
#                           total_size = rep(mean(total_size, na.rm = TRUE), 2),
#                           bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2),
#                           agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
#                           border = rep(mean(border, na.rm = TRUE), 2),
#                           construction = rep(mean(construction, na.rm = TRUE), 2),
#                           pollution = rep(mean(pollution, na.rm = TRUE), 2),
#                           shipping = rep(mean(shipping, na.rm = TRUE), 2),
#                           fish = rep(mean(fish, na.rm = TRUE), 2),
#                           power = rep(mean(power, na.rm = TRUE), 2),
#                           river = rep(mean(river, na.rm = TRUE), 2)
#                )
# )
# 
# agr_df_1 <- with(dyadicdat,
#                data.frame(salience = rep(summary(dyadicdat$salience)[5], 2),
#                           symmetry = c(summary(symmetry)[1], summary(symmetry)[6]),
#                           border_length = rep(mean(border_length, na.rm = TRUE), 2),
#                           total_size = rep(mean(total_size, na.rm = TRUE), 2),
#                           bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2),
#                           agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
#                           border = rep(mean(border, na.rm = TRUE), 2),
#                           construction = rep(mean(construction, na.rm = TRUE), 2),
#                           pollution = rep(mean(pollution, na.rm = TRUE), 2),
#                           shipping = rep(mean(shipping, na.rm = TRUE), 2),
#                           fish = rep(mean(fish, na.rm = TRUE), 2),
#                           power = rep(mean(power, na.rm = TRUE), 2),
#                           river = rep(mean(river, na.rm = TRUE), 2)
#                )
# )
# 
# com_df <- with(com_dat,
#                data.frame(interdependence =  rep(mean(interdependence, na.rm = TRUE), 2),
#                           salience = rep(mean(salience, na.rm = TRUE), 2),
#                           symmetry = c(summary(symmetry)[1], summary(symmetry)[6]),
#                           border_length = rep(mean(border_length, na.rm = TRUE), 2),
#                           total_size = rep(mean(total_size, na.rm = TRUE), 2),
#                           bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2),
#                           agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
#                           commission_cum = rep(mean(enforcement_cum, na.rm = TRUE), 2)
#                )
# )
# enf_df <- with(enf_dat,
#                data.frame(interdependence =  rep(mean(interdependence, na.rm = TRUE), 2),
#                           salience = rep(mean(salience, na.rm = TRUE), 2),
#                           symmetry = c(summary(symmetry)[1], summary(symmetry)[6]),
#                           border_length = rep(mean(border_length, na.rm = TRUE), 2),
#                           total_size = rep(mean(total_size, na.rm = TRUE), 2),
#                           bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2),
#                           agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
#                           enforcement_cum = rep(mean(enforcement_cum, na.rm = TRUE), 2)
#                )
# )
# mon_df <- with(mon_dat,
#                data.frame(interdependence =  rep(mean(interdependence, na.rm = TRUE), 2),
#                           salience = rep(mean(salience, na.rm = TRUE), 2),
#                           symmetry = c(summary(symmetry)[1], summary(symmetry)[6]),
#                           border_length = rep(mean(border_length, na.rm = TRUE), 2),
#                           total_size = rep(mean(total_size, na.rm = TRUE), 2),
#                           bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2),
#                           agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
#                           monitoring_cum = rep(mean(monitoring_cum, na.rm = TRUE), 2)
#                )
# )
# con_df <- with(con_dat,
#                data.frame(interdependence =  rep(mean(interdependence, na.rm = TRUE), 2),
#                           salience = rep(mean(salience, na.rm = TRUE), 2),
#                           symmetry = c(summary(symmetry)[1], summary(symmetry)[6]),
#                           border_length = rep(mean(border_length, na.rm = TRUE), 2),
#                           total_size = rep(mean(total_size, na.rm = TRUE), 2),
#                           bi_lingue = rep(mean(bi_lingue, na.rm = TRUE), 2),
#                           agr_cum = rep(mean(agr_cum, na.rm = TRUE), 2),
#                           conflict_cum = rep(mean(conflict_cum, na.rm = TRUE), 2)
#                )
# )
# 
# fit_agr <- survfit(res_agreement_1, newdata = agr_df_1)
# fit_com <- survfit(res_commission, newdata = com_df)
# fit_mon <- survfit(res_monitoring, newdata = mon_df)
# fit_con <- survfit(res_conflict, newdata = con_df)
# 
# par(mfrow=c(1,3),mar=c(3, 2.5, 2, 2))
# ggsurvplot(fit_agr, 
#            conf.int = TRUE, 
#            legend.labs=c("1^st quartile", "3^rd qaurtile"),
#            xlab = "Time in days", 
#            break.time.by = 4000, 
#            ggtheme = theme_light(), # customize plot and risk table with a theme.
#            # risk.table = "abs_pct",  # absolute number and percentage at risk.
#            # risk.table.y.text.col = T,# colour risk table text annotations.
#            # risk.table.y.text = FALSE,# show bars instead of names in text annotations
#            # in legend of risk table.
#            # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
#            surv.median.line = "hv",  # add the median survival pointer.
#            # xlim = c(0, 2000),
#            # ggtheme = theme_minimal(), 
#            data = agr_df)
# ggsurvplot(fit_mon, conf.int = TRUE, legend.labs=c("1^st quartile", "3^rd qaurtile"),
#            ggtheme = theme_minimal(),  data = mon_df)
# ggsurvplot(fit_con, conf.int = TRUE, legend.labs=c("1^st quartile", "3^rd qaurtile"),
#            ggtheme = theme_minimal(),  data = con_df)
# ggsurvplot(fit_com, conf.int = TRUE, legend.labs=c("1^st quartile", "3^rd qaurtile"),
#            ggtheme = theme_minimal(),  data = com_df)

