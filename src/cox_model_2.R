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
# (6) assumptions and rectifications

######################################
# (1) treaty design data
######################################

# change directory here to replicate 
setwd("/Volumes/Transcend/Uni/Zivi Eawag/Script/Model/ModelInput")

load("/Volumes/Transcend/Uni/Zivi Eawag/Script/Model/ModelInput/dyadicdat.RData")

dyadicdat <- dyadicdat %>% filter(bilateral == 1)

dyadicdat %>% 
  group_by(ID_SM) %>% 
  filter(bilateral == 1) %>% 
  slice(1)

# # subset to those issues included
dyadicdat <- dyadicdat %>%
  filter(pollution == 1 | fish == 1 | shipping == 1 | construction == 1 )

# create DV for agreements
dyadicdat <- dyadicdat %>%
  mutate(treaty_yes = ifelse(!is.na(ID_SM), 1, 0))

dyadicdat <- dyadicdat %>%
  mutate(year = as.numeric(format(as.Date(time, origin = "1970-01-01"),'%Y')),
         date = as.numeric(format(as.Date(time, origin = "1970-01-01"),'%Y%m%d')))

##########
# (2) risk set
##########

# risk set rivers
load("cantons_border.RData")
cont_adj <- cont_adj[order(rownames(cont_adj)), order(colnames(cont_adj))]
cont_adj[cont_adj >0] <- 1
dyads <- as.data.frame(get.edgelist(graph_from_adjacency_matrix(cont_adj, mode = "undirected")))

# make sure the dyads are ordered by name, these are the cantons in the risk set. Only contiguous cantons can sign agreements
colnames(dyads) <- c("canton1", "canton2")
dyads <- dyads %>%
  mutate(canton1 = as.character(canton1),
         canton2 = as.character(canton2),
         cantons = ifelse(canton1 < canton2, paste(canton1, canton2), paste(canton2, canton1))) %>%
  distinct(cantons, .keep_all = T)

############
# (3) variables
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
         monitoring_cum = cumsum(mon_cant) -1,
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
library(TraMineRextras)
seq <- dyadicdat %>%
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
  seqecreate(., id = id, timestamp = time, event = event) 

seq %>%
  seqpcplot(.,
            # alphabet = alph,
            filter = list(type = "function",
                          value = "cumfreq",
                          level = .6),
            order.align = "last",
            # ltype = "non-embeddable",
            cex = 1.5, lwd = .4,
            lcourse = "downwards")

seq_sub <- seqefsub(seq, pmin.support = 0.05)

ggplot(seq_sub) + 
  geom_bar(mapping = aes(y=Support, fill = Period) ,width=.3,
           stat="identity", position = "dodge") +
  theme_classic() +
  theme(
    # text = element_text(size=15),
    # axis.text=element_text(size=15, angle = 90),
    # axis.ticks.length.x = unit(.2, "cm"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
    # panel.grid.major.y = element_line(colour = "grey40"),
    # panel.grid.major.x = element_blank(),
    # panel.spacing = unit(3, "lines"),
    # strip.background = element_rect(color = "white"),
    # strip.text.x = element_text(size = 15),
    # strip.placement = "inside",
    # legend.position = "none"
  ) 

seq <- as.data.frame(seq_sub$data$Support)
seq$seq_sub <- as.character(seq_sub$subseq)
colnames(seq) <- c("support", "subseq")

# or, order by prevalence:
seq$subseq <- factor(seq$subseq, levels = seq$subseq[rev(order(seq$support))])

ggplot(seq[1:20,], aes(x = subseq, y= support)) +
  geom_bar( position="dodge",stat="identity", width=.5) +
  # coord_flip() +
  ylab("Proportion of Subsequences") +   xlab("Sub Sequence") +
  theme_classic() +
  theme(
    text = element_text(size=15),
    axis.text=element_text(size=15, angle = 90),
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

# plot(seq_sub[1:15], col = "blue", ylab = "Frequency", xlab = "Subsequences", cex = 1.5)

# transform to event sequences
# https://stackoverflow.com/questions/28965909/plot-event-sequences-event-sequences-clustering and 
# Studer, M., Müller, N.S., Ritschard, G. & Gabadinho, A. (2010), "Classer, discriminer et visualiser des séquences d'événements", In Extraction et gestion des connaissances (EGC 2010), Revue des nouvelles technologies de l'information RNTI. Vol. E-19, pp. 37-48.
# Ritschard, G., Bürgin, R., and Studer, M. (2014), "Exploratory Mining of Life Event Histories", In McArdle, J.J. & Ritschard, G. (eds) Contemporary Issues in Exploratory Data Mining in the Behavioral Sciences. Series: Quantitative Methodology, pp. 221-253. New York: Routledge.

# idcost <- 1:3
# dd <- seqedist(seq, idcost=idcost, vparam=.1)
# 
# library(cluster)
# clusterward1 <- agnes(dd, diss = TRUE, method = "ward")
# plot(clusterward1, which.plot = 2)
# cl1.4 <- cutree(clusterward1, k = 4)
# cl1.4fac <- factor(cl1.4, labels = paste("Type", 1:4))
# 
time_to_event <- function(dyadicdat, var1, var2){
  library( data.table )
  DT <- as.data.table(dyadicdat)
  #new clumn-name will be set autmatically based on vars
  colname <- paste("t", var1, "from",  var2, sep="_")
  #melt to long
  DT.melt <- melt(DT, id.vars = c("cantons", "time"), measure.vars = c(var2, var1))
  #only keep 1's
  DT.melt <- DT.melt[ value == 1, ]
  #set keys
  setkey( DT.melt, cantons, time )
  #get time of previous type2 for all rows with type2
  temp <- DT.melt[ variable == var1 & value == 1, ][ DT.melt,
                                                     (colname) := {
                                                       #create on-the-fly subset
                                                       val = DT.melt[ cantons == i.cantons & value == 1 & variable == var2 & time < i.time, ]
                                                       list( min( i.time - val$time ) )
                                                     }, by = .EACHI ][]

  temp[ is.infinite( get(colname) ), (colname) := NA ][]

  # join back to original
  # use eval+parse to keep the colname variable
  expr = paste0("DT[ temp, (colname) := i.", colname, ", on = .(cantons, time)]")
  eval(parse(text=expr))

  dyadicdat <- as.data.frame(DT)
  return(dyadicdat)
}

plot_dat <- dyadicdat[dyadicdat$year >= 1980 &
                        dyadicdat$year != 0,]
# apply this function
plot_dat <- time_to_event(plot_dat, "fish", "treaty_yes")
plot_dat <- time_to_event(plot_dat, "construction", "treaty_yes")
plot_dat <- time_to_event(plot_dat, "pollution", "treaty_yes")
plot_dat <- time_to_event(plot_dat, "shipping", "treaty_yes")

plot_dat <- time_to_event(plot_dat, "fish", "commission")
plot_dat <- time_to_event(plot_dat, "construction", "commission")
plot_dat <- time_to_event(plot_dat, "pollution", "commission")
plot_dat <- time_to_event(plot_dat, "shipping", "commission")

plot_dat <- time_to_event(plot_dat, "fish", "conflict")
plot_dat <- time_to_event(plot_dat, "construction", "conflict")
plot_dat <- time_to_event(plot_dat, "pollution", "conflict")
plot_dat <- time_to_event(plot_dat, "shipping", "conflict")

plot_dat <- time_to_event(plot_dat, "conflict", "commission")
plot_dat <- time_to_event(plot_dat, "conflict", "monitoring")

plot_dat <- time_to_event(plot_dat, "commission", "commission")
plot_dat <- time_to_event(plot_dat, "monitoring", "monitoring")
plot_dat <- time_to_event(plot_dat, "conflict", "conflict")

plot_dat %>%
  select(cantons, time, commission, monitoring, conflict,
         t_fish_from_treaty_yes, t_construction_from_treaty_yes,
         t_pollution_from_treaty_yes, t_shipping_from_treaty_yes
         ) %>%
  as.data.frame()

p <- plot_dat %>%
  reshape2::melt(., id.vars = "time", measure.vars = c("t_fish_from_treaty_yes", "t_construction_from_treaty_yes", "t_pollution_from_treaty_yes", "t_shipping_from_treaty_yes",
                                                       "t_fish_from_commission", "t_construction_from_commission", "t_pollution_from_commission", "t_shipping_from_commission",
                                                       "t_fish_from_conflict", "t_construction_from_conflict", "t_pollution_from_conflict", "t_shipping_from_conflict")) %>%
  dplyr::filter(!is.na(value)) 

p$variable2[p$variable == "t_commission_from_commission" | p$variable == "t_commission_from_monitoring" | p$variable == "t_commission_from_conflict"] <- "commission"
p$variable2[p$variable == "t_monitoring_from_commission" | p$variable == "t_monitoring_from_monitoring" | p$variable == "t_monitoring_from_conflict"] <- "monitoring"
p$variable2[p$variable == "t_conflict_from_commission" | p$variable == "t_conflict_from_monitoring" | p$variable == "t_conflict_from_conflict"] <- "conflict"

ggplot(p, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot(aes(x=variable, y=value, fill=variable)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_minimal() 

## correlation between mechanisms
chi <- matrix(nrow = 6, ncol = 3)

colnames(chi) <- c("commission", "monitoring",  "conflict")
rownames(chi) <- c("commission", "", "monitoring", "",  "conflict", "")

chi[1,2] <- round(chisq.test(table(plot_dat$commission, plot_dat$monitoring))$statistic, 3)
chi[2,2] <- paste0( "(", round(chisq.test(table(plot_dat$commission, plot_dat$monitoring))$p.value, 3), ")")
chi[1,3] <- round(chisq.test(table(plot_dat$commission, plot_dat$conflict))$statistic, 3)
chi[2,3] <- paste0( "(", round(chisq.test(table(plot_dat$commission, plot_dat$conflict))$p.value, 3), ")")
chi[3,3] <- round(chisq.test(table(plot_dat$conflict, plot_dat$monitoring))$statistic, 3)
chi[4,3] <- paste0( "(", round(chisq.test(table(plot_dat$commission, plot_dat$monitoring))$p.value, 3), ")")
chi[is.na(chi)] <- "-"

latex(chi, file="")

## simple descriptives of co-appearance
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


colnames(co_app) <- c("commission", "monitoring", "conflict")
rownames(co_app) <- c("pollution", "shipping", "fish", "construction")

latex(co_app, file="")

# how many agreeement per issue?
sum(plot_dat$pollution)
sum(plot_dat$shipping)
sum(plot_dat$fish)
sum(plot_dat$construction)

sum(plot_dat$commission)
sum(plot_dat$monitoring)
sum(plot_dat$conflict)

sum(plot_dat$treaty_yes)

## plot 2
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

# plot the number of new agreements + design principles
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
                                              "conflict" = "Conflict Resolution (yes)"))) +
  theme_classic() +
  # coord_flip() +
  theme(
    text = element_text(size=15),
    axis.text=element_text(size=15, angle = 90),
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

# plot the number of new agreements by issue
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
    text = element_text(size=15),
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

# ## plot 2
# dyadicdat %>%
#   group_by(ID_SM) %>%
#   slice(1) %>%
#   as_tibble() %>%
#   # filter(year >= 1980) %>%
#   mutate(date = as.Date(time, origin = "1970-01-01"),
#          year = as.numeric(format(date,'%Y')),
#          agreement_yes = ifelse(!is.na(ID_SM), 1, 0),
#          Period = NA,
#          Period = replace(Period,  year< 1930, "<1930"), 
#          Period = replace(Period, year>= 1930 & year< 1940, "1930-1939"), 
#          Period = replace(Period, year>= 1940 & year< 1950, "1940-1949"), 
#          Period = replace(Period, year>= 1950 & year< 1960, "1950-1959"), 
#          Period = replace(Period, year>= 1960 & year< 1970, "1960-1969"), 
#          Period = replace(Period, year>= 1970 & year< 1980, "1970-1979"), 
#          Period = replace(Period, year>= 1980 & year< 1990, "1980-1989"), 
#          Period = replace(Period, year>= 1990 & year< 2000, "1990-1999"),
#          Period = replace(Period, year>= 2000 & year< 2010, "2000-2009"),
#          Period = replace(Period, year>= 2010 & year<= 2020, "2010-2020")) %>%
#   reshape2::melt(., id.vars = "Period", measure.vars = c("agreement_yes", "commission", "monitoring", "conflict")) %>%
#   group_by(Period, variable) %>%
#   summarise(Frequency = sum(value)) %>%
#   ggplot() + 
#   geom_bar(aes(y = Frequency, x = Period, width=.7, fill = variable, group = variable), # color = variable, 
#            width=.7,
#            stat="identity", position = position_dodge(width=.8)) +
#   scale_fill_manual(values=c("#bdbdbd","#969696", "#525252", "#252525" ), name = "Condition") + 
#   # facet_wrap(~variable) +
#   theme_classic() +
#   # coord_flip() +
#   theme(
#     axis.text=element_text(size=10, angle = 90),
#     panel.grid.major.y = element_line(colour = "grey80"),
#     panel.grid.major.x = element_blank(),
#     panel.spacing = unit(1, "lines")) 
# dev.off()
# 
# devtools::install_github("davedgd/ggNestedBarChart")
# library(ggNestedBarChart)
# 
# p1 <- dyadicdat %>%
#   group_by(ID_SM) %>%
#   slice(1) %>%
#   as_tibble() %>%
#   # filter(year >= 1980) %>%
#   mutate(date = as.Date(time, origin = "1970-01-01"),
#          year = as.numeric(format(date,'%Y')),
#          agreement_yes = ifelse(!is.na(ID_SM), 1, 0),
#          Period = NA,
#          Period = replace(Period,  year< 1930, "<1930"), 
#          Period = replace(Period, year>= 1930 & year< 1940, "1930-1939"), 
#          Period = replace(Period, year>= 1940 & year< 1950, "1940-1949"), 
#          Period = replace(Period, year>= 1950 & year< 1960, "1950-1959"), 
#          Period = replace(Period, year>= 1960 & year< 1970, "1960-1969"), 
#          Period = replace(Period, year>= 1970 & year< 1980, "1970-1979"), 
#          Period = replace(Period, year>= 1980 & year< 1990, "1980-1989"), 
#          Period = replace(Period, year>= 1990 & year< 2000, "1990-1999"),
#          Period = replace(Period, year>= 2000 & year< 2010, "2000-2009"),
#          Period = replace(Period, year>= 2010 & year<= 2020, "2010-2020")) %>%
#   reshape2::melt(., id.vars = "Period", measure.vars = c("agreement_yes", "commission", "monitoring", "conflict")) %>%
#   group_by(Period, variable) %>%
#   summarise(Frequency = sum(value)) %>%
#   ggplot(., aes(y = Frequency, x = Period, width=.7, fill = variable, group = variable), stat = "identity") +
#   geom_bar(stat = "identity") +
#   # facet_wrap(~variable) +  
#   scale_fill_manual(values=c("#bdbdbd","#969696", "#525252", "#252525" ), name = "Condition") + 
#   facet_wrap(vars(variable, Period), strip.position = "top", scales = "free_x", nrow = 1) +
#   theme_bw(base_size = 13) +
#   theme(panel.spacing = unit(0, "lines"),
#         strip.background = element_rect(color = "black", size = 0, fill = "grey92"),
#         strip.placement = "outside",
#         axis.text=element_text(size=10, angle = 90),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 0),
#         panel.background = element_rect(fill = "white"),
#         legend.position = "none") +
#   scale_y_continuous(expand = expand_scale(mult = c(0, .1))) 
#   # geom_text(aes(label = paste0(Frequency, "%")), position = position_stack(0.5), color = "white", fontface = "bold")
# 
# p1
# 
# ggNestedBarChart(p1)

## standard summary table
library(summarytools)
sum_dat <- dyadicdat[dyadicdat$year >= 1980,] %>%  ungroup() %>%
  dplyr::select(treaty_yes, commission, monitoring, conflict, agr_cum, commission_cum, monitoring_cum, conflict_cum, pollution, shipping, fish, construction, 
                bi_lingue, symmetry, border_length, total_size) 

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
                       "Pollution",
                       "Shipping",
                       "Fish",
                       "Construction",
                       "Bi-lingue",
                       "Symmetry",
                       "Border Length [km]",
                       "Size  [km$^2$]")
latex(sum_tab, file="")

# are the data set up correctly?
dyadicdat %>%
  arrange(cantons, time) %>%
  dplyr::select(ID_SM, cantons, time, year, date, commission, monitoring, conflict, commission_cum, monitoring_cum, conflict_cum, agr_cum, treaty_yes) %>%
  as.data.frame()


dyadicdat <- dyadicdat %>%
  mutate(any_ext_yes = ifelse(mon_comm == 1 | conf_body == 1 | conf_fed == 1 | commission == 1, 1, 0),
         any_mech_yes = ifelse(commission == 1 | monitoring == 1 | conflict == 1, 1, 0))

############
# (5) analysis of repeated events
############


res_agreement_main_1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_agreement_main_2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_any_mech_main_1 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_any_mech_main_2 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                           +	pollution	+ shipping +	fish + construction
                           + symmetry + bi_lingue,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_main_1 <-   coxph(Surv(time, commission) ~  
                            + agr_cum + commission_cum + monitoring_cum + conflict_cum,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_main_2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish + construction
                          + symmetry + bi_lingue  
                          ,
                          cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_main_1 <-   coxph(Surv(time, mon_cant) ~ 
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          ,
                          cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_main_2 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping +	fish + construction
                          + symmetry + bi_lingue  
                          ,
                          cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_main_1 <-     coxph(Surv(time, conflict) ~    
                          + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          ,
                          cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_main_2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                          +	pollution	+ shipping  + construction
                          + symmetry + bi_lingue  
                          ,
                          cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

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
          # star.char = c("\\cdot", "*", "**", "***"),
          # star.cutoffs = c(.1, .05, .01, .001),
          style = "ajs",
          omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("AIC", as.character(round(unlist(lapply(model_list, function(x) AIC(x))), 3))),
                           c("BIC", as.character(round(unlist(lapply(model_list, function(x) BIC(x))), 3))),
                           c("Events", as.character(unlist(lapply(model_list, function(x) (x)$nevent))))
                           )
)


res_agreement_1.1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum ,
                            cluster = dyad_id, data = dyadicdat)
res_agreement_2.1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue,
                            cluster = dyad_id, data = dyadicdat)
res_any_mech_1.1 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum,
                           cluster = dyad_id, data = dyadicdat)
res_any_mech_2.1 <-    coxph(Surv(time, any_mech_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                           +	pollution	+ shipping +	fish + construction
                           + symmetry + bi_lingue,
                           cluster = dyad_id, data = dyadicdat)
res_commission_1.1 <-   coxph(Surv(time, commission) ~  
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum,
                            cluster = dyad_id, data = dyadicdat)
res_commission_2.1 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_monitoring_1.1 <-   coxph(Surv(time, mon_cant) ~ 
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_monitoring_2.1 <-   coxph(Surv(time, mon_cant) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_conflict_1.1 <-     coxph(Surv(time, conflict) ~    
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat)
res_conflict_2.1 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping  + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat)

model_list <- list(res_agreement_1, res_agreement_2, res_agreement_1.1, res_agreement_2.1, 
                   res_any_mech_1, res_any_mech_2, res_any_mech_1.1, res_any_mech_2.1,
                   res_commission_1, res_commission_2, res_commission_1.1, res_commission_2.1, 
                   res_monitoring_1, res_monitoring_2, res_monitoring_1.1, res_monitoring_2.1, 
                   res_conflict_1, res_conflict_2, res_conflict_1.1, res_conflict_2.1)

library(stargazer)
stargazer(model_list,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements after 1980-01-01",
          # column.separate = c(3, 1),
          # dep.var.labels=c("Agreement (yes)", "Commission (yes)", "Monitoring (yes)", "Conflict (yes)"),
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

results <- function(res_agreement_1){
  names <- names(summary(res_agreement_1)$coefficients[, 2])
  yhat  <- summary(res_agreement_1)$coefficients[, 2]
  se    <- summary(res_agreement_1)$coefficients[, 4]
  z <- qnorm(.975)
  lower <- yhat - z*se
  upper <- yhat + z*se
  index <- 1:length(names)
  
  df <- data.frame(names, yhat, se, lower, upper, index)
  
  return(df)
}

results_list <- lapply(model_list, function(x) results(x))
results_list <- lapply(seq_along(results_list), function(x) cbind(results_list[[x]], unique.id=x))
results <- do.call("rbind", results_list)
results$dv[results$unique.id %in% c(1:4)] <- "Agreement"
results$dv[results$unique.id %in% c(5:8)] <- "Any Mechanims"
results$dv[results$unique.id %in% c(9:12)] <- "Commission"
results$dv[results$unique.id %in% c(13:16)] <- "Montoring"
results$dv[results$unique.id %in% c(17:20)] <- "Conflict"

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
  scale_shape_manual(name = "model",
                     labels = c("Main Effects", "Full Model"),
                     values = c(1,2)) +
  scale_colour_manual(name = "shape",
                      labels = c("1980-2020", "1841-2020"),
                      values = c("grey","black")) +
  # scale_shape_manual(values = c(16,17,16,17, 16,17)) +
  theme(text=element_text(size=15, color="black"),
        strip.text.x = element_text(size = 15),
        panel.spacing = unit(3, "lines"),
        axis.text.x = element_text(size=15, vjust=0.5, color = 'black'),  # x-axis labels
        axis.text.y = element_text(size=15, vjust=0.5, color = 'black'),  # y-axis labels
        axis.title.x = element_text(size=17.5, vjust=0.1),                # x-title justification  
        axis.title.y = element_text(size=17.5, vjust=1.5),                 # y-title justification
        legend.title = element_blank(),
        legend.position = "bottom",       
        legend.text = element_text(size = 14.5)
        ) +
  facet_wrap(~dv, ncol=5,
             labeller = labeller(variable = c("treaty_yes" = "Agreement (yes)", 
                                              "commission" = "Commission (yes)",  
                                              "monitoring" = "Monitoring (yes)", 
                                              "conflict" = "Conflict Resolution (yes)"))) 


res_agreement_1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_agreement_2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_agreement_3 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry 
                            + total_size + border_length
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_1 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_commission_3 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry  + total_size + border_length
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_1 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_2 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_monitoring_3 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry  + total_size + border_length
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_1 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_conflict_3 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping  + construction
                            + symmetry  + total_size + border_length
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

model_list <- list(res_agreement_1, res_agreement_2, res_agreement_3, 
                   res_commission_1, res_commission_2, res_commission_3, 
                   res_monitoring_1, res_monitoring_2, res_monitoring_3, 
                   res_conflict_1, res_conflict_2, res_conflict_3)

plot(summary(res_agreement_1))

stargazer(model_list,
          type = "latex",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, bilateral agreements after 1980-01-01",
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
                             "Total Size [km$^2$]",
                             "Border Length [km]"),
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

cox.zph(res_agreement_2.1)
cox.zph(res_any_mech_2.1)
cox.zph(res_commission_2.1)
cox.zph(res_monitoring_2.1)
cox.zph(res_conflict_2.1)

res_agreement_2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum +  tt(agr_cum) + commission_cum + tt(commission_cum) + monitoring_cum  + tt(monitoring_cum) + conflict_cum + tt(conflict_cum)
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data =  dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], tt = function(x, t, ...) x * log(t+20))
res_agreement_3 <-    coxph(Surv(time, treaty_yes) ~ agr_cum +  tt(agr_cum) + commission_cum + tt(commission_cum) + monitoring_cum  + tt(monitoring_cum) + conflict_cum + tt(conflict_cum)
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,], tt = function(x, t, ...) x * log(t+20))
res_commission_1 <-   coxph(Surv(time, commission) ~  
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_commission_2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_commission_3 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_monitoring_1 <-   coxph(Surv(time, monitoring) ~ 
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_monitoring_2 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_monitoring_3 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_conflict_1 <-     coxph(Surv(time, conflict) ~    
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_conflict_2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping + construction
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))
res_conflict_3 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping  + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat, tt = function(x, t, ...) x * log(t+20))

model_list <- list(res_agreement_1, res_agreement_2, res_agreement_3, 
                   res_commission_1, res_commission_2, res_commission_3, 
                   res_monitoring_1, res_monitoring_2, res_monitoring_3, 
                   res_conflict_1, res_conflict_2, res_conflict_3)

stargazer(model_list,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, all bilateral agreements",
          # column.labels   = c("Mechanisms", "Agreement"),
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

res_agreement_1 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_agreement_2 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_agreement_3 <-    coxph(Surv(time, treaty_yes) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_commission_1 <-   coxph(Surv(time, commission) ~  
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_commission_2 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_commission_3 <-   coxph(Surv(time, commission) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_monitoring_1 <-   coxph(Surv(time, monitoring) ~ 
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_monitoring_2 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_monitoring_3 <-   coxph(Surv(time, monitoring) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping +	fish + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_conflict_1 <-     coxph(Surv(time, conflict) ~    
                              + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_conflict_2 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping + construction
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])
res_conflict_3 <-     coxph(Surv(time, conflict) ~ agr_cum + commission_cum + monitoring_cum + conflict_cum 
                            +	pollution	+ shipping  + construction
                            + symmetry + bi_lingue  
                            ,
                            cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1945 | dyadicdat$year == 0,])

model_list <- list(res_agreement_1, res_agreement_2, res_agreement_3, 
                   res_commission_1, res_commission_2, res_commission_3, 
                   res_monitoring_1, res_monitoring_2, res_monitoring_3, 
                   res_conflict_1, res_conflict_2, res_conflict_3)

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


res_ext_1  <-   coxph(Surv(time, any_ext_yes) ~ 
                    agr_cum + commission_cum + monitoring_cum + conflict_cum 
                    ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_ext_2  <-   coxph(Surv(time, any_ext_yes) ~ 
                    agr_cum + commission_cum + monitoring_cum + conflict_cum 
                    +	pollution	+ shipping + fish + construction
                    ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])
res_ext_3  <-   coxph(Surv(time, any_ext_yes) ~ 
                    agr_cum + commission_cum + monitoring_cum + conflict_cum 
                    +	pollution	+ shipping + fish + construction
                    + symmetry + bi_lingue  
                    ,
                    cluster = dyad_id, data = dyadicdat[dyadicdat$year >= 1980 | dyadicdat$year == 0,])

model_list <- list(res_agr_1, res_agr_2, res_agr_3, res_mech_1, res_mech_2, res_mech_3, res_ext_1, res_ext_2, res_ext_3)

stargazer(model_list,
          type = "text",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, all bilateral agreements after 1980-01-01",
          column.separate = c(3, 1),
          dep.var.labels=c("Agreement (yes)", "Any Mechanism (yes)", "Any external Mechanism (yes)"),
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

res_agr_1 <-    coxph(Surv(time, treaty_yes) ~    
                        + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_agr_2 <-    coxph(Surv(time, treaty_yes) ~    
                        + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      +	pollution	+ shipping + construction
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_agr_3 <-    coxph(Surv(time, treaty_yes) ~    
                        + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      +	pollution	+ shipping + construction
                      + salience + symmetry + bi_lingue  
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_mech_1 <-   coxph(Surv(time, any_mech_yes) ~                    
                        + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_mech_2 <-   coxph(Surv(time, any_mech_yes) ~                    
                        + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_mech_3 <-   coxph(Surv(time, any_mech_yes) ~                    
                        + agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      +	pollution	+ shipping + construction
                      + salience + symmetry + bi_lingue  
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_ext_1  <-   coxph(Surv(time, any_ext_yes) ~ 
                        agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_ext_2  <-   coxph(Surv(time, any_ext_yes) ~ 
                        agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      +	pollution	+ shipping  + construction
                      ,
                      cluster = dyad_id, data = dyadicdat)
res_ext_3  <-   coxph(Surv(time, any_ext_yes) ~ 
                        agr_cum + commission_cum + monitoring_cum + conflict_cum 
                      +	pollution	+ shipping  + construction
                      + salience + symmetry + bi_lingue  
                      ,
                      cluster = dyad_id, data = dyadicdat)

model_list <- list(res_agr_1, res_agr_2, res_agr_3, res_mech_1, res_mech_2, res_mech_3, res_ext_1, res_ext_2, res_ext_3)

stargazer(model_list,
          type = "latex",
          booktaps = T,
          title = "Cox Proportional Hazard Models: risk set contiguity through land, all bilateral agreements",
          column.separate = c(3, 1),
          dep.var.labels=c("Agreement (yes)", "Any Mechanism (yes)", "Any external Mechanism (yes)"),
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

