# install.packages('tsna',repos='http://cran.us.r-project.org', 
#                  dependencies=TRUE)
# install.packages('ndtv',repos='http://cran.us.r-project.org', 
#                  dependencies=TRUE)
library(tsna)
library(ndtv) # also loads animation and networkDynamic
library(igraph)
library(dplyr)

# load dyadic data
load("/Volumes/Transcend/Uni/Zivi Eawag/Script/Konkordats/ModelInput/dyadicdat.RData")

dyadicdat$date <- as.Date(dyadicdat$time, origin = "1970-01-01")

dyadicdat <- dyadicdat %>% arrange(date)
dyadicdat$canton1 <- as.factor(dyadicdat$canton1)
dyadicdat$canton2 <- as.factor(dyadicdat$canton2)
dyadicdat$ID_SM <- as.factor(dyadicdat$ID_SM)

library(forcats) # understand why changing levels changes values of factors, too stupid to understand just now
dyadicdat$canton1 <- fct_expand(dyadicdat$canton1, levels = sort(unique(c(as.character(dyadicdat$canton1), as.character(dyadicdat$canton2)))))
dyadicdat$canton2 <- fct_expand(dyadicdat$canton2, levels = sort(unique(c(as.character(dyadicdat$canton1), as.character(dyadicdat$canton2)))))

dyadicdat <- dyadicdat %>% filter(bilateral == 1)
  
lst <- list()

for (i in as.character(unique(dyadicdat$date))){
 temp_mat1 <- as.matrix(table(dyadicdat$canton1[dyadicdat$date <= i], dyadicdat$ID_SM[dyadicdat$date <= i]))
 temp_mat2 <- as.matrix(table(dyadicdat$canton2[dyadicdat$date <= i], dyadicdat$ID_SM[dyadicdat$date <= i]))
 temp_mat1 <- temp_mat1[sort(rownames(temp_mat1)), sort(colnames(temp_mat1))]
 temp_mat2 <- temp_mat2[sort(rownames(temp_mat2)), sort(colnames(temp_mat2))]
 temp_mat <- temp_mat1 + temp_mat2
 lst[[i]] <- temp_mat %*% t(temp_mat)
 diag(lst[[i]]) <- 0
 lst[[i]] <- as.network(lst[[i]], 
                        vertex.attrnames = rownames(lst[[i]]),
                        directed = F,
                        ignore.eval = F,
                        names.eval='weight')
}

summary( lst[[i]])

dyn <- networkDynamic(network.list = lst, 
                      vertex.TEA.names=TRUE,
                      create.TEAs=TRUE, 
                      edge.TEA.names=c('weight','type'))

get.edge.attribute.active(dyn,'weight',at=2)

render.d3movie(dyn, 
               displaylabels=TRUE,
               edge.lwd='weight',
               output.mode = 'htmlWidget')

timeline(dyn)

proximity.timeline(dyn, 
                   default.dist=6,
                   mode='sammon',
                   labels.at=17,
                   vertex.cex=4)


