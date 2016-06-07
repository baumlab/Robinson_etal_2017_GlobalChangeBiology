
### Code accompanying:

## Robinson et al.
## Fishing degrades size structure of coral reef fish communities
## In review at Global Change Biology - May 2016
## Preprint posted on PeerJ - June 2016


## Script to build predictive models for size spectra and biomass estimates
## Calculates standardised t-values and model tables for top-ranked model set
# Following recommendations of Cade BS (2015) Model averaging and muddled multimodel inference. Ecology, 96, 2370–2382.

### Load required packages
require(mgcv)
require(MuMIn)
require(lme4)
require(lmerTest)
options(na.action = "na.fail") 

## setwd and load dataframe with spectrum slope estimates, biomass estimates, and scaled explanatory covariates
setwd("git_repos/fishing-reefs-spectra")
isd<-read.csv("Analyses/results/CREP_ISD_biomass_nosharks.csv")

## set optimizer to be used in mixed models
ctrl <- lmeControl(opt='optim')

### Estimating R2 with this function:: 
## from Jarret Byrnes: http://stats.stackexchange.com/questions/111150/calculating-r2-in-mixed-models-using-nakagawa-schielzeths-2013-r2glmm-me
#Overall model R2
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

#----------------------------------------------------------------------------
############### BUILD MODELS ########################


######### Size spectrum models 

### Distance to market
ISD_dist_formula<-formula(exponent ~ logdistance   + isltype_at.low + isltype_at.high +
								 min_SST + prod + logbathy + complexity + logreefarea75 + 
							  loglandarea +
						  (1|OBS_YEAR))

isd.dist<-lmer(ISD_dist_formula, isd, REML=FALSE)
summary(isd.dist)
## check for collinearity
vif.mer(isd.dist)


### Human population density

ISD_pop_formula<-formula(exponent ~ logpoparea   + isltype_at.low + isltype_at.high + 
								 + min_SST + prod + logbathy + complexity + logreefarea75 + 
							  loglandarea   +
						    (1|OBS_YEAR), isd)

isd.pop<-lmer(ISD_pop_formula, isd, REML=FALSE)
summary(isd.pop)
## check for collinearity
vif.mer(isd.pop)	

######### Fish biomass models 

### Distance to market
biom_dist<-glmer(biomass ~ logdistance + isltype_at.high + isltype_at.low + prod + min_SST + complexity +
              logbathy + loglandarea + logreefarea75 +
                 (1|OBS_YEAR),na.action = "na.fail",
                  data=isd,family=Gamma(link=log))

summary(biom_dist)

### Human population density
biom_pop<-glmer(biomass ~ logpoparea  + isltype_at.high + isltype_at.low + prod + min_SST + complexity +
              logbathy + loglandarea + logreefarea75 +
                 (1|OBS_YEAR),na.action = "na.fail",
                  data=isd,family=Gamma(link=log))

summary(biom_pop)


#----------------------------------------------------------------------------
############### Compile model candidate set ranked by AICc ###############
############### Compile model sets using MuMIn 'dredge' function ###############


#----------------------------------------------------------------------------
############### Distance to market  ~ ISD
dred.dist<-dredge(isd.dist,evaluate=TRUE, rank="AICc", trace =2,
			subset=(isltype_at.high && isltype_at.low))


#Subset for models within deltaAICc 7 to determine top model set 
dist.isd.top.models <- get.models(dred.dist, subset = delta < 7)
model.sel(dist.isd.top.models)

## Now estimate weighted t-statistics

#Model average absolute value of t-statistics for measure of variable importance
dist.isd.imp <- as.data.frame(matrix(NA, nrow = (ncol(dred.dist) - 6), ncol = 4))
colnames(dist.isd.imp) <- c("Var", "RI.t.abs", "RI.t.ratio",  "var.t")

#Loop through all the variables in the model
for (i in 1:(ncol(dred.dist) - 6)){
  var.temp <- colnames(dred.dist)[i+1]
  dist.isd.imp[i,"Var"] <- var.temp
  
  #Subset for only models for a given predictor, don't recalulate the model weights or deltaAICcs
  dred.dist_x <- dred.dist[i = !is.na(dred.dist[,paste(var.temp)]), j = 1:ncol(dred.dist), recalc.weights =FALSE, recalc.delta = FALSE]
  
  #Pull out the t-statistics for the given variables for each model in which it appears
  #Put all values in a data frame
  t.values <- as.data.frame(matrix(NA, nrow = nrow(dred.dist_x), ncol = 4))
  colnames(t.values) <- c('varx',  "model.wt", 'imp.t.abs', "imp.t.ratio")
  for (l in 1:nrow(dred.dist_x)){
    #put in variable name
    t.values[l,"varx"] <- var.temp
    
    #Pull out one model for the set that includes the variable of interest
    temp_mod <- get.models(dred.dist, which(!is.na(dred.dist[,paste(var.temp)]))[l])[[1]]
    
    #Pull out the model weight for this model
    wt.temp <- dred.dist$weight[which(!is.na(dred.dist[,paste(var.temp)]))][l]
    t.values[l,"model.wt"] <- wt.temp
    
    #abs value of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])
    t.values[l, "imp.t.abs"] <-  RI.x.temp 
    
    #ratios of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])/max(abs(coef(summary(temp_mod))[-1,3]))
    t.values[l, "imp.t.ratio"] <-  RI.x.temp 
  }
  #AICc weighted average of t-statistics for variable of interest
  avgRI.x.abs <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[2]*x[1]))
  dist.isd.imp[i,"RI.t.abs"] <- avgRI.x.abs
  
   #AICc weighted average ratio of t-statistics for variable of interest
  avgRI.x <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[3]*x[1]))
  dist.isd.imp[i,"RI.t.ratio"] <- avgRI.x

  ### variance
  varT.x<-sum(t.values[,2]*((t.values[,3]-avgRI.x.abs)^2))
  dist.isd.imp[i, "var.t"] <- varT.x
}


## get r-squareds for ISD _ distance models
top.dist<-get.models(dred.dist, subset=delta<7)
dist.isd.top.models<-data.frame(subset(dred.dist, delta<7))	
models<-rownames(dist.isd.top.models)
r2.isd.dist<-data.frame(ncol=1, models)

for (i in 1:length(top.dist)){

r2.isd.dist$r2[i]<-as.numeric(r.squaredGLMM(top.dist[[models[i]]])[1])
}

dist.isd.top.models$R2<-r2.isd.dist$r2



#----------------------------------------------------------------------------
############### Human population  ~ ISD


ptm<-proc.time()
dred.pop<-dredge(isd.pop,evaluate=TRUE, rank="AICc", trace =2,  subset = (isltype_at.high && isltype_at.low))
proc.time() - ptm

#Subset for models within deltaAICc 7 to determine top model set 
pop.isd.top.models <- get.models(dred.pop, subset = delta < 7)
model.sel(pop.isd.top.models)

#Model average absolute value of t-statistics for measure of variable importance
pop.isd.imp <- as.data.frame(matrix(NA, nrow = (ncol(dred.pop) - 6), ncol = 3))
colnames(pop.isd.imp) <- c("Var", "RI.t.abs", "RI.t.ratio")

#Loop through all the variables in the model
pop.isd.imp <- as.data.frame(matrix(NA, nrow = (ncol(dred.pop) - 6), ncol = 4))
colnames(pop.isd.imp) <- c("Var", "RI.t.abs", "RI.t.ratio",  "var.t")

#Loop through all the variables in the model
for (i in 1:(ncol(dred.pop) - 6)){
  var.temp <- colnames(dred.pop)[i+1]
  pop.isd.imp[i,"Var"] <- var.temp
  
  #Subset for only models for a given predictor, don't recalulate the model weights or deltaAICcs
  dred.pop_x <- dred.pop[i = !is.na(dred.pop[,paste(var.temp)]), j = 1:ncol(dred.pop), recalc.weights =FALSE, recalc.delta = FALSE]
  
  #Pull out the t-statistics for the given variables for each model in which it appears
  #Put all values in a data frame
  t.values <- as.data.frame(matrix(NA, nrow = nrow(dred.pop_x), ncol = 4))
  colnames(t.values) <- c('varx',  "model.wt", 'imp.t.abs', "imp.t.ratio")
  for (l in 1:nrow(dred.pop_x)){
    #put in variable name
    t.values[l,"varx"] <- var.temp
    
    #Pull out one model for the set that includes the variable of interest
    temp_mod <- get.models(dred.pop, which(!is.na(dred.pop[,paste(var.temp)]))[l])[[1]]
    
    #Pull out the model weight for this model
    wt.temp <- dred.pop$weight[which(!is.na(dred.pop[,paste(var.temp)]))][l]
    t.values[l,"model.wt"] <- wt.temp
    
    #abs value of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])
    t.values[l, "imp.t.abs"] <-  RI.x.temp 
    
    #ratios of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])/max(abs(coef(summary(temp_mod))[-1,3]))
    t.values[l, "imp.t.ratio"] <-  RI.x.temp 
  }
  #AICc weighted average of t-statistics for variable of interest
  avgRI.x.abs <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[2]*x[1]))
  pop.isd.imp[i,"RI.t.abs"] <- avgRI.x.abs
  
  #AICc weighted average ratio of t-statistics for variable of interest
  avgRI.x <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[3]*x[1]))
  pop.isd.imp[i,"RI.t.ratio"] <- avgRI.x

  
  ### variance
  varT.x<-sum(t.values[,2]*((t.values[,3]-avgRI.x.abs)^2))
  pop.isd.imp[i, "var.t"] <- varT.x
}


## get r-squareds for ISD _ pop models
top.pop<-get.models(dred.pop, subset=delta<7)
pop.isd.top.models<-data.frame(subset(dred.pop, delta<7))	
models<-rownames(pop.isd.top.models)

r2.isd.pop<-data.frame(ncol=1, models)

for (i in 1:length(top.pop)){

r2.isd.pop$r2[i]<-as.numeric(r.squaredGLMM(top.pop[[models[i]]])[1])
}

pop.isd.top.models$R2<-r2.isd.pop$r2

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Now for biomass models
# Compile candidate model set

ptm<-proc.time()
dred.dist<-dredge(biom_dist,evaluate=TRUE, rank="AICc", trace =2, subset=(isltype_at.high && isltype_at.low))
proc.time() - ptm

ptm<-proc.time()
dred.pop<-dredge(biom_pop,evaluate=TRUE, rank="AICc", trace =2, subset=(isltype_at.high && isltype_at.low))
proc.time() - ptm


#----------------------------------------------------------------------------
############### Distance to market  ~ biomass t-values

#Subset for models within deltaAICc 7 to determine top model set 
dist.biom.top.models <- get.models(dred.dist, subset = delta < 7)
model.sel(dist.biom.top.models)

#Model average absolute value of t-statistics for measure of variable importance
dist.biom.imp <- as.data.frame(matrix(NA, nrow = (ncol(dred.dist) - 6), ncol = 3))
colnames(dist.biom.imp) <- c("Var", "RI.t.abs", "RI.t.ratio")

#Loop through all the variables in the model
dist.biom.imp <- as.data.frame(matrix(NA, nrow = (ncol(dred.dist) - 6), ncol = 4))
colnames(dist.biom.imp) <- c("Var", "RI.t.abs", "RI.t.ratio",  "var.t")

#Loop through all the variables in the model
for (i in 1:(ncol(dred.dist) - 6)){
  var.temp <- colnames(dred.dist)[i+1]
  dist.biom.imp[i,"Var"] <- var.temp
  
  #Subset for only models for a given predictor, don't recalulate the model weights or deltaAICcs
  dred.dist_x <- dred.dist[i = !is.na(dred.dist[,paste(var.temp)]), j = 1:ncol(dred.dist), recalc.weights =FALSE, recalc.delta = FALSE]
  
  #Pull out the t-statistics for the given variables for each model in which it appears
  #Put all values in a data frame
  t.values <- as.data.frame(matrix(NA, nrow = nrow(dred.dist_x), ncol = 4))
  colnames(t.values) <- c('varx',  "model.wt", 'imp.t.abs', "imp.t.ratio")
  for (l in 1:nrow(dred.dist_x)){
    #put in variable name
    t.values[l,"varx"] <- var.temp
    
    #Pull out one model for the set that includes the variable of interest
    temp_mod <- get.models(dred.dist, which(!is.na(dred.dist[,paste(var.temp)]))[l])[[1]]
    
    #Pull out the model weight for this model
    wt.temp <- dred.dist$weight[which(!is.na(dred.dist[,paste(var.temp)]))][l]
    t.values[l,"model.wt"] <- wt.temp
    
    #abs value of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])
    t.values[l, "imp.t.abs"] <-  RI.x.temp 
    
    #ratios of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])/max(abs(coef(summary(temp_mod))[-1,3]))
    t.values[l, "imp.t.ratio"] <-  RI.x.temp 
  }
  #AICc weighted average of t-statistics for variable of interest
  avgRI.x.abs <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[2]*x[1]))
  dist.biom.imp[i,"RI.t.abs"] <- avgRI.x.abs
  
  #AICc weighted average ratio of t-statistics for variable of interest
  avgRI.x <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[3]*x[1]))
  dist.biom.imp[i,"RI.t.ratio"] <- avgRI.x

  
  ### variance
  varT.x<-sum(t.values[,2]*((t.values[,3]-avgRI.x.abs)^2))
  dist.biom.imp[i, "var.t"] <- varT.x
}

## get r-squareds for biomass _ distance models

top.dist<-get.models(dred.dist, subset=delta<7)
	dist.biom.top.models<-data.frame(subset(dred.dist, delta<7))	
models<-rownames(dist.biom.top.models)

r2.biom.dist<-data.frame(ncol=1, models)

for (i in 1:length(top.dist)){

r2.biom.dist$r2[i]<-r2.corr.mer(top.dist[[models[i]]])

}

dist.biom.top.models$R2<-r2.biom.dist$r2


#----------------------------------------------------------------------------
############### Population density  ~ biomass t-values

#Subset for models within deltaAICc 7 to determine top model set 
pop.biom.top.models <- get.models(dred.pop, subset = delta < 7)
model.sel(pop.biom.top.models)


#Model average absolute value of t-statistics for measure of variable importance

#Loop through all the variables in the model
pop.biom.imp <- as.data.frame(matrix(NA, nrow = (ncol(dred.pop) - 6), ncol = 4))
colnames(pop.biom.imp) <- c("Var", "RI.t.abs", "RI.t.ratio",  "var.t")

#Loop through all the variables in the model
for (i in 1:(ncol(dred.pop) - 6)){
  var.temp <- colnames(dred.pop)[i+1]
  pop.biom.imp[i,"Var"] <- var.temp
  
  #Subset for only models for a given predictor, don't recalulate the model weights or deltaAICcs
  dred.pop_x <- dred.pop[i = !is.na(dred.pop[,paste(var.temp)]), j = 1:ncol(dred.pop), recalc.weights =FALSE, recalc.delta = FALSE]
  
  #Pull out the t-statistics for the given variables for each model in which it appears
  #Put all values in a data frame
  t.values <- as.data.frame(matrix(NA, nrow = nrow(dred.pop_x), ncol = 4))
  colnames(t.values) <- c('varx',  "model.wt", 'imp.t.abs', "imp.t.ratio")
  for (l in 1:nrow(dred.pop_x)){
    #put in variable name
    t.values[l,"varx"] <- var.temp
    
    #Pull out one model for the set that includes the variable of interest
    temp_mod <- get.models(dred.pop, which(!is.na(dred.pop[,paste(var.temp)]))[l])[[1]]
    
    #Pull out the model weight for this model
    wt.temp <- dred.pop$weight[which(!is.na(dred.pop[,paste(var.temp)]))][l]
    t.values[l,"model.wt"] <- wt.temp
    
    #abs value of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])
    t.values[l, "imp.t.abs"] <-  RI.x.temp 
    
    #ratios of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])/max(abs(coef(summary(temp_mod))[-1,3]))
    t.values[l, "imp.t.ratio"] <-  RI.x.temp 
  }
  #AICc weighted average of t-statistics for variable of interest
  avgRI.x.abs <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[2]*x[1]))
  pop.biom.imp[i,"RI.t.abs"] <- avgRI.x.abs
  
  #AICc weighted average ratio of t-statistics for variable of interest
  avgRI.x <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[3]*x[1]))
  pop.biom.imp[i,"RI.t.ratio"] <- avgRI.x

  
  ### variance
  varT.x<-sum(t.values[,2]*((t.values[,3]-avgRI.x.abs)^2))
  pop.biom.imp[i, "var.t"] <- varT.x
}


## get r-squareds for biomass _ pop models
top.pop<-get.models(dred.pop, subset=delta<7)
pop.biom.top.models<-data.frame(subset(dred.pop, delta<7))	
models<-rownames(pop.biom.top.models)
head(pop.biom.top.models)

r2.biom.pop<-data.frame(ncol=1, models)

for (i in 1:length(top.pop)){

r2.biom.pop$r2[i]<-r2.corr.mer(top.pop[[models[i]]])

}
pop.biom.top.models$R2<-r2.biom.pop$r2


#------------------------------------------------------------
############## Save tvalues and top-ranked models as csv files #################

write.csv(dist.isd.imp, file="CREP_only_ISD_distance_tvalues.csv")
write.csv(dist.isd.top.models, file="CREP_only_ISD_distance_topAIC.csv")

write.csv(pop.isd.imp, file="CREP_only_ISD_population_tvalues.csv")
write.csv(pop.isd.top.models, file="CREP_only_ISD_population_topAIC.csv")

write.csv(dist.biom.imp, file="CREP_only_biomass_distance_tvalues.csv")
write.csv(dist.biom.top.models, file="CREP_only_biomass_distance_topAIC.csv")

write.csv(pop.biom.imp, file="CREP_only_biomass_population_tvalues.csv")
write.csv(pop.biom.top.models, file="CREP_only_biomass_population_topAIC.csv")