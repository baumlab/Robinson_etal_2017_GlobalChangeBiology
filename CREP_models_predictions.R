
### Code accompanying:

## Robinson et al.
## Fishing degrades size structure of coral reef fish communities
## In review at Global Change Biology - May 2016
## Preprint posted on PeerJ - June 2016


## Script to make predictions using candidate model set
## Predictions are weighted by AIC weights, and made across the range of the focal predictor (spectrum slope or biomass)
## while holding all other covarites to their mean
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

## read in standardised t-values 
isd.dist.t<-read.csv("Analyses-final/model_tables/no-sharks/CRED_only_noshark_isd_distance_tvalues.csv")
isd.pop.t<-read.csv("Analyses-final/model_tables/no-sharks/CRED_only_noshark_isd_population_tvalues.csv")
biom.dist.t<-read.csv("Analyses-final/model_tables/no-sharks/CRED_only_noshark_biomass_distance_tvalues.csv")
biom.pop.t<-read.csv("Analyses-final/model_tables/no-sharks/CRED_only_noshark_biomass_population_tvalues.csv")


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
############### Compile model sets using MuMIn 'dredge' function ###############

#----------------------------------------------------------------------------
# Distance to market - ISD
dred.dist<-dredge(isd.dist,evaluate=TRUE, rank="AICc", trace =2,
			subset=(isltype_at.high && isltype_at.low))


#Subset for models within deltaAICc 7 to determine top model set 
dist.isd.top.models <- get.models(dred.dist, subset = delta < 7)
model.sel(dist.isd.top.models)

len.dist.isd<-length(dist.isd.top.models)
#recalc model weights for the top model set
top.weights.dist.isd <- dred.dist$weight[1:len.dist.isd]/sum(dred.dist$weight[1:len.dist.isd])
sum(top.weights.dist.isd)


dist.isd.preds<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(dist.isd.preds)<-isd.dist.t$Var
dist.isd.preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(dist.isd.preds.var)<-isd.dist.t$Var


## looping through all variables
for (y in 1:(nrow(isd.dist.t))){
  #Pull out the variable of interest
  var.temp <-isd.dist.t$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at it's mean
    newdata <- expand.grid(logdistance = 0, isltype_at.high = 0, isltype_at.low = 0,
                           prod = 0,min_SST = 0 ,  complexity = 0,logbathy=0,
                           loglandarea= 0, logreefarea75 = 0, test.var = seq(from = min(isd[,paste(var.temp)]), 
                                                                                   to = max(isd[,paste(var.temp)]), 
                                                                                   length.out = 100))
    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]

    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len.dist.isd , ncol = nrow(newdata))
    raw.pred <-  matrix(NA, nrow = len.dist.isd , ncol = nrow(newdata))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdata))
    for (j in 1:len.dist.isd){
      temp.mod <- get.models(dred.dist, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdata, re.form = NA) # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred* top.weights.dist.isd[j] # Weight the prediction by the corresponding model weight
      wt.pred[j,] <- pred.wt.temp
    }
    
    #Calculated the weighted mean prediction
    sum.pred.wt <- colSums(wt.pred)
    # Caculate weighted sample variance
    
    # calculate the weighted sample variance for each model
    for (i in 1:length(pred)){
    var.wt.pred[i]<-sum(top.weights.dist.isd*((raw.pred[,i]-sum.pred.wt[i])^2))
      }
    ## add to prediction dataframe
    dist.isd.preds[,y]<-sum.pred.wt
    dist.isd.preds.var[,y]<-as.vector(var.wt.pred)
}

write.csv(dist.isd.preds, file="weighted_preds_CRED_isd_dist.csv")
write.csv(dist.isd.preds.var, file="weighted_variance_CRED_isd_dist.csv")


#----------------------------------------------------------------------------
# Population density - ISD

### 25 minutes..........
ptm<-proc.time()
dred.pop<-dredge(isd.pop,evaluate=TRUE, rank="AICc", trace =2,  subset = (isltype_at.high && isltype_at.low))
proc.time() - ptm

#Subset for models within deltaAICc 7 to determine top model set 
pop.isd.top.models <- get.models(dred.pop, subset = delta < 7)
model.sel(pop.isd.top.models)



len.pop.isd<-length(pop.isd.top.models)
#recalc model weights for the top model set
top.weights.pop.isd <- dred.pop$weight[1:len.pop.isd]/sum(dred.pop$weight[1:len.pop.isd])
sum(top.weights.pop.isd)


pop.isd.preds<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(pop.isd.preds)<-isd.pop.t$Var
pop.isd.preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(pop.isd.preds.var)<-isd.pop.t$Var
## looping through all variables
for (y in 1:(nrow(isd.pop.t))){
  #Pull out the variable of interest
  var.temp <-isd.pop.t$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at it's mean
    newdata <- expand.grid(logpoparea = 0, isltype_at.high = 0, isltype_at.low = 0,
                           prod = 0,min_SST = 0 ,  complexity = 0,logbathy=0,
                           loglandarea= 0, logreefarea75 = 0, test.var = seq(from = min(isd[,paste(var.temp)]), 
                                                                                   to = max(isd[,paste(var.temp)]), 
                                                                                   length.out = 100))
    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]

    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len.pop.isd , ncol = nrow(newdata))
    raw.pred <-  matrix(NA, nrow = len.pop.isd , ncol = nrow(newdata))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdata))
    for (j in 1:len.pop.isd){
      temp.mod <- get.models(dred.pop, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdata, re.form = NA) # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred* top.weights.pop.isd[j] # Weight the prediction by the corresponding model weight
      wt.pred[j,] <- pred.wt.temp
    }
    
    #Calculated the weighted mean prediction
    sum.pred.wt <- colSums(wt.pred)
    # Caculate weighted sample variance
    
       # calculate the weighted sample variance for each model
    for (i in 1:length(pred)){
    var.wt.pred[i]<-sum(top.weights.pop.isd*((raw.pred[,i]-sum.pred.wt[i])^2))
      }
    ## add to prediction dataframe
    pop.isd.preds[,y]<-sum.pred.wt
    pop.isd.preds.var[,y]<-as.vector(var.wt.pred)
}

write.csv(pop.isd.preds, file="weighted_preds_CRED_isd_pop.csv")
write.csv(pop.isd.preds.var, file="weighted_variance_CRED_isd_pop.csv")



#----------------------------------------------------------------------------
# Distance to market - biomass

ptm<-proc.time()
dred.dist<-dredge(biom_dist,evaluate=TRUE, rank="AICc", trace =2, subset=(isltype_at.high && isltype_at.low))
proc.time() - ptm

#Subset for models within deltaAICc 7 to determine top model set 
dist.biom.top.models <- get.models(dred.dist, subset = delta < 7)
model.sel(dist.biom.top.models)



len.dist.biom<-length(dist.biom.top.models)
#recalc model weights for the top model set
top.weights.dist.biom <- dred.dist$weight[1:len.dist.biom]/sum(dred.dist$weight[1:len.dist.biom])
sum(top.weights.dist.biom)


dist.biom.preds<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(dist.biom.preds)<-biom.dist.t$Var
dist.biom.preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(dist.biom.preds.var)<-biom.dist.t$Var
## looping through all variables
for (y in 1:(nrow(biom.dist.t))){
  #Pull out the variable of interest
  var.temp <-biom.dist.t$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at it's mean
    newdata <- expand.grid(logdistance = 0, isltype_at.high = 0, isltype_at.low = 0,
                           prod = 0,min_SST = 0 ,  complexity = 0,logbathy=0,
                           loglandarea= 0, logreefarea75 = 0, test.var = seq(from = min(isd[,paste(var.temp)]), 
                                                                                   to = max(isd[,paste(var.temp)]), 
                                                                                   length.out = 100))
    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]

    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len.dist.biom , ncol = nrow(newdata))
    raw.pred <-  matrix(NA, nrow = len.dist.biom , ncol = nrow(newdata))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdata))
    for (j in 1:len.dist.biom){
      temp.mod <- get.models(dred.dist, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdata,type="response", re.form = NA) # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred* top.weights.dist.biom[j] # Weight the prediction by the corresponding model weight
      wt.pred[j,] <- pred.wt.temp
    }
    
    #Calculated the weighted mean prediction
    sum.pred.wt <- colSums(wt.pred)

    
       # calculate the weighted sample variance for each model
    for (i in 1:length(pred)){
    var.wt.pred[i]<-sum(top.weights.dist.biom*((raw.pred[,i]-sum.pred.wt[i])^2))
      }
    ## add to prediction dataframe
    dist.biom.preds[,y]<-sum.pred.wt
    dist.biom.preds.var[,y]<-as.vector(var.wt.pred)
}

write.csv(dist.biom.preds, file="weighted_preds_CRED_biom_dist.csv")
write.csv(dist.biom.preds.var, file="weighted_variance_CRED_biom_dist.csv")


#----------------------------------------------------------------------------
# Population density - biomass
ptm<-proc.time()
dred.pop<-dredge(biom_pop,evaluate=TRUE, rank="AICc", trace =2, subset=(isltype_at.high && isltype_at.low))
proc.time() - ptm


#Subset for models within deltaAICc 7 to determine top model set 
pop.biom.top.models <- get.models(dred.pop, subset = delta < 7)
model.sel(pop.biom.top.models)


len.pop.biom<-length(pop.biom.top.models)
#recalc model weights for the top model set
top.weights.pop.biom <- dred.pop$weight[1:len.pop.biom]/sum(dred.pop$weight[1:len.pop.biom])
sum(top.weights.pop.biom)


pop.biom.preds<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(pop.biom.preds)<-biom.pop.t$Var
pop.biom.preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=9))
colnames(pop.biom.preds.var)<-biom.pop.t$Var
## looping through all variables
for (y in 1:(nrow(biom.pop.t))){
  #Pull out the variable of interest
  var.temp <-biom.pop.t$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at it's mean
    newdata <- expand.grid(logpoparea = 0, isltype_at.high = 0, isltype_at.low = 0,
                           prod = 0,min_SST = 0 ,  complexity = 0,logbathy=0,
                           loglandarea= 0, logreefarea75 = 0, test.var = seq(from = min(isd[,paste(var.temp)]), 
                                                                                   to = max(isd[,paste(var.temp)]), 
                                                                                   length.out = 100))
    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]

    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len.pop.biom , ncol = nrow(newdata))
    raw.pred <-  matrix(NA, nrow = len.pop.biom , ncol = nrow(newdata))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdata))
    for (j in 1:len.pop.biom){
      temp.mod <- get.models(dred.pop, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdata,type="response", re.form = NA) # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred* top.weights.pop.biom[j] # Weight the prediction by the corresponding model weight
      wt.pred[j,] <- pred.wt.temp
    }
    
    #Calculated the weighted mean prediction
    sum.pred.wt <- colSums(wt.pred)

    
       # calculate the weighted sample variance for each model
    for (i in 1:length(pred)){
    var.wt.pred[i]<-sum(top.weights.pop.biom*((raw.pred[,i]-sum.pred.wt[i])^2))
      }
    ## add to prediction dataframe
    pop.biom.preds[,y]<-sum.pred.wt
    pop.biom.preds.var[,y]<-as.vector(var.wt.pred)
}

write.csv(pop.biom.preds, file="weighted_preds_CRED_biom_pop.csv")
write.csv(pop.biom.preds.var, file="weighted_variance_CRED_biom_pop.csv")



