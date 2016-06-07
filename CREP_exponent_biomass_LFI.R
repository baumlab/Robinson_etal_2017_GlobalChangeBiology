
### Code accompanying:

## Robinson et al.
## Fishing degrades size structure of coral reef fish communities
## In review at Global Change Biology - May 2016
## Preprint posted on PeerJ - June 2016


## Script to make predictions using candidate model set

### Load required packages
require(mgcv)
require(MuMIn)
require(lme4)
require(lmerTest)
options(na.action = "na.fail") 

## setwd and load dataframe with spectrum slope estimates, biomass estimates,LFI, and scaled explanatory covariates
setwd("git_repos/fishing-reefs-spectra")
isd<-read.csv("Analyses/results/CREP_ISD_biomass_nosharks.csv")


#----------------------------------------------------------------------------
# Exponent ~ biomass models 

biom.remote.mod<-lmer(exponent ~ biomass + (1 | OBS_YEAR), isd[isd$STATE=="Remote",], REML=TRUE)
summary(biom.remote.mod)
r.squaredGLMM(biom.remote.mod)

# examine model fitting
plot(resid(biom.remote.mod))
hist(resid(biom.remote.mod))
plot(resid(biom.remote.mod), isd$biomass[isd$STATE=="Remote"])

# get slopes + intercepts
rem.int<-coef(summary(biom.remote.mod))[ , "Estimate"][1]
rem.slope<-coef(summary(biom.remote.mod))[ , "Estimate"][2]
rem.ci.int<-confint(biom.remote.mod)[3,]
rem.ci.slope<-confint(biom.remote.mod)[4,]

#----------------------------------------------------------------------------
## repeat for populated islands 
biom.populated.mod<-lmer(exponent ~ biomass + (1 | OBS_YEAR), isd[isd$STATE=="Populated",], REML=TRUE)
summary(biom.populated.mod)
r.squaredGLMM(biom.populated.mod)

# examine model fitting
plot(resid(biom.populated.mod))
hist(resid(biom.populated.mod))
plot(resid(biom.populated.mod), isd$biomass[isd$STATE=="Populated"])

# get slopes + intercepts
pop.int<-coef(summary(biom.populated.mod))[ , "Estimate"][1]
pop.slope<-coef(summary(biom.populated.mod))[ , "Estimate"][2]
pop.ci.int<-confint(biom.populated.mod)[3,]
pop.ci.slope<-confint(biom.populated.mod)[4,]

#----------------------------------------------------------------------------
### Supplementary analysis: exponent ~ LFI relationship

# Populated islands
lfi.pop.mod<-lmer(exponent ~ LFI + (1 | OBS_YEAR), isd[isd$STATE=="Populated",], REML=TRUE)
summary(lfi.pop.mod)
r.squaredGLMM(lfi.pop.mod)

# examine model fitting
plot(resid(lfi.pop.mod))
hist(resid(lfi.pop.mod))
plot(resid(lfi.pop.mod), isd$lfi[isd$STATE=="Populated"])

# get slopes + intercepts
lfi.pop.int<-coef(summary(lfi.pop.mod))[ , "Estimate"][1]
lfi.pop.slope<-coef(summary(lfi.pop.mod))[ , "Estimate"][2]
lfi.pop.ci.int<-confint(lfi.pop.mod)[3,]
lfi.pop.ci.slope<-confint(lfi.pop.mod)[4,]

## Remote Islands
lfi.remote.mod<-lmer(exponent ~ LFI + (1 | OBS_YEAR), isd[isd$STATE=="Remote",], REML=TRUE)
summary(lfi.remote.mod)
r.squaredGLMM(lfi.remote.mod)

# examine model fitting
plot(resid(lfi.remote.mod))
hist(resid(lfi.remote.mod))
plot(resid(lfi.remote.mod), isd$lfi[isd$STATE=="Populated"])

# get slopes + intercepts
lfi.remote.int<-coef(summary(lfi.remote.mod))[ , "Estimate"][1]
lfi.remote.slope<-coef(summary(lfi.remote.mod))[ , "Estimate"][2]
lfi.remote.ci.int<-confint(lfi.remote.mod)[3,]
lfi.remote.ci.slope<-confint(lfi.remote.mod)[4,]
