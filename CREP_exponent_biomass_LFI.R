
### Code accompanying:

## Robinson et al.
## Fishing degrades size structure of coral reef fish communities
## In review at Global Change Biology - May 2016
## Preprint posted on PeerJ - June 2016


## Script to make predictions using candidate model set

### Load required packages
require(mgcv)
require(MuMIn)
require(plotrix)
require(lme4)
require(lmerTest)
require(RColorBrewer)
options(na.action = "na.fail") 

## setwd and load dataframe with spectrum slope estimates, biomass estimates,LFI, and scaled explanatory covariates
setwd("git_repos/fishing-reefs-spectra")
isd<-read.csv("Analyses/results/CREP_ISD_biomass_nosharks.csv")


#----------------------------------------------------------------------------
##########################################
### exponent ~ biomass #######################
##########################################

biom.mod<-lmer(exponent ~ biomass*STATE + (1 | OBS_YEAR), isd, REML=TRUE)
summary(biom.mod)
r.squaredGLMM(biom.mod)
plot(resid(biom.mod))
hist(resid(biom.mod))
plot(resid(biom.mod), isd$biomass)


# get slopes + intercepts
x <- seq(min(isd$biomass), max(isd$biomass), length.out = 100)
newdat <- data.frame(expand.grid(biomass=x, STATE=c("Remote","Populated"), exponent=0))

mm <- model.matrix(terms(biom.mod),newdat)
newdat$exponent <- predict(biom.mod,newdat,re.form=NA)


pvar1 <- diag(mm %*% tcrossprod(vcov(biom.mod),mm))
tvar1 <- pvar1+VarCorr(biom.mod)$OBS_YEAR[1]  ## must be adapted for more complex models
cmult <- 2 ## could use 1.96
newdat <- data.frame(
    newdat
    # fixed effects uncertainty, min (plo) and max (phi)
    , plo = newdat$exponent-cmult*sqrt(pvar1)
    , phi = newdat$exponent+cmult*sqrt(pvar1)
    # fixed and random effects uncertainty, min (tlo) and max (thi)
    , tlo = newdat$exponent-cmult*sqrt(tvar1)
    , thi = newdat$exponent+cmult*sqrt(tvar1)
)


######################## FIGURE 4 ##########################################

cols<-c(brewer.pal(6,"Set1"), '#ffffb2', '#bd0026')
### estimate number of sites per country, and assign colours based on population density
predict<-read.csv("predictors_merged_nSPC_SPC_islandlevel.csv")
isd$pop.plot<-log10(predict$distance[match( isd$ISLAND, predict$Site)])
labels<-aggregate(ISLAND ~ lat + lon + REGION + STATE+ pop.plot, isd, unique)

## create vector for colours in map
col_seq<-data.frame(seq(min((labels$pop.plot)), max((labels$pop.plot))+0.01, 0.01))
col_seq$cols<-colorRampPalette(c(cols[8], cols[6]))(dim(col_seq)[1])
colnames(col_seq)<-c("Pop", "col")
col_seq$Pop<-round(col_seq$Pop, 2)
isd$cols<-col_seq$col[match(round(isd$pop.plot,2), col_seq$Pop)]



popul<-isd[isd$STATE=="Populated",]
remote<-isd[isd$STATE=="Remote",]
par(mar=c(5,5,4,2), mgp=c(3,1,0), cex.axis=1.1, mfrow=c(1,1), tck=0.03, cex.axis=1.4)
plot(popul$biomass, popul$exponent, xlim=c(0,2000), ylim=c(-2,-1), col="transparent",
						 axes=FALSE, xlab="", ylab="")
axis(1, at =seq(0, 2000, 500))
axis(2, at=seq(-2, -1, 0.2))
clip(min(popul$biomass)-50, max(popul$biomass)+50, min(popul$exponent)-0.1, max(popul$exponent)+0.35)
points(popul$biomass, popul$exponent, pch=21, bg=popul$col, cex=1.8)
pop.dat<-newdat[newdat$STATE=="Populated",]
lines(pop.dat$biomass, pop.dat$exponent, col=col_seq$col[50], lwd=2.5)
lines(pop.dat$biomass, pop.dat$plo, col=col_seq$col[50], lwd=2, lty=2)
lines(pop.dat$biomass, pop.dat$phi, col=col_seq$col[50], lwd=2, lty=2)

## remote plots
clip(min(remote$biomass)-50, max(remote$biomass)+50, min(remote$exponent)-0.1, max(remote$exponent)+0.35)
points(remote$biomass, remote$exponent, pch=22, bg=remote$cols, cex=1.8)
rem.dat<-newdat[newdat$STATE=="Remote",]
lines(rem.dat$biomass, rem.dat$exponent, col=col_seq$col[220], lwd=2.5)
lines(rem.dat$biomass, rem.dat$plo, col=col_seq$col[220], lwd=2, lty=2)
lines(rem.dat$biomass, rem.dat$phi, col=col_seq$col[220], lwd=2, lty=2)

mtext(1, text=expression(paste("Reef fish biomass (kg ha"^"-1", ")")), line=3, cex=1.5)
mtext("Size spectrum exponent", side=2, line=3, cex=1.5)

par(xpd=FALSE)
	color.legend(1900, -1.5, 2000, -1.1, legend=seq(0, 3, 1), rect.col=(col_seq$col),align="lt", cex=1.2, gradient="y")
	text(1850, -1.3, expression(paste("Log"[10], "proximity to market")), srt=90, cex=1.2)




#----------------------------------------------------------------------------
### Supplementary analysis: exponent ~ LFI relationship


##########################################
### exponent ~ LFI #######################
##########################################

lfi.mod<-lmer(exponent ~ LFI*STATE + (1 | OBS_YEAR), isd, REML=TRUE)
summary(lfi.mod)
r.squaredGLMM(lfi.mod)
plot(resid(lfi.mod))
hist(resid(lfi.mod))
plot(resid(lfi.mod), isd$lfi)

# # get slopes + intercepts

x <- seq(min(isd$LFI), max(isd$LFI), length.out = 100)
lfidat <- data.frame(expand.grid(LFI=x, STATE=c("Remote","Populated"), exponent=0))

mm <- model.matrix(terms(lfi.mod),lfidat)
lfidat$exponent <- predict(lfi.mod,lfidat,re.form=NA)


pvar1 <- diag(mm %*% tcrossprod(vcov(lfi.mod),mm))
tvar1 <- pvar1+VarCorr(lfi.mod)$OBS_YEAR[1]  ## must be adapted for more complex models
cmult <- 2 ## could use 1.96
lfidat <- data.frame(
    lfidat
    # fixed effects uncertainty, min (plo) and max (phi)
    , plo = lfidat$exponent-cmult*sqrt(pvar1)
    , phi = lfidat$exponent+cmult*sqrt(pvar1)
    # fixed and random effects uncertainty, min (tlo) and max (thi)
    , tlo = lfidat$exponent-cmult*sqrt(tvar1)
    , thi = lfidat$exponent+cmult*sqrt(tvar1)
)



######################## FIGURE S1 ##########################################


par(mar=c(5,5,4,2), mgp=c(3,1,0), cex.axis=1.4, mfrow=c(1,1), tck=0.03)
plot(isd$exponent ~ isd$LFI, ylim=c(-2,-1), xlim=c(0.05,0.8), col="transparent",
						 axes=FALSE, xlab="", ylab="")
axis(2, at =seq(-2, -1, 0.1))
axis(1, at=seq(0, 1, 0.2))
clip(min(popul$LFI)-0.1, max(popul$LFI)+0.1, min(popul$exponent)-0.1,max(popul$exponent)+0.1)
pop.lfi<-lfidat[lfidat$STATE=="Populated",]
lines(pop.lfi$LFI, pop.lfi$exponent, col=col_seq$col[50], lwd=2.5)
lines(pop.lfi$LFI, pop.lfi$plo, col=col_seq$col[50], lwd=2, lty=2)
lines(pop.lfi$LFI, pop.lfi$phi, col=col_seq$col[50], lwd=2, lty=2)
points(popul$LFI, popul$exponent, pch=22, bg=popul$cols, cex=1.8)

clip(min(remote$LFI)-0.1, max(remote$LFI)+0.1, min(remote$exponent)-0.1,max(remote$exponent)+0.1)
rem.lfi<-lfidat[lfidat$STATE=="Remote",]
lines(rem.lfi$LFI, rem.lfi$exponent, col=col_seq$col[220], lwd=2.5)
lines(rem.lfi$LFI, rem.lfi$plo, col=col_seq$col[220], lwd=2, lty=2)
lines(rem.lfi$LFI, rem.lfi$phi, col=col_seq$col[220], lwd=2, lty=2)
points(remote$LFI, remote$exponent, pch=22, bg=remote$cols, cex=1.8)


mtext("Size spectrum exponent", side=2, line=3, cex=1.5)
mtext("Large fish indicator", side=1, line=3, cex=1.5)

require(plotrix)
par(xpd=FALSE)
	color.legend(0.8, -1.5, 0.82, -1.1, legend=seq(0, 3, 1), rect.col=(col_seq$col),align="lt", cex=1.2, gradient="y")
	text(0.78, -1.3, expression(paste("Log"[10], "proximity to market")), srt=90, cex=1.2)



