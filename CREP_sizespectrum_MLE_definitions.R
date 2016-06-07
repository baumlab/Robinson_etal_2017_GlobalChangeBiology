#### Code accompanying:

## Robinson et al.
## Fishing degrades size structure of coral reef fish communities
## In review at Global Change Biology - May 2016
## Preprint posted on PeerJ - June 2016


## Script to estimate size spectra from UVC reef fish data


#----------------------------------------------------------------------------
### R function definitions ###


# Negative log-likelihood function for bounded power law
#  distribution (PLB model). From equation (A.23) in Edwards (2011),
#  Ecology 92:1247-1257.


powlawbound = function(mu)
 

  {              #  Negative log-likelihood is:
  if(mu == 1)    # See equation (A.25) in Ecology paper.
    { n * log( log( b/a)) + sumlogx}
  else
    {  -n * log( (mu - 1) / (a^(1 - mu) - b^(1 - mu))) + mu * sumlogx}
                 # This works for mu<1 as well as mu>1
  }


### convert abundance ~ size dataframe into vector of body sizes
dataorder<-function(x, catch, size) data.frame(rep(x[,size], times=x[,catch] ))

## calculate MLE parameters from abundance ~ size data
mle_params<-function(df, size, catch){
    #df$mass[df$mass==0]<-0.001
  
df<-dataorder(df, size=size, catch=catch)
rawvalsall<<-df[,1] 
x <<-df[,1]
a <<- min(rawvalsall)
b <<- max(rawvalsall)
n <<- length(rawvalsall)

sumlogx <<- sum(log(rawvalsall)) # for powlawbound function
sumx<<-sum(rawvalsall)
}

## numerically optimise powlawbound to find MLE of ISD exponent (i.e. size spectrum slope)

plb_mlefunc<- function(df, size, catch) {
  mle_params(df, size, catch)

outpowlaw<- nlm(powlawbound, mu.start2vec[1]) 

mumle<- outpowlaw$estimate
Cpowlaw<- (mumle - 1) / (a^(1-mumle) - b^(1-mumle)) 
negloglikmlepowlaw<- outpowlaw$minimum

AICpowlaw<- 2 * negloglikmlepowlaw + 2*numpar

# Now for 95% CI of MLE:    
muvarynegloglik<- vector()       # negative log lik for each mu value
  for(i in 1:length(muvec))
    {
      muvarynegloglik[i] = powlawbound(muvec[i])
    }
critval1<- negloglikmlepowlaw + qchisq(0.95,1)/2
                    # 1 dof. Hilborn and Mangel (1997) p163.
muin95<- muvec[ muvarynegloglik < critval1] 
                    # mu values definitely in 95% confidence interval
minmuin95<- min(muin95)
maxmuin95<- max(muin95)

coefs<-data.frame(mumle, Cpowlaw, minmuin95, maxmuin95, AICpowlaw)
colnames(coefs)<-c("exponent_num", "constant", "CImin", "CImax", "AIC")
coefs
}

### End R function definitions ###
#----------------------------------------------------------------------------
### R parameter definitions ###

numpar = 3    # num of pars for bounded model

#  They are the starting value for optimisation
#  and the start, end, and increments to use
#  to generate values to calculate 95%
#  confidence intervals.
muvecstartvec    = vector()
muvecendvec      = vector()
muvecincvec      = vector()
mu.start2vec     = vector()

# Initial starting values for finding MLE using nlm. 
mu.start2vec[1] = c(0.5) 
# range for mu bounded.
muvecstartvec[1]  = c(0.5)
muvecendvec[1]  = c(2)
muvecincvec[1]  = c(0.001)

#  Range of mu to try for PLB model, can't be <=1
muvec = seq(muvecstartvec[1], muvecendvec[1], muvecincvec[1])

### End of R parameter definitions ###

#----------------------------------------------------------------------------
####  size spectra estimation #####

# load in UVC data. Counts (= "number") of body sizes (="MASS") of reef fish > 10g
load("CREP_fishsizes.Rdata")

## size spectrum slope for each island by year
isd<-ddply(CREP, .(ISLAND, OBS_YEAR), plb_mlefunc, size="MASS", catch="number")

## Size spectrum slope = -exponent_num from function outputs. ##

write.csv(isd, file="CREP_ISD_estimates.csv")
