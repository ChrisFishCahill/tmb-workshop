#Cahill 5 November 2019

setwd( "C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/Week 3" )
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
library(ggplot2)
library(lme4)
library(TMB)

#---------------------------------------
# Simulate data
#Normal mixed effects model, i.e., an LMM
#y_i ~ N(mu_i, SD^2)
#where mu_i = beta0 + eps_g[group]
# eps_g ~ N(0, SD_group^2)
#---------------------------------------

set.seed(13) 

# Simulate predictors
Ngroup = 25
group_i = rep(1:Ngroup, each=20)
eps_g = rnorm(length(unique(group_i)), mean=0, sd=1) #group devs/raneffs/error terms
beta0 = 0 #global mean  

#Simulate the response data:
y_i = eps_g[group_i] + beta0 + rnorm(length(group_i), mean=0, sd=1)

data <- data.frame("y_i" = y_i, "group_i"=group_i)

ggplot(data,aes(y=y_i, x=group_i))+
  geom_point()+
  geom_hline(yintercept=0, linetype=2) + 
  ylab("Value") + xlab("Group")

#Fit a linear model with each intercept as a fixed effect 
#I think Gelman and Hill call this "Complete pooling"

lm <- lm(data$y_i~as.factor(data$group_i))
coeffs <- c(coef(lm)[1], coef(lm)[1] + coef(lm)[2:Ngroup]) 
data$pooled <- coeffs[data$group_i]
data$truth <- compare$True[data$group_i]

#run in R using lme4 a la Boker and Bates
#Uses adaptive Gaussian-Hermite Quadrature
Use_REML = FALSE
Lme = lmer( y_i ~ 1|factor(group_i), REML=Use_REML)
summary(Lme)

#now in TMB
# Compile model
Version = "linear_mixed_model"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "Ngroups"=length(unique(group_i)), "g_i"=group_i-1, "y_i"=y_i) 
Parameters = list( "beta0"=4, "log_SD"=2, "log_SDG"=2, "eps_g"=rep(0,Data$Ngroups) )
Parameters
Random = c("eps_g")
if( Use_REML==TRUE ) Random = union( Random, "beta0")

# Build object
dyn.load( dynlib("linear_mixed_model") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
SD = sdreport( Obj )

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

# Get reporting and SEs
Report = Obj$report()
ParHat = as.list( Opt$SD, "Estimate" )

ParHat

#Compare some stuff
# Global mean
c( fixef(Lme), ParHat$beta0, mean(y_i) )

# Random effects
compare = data.frame( "True"=eps_g, ranef(Lme)[['factor(group_i)']], "TMB"=ParHat$eps_g )
colnames(compare)[2] <- "Lme4"
compare

#So, gauss-hermite quadrature and Laplace are identical out to ~7 decimal places

ggplot(compare, aes(x=Lme4, y=TMB)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype=2)

# Variances
summary(Lme)
unlist( Report[c("SDG","SD")] )

#Cash money 
mu_preds <- ParHat$beta0 + eps_g #These are the best estimates of the group means
data$mle = mu_preds[data$group_i]

#Plot everything
ggplot(data,aes(y=y_i, x=group_i))+
  geom_point()+
  geom_point(mapping=aes(y=mle, x=group_i, size=1.3), colour="blue") + 
  geom_point(mapping=aes(y=pooled, x=group_i, size=1.3), colour="darkorange1") + 
  geom_point(mapping=aes(y=truth, x=group_i, size=1.3), colour="grey") + 
  geom_abline(intercept = ParHat$beta0, slope=0, linetype=2, size=1) + 
  geom_abline(intercept = beta0, slope=0, linetype=2, size=1, colour="grey") + 
  theme(legend.position = "none") 

#blue is the mixed model, orange is the linear regression, grey is truth
#Now let's do the REML version
Random 
Use_REML = TRUE
if( Use_REML==TRUE ) Random = union( Random, "beta0")
Random #so, basically just telling TMB to integrate more stuff

Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)

# Optimize
Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
SD = sdreport( Obj )

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

ParHat = as.list( Opt$SD, "Estimate" )
ParHat
SEHat  = as.list( Opt$SD, "Std. Error" )
SEHat

#--------------------------------------------
#Now let's do it for an overdispersed Poisson model, i.e., a GLMM
#
#Suppose we are interested in counting Kentucky Jaguar Worms (obviously in Kentucky)
#Let's say that we conduct worm counts at a bunch of sites
#Let's further say that we are interested in accounting for observation error because
#Jaguar worms are elusive beasts and the technicians like Booker's bourbon
#
#Let's make a model where
# y_i ~ Poisson(lambda_i)
# where lambda_i = beta0 + eps_g + eps_obs
# & eps_g ~ N(0, SD_g^2)
# & eps_obs ~ N(0, SD_obs^2)

#Let's say there are 10 observations for each of 10 sites because I make the rules
#--------------------------------------------


