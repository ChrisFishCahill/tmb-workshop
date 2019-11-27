#Simulate some fake data, analyze it, make lots of plots
require(TMB)
require(ggplot2)
require(tidyverse)
require(TMBhelper)
require(mgcv)

Nyears = 50 #number of years of data
B_init = 20 #initial population size
mu_lambda = 1.05 #Mean annual population growth rate
SigP = 0.15 #Process (temporal) variation in growth rate lambda
SigO = 7 #Variance of the observation term 

y_obs = B = Nyears
B[1] = B_init #initialize  the N vector

#Simulate lambdas with process error
set.seed(42)
lambda_t = rnorm(Nyears-1, mu_lambda, SigP)

#Simulate the exponential growth model:
for(t in 1:(Nyears-1)){
  B[t+1] = B[t]*lambda_t[t]
}

#Generate the observed data | True Biomass
for(t in 1:Nyears){
  y_obs[t] = rnorm(1, B[t], SigO) 
}

#Fit a smoother across time (not the same model, but important alternative)
#could also do a loess, but the Galpern posse seems to have a thing for Simon's wiggle-lines
Gam = gam( y_obs ~ s(I(1:Nyears)) )
ypred_t = predict(Gam, se=TRUE)

# Compile model
setwd( "C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 4" )
Version = "state_space_exponential"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "Nyears"=Nyears, "Y_obs_t"=y_obs )
Parameters = list( "B0"=0, "log_sigmaP"=1, "log_sigmaO"=1, "mu_lambda"=1, "lambda_t"=rep(0,Nyears) )
Random = c("lambda_t")

Use_REML = TRUE
if( Use_REML==TRUE ) Random = union( Random, c("B0","mu_lambda") )

# Build object
dyn.load( dynlib("state_space_exponential") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt = fit_tmb( obj=Obj, newtonsteps=1 )
Opt

SD = sdreport( Obj )

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

# Get reporting and SEs
Report = Obj$report()
ParHat = as.list( Opt$SD, "Estimate" )
ParHat[["biomass_t"]] = SD$value[names(SD$value)=="biomass_t"]

ParHat

#Plot stuff
plot.data = data.frame(true_biomass=B, true_lambda=c(mu_lambda, lambda_t), 
                       y_obs = y_obs, gam_lower = ypred_t$fit - ypred_t$se.fit*1.96,gam_mu = ypred_t$fit,  
                       gam_upper = ypred_t$se.fit*1.96+ypred_t$fit, 
                       tmb_lower = ParHat$biomass_t - SD$sd[names(SD$value)=="biomass_t"]*1.96 , tmb_mu = ParHat$biomass_t,  
                       tmb_upper = SD$sd[names(SD$value)=="biomass_t"]*1.96 + ParHat$biomass_t, 
                       tmb_lambda_lower = ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96,  
                       tmb_lambda_mu = ParHat$lambda_t, 
                       tmb_lambda_upper = ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96
                       )
plot.data

plot(0, 0, ylim = c(min(plot.data$tmb_lower), max(plot.data$tmb_upper)), 
     xlim = c(0.5, Nyears), ylab = "Population
     Biomass (tons)", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:Nyears, Nyears:1), y = c(plot.data$tmb_lower, plot.data$tmb_upper[Nyears:1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
points(plot.data$y_obs, type="p", pch=16, col="black")
lines(plot.data$true_biomass, type="l", col="darkorange")# truth
lines(plot.data$tmb_mu, type="l", col="blue") #estimated
lines(plot.data$gam_mu, type="l", col="red") #estimated gam

legend(x = 30, y = 45, legend = c("True", "Observed", "Estimated (TMB)", 
                                                       "Estimated (Gam)"),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("darkorange","black", "blue", "red"),
       bty = "n", cex = 1)

#Here's a thing gam cannot give you (to my knowledge):
plot(0, 0, ylim = c(min(plot.data$tmb_lambda_lower), max(plot.data$tmb_lambda_upper)), 
     xlim = c(0.5, Nyears), ylab = "Lambda", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:Nyears, Nyears:1), y = c(plot.data$tmb_lambda_lower, plot.data$tmb_lambda_upper[Nyears:1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
lines(ParHat$lambda_t)
lines(plot.data$true_lambda, col="darkorange")

#---------------------------------------------------------
#Do it again but predict left out years this time
dyn.unload( dynlib("state_space_exponential") )

Data = list( "Nyears"=Nyears, "Y_obs_t"=y_obs )
Data$Y_obs_t[47:50]=NA
Parameters = list( "B0"=20, "log_sigmaP"=1, "log_sigmaO"=1, "mu_lambda"=1, "lambda_t"=rep(0,Nyears) )

compile( paste0(Version,".cpp") )
dyn.load( dynlib("state_space_exponential") )

Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt = fit_tmb( obj=Obj, newtonsteps=1 )
Opt

SD = sdreport( Obj )
SD

ParHat = as.list( Opt$SD, "Estimate" )
ParHat$lambda_t
ParHat[["biomass_t"]] = SD$value[names(SD$value)=="biomass_t"]
Report = Obj$report()

plot.data = data.frame(true_biomass=B, true_lambda=c(mu_lambda, lambda_t), 
                       y_obs = y_obs, gam_lower = ypred_t$fit - ypred_t$se.fit*1.96,gam_mu = ypred_t$fit,  
                       gam_upper = ypred_t$se.fit*1.96+ypred_t$fit, 
                       tmb_lower = ParHat$biomass_t - SD$sd[names(SD$value)=="biomass_t"]*1.96 , tmb_mu = ParHat$biomass_t,  
                       tmb_upper = SD$sd[names(SD$value)=="biomass_t"]*1.96 + ParHat$biomass_t, 
                       tmb_lambda_lower = ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96,  
                       tmb_lambda_mu = ParHat$lambda_t, 
                       tmb_lambda_upper = ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96
)


plot(0, 0, ylim = c(min(plot.data$tmb_lower), max(plot.data$tmb_upper)), 
     xlim = c(0.5, Nyears), ylab = "Population
     Biomass (tons)", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:Nyears, Nyears:1), y = c(plot.data$tmb_lower, plot.data$tmb_upper[Nyears:1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
points(plot.data$y_obs, type="p", pch=16, col="black")
lines(plot.data$true_biomass, type="l", col="darkorange")# truth
lines(plot.data$tmb_mu, type="l", col="blue") #estimated

legend(x = 30, y = 25, legend = c("True", "Observed", "Estimated (TMB)"),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("darkorange","black", "blue"),
       bty = "n", cex = 1)

#How could we get these confidence intervals to stay positive? 
#i.e., why did chris do something dumb that he only now realized?
#---------------------------------------------------------
#Redo it with 
#SigP = SigO = 1; B0 = 20;
#for a 'fun' learning experience
#---------------------------------------------------------