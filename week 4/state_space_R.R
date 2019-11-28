#Time-series / state-space models in TMB
#Cahill 27 Nov 2019 
rm(list=ls(all=TRUE))
#Simulate some fake data, analyze it, make lots of plots
require(TMB)
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
gam.check(Gam)

# Compile model
setwd( "C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 4" )
Version = "state_space_exponential"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "Nyears"=Nyears, "Y_obs_t"=y_obs )
Parameters = list( "logB0"=0, "log_sigmaP"=1, "log_sigmaO"=1, "mu_lambda"=1, "lambda_t"=rep(0,Nyears-1) )
Random = c("lambda_t")

Use_REML = TRUE
if( Use_REML==TRUE ) Random = union( Random, c("logB0","mu_lambda") )

# Build object
dyn.load( dynlib("state_space_exponential") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt = fit_tmb( obj=Obj, newtonsteps=1 )
Opt

SD = sdreport( Obj )
#SEHat  = as.list( Opt$SD, "Std. Error" )

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

# Get reporting and SEs
Report = Obj$report()
ParHat = as.list( Opt$SD, "Estimate" )
ParHat[["biomass_t"]] = SD$value[names(SD$value)=="biomass_t"]

ParHat

#Plot stuff
plot.data = data.frame(true_biomass=B, true_lambda=c(lambda_t, NA), 
                       y_obs = y_obs, gam_lower = ypred_t$fit - ypred_t$se.fit*1.96,gam_mu = ypred_t$fit,  
                       gam_upper = ypred_t$se.fit*1.96+ypred_t$fit, 
                       tmb_lower = ParHat$biomass_t - SD$sd[names(SD$value)=="biomass_t"]*1.96 , tmb_mu = ParHat$biomass_t,  
                       tmb_upper = SD$sd[names(SD$value)=="biomass_t"]*1.96 + ParHat$biomass_t, 
                       tmb_lambda_lower = c(ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96, NA),  
                       tmb_lambda_mu = c(ParHat$lambda_t, NA), 
                       tmb_lambda_upper = c(ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96, NA)
                       )
plot.data
png( file="results.png", width=8, height=5, res=800, units="in" )
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

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
dev.off()

#Here's a thing gam cannot give you (to my knowledge):

#Plot stuff
plot.data = data.frame(true_lambda=lambda_t,
                       tmb_lambda_lower = ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96, 
                       tmb_lambda_mu = ParHat$lambda_t, 
                       tmb_lambda_upper = ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96
)

png( file="lambda.png", width=8, height=5, res=800, units="in" )
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

plot(0, 0, ylim = c(0.5, 1.5), 
     xlim = c(0.5, Nyears-1), ylab = "Lambda", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:nrow(plot.data), (nrow(plot.data):1)), 
        y = c(plot.data$tmb_lambda_lower, plot.data$tmb_lambda_upper[(Nyears-1):1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
lines(ParHat$lambda_t)
lines(plot.data$true_lambda, pch=16, type="b", col="darkorange")
abline(h=1.0, col="blue", lty=2)
dev.off()

#---------------------------------------------------------
#Do it again but predict left out years this time
dyn.unload( dynlib("state_space_exponential") )

Data = list( "Nyears"=Nyears, "Y_obs_t"=y_obs )
Data$Y_obs_t[47:50]=NA
Parameters = list( "logB0"=0, "log_sigmaP"=1, "log_sigmaO"=1, "mu_lambda"=1, "lambda_t"=rep(0,Nyears-1) )

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

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

ParHat = as.list( Opt$SD, "Estimate" )
ParHat$lambda_t
ParHat[["biomass_t"]] = SD$value[names(SD$value)=="biomass_t"]
Report = Obj$report()

plot.data = data.frame(true_biomass=B, true_lambda=c(lambda_t, NA), 
                       y_obs = y_obs, gam_lower = ypred_t$fit - ypred_t$se.fit*1.96,gam_mu = ypred_t$fit,  
                       gam_upper = ypred_t$se.fit*1.96+ypred_t$fit, 
                       tmb_lower = ParHat$biomass_t - SD$sd[names(SD$value)=="biomass_t"]*1.96 , tmb_mu = ParHat$biomass_t,  
                       tmb_upper = SD$sd[names(SD$value)=="biomass_t"]*1.96 + ParHat$biomass_t, 
                       tmb_lambda_lower = c(ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96, NA),  
                       tmb_lambda_mu = c(ParHat$lambda_t, NA), 
                       tmb_lambda_upper = c(ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96, NA), 
                       year=1:Nyears
)

png( file="predict.png", width=8, height=5, res=800, units="in" )
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

plot(0, 0, ylim = c(min(plot.data$tmb_lower), max(plot.data$tmb_upper)), 
     xlim = c(0.5, Nyears), ylab = "Population
     Biomass (tons)", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:Nyears, Nyears:1), y = c(plot.data$tmb_lower, plot.data$tmb_upper[Nyears:1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
points(plot.data$y_obs, type="p", pch=16, col="black")
lines(plot.data$true_biomass, type="l", col="darkorange")# truth
lines(plot.data$tmb_mu, type="l", col="blue") #estimated
points(plot.data$y_obs[47:50]~plot.data$year[47:50], col="red", pch=16)

legend(x = 30, y = 25, legend = c("True", "Observed", "Estimated (TMB)", "Predicted"),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("darkorange","black", "blue", "red"),
       bty = "n", cex = 1)

dev.off()
#---------------------------------------------------------
#Redo it with 
#SigP = SigO = 1; B0 = 20;
#for a 'fun' learning experience
#i.e., see Auger-methe et al. 2016
#---------------------------------------------------------