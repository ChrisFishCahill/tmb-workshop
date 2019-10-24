#data:
adults    <- c(2.4,3.2,3.9,5.7,6,7.4,8.2,10,10.1,10.4,11.3,12.8,18,24)
juveniles <- c(11.6,7.1,14.3,19.1,12.4,19.7,31.5,18.5,22.1,26.9,19.2,21,18.1,26.8)
temps     <- c(12.00,15.04,13.45,12.03,13.00,13.48,9.24,15.33,12.08,9.50,17.10,12.16,22.00,13.64)
data      <- data.frame(adults=adults, juveniles=juveniles, temps=temps)
Nobs      <- nrow(data)
plot(data)

library(TMB)
setwd("C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 1/week 1 soln")
compile( "hwk1.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("hwk1") )
Params = list("parm_vec"=rep(1, 4))

Data = list( "Y_i"=data$juveniles, "A_i"=data$adults, "T_i"=data$temps )
Obj = MakeADFun( data=Data, parameters=Params, DLL="hwk1")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )

opt_tmb = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )

opt_tmb$par # estimated parameters

SD = sdreport( Obj ) # standard errors
SD

#Some checks 
final_gradient = Obj$gr( opt_tmb$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not (converged) today, Satan!")

opt_tmb$objective #34.28401 = nll

AIC = 2*length(opt_tmb$par) + 2*opt_tmb$objective #2*n + 2*nll
AICc = AIC + ( 2*length(opt_tmb$par)^2 + 2*length(opt_tmb$par) ) / (nrow(data) - length(opt_tmb$par) - 1)

AIC  
AICc 
opt = TMBhelper::Optimize(obj=Obj, getsd=T, newtonsteps=1)

#-------------------------------
#Simulation experiment in R
#-------------------------------

#Declare parameters from the previous fit
alpha = opt_tmb$par[1]
beta = opt_tmb$par[2]
c_temp = opt_tmb$par[3]
sd = exp(opt_tmb$par[4])

Nobs = nrow(data)
Nsim = 100
SimResults = matrix(NA, nrow=Nsim, ncol=length(opt_tmb$par))
colnames(SimResults) = c("alpha", "beta", "c", "sd")

dyn.load( dynlib("hwk1") )
Params = list("parm_vec"=opt_tmb$par)

#Calculate the systematic component of the model | mles from fit: 
Juveniles_det_i = (alpha*data$adults) / ((1 + beta*data$adults)) + c_temp*data$temps

set.seed(3)
replicate=1

while( replicate < (Nsim+1) ){
  #Simulate data
  Juveniles_sim_i = rnorm(Nobs, Juveniles_det_i, sd)
  
  #Create tagged data list and AD function
  Data = list( "Y_i"=Juveniles_sim_i, "A_i"=data$adults, "T_i"=data$temps )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="hwk1")
  
  #Optimize and calculate SEs
  opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
                 
  #Convergence checks 
  SD = sdreport( Obj )
  final_gradient = Obj$gr( opt$par )
  
  #If estimation behaves, record results and advance loop
  if( all(abs(final_gradient)<0.0001) | SD$pdHess==TRUE ) {  
    SimResults[replicate, ] = opt$par
    replicate = replicate + 1 
  }
}

par(mfrow=c(2,2))
hist(SimResults[,1], main="alpha")
abline(v=alpha, lty=2, lwd=3, col="steelblue")

hist(SimResults[,2], main="beta")
abline(v=beta, lty=2, lwd=3, col="steelblue")

hist(SimResults[,3], main="c")
abline(v=c_temp, lty=2, lwd=3, col="steelblue")

hist(SimResults[,4], main="logsd")
abline(v=log(sd), lty=2, lwd=3, col="steelblue")

#-----------------------------
#Simulation experiment with TMB
#-----------------------------

SimResults2 <- replicate(Nsim, {
  simdata <- Obj$simulate(par=Obj$par, complete=TRUE)
  obj2 <- MakeADFun(simdata, Params, DLL="hwk1")
  nlminb(obj2$par, obj2$fn, obj2$gr)$par
})

par(mfrow=c(2,2))
hist(SimResults2[1,], main="alpha")
abline(v=alpha, lty=2, lwd=3, col="steelblue")

hist(SimResults2[2,], main="beta")
abline(v=beta, lty=2, lwd=3, col="steelblue")

hist(SimResults2[3,], main="c")
abline(v=c_temp, lty=2, lwd=3, col="steelblue")

hist(SimResults2[4,], main="logsd")
abline(v=log(sd), lty=2, lwd=3, col="steelblue")
