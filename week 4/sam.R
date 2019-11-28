require(TMB)
require(TMBhelper)
nt = 100
x0 = 3
sigmaP = 0.2
sigmaO = 0.2
alpha = 0.04

# Simulate predictors
x_t = y_t = rep(NA, nt)
x_t[1] = x0

#latent process
for( t in 2:nt ){
  x_t[t] = x_t[t-1] + rnorm( 1, mean=alpha, sd=sigmaP )
}

#observation model
for( t in 1:nt ){
  y_t[t] = x_t[t] + rnorm( 1, mean=0, sd=sigmaO )
}

plot(y_t ~ I(1:nt)) #observations = y_t
lines(x_t ~ I(1:nt), col="blue", lty=3) #truth = x_t

setwd( "C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 4" )
Version = "kalman"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "nt"=nt, "y_t"=y_t )
Parameters = list( "x0"=0, "log_sigmaP"=1, "log_sigmaO"=1, "alpha"=1, "x_t"=rep(0,nt) )
Random = c("x_t")

Use_REML = TRUE
if( Use_REML==TRUE ) Random = union( Random, c("x0","alpha") )

# Build object
dyn.load( dynlib("kalman") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt = fit_tmb( obj=Obj, newtonsteps=1 )
Opt

SD = sdreport( Obj )
