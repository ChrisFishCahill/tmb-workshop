#Cahill 23 Oct 2019
#Simulate data for Peninsular homing clams 
#
#Model form:
#y_i ~ Poisson(lambda_i)
#E(y_i) = lambda_i
#Var(y_i) = lambda_i 
#log(lambda_i) = Beta0 + Beta1*x1 + Beta2*x2 
# 
#Beta1 is some measure of habitat quality 
#Beta2 is temperature

library( TMB )
library( ggplot2 )

#Pick some parameters
Nobs = 1000

Beta0 = 4.9
Beta1 = 0.75
Beta2 = -0.3
set.seed(1)
x1 = rnorm(Nobs, mean=0, sd=1) 
x2 = rnorm(Nobs, mean=0, sd=1) 

#define the linear predictor
lambda_i = Beta0 + Beta1*x1 + Beta2*x2^2
y_i = rpois(Nobs, lambda_i)

data = data.frame("Y_i"=y_i, "x1"=x1, "x2"=x2)

plot(data) 
plot(data$Y_i~data$x2^2)
glm.fit = glm(data$Y_i~data$x1 + data$x2^2, family="poisson")
summary(glm.fit)
#------------

#Compile and load the .cpp 

setwd("C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 2")
compile( "poisson.cpp" )
dyn.load( dynlib("poisson") )

# Step 2 -- build inputs and object
Data = list( "y_i"=data$Y_i, "X_ij"= model.matrix(~ 1 + data$x1 + data$x2^2 ))
Params = list("b_j"=rep(0, ncol(Data$X_ij)))
Obj = MakeADFun( data=Data, parameters=Params, DLL="poisson")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )

opt = TMBhelper::Optimize(obj=Obj, getsd=T, newtonsteps=1)

SD = sdreport( Obj )
final_gradient = Obj$gr( opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

#Extract the intercept and SE
ParHat = as.list( opt$SD, "Estimate" )
SEHat  = as.list( opt$SD, "Std. Error" )

opt$objective #nll
logLik(glm.fit) #ll from glm 

#-----------------

