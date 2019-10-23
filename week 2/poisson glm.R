#Cahill 23 Oct 2019
#Simulate data for Gorbachev's Holy Goose 
#
#Model form:
#y_i ~ Poisson(lambda_i)
#E(y_i) = lambda_i
#Var(y_i) = lambda_i 
#log(lambda_i) = Beta0 + Beta1*x1 + Beta2*x2 + Beta3*x3
# 

library( TMB )
library( ggplot2 )

#Pick some parameters
Nobs = 1000

mean.lambda = 10  #mean expected count
Beta0 = log(mean.lambda) #mean expected count log space
Beta1 = -2 #effect of elevation
Beta2 = 2 #effect of forest cover
Beta3 = 1 #interaction effect of elevation and forest

set.seed(3)
x1 = runif(Nobs, -1, 1) # elevation
x2 = runif(Nobs, -1, 1) # cover

#define the linear predictor
log.lambda = Beta0 + Beta1*x1 + Beta2*x2 + Beta3*x1*x2 

lambda = exp(log.lambda)

y_i = rpois(Nobs, lambda)
  
plot(y_i)

data = data.frame("Y_i"=y_i, "x1"=x1, "x2"=x2)

plot(data) 
glm.fit = glm(data$Y_i~data$x1 + data$x2 + data$x1:data$x2 , family="poisson")
summary(glm.fit)

#------------

#Compile and load the .cpp 

setwd("C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 2")
compile( "poisson.cpp" )
dyn.load( dynlib("poisson") )

#Check out model.matrix()
head(model.matrix(~ 1 + data$x1 + data$x2 + data$x1:data$x2))

# Step 2 -- build inputs and object
Data = list( "y_i"=data$Y_i, "X_ij"= model.matrix(~ 1 + data$x1 + data$x2 + data$x1:data$x2))
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

opt
#-----------------

