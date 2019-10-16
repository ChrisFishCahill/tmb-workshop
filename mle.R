#Cahill 16 Oct 2019
#Simulate some fake data and estimate it using several methods
#We'll do this for a trivial example even though this model has analytical solutions
#
#Model form:
#Y_i = Intercept + Slope*Predictor + error_i where error_i ~ N(0, sd)
#

library( TMB )
library( ggplot2 )

#Pick some parameters
Nobs = 100
Beta0 = 3.7
Beta1 = 0.5
LogSD = 0.2

#Simulate some covariates and then observed data:
set.seed(1)
X1 = rnorm(Nobs, mean=0, sd=1) 
y_i = mu_i = Beta0 + Beta1*X1 + rnorm(Nobs, mean=0, sd=exp(LogSD))

data = data.frame("Y_i"=y_i, "x1"=X1)

#Plot it so we know what's going on
p <- ggplot(data = data,
            aes(x = x1,
                y = Y_i)) +
     geom_point() + 
     geom_smooth(method = "lm")

p
#----------------------------------
#Write a NLL function to be optimized in R:

parms = c(4.2, 0.2, 0.3) #starting values for Beta0, Beta1, LogSD

#Write the objective function
getNegLogLike <- function(parms){
  Beta0 <- parms[1] #unpack the parameters
  Beta1 <- parms[2]
  SD <- exp(parms[3]) 
  Nobs  <- nrow(data)
  #Note that data is a global variable so can call it internally in this f(x)
  Y_pred_i   <-  rep(NA, Nobs)
  Y_pred_i   <-  Beta0 + Beta1*data$x1 #Create a vector of predictions
  NegLogLike   <-  -sum(log(dnorm(data$Y_i, Y_pred_i, sd=SD))) #calculate the normal dist neg log like
  
  #Other ways to code the NegLL:
  #NegLogLike2  <-  -sum(dnorm(data$Y_i, Y_pred_i, sd=SD, log=T)) #neg log like another way
  #NegLogLike  <-  -sum(log( (1 / sqrt(2*pi*SD^2)) * exp( - ((data$Y_i-Y_pred_i)^2/(2*SD^2)) ) )) #nll nasty way
  
  return(NegLogLike)
}

#Does it give a value when we call it?
getNegLogLike(parms)

#Fit 1: An acid-test for all other calculated log likelihoods:
logLik(lm(data$Y_i ~ data$x1))

#Fit 2--calculate the neg log like via nlminb
?nlminb
opt_nlminb <- nlminb(parms, getNegLogLike)
opt_nlminb$par
opt_nlminb$objective #neg log like

#Fit 3--calculate neg log like via optim
opt_optim <- optim(parms, getNegLogLike, hessian=T)
opt_optim$par
opt_optim$value #neg log like

#Extract the SE and calcualte 95% Confidence Intervals:
#------------------------------------------------------------#
#Time for some math--This is here for extra information
#to help show how one gets standard errors for the parameters

#The hessian matrix is the second derivative of the 
#Objective function, so if the objective function is 
#negative Log Likelihood the hessian is the observed 
#"Fisher information."  The inverse of the hessian is 
#thus an asymptotic estimate of the variance covariance matrix of
#the parameters
#------------------------------------------------------------#

varcov <- solve(opt_optim$hessian) #solve hessian for variance-covariance matrix
se <- sqrt(diag(varcov))
se

#Plot it
NumberOfBetas <-  length(parms)
Combined <- rbind(opt_optim$par[1] + se[1]%o%qnorm(c(0.025,0.50, 0.975)), #+ / - 1.96*SE = 95% confidence interval
                  opt_optim$par[2] + se[2]%o%qnorm(c(0.025,0.50, 0.975)),
                  opt_optim$par[3] + se[3]%o%qnorm(c(0.025,0.50, 0.975)))
                  
colnames(Combined) <- c("lower95", "MLE", "upper95")
Combined <- as.data.frame(Combined)

Combined$WhichVariable <- factor(c("Beta0", "Beta1", "LogSD"))

Combined$Truth <- c(Beta0, Beta1, LogSD)

p <- ggplot()
p <- p + geom_point(data = Combined,
                    aes(x = WhichVariable,
                        y = MLE))

p <- p + geom_errorbar(data = Combined,
                       aes(x = WhichVariable,
                           ymax = upper95,
                           ymin = lower95),
                       width=0.2)

p <- p + xlab("Parameter") + ylab("Value")
p <- p + theme(text = element_text(size = 15)) + theme(axis.title.x = element_blank())

p <- p + theme(legend.position="none")

p <- p + geom_jitter(data = Combined,
               aes(x = WhichVariable,
                   y = Truth), 
               color='darkblue', shape = 2, width = 0.05, height = 0)
p

#-------------------------------------
#Now let's do it in TMB

# Step 1 -- make and compile template file
setwd("C:/Users/Chris Cahill/Documents/GitHub/tmb-workshop/week 1")
compile( "linear_model_v1.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v1") )
Params = list("Beta0"=0, "Beta1"=0, "logSD"=0)
Data = list( "Y_i"=data$Y_i, "x1_i"=data$x1 )
Obj = MakeADFun( data=Data, parameters=Params, DLL="linear_model_v1")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )
opt_tmb = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
opt_tmb$diagnostics = data.frame( "name"=names(Obj$par), "Est"=opt_tmb$par, "final_gradient"=as.vector(Obj$gr(opt_tmb$par)))

opt_tmb$par # estimated parameters

SD = sdreport( Obj ) # standard errors
SD

#Some checks 
final_gradient = Obj$gr( opt_tmb$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not (converged) today, Satan!")

SD$value #Prediction for each data point
SD$sd #SEs for each data point

#End end end
