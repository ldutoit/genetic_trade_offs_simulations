#### Two functions are here, the visualisation and the numerical_predictions
#### Those functions have no optional arguments and output in the global environment. 
#### As such, they are not meant to be used outside of the core_simulations.R
###Numerical predictions, employed in core_simulations.R
num_pred<-function(){
  for(i in 1:1000){
    theta <<- acos(corr.sel.v2[i])
    
    #Predicted mean and variance for s = (s.f + s.m)/2
    mean.s <<- -0.5*n*m^2
    var.s <<- 0.5*m^2*((1 + cos(theta))*z^2 + n*m^2)
    
    DFE.all = function(s){
      dnorm(s, mean.s, sqrt(var.s))
    }
    
    #proportion of mutations that are beneficial
    Pr.ben <<- integrate(DFE.all, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for beneficial mutations
    SA.new = function(s){
      (1 - erf(s/(m*z*sqrt(1 - cos(theta)))))*dnorm(s, mean.s, sqrt(var.s))/Pr.ben
    }
    
    Pr.SA.new[i] <<- integrate(SA.new, 0, Inf, stop.on.error = FALSE)$value
    
    #probability that a random beneficial mutation is established (constant N)
    integrand = function(s){
      2*s*dnorm(s, mean.s, sqrt(var.s))/Pr.ben
    }
    
    Pr.establishment <<- integrate(integrand, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for established mutations
    SA.est = function(s){
      (1 - erf(s/(m*z*sqrt(1 - cos(theta)))))*2*s*dnorm(s, mean.s, sqrt(var.s))/(Pr.establishment*Pr.ben)
    }
    
    Pr.SA.est[i] <<- integrate(SA.est, 0, Inf, stop.on.error = FALSE)$value
    
    #probability that a random beneficial mutation is established (declining N)
    integrand = function(s){
      2*((1 + s)*R.decline - 1)*dnorm(s, mean.s, sqrt(var.s))
    }
    
    Pr.est.decline <<- integrate(integrand, lower = s.min, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for established mutations
    SA.est.decline = function(s){
      (1 - erf(s/(m*z*sqrt(1 - cos(theta)))))*2*((1 + s)*R.decline - 1)*dnorm(s, mean.s, sqrt(var.s))/(Pr.est.decline)
    }
    
    Pr.SA.decline[i] <<- integrate(SA.est.decline, s.min, Inf, stop.on.error = FALSE)$value
    
    #probability that a random beneficial mutation is established (expanding N)
    integrand = function(s){
      2*((1 + s)*R.expand - 1)*dnorm(s, mean.s, sqrt(var.s))
    }
    
    Pr.est.expand <<- integrate(integrand, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for established mutations
    SA.est.expand = function(s){
      (1 - erf(s/(m*z*sqrt(1 - cos(theta)))))*2*((1 + s)*R.expand - 1)*dnorm(s, mean.s, sqrt(var.s))/(Pr.est.expand)
    }
    
    Pr.SA.expand[i] <<- integrate(SA.est.expand, 0, Inf, stop.on.error = FALSE)$value
  }
}  

#### Visualisations
visualise<-function(){
  plot(corr.sel.v2, Pr.SA.new, xlab = "Alignment of selection", ylab = "Proportion trade-off", type = "l", ylim = c(0, 1), lwd = 3,main="n=50, z=1")
  lines(corr.sel.v2, Pr.SA.decline, col = "palegreen", lwd = 3)
  points(x=dec$corr.sel,y=dec$Pr_SA,pch=19,col="palegreen")
  lines(corr.sel.v2, Pr.SA.est, col = "green3", lwd = 3)
  points(x=ct$corr.sel,y=ct$Pr_SA,pch=19,col="green3")
  lines(corr.sel.v2, Pr.SA.expand, col="darkgreen", lwd = 3)
  points(x=inc$corr.sel,y=inc$Pr_SA,pch=19,col="darkgreen")
  points(x=ct$corr.sel,y=ct$Pr_SA_new,pch=19)
}

