# Code to produce supplementary figure with 
#
# Panel (A): Schematic illustration of sampling scheme in 2-D
# Panel (B): Plot overlaying numerical predictions and
#            simulation results for both symmetric and 
#            asymmetric displacements from the optima,
#            ranging from cos(theta_sel) = 0 to 1. 


rm(list=ls())
############################################
### Package Dependencies
library(plotrix)

############################################
## Some simple dependency functions 

proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

z1.len  <-  function(theta) {
  2*sqrt(1 - sin(theta)^2)
}

#' Gaussian Error Function 
erf <- function(x) {
    2*pnorm(x * sqrt(2)) - 1
}




############################################
# Approximations for completely aligned 
# orientations. Eqs. S16 - S18

xBar  <-  function(n=50, m=0.05, zBar) {
  (n*m) / (2*sqrt(zBar^2 + (1/2)*n*m^2))
}
 
x.min <-  function(n=50, m=0.05, z.min) {
  (n*m) / (2*sqrt(z.min^2 + (1/2)*n*m^2))
}

Pr.tradeoff.Aligned.New  <-  function(xBar, x.min) {
  (erf(x.min/sqrt(2)) - erf(xBar/sqrt(2))) / (1 - erf(xBar/sqrt(2)) )
}

Pr.Tradeoff.Aligned.R1  <-  function(n=50, xBar, x.min) {

  (   sqrt(2)/sqrt(pi) * (exp(- xBar^2/2) - exp(- x.min^2/2)) + xBar*(erf(xBar/sqrt(2)) -  erf(x.min/sqrt(2))) ) / 
    ( sqrt(2)/sqrt(pi) * exp(- xBar^2/2) - xBar * (1 - erf(xBar/sqrt(2))) ) 
}


Pr.Tradeoff.Aligned.expand  <-  function(R, zBar, xBar, x.min) {
    (  (R*(2*zBar^2 * xBar*sqrt(2))/(n*sqrt(pi))) * (exp(-xBar^2/2 ) - exp(-x.min^2/2) )) / 
    (  (R*(2*zBar^2 * xBar*sqrt(2))/(n*sqrt(pi))) * exp(-xBar^2/2) + ( (R - 1) - R*2*zBar^2 * xBar^2/n)*(1 - erf(xBar/sqrt(2))) ) +
    ( ((R - 1) - R*2*zBar^2 * xBar^2/n) * (erf(x.min/sqrt(2)) - erf(xBar/sqrt(2))) ) / 
    (  (R*(2*zBar^2 * xBar*sqrt(2))/(n*sqrt(pi))) * exp(-xBar^2/2) + ( (R - 1) - R*2*zBar^2 * xBar^2/n)*(1 - erf(xBar/sqrt(2))) )
}


#  !!! THIS FUNCTION IS MISBEHAVING. NOT SURE IF THERE IS A TYPO OR WHAT !!!
Pr.Tradeoff.Aligned.decline  <-  function(R, zBar, xBar, x.min) {
  psi  <-  ((1-R)*n) / ((2*(zBar^2)*R*xBar) )

  ( (R*2*zBar^2*xBar*sqrt(2))/(n*sqrt(pi)) * (exp(-(xBar + psi)^2/2) - exp(-x.min^2/2) )) / 
  ( (R*2*zBar^2*xBar*sqrt(2))/(n*sqrt(pi)) *  exp(-(xBar + psi)^2/2) + ( (R - 1) - (R*2*zBar^2 * xBar^2)/n)*(1 - erf((xBar + psi)/sqrt(2)) ) ) +
    (( (R - 1) - (R*2*zBar^2 * xBar^2)/n)  * (erf(x.min/sqrt(2)) - erf((xBar + psi)/sqrt(2)) ) ) / 
  ( (R*2*zBar^2*xBar*sqrt(2))/(n*sqrt(pi)) *  exp(-(xBar + psi)^2/2) + ( (R - 1) - (R*2*zBar^2 * xBar^2)/n)*(1 - erf((xBar + psi)/sqrt(2)) ) )

}






############################################
## Code for numerical predictions
## and simulations
##
## - Code for symmetric cases are from Ludo
## - Code for asymmetric cases are mine, but 
##   are only slightly edited versions of Ludo's
##   functions to get the asymmetric results
##


# Numerical prediction lines for SYMMETRIC DISPLACEMENTS
# this is a contained function, need to define parameters before calling
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

# Numerical prediction lines for ASYMMETRIC DISPLACEMENTS
# this is a contained function, need to define parameters before calling
num_pred_asym  <-  function(){
  for(i in 501:1000){
    theta <<-  acos(corr.sel.v2[i])
    z1    <<-  2*sqrt(1 - sin(theta)^2)
    z2    <<-  1

    #Predicted mean and variance for s = (s.f + s.m)/2
    mean.s <<- -0.5*n*m^2
    var.s <<- 0.5*m^2*(0.5*z1^2 + 0.5*z2^2 + z1*z2*cos(theta) + n*m^2)
    
    DFE.all = function(s){
      dnorm(s, mean.s, sqrt(var.s))
    }
    
    #proportion of mutations that are beneficial
    Pr.ben <<- integrate(DFE.all, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for beneficial mutations
    Afun  =  function(s, n, m,z1, z2){
      (2*s +( (z1^2 - z2^2)*(2*s + n*m^2) ) / (z1^2 + z2^2 + 2*z1*z2*cos(theta) + 2*n*m^2)) / ( m*sqrt(z1^2 + z2^2 - 2*z1*z2*cos(theta) - (z1^2 - z2^2)^2 / (z1^2 + z2^2 +2*z1*z2*cos(theta) + 2*n*m^2)) )
    }
    Bfun  =  function(s, n, m,z1, z2){
      (2*s - ( (z1^2 - z2^2)*(2*s + n*m^2) ) / (z1^2 + z2^2 +2*z1*z2*cos(theta) + 2*n*m^2)) / ( m*sqrt(z1^2 + z2^2 - 2*z1*z2*cos(theta) - (z1^2 - z2^2)^2 / (z1^2 + z2^2 +2*z1*z2*cos(theta) + 2*n*m^2)) )
    }

    SA.new = function(s){
      (1 - (1/2)*(erf(Afun(s, n, m, z1, z2)/sqrt(2)) + erf(Bfun(s, n, m,z1, z2)/sqrt(2))))*dnorm(s, mean.s, sqrt(var.s))/Pr.ben
    }
    
    Pr.SA.new[i] <<- integrate(SA.new, 0, Inf, stop.on.error = FALSE)$value
    
    #probability that a random beneficial mutation is established (constant N)
    integrand = function(s){
      2*s*dnorm(s, mean.s, sqrt(var.s))/Pr.ben
    }
    
    Pr.establishment <<- integrate(integrand, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for established mutations
    SA.est = function(s){
      (1 - (1/2)*(erf(Afun(s, n, m,z1, z2)/sqrt(2)) + erf(Bfun(s, n, m,z1, z2)/sqrt(2))))*2*s*dnorm(s, mean.s, sqrt(var.s))/(Pr.establishment*Pr.ben)
    }
    
    Pr.SA.est[i] <<- integrate(SA.est, 0, Inf, stop.on.error = FALSE)$value
    
    #probability that a random beneficial mutation is established (declining N)
    integrand = function(s){
      2*((1 + s)*R.decline - 1)*dnorm(s, mean.s, sqrt(var.s))
    }
    
    Pr.est.decline <<- integrate(integrand, lower = s.min, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for established mutations
    SA.est.decline = function(s){
      (1 - (1/2)*(erf(Afun(s, n, m,z1, z2)/sqrt(2)) + erf(Bfun(s, n, m,z1, z2)/sqrt(2))))*2*((1 + s)*R.decline - 1)*dnorm(s, mean.s, sqrt(var.s))/(Pr.est.decline)
    }
    
    Pr.SA.decline[i] <<- integrate(SA.est.decline, s.min, Inf, stop.on.error = FALSE)$value
    
    #probability that a random beneficial mutation is established (expanding N)
    integrand = function(s){
      2*((1 + s)*R.expand - 1)*dnorm(s, mean.s, sqrt(var.s))
    }
    
    Pr.est.expand <<- integrate(integrand, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for established mutations
    SA.est.expand = function(s){
      (1 - (1/2)*(erf(Afun(s, n, m,z1, z2)/sqrt(2)) + erf(Bfun(s, n, m,z1, z2)/sqrt(2))))*2*((1 + s)*R.expand - 1)*dnorm(s, mean.s, sqrt(var.s))/(Pr.est.expand)
    }
    
    Pr.SA.expand[i] <<- integrate(SA.est.expand, 0, Inf, stop.on.error = FALSE)$value
  }
}  

##########################
# Symmetric displacements
# Recreate main results from Fig. 3b
fixSA.netbenefit.sym <- function(n, n.ben.tofix, max.mutations, z, m, R, output.dir, file.name, visualize=TRUE){
  
  ## Defining more parameters:
  
  # 1. Position of optima, male and female:
  z.f = z #female distance
  z.m = z #male distance
  O = rep(0, n) #location of optimum
  
  # 2. Population size:
  N.min= 10^3
  N.max= 10^7
  N.ct= 10^5
  if(R<1){
    N0=N.max
    N.t = function(t){
      max(N0*R^t, N.min)}
  } else if(R>1){
    N0=N.min
    N.t = function(t){
      min(N0*R^t, N.max)}
  } else{
    N0=N.ct
    N.t=function(t){
      N.ct}
  }
  
  # 3. Create table with information about simulation parameters:
  myparameters<-data.frame("n"=n, "m"=m, "R"=R) 
  
  # 4. Define output file name, if not given:
  get.name=c()
  if(file.name=="all"){
    for(i in seq_along(colnames(myparameters))){
      get.name[i]<-paste(colnames(myparameters)[i], myparameters[,i], sep="")
    }}else{
      for(i in seq_along(file.name)){
        get.name[i]<-paste(file.name[i], x[,file.name[i]], sep="")
      }
    }
  thefile=paste(get.name, collapse="_")
  
  # 5. Make sure directory name ends with "/":
  if(substr(output.dir, nchar(output.dir), nchar(output.dir)) %in% "/"){
    thedir<-output.dir
  }else{
    thedir<-paste(output.dir, "/", sep="")
  }
  
  #Hold results of all simulations
  myresults<-list()
  
  ## 6. Loop across values of corr.sel:
  all.cs=seq(0, 1, 0.2) # create vector with values of corr.sel from 0 to 1, with increments of 0.2
  for(j in seq_along(all.cs)){
    
    # Position of phenotypes, given angle of selection
    corr.sel=all.cs[j]
    theta = acos(corr.sel) #angle between selection vectors
    A.f = c(z.f, rep(0, n-1)) #female wild type phenotype
    A.m = c(z.m*cos(theta), z.m*sin(theta), rep(0, n-2))  #male wild type phenotype
      
    ## Mutational information
    ben = SA = SA.ben = SA.ben.fix = n.ben.fix = 0
    
    # Number of mutations to start with
    no.mutations = 0
    
    # Create vectors to hold selection coefficients of simulated mutations, so that we can later estimate the fraction of trade-offs among fixed mutations:
    s.f = vector() # female selection coefficients
    s.m = vector() # male selection coefficients
    fixations <- vector() # fixed or lost
    
    #print(paste("Starting simulations", j, "out of", length(all.cs), "- for corr.sel=", corr.sel))
    
    # While we still haven't reached the target, keep simulating new mutations
    while(n.ben.fix < n.ben.tofix & no.mutations < max.mutations){
      
      no.mutations = no.mutations + 1
      x = rnorm(n) #vector of standard normal random variables
      r = m*sqrt(rchisq(1, n, ncp = 0)) # Magnitude of mutations
      M = r*x/sqrt(sum(x^2)) #effects on the set of traits
      
      #Female and male distances from optimum
      z.mut.f = sqrt(sum((A.f + M - O)^2))
      z.mut.m = sqrt(sum((A.m + M - O)^2))
      
      #heterozygous and homozygous fitness effects, per sex
      s.f[no.mutations] = exp(-0.5*z.mut.f^2)/exp(-0.5*z.f^2) - 1
      s.m[no.mutations] = exp(-0.5*z.mut.m^2)/exp(-0.5*z.m^2) - 1

      
      #Sex-specific fitnesses of each haploid genotype
      wa.f <- 1 + s.f[no.mutations]
      wa.m <- 1 + s.m[no.mutations]
      wA.f <- 1
      wA.m <- 1
      
      #Internal equilibrium
      int.eq <- -((s.f[no.mutations]+s.m[no.mutations])/(2*s.f[no.mutations]*s.m[no.mutations]))
      #Replace invalid internal equilibria with 1
      int.eq[int.eq<0 | int.eq>1] <- 1
      
      #Initial frequency of a allele in the overall population
      p.0 <- 1/N0
      p <- p.0
      time = 0
      
      #If the allele is net-beneficial, evolve it to fixation or loss
#      if(R > 1) {
#        netBenCrit  <- (s.f[no.mutations]+s.m[no.mutations])/2 > (1 - R)/R  
#      } else{netBenCrit  <-  (wa.f+wa.m)/2 > 1 }
      netBenCrit  <-  (wa.f+wa.m)/2 > 1 
      if(netBenCrit){

        while (p>0 & p<1 & p<int.eq){  # While the mutation is not fixed
          #Pop size:
          N = round(N.t(time))
          
          #Sex-specific allele frequencies after selection    
          p.f <- (p*wa.f) / ((p*wa.f)+((1-p)*wA.f))
          p.m <- (p*wa.m) / ((p*wa.m)+((1-p)*wA.m))
          #Wright-Fisher sampler:
          nomin.1<-rbinom(1,round(N/2),prob=p.f)
          nomin.2<-rbinom(1,round(N/2),prob=p.m)
          nomin=nomin.1+nomin.2
          p <- nomin/N
          
          time=time+1
        }
        fixations[no.mutations] <- p
        #If fixation, record it and add it to fixation counter
        if(p>=int.eq)
        {n.ben.fix = n.ben.fix + 1 # count to n.ben.tofix
        } 
      } 
    }
    ## Compute results once target is reached:
    #Beneficial alleles
    ben = sum((s.f + s.m)/2 > 0)
    #SA alleles
    SA = sum(s.f*s.m < 0)
    #SA and beneficial alleles
    SA.ben =  sum(s.f*s.m<0 & (s.f + s.m)/2 > 0)
    #Beneficial alleles that fix
    ben.fix = sum( (s.f + s.m)/2 > 0 & fixations>0)
    #SA beneficial alleles that fix
    SA.ben.fix = sum( s.f*s.m<0 & (s.f + s.m)/2 > 0 & fixations>0)
    print(paste("Total mut.=",no.mutations,", total ben mut.=",ben,", N ben. fix.=",n.ben.fix,", Pr(SA)=",SA.ben.fix/ben.fix))
    
    # close the while loop
    #print(paste("Finished simulations", j, "out of", length(all.cs), "- for corr.sel=", corr.sel))
    myresults[[length(myresults)+1]]<-cbind(myparameters, data.frame("z.m"=z.m,"z.f"=z.f,"theta"=theta, "corr.sel"=corr.sel, "total.mut"=no.mutations, "netben.mut"=ben,"SA.mut"=SA,"netben.SA.mut"=SA.ben,"netben.fix"=ben.fix,"SA.ben.fix"=SA.ben.fix)) 
  }
  
  #Save outputs:
  write.csv(do.call(rbind, myresults), file=paste(thedir, "Symmetry_table_", thefile, ".csv", sep=""),quote=F, row.names = F)
  
  print(paste("output written in:",thedir,thefile,".csv",sep=""))
  
}

# Rerun main sims for R > 1 with new criteria for "net beneficial"
#fixSA.netbenefit.sym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7, z=1, m=0.05, R = 1.02, output.dir=dir, file.name="all")


## Parameters
# The following parameters are fixed within the function:
# - N.min: minimum population size, starting value for expanding populations and floor values for declining populations (set to 10^3)
# - N.max: maximum population size, ceiling value for expanding populations and starting values for declining populations  (set to 10^7)
# - N.ct: population size for a stable population (set to 10^5)
# - corr.sel: between sex correlation for selection (set to run for values of 0 to 1, with increments of 0.2)

# The function uses the following arguments:
# - n: number of traits or dimensionality (default is 50)
# - n.ben.tofix: the number of net beneficial to fix before the simulation is considered complete (default is 10^3)
# -  max.mutations: the maximum number of mutations that can be simulated before the simulations stop to avoid near infinite sims ( default id 10^9)
# - m: standard deviation of the mutation effects per trait (scaled; default is 0.05)
# - R: rate of population growth (default is 1, indicating constant population size)
# - N0: population size at time 0 (default is 10^5)
# - output.dir: path where to save output files (default is working directory)
# - file.name: argument to define which parameter values should be shown in the file name of the output files. Options are a vector of minimum length 1 (e.g. file.name="R", or file.name=c("R", "z")), where "all" creates a file name indicating n, m, z, R and N0 (default is "all")

fixSA.netbenefit.asym <- function(n, n.ben.tofix, max.mutations,  m, R, output.dir, file.name, visualize=TRUE){
  
  ## Defining more parameters:
  
  # 1. Position of optima, male and female:
  O = rep(0, n) #location of optimum
  
  # 2. Population size:
  N.min= 10^3
  N.max= 10^7
  N.ct= 10^5
  if(R<1){
    N0=N.max
    N.t = function(t){
      max(N0*R^t, N.min)}
  } else if(R>1){
    N0=N.min
    N.t = function(t){
      min(N0*R^t, N.max)}
  } else{
    N0=N.ct
    N.t=function(t){
      N.ct}
  }
  
  # 3. Create table with information about simulation parameters:
  myparameters<-data.frame("n"=n, "m"=m, "R"=R) 
  
  # 4. Define output file name, if not given:
  get.name=c()
  if(file.name=="all"){
    for(i in seq_along(colnames(myparameters))){
      get.name[i]<-paste(colnames(myparameters)[i], myparameters[,i], sep="")
    }}else{
      for(i in seq_along(file.name)){
        get.name[i]<-paste(file.name[i], x[,file.name[i]], sep="")
      }
    }
  thefile=paste(get.name, collapse="_")
  
  # 5. Make sure directory name ends with "/":
  if(substr(output.dir, nchar(output.dir), nchar(output.dir)) %in% "/"){
    thedir<-output.dir
  }else{
    thedir<-paste(output.dir, "/", sep="")
  }
  
  #Hold results of all simulations
  myresults<-list()
  
  ## 6. Loop across values of corr.sel:
  #all.cs=seq(0, 1, 0.2) # create vector with values of corr.sel from 0 to 1, with increments of 0.2
  for(theta in seq(pi/3, 0,length=10)){
    z.f = 1
    z.m = 2*sqrt(1 - (sin(theta))**2)##COLIN
    # Position of phenotypes, given angle of selection
    #corr.sel=all.cs[j]
    #theta = acos(corr.sel) #angle between selection vectors
    A.f = c(1, rep(0, n-1)) #female wild type phenotype
    A.m = c(z.m*cos(theta), z.m*sin(theta), rep(0, n-2))  #male wild type phenotype
    
    ## Mutational information
    ben = SA = SA.ben = SA.ben.fix = n.ben.fix = 0
    
    # Number of mutations to start with
    no.mutations = 0
    
    # Create vectors to hold selection coefficients of simulated mutations, so that we can later estimate the fraction of trade-offs among fixed mutations:
    s.f = vector() # female selection coefficients
    s.m = vector() # male selection coefficients
    fixations <- vector() # fixed or lost
    
    #print(paste("Starting simulations", j, "out of", length(all.cs), "- for corr.sel=", corr.sel))
    
    # While we still haven't reached the target, keep simulating new mutations
    while(n.ben.fix < n.ben.tofix & no.mutations < max.mutations){
      
      no.mutations = no.mutations + 1
      x = rnorm(n) #vector of standard normal random variables
      r = m*sqrt(rchisq(1, n, ncp = 0)) # Magnitude of mutations
      M = r*x/sqrt(sum(x^2)) #effects on the set of traits
      
      #Female and male distances from optimum
      z.mut.f = sqrt(sum((A.f + M - O)^2))
      z.mut.m = sqrt(sum((A.m + M - O)^2))
      
      #heterozygous and homozygous fitness effects, per sex
      s.f[no.mutations] = exp(-0.5*z.mut.f^2)/exp(-0.5*z.f^2) - 1
      s.m[no.mutations] = exp(-0.5*z.mut.m^2)/exp(-0.5*z.m^2) - 1

      
      #Sex-specific fitnesses of each haploid genotype
      wa.f <- 1 + s.f[no.mutations]
      wa.m <- 1 + s.m[no.mutations]
      wA.f <- 1
      wA.m <- 1
      
      #Internal equilibrium
      int.eq <- -((s.f[no.mutations]+s.m[no.mutations])/(2*s.f[no.mutations]*s.m[no.mutations]))
      #Replace invalid internal equilibria with 1
      int.eq[int.eq<0 | int.eq>1] <- 1
      
      #Initial frequency of a allele in the overall population
      p.0 <- 1/N0
      p <- p.0
      time = 0
      
      #If the allele is net-beneficial, evolve it to fixation or loss
#      if(R > 1) {
#        netBenCrit  <- (s.f[no.mutations]+s.m[no.mutations])/2 > (1 - R)/R  
#      } else{netBenCrit  <-  (wa.f+wa.m)/2 > 1 }
      netBenCrit  <-  (wa.f+wa.m)/2 > 1 

      if(netBenCrit){

        while (p>0 & p<1 & p<int.eq){  # While the mutation is not fixed
          #Pop size:
          N = round(N.t(time))
          
          #Sex-specific allele frequencies after selection    
          p.f <- (p*wa.f) / ((p*wa.f)+((1-p)*wA.f))
          p.m <- (p*wa.m) / ((p*wa.m)+((1-p)*wA.m))
          #Wright-Fisher sampler:
          nomin.1<-rbinom(1,round(N/2),prob=p.f)
          nomin.2<-rbinom(1,round(N/2),prob=p.m)
          nomin=nomin.1+nomin.2
          p <- nomin/N
          
          time=time+1
        }
        fixations[no.mutations] <- p
        #If fixation, record it and add it to fixation counter
        if(p>=int.eq)
        {n.ben.fix = n.ben.fix + 1 # count to n.ben.tofix
        } 
      } 
    }
    ## Compute results once target is reached:
    #Beneficial alleles
    ben = sum((s.f + s.m)/2 > 0)
    #SA alleles
    SA = sum(s.f*s.m < 0)
    #SA and beneficial alleles
    SA.ben =  sum(s.f*s.m<0 & (s.f + s.m)/2 > 0)
    #Beneficial alleles that fix
    ben.fix = sum( (s.f + s.m)/2 > 0 & fixations>0)
    #SA beneficial alleles that fix
    SA.ben.fix = sum( s.f*s.m<0 & (s.f + s.m)/2 > 0 & fixations>0)
    print(paste("Total mut.=",no.mutations,", total ben mut.=",ben,", N ben. fix.=",n.ben.fix,", Pr(SA)=",SA.ben.fix/ben.fix))
    
    # close the while loop
    #print(paste("Finished simulations", j, "out of", length(all.cs), "- for corr.sel=", corr.sel))
    myresults[[length(myresults)+1]]<-cbind(myparameters, data.frame("z.m"=z.m,"z.f"=z.f,"theta"=theta, "total.mut"=no.mutations, "netben.mut"=ben,"SA.mut"=SA,"netben.SA.mut"=SA.ben,"netben.fix"=ben.fix,"SA.ben.fix"=SA.ben.fix)) 
  }
  
  #Save outputs:
  write.csv(do.call(rbind, myresults), file=paste(thedir, "assymetry_table_", thefile, ".csv", sep=""),quote=F, row.names = F)
  
  print(paste("output written in:",thedir,thefile,".csv",sep=""))
  
}



