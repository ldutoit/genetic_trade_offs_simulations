### This file defines a function that runs the Wright-Fisher distribution of fitness effects with sexes. It focuses on the assymetric case, and does the plotting as well.
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
    var.s  <<- 0.5*m^2*((z1^2 + z2^2)/2 + z1*z2*cos(theta) + n*m^2)
    
    DFE.all = function(s){
      dnorm(s, mean.s, sqrt(var.s))
    }
    
    #proportion of mutations that are beneficial
    Pr.ben <<- integrate(DFE.all, lower = 0, upper = Inf, stop.on.error = FALSE)$value
    
    #Pr(SA) for beneficial mutations
    Afun  =  function(s, n, m,z1, z2){
      (2*s +( (z1^2 - z2^2)*(2*s + n*m^2) ) / (z1^2 + z2^2 + 2*z1*z2*cos(theta) + 2*n*m^2)) / ( m*sqrt(z1^2 + z2^2 - 2*z1*z2*cos(theta) - ( (z1^2 - z2^2) ) / (z1^2 + z2^2 +2*z1*z2*cos(theta) + 2*n*m^2)) )
    }
    Bfun  =  function(s, n, m,z1, z2){
      (2*s - ( (z1^2 - z2^2)*(2*s + n*m^2) ) / (z1^2 + z2^2 +2*z1*z2*cos(theta) + 2*n*m^2)) / ( m*sqrt(z1^2 + z2^2 - 2*z1*z2*cos(theta) - ( (z1^2 - z2^2) ) / (z1^2 + z2^2 +2*z1*z2*cos(theta) + 2*n*m^2)) )
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
#num_pred_asym()

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
      if(R > 1) {
        netBenCrit  <- (s.f[no.mutations]+s.m[no.mutations])/2 > (1 - R)/R  
      } else{netBenCrit  <-  (wa.f+wa.m)/2 > 1 }

      if(netBenCrit){
#      if((wa.f+wa.m)/2>1){
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
fixSA.netbenefit.sym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7, z=1, m=0.05, R = 1.02, output.dir=dir, file.name="all")


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
      if(R > 1) {
        netBenCrit  <- (s.f[no.mutations]+s.m[no.mutations])/2 > (1 - R)/R  
      } else{netBenCrit  <-  (wa.f+wa.m)/2 > 1 }

      if(netBenCrit){
#      if((wa.f+wa.m)/2>1){
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


###SIMULATIONS

dir = "example_results/"
fixSA.netbenefit.asym(n=2,  n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1, output.dir=dir, file.name="all")
fixSA.netbenefit.asym(n=10, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1, output.dir=dir, file.name="all")
fixSA.netbenefit.asym(n=20, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1, output.dir=dir, file.name="all")
fixSA.netbenefit.asym(n=30, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1, output.dir=dir, file.name="all")
fixSA.netbenefit.asym(n=40, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1, output.dir=dir, file.name="all")
fixSA.netbenefit.asym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1, output.dir=dir, file.name="all")

fixSA.netbenefit.asym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1.02, output.dir=dir, file.name="all")
fixSA.netbenefit.asym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 0.98, output.dir=dir, file.name="all")


###VISUALISATIONS
tradeoffs_n2            <- read.csv(paste0(dir,"/asymetry_table_n2_m0.05_R1.csv"), h=T)
tradeoffs_n10           <- read.csv(paste0(dir,"/asymetry_table_n10_m0.05_R1.csv"),h=T)
tradeoffs_n20           <- read.csv(paste0(dir,"/asymetry_table_n20_m0.05_R1.csv"),h=T)
tradeoffs_n30           <- read.csv(paste0(dir,"/asymetry_table_n30_m0.05_R1.csv"),h=T)
tradeoffs_n40           <- read.csv(paste0(dir,"/asymetry_table_n40_m0.05_R1.csv"),h=T)
tradeoffs_n50           <- read.csv(paste0(dir,"/asymetry_table_n50_m0.05_R1.csv"),h=T)

# Stable Populations, R = 1
tradeoffs_n2$Pr_SA     <- tradeoffs_n2$SA.ben.fix/tradeoffs_n2$netben.fix
tradeoffs_n2$Pr_SA_new <- tradeoffs_n2$netben.SA.mut/tradeoffs_n2$netben.mut

tradeoffs_n10$Pr_SA     <- tradeoffs_n10$SA.ben.fix/tradeoffs_n10$netben.fix
tradeoffs_n10$Pr_SA_new <- tradeoffs_n10$netben.SA.mut/tradeoffs_n10$netben.mut

tradeoffs_n20$Pr_SA     <- tradeoffs_n20$SA.ben.fix/tradeoffs_n20$netben.fix
tradeoffs_n20$Pr_SA_new <- tradeoffs_n20$netben.SA.mut/tradeoffs_n20$netben.mut

tradeoffs_n30$Pr_SA     <- tradeoffs_n30$SA.ben.fix/tradeoffs_n30$netben.fix
tradeoffs_n30$Pr_SA_new <- tradeoffs_n30$netben.SA.mut/tradeoffs_n30$netben.mut

tradeoffs_n40$Pr_SA     <- tradeoffs_n40$SA.ben.fix/tradeoffs_n40$netben.fix
tradeoffs_n40$Pr_SA_new <- tradeoffs_n40$netben.SA.mut/tradeoffs_n40$netben.mut

tradeoffs_n50$Pr_SA     <- tradeoffs_n50$SA.ben.fix/tradeoffs_n50$netben.fix
tradeoffs_n50$Pr_SA_new <- tradeoffs_n50$netben.SA.mut/tradeoffs_n50$netben.mut

# Expanding Populations, R = 1.02
tradeoffs_R1.02_n50           <- read.csv(paste0(dir,"/assymetry_table_n50_m0.05_R1.02.csv"),h=T)
tradeoffs_R1.02_n50$Pr_SA  <-  tradeoffs_R1.02_n50$SA.ben.fix/tradeoffs_R1.02_n50$netben.fix

# Declining Populations, R = 0.98
tradeoffs_R0.98_n50           <- read.csv(paste0(dir,"/assymetry_table_n50_m0.05_R0.98.csv"),h=T)
tradeoffs_R0.98_n50$Pr_SA  <-  tradeoffs_R0.98_n50$SA.ben.fix/tradeoffs_R0.98_n50$netben.fix


#Plot theta
pdf(file = "./plot_assymetry_n_2-50.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches
# Proportion of net beneficial mutations that exhibit trade-offs
par(mfrow=c(1,2))
plot(x=cos(tradeoffs_n50$theta),y=tradeoffs_n50$Pr_SA_new,pch=19,col="black", xlab = expression(paste("Alignment of selection: ", cos(italic(theta[sel])))), ylab = "Proportion trade-off", main = "New mutations", type = "p", ylim = c(0, 0.6), lwd = 1)
points(x=cos(tradeoffs_n40$theta),y=tradeoffs_n40$Pr_SA_new,pch=21,col="black", bg='grey90', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n30$theta),y=tradeoffs_n30$Pr_SA_new,pch=21,col="black", bg='grey80', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n20$theta),y=tradeoffs_n20$Pr_SA_new,pch=21,col="black", bg='grey70', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n10$theta),y=tradeoffs_n10$Pr_SA_new,pch=21,col="black", bg='grey60', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n2$theta),y=tradeoffs_n2$Pr_SA_new,pch=21,col="black", bg='grey50', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)

# Proportion of fixed mutations that are also SA
plot(x=cos(tradeoffs_n50$theta),y=tradeoffs_n50$Pr_SA,pch=19,col="black", xlab = expression(paste("Alignment of selection: ", cos(italic(theta[sel])))), ylab = "", main = "Established mutations", type = "p", ylim = c(0, 0.6), lwd = 1)
points(x=cos(tradeoffs_n40$theta),y=tradeoffs_n40$Pr_SA,pch=21,col="black", bg='grey90', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n30$theta),y=tradeoffs_n30$Pr_SA,pch=21,col="black", bg='grey80', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n20$theta),y=tradeoffs_n20$Pr_SA,pch=21,col="black", bg='grey70', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n10$theta),y=tradeoffs_n10$Pr_SA,pch=21,col="black", bg='grey60', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)
points(x=cos(tradeoffs_n2$theta),y=tradeoffs_n2$Pr_SA,pch=21,col="black", bg='grey50', xlab = "z.m (z.f=1)", ylab = "Proportion trade-off", type = "p", ylim = c(0, 1), lwd = 1)

#Legend
usr  <-  par('usr')
legend( x       =  usr[2],
        y       =  usr[4],
        legend  =  c(
                     expression(paste(n==50)),
                     expression(paste(n==40)),
                     expression(paste(n==30)),
                     expression(paste(n==20)),
                     expression(paste(n==10)),
                     expression(paste(n==2))),
         col     =  'black',
         pch     =  c(21,21,21,21,21,21),
         pt.bg   =  c('black',
                      'grey90',
                      'grey80',
                      'grey70',
                      'grey60',
                      'grey50'),
         cex     =  0.65,
         pt.cex  =  0.75,
         xjust   =  1,
         yjust   =  1,
         bty     =  'n',
         border  =  NA)

dev.off()

## Overlay over demographic scenarios

dec <- read.csv(paste0(dir,"table_n50_m0.05_z1_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
#inc <- read.csv(paste0(dir,"table_n50_m0.05_z1_R1.02.csv"))
inc <- read.csv(paste0(dir,"Symmetry_table_n50_m0.05_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n50_m0.05_z1_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

# Numerical predictions 
#Parameters
n = 50 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 1 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

#num_pred() # outputs, Pr.SA.decline,Pr.SA.est, Pr.SA.expand

###########################
### Make the overlay plot
###########################
pdf(file = "./plot_assymetry_overlay_n_2-50.pdf",   # Directory + filename
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# make plot
par(omi=c(0.25, 0.5, 0.25, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
plot(NA, xlab = " ", ylab = " ", type = "l", ylim = c(0, 0.8),xlim=c(0,1), lwd = 3, main=" ",col="white")

# Plot Analytic approximations for Symmetric displacements
num_pred()
lines(corr.sel.v2, Pr.SA.new, col = "black", lwd = 3)
lines(corr.sel.v2, Pr.SA.decline, col = "palegreen", lwd = 3)
lines(corr.sel.v2, Pr.SA.est, col = "green3", lwd = 3)
lines(corr.sel.v2, Pr.SA.expand, col="darkgreen", lwd = 3)

# Plot main results for Symmetric displacements
points(x=ct$corr.sel,y=ct$Pr_SA_new, pch=21, col=adjustcolor("black",alpha=1),      bg=adjustcolor("black",alpha=0.8))
points(x=inc$corr.sel,y=inc$Pr_SA,   pch=21, col=adjustcolor("darkgreen",alpha=1),  bg=adjustcolor("darkgreen",alpha=0.8))
points(x=ct$corr.sel,y=ct$Pr_SA,     pch=21, col=adjustcolor("green",alpha=1),      bg=adjustcolor("green",alpha=0.8))
points(x=dec$corr.sel,y=dec$Pr_SA,   pch=21, col=adjustcolor("lightgreen",alpha=1), bg=adjustcolor("lightgreen",alpha=0.8))

# Plot Analytic approximations for Asymmetric displacements
num_pred_asym()
lines(corr.sel.v2[501:1000], Pr.SA.new[501:1000], col = "black", lwd = 3, lty=3)
lines(corr.sel.v2[501:1000], Pr.SA.expand[501:1000], col="darkgreen", lwd = 3, lty=3)
lines(corr.sel.v2[501:1000], Pr.SA.est[501:1000], col = "green3", lwd = 3, lty=3)
lines(corr.sel.v2[501:1000], Pr.SA.decline[501:1000], col = "palegreen", lwd = 3, lty=3)

# Plot supp results for Asymmetric displacements
points(cos(tradeoffs_n50$theta)[-9],tradeoffs_n50$Pr_SA_new[-9],         pch=23, col=adjustcolor("black",alpha=1),     bg=adjustcolor("black",alpha=0.5))
points(cos(tradeoffs_R1.02_n50$theta)[-9],tradeoffs_R1.02_n50$Pr_SA[-9], pch=23, col=adjustcolor("darkgreen",alpha=1),  bg=adjustcolor("darkgreen",alpha=0.5))
points(cos(tradeoffs_n50$theta)[-9],tradeoffs_n50$Pr_SA[-9],             pch=23, col=adjustcolor("green",alpha=1), bg=adjustcolor("green",alpha=0.5))
points(cos(tradeoffs_R0.98_n50$theta)[-9],tradeoffs_R0.98_n50$Pr_SA[-9], pch=23, col=adjustcolor("lightgreen",alpha=1), bg=adjustcolor("lightgreen",alpha=0.5))


proportionalLabel(0.5, 1.075, expression(paste("Symmetric vs. Asymmetric Displacements, ", italic(n)==50)),cex=1.4, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel(-0.15, 0.5, expression(paste("Proportion trade-off")),cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.15, expression(paste("Alignment of selection ", (cos(theta[sel])))),cex=1.2, adj=c(0.5, 0.5), xpd=NA)

#Legend
usr  <-  par('usr')
legend( x       =  usr[2]*0.675,
        y       =  usr[4],
        legend  =  c(
          expression("New mutations"),
          expression("R > 1"),
          expression("R = 1"),
          expression("R < 1")),
        col     =   c('black',
                      'darkgreen',
                      'green',
                      'lightgreen'),
        pch     =  21,
        pt.bg   =  c(adjustcolor("black",alpha=0.8),
                     adjustcolor("darkgreen",alpha=0.8),
                     adjustcolor("green",alpha=0.8),
                     adjustcolor("lightgreen",alpha=0.8)),
        cex     =  1,
        pt.cex  =  1,
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA)
legend( x       =  usr[2]*0.76,
        y       =  usr[4],
        legend  =  c(
          expression(" "),
          expression(" "),
          expression(" "),
          expression(" ")),
        col     =   c('black',
                      'darkgreen',
                      'green',
                      'lightgreen'),
        pch     =  23,
        pt.bg   =  c(adjustcolor("black",alpha=0.5),
                     adjustcolor("darkgreen",alpha=0.5),
                     adjustcolor("green",alpha=0.5),
                     adjustcolor("lightgreen",alpha=0.5)),
        cex     =  1,
        pt.cex  =  1,
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA)
#proportionalLabel(0.74, 0.925, "}", adj=c(0, 1), cex=4)
proportionalLabel(0.79, 0.915, "Asymmetric", adj=c(0, 1), cex=1)
proportionalLabel(0.77, 0.875, "Displacements", adj=c(0, 1), cex=1)
dev.off()