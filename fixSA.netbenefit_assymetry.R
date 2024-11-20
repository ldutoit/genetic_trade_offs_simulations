### This file defines a function that runs the Wright-Fisher distribution of fitness effects with sexes.

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
# - z: distance from phenotype to the optimum (assumed to be the same in each sex; default is 1)
# - R: rate of population growth (default is 1, indicating constant population size)
# - N0: population size at time 0 (default is 10^5)
# - output.dir: path where to save output files (default is working directory)
# - file.name: argument to define which parameter values should be shown in the file name of the output files. Options are a vector of minimum length 1 (e.g. file.name="R", or file.name=c("R", "z")), where "all" creates a file name indicating n, m, z, R and N0 (default is "all")
# - visualize: TRUE/FALSE argument defining if plots (sm~sf and s.diff~s.mean) should be outputted or not (default is TRUE)

fixSA.netbenefit<-function(n, n.ben.tofix, max.mutations,  m, z, R, output.dir, file.name, visualize=TRUE){
  
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
  myparameters<-data.frame("n"=n, "m"=m, "z"=z, "R"=R) 
  
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
    
    print(paste("Starting simulations", j, "out of", length(all.cs), "- for corr.sel=", corr.sel))
      
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
      if((wa.f+wa.m)/2>1){  
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
    print(paste("Finished simulations", j, "out of", length(all.cs), "- for corr.sel=", corr.sel))
    myresults[[j]]<-cbind(myparameters, data.frame("corr.sel"=corr.sel, "total.mut"=no.mutations, "netben.mut"=ben,"SA.mut"=SA,"netben.SA.mut"=SA.ben,"netben.fix"=ben.fix,"SA.ben.fix"=SA.ben.fix)) 
  }
  
  #Save outputs:
  write.csv(do.call(rbind, myresults), file=paste(thedir, "table_", thefile, ".csv", sep=""),quote=F, row.names = F)
  
  print(paste("output written in:",thedir,thefile,".csv",sep=""))

}

