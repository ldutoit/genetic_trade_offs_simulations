### This file runs the simulations analyses for 
### supplementary figure 5, including the assymetric case.

################################################
## Wright - Fisher simulations
################################################

source('./suppFig-functions.R')

################################################
## Run simulations
################################################

dir<-"simResults/" # output directory, created if it does not exist. Finish by /


#create output directory if needed
if (!dir.exists(dir)) {
  dir.create(dir)
  message(paste(dir,"Directory created."))
} else {
  warning("Directory already exists.")
}


# Run the simulations for z.1 = z.2

 fixSA.netbenefit.sym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7, z=1, m=0.05, R = 1,    output.dir=dir, file.name="all")
 fixSA.netbenefit.sym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7, z=1, m=0.05, R = 1.02, output.dir=dir, file.name="all")
 fixSA.netbenefit.sym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7, z=1, m=0.05, R = 0.98, output.dir=dir, file.name="all")
 #Run the simulations for z.1 =/= z.2
 fixSA.netbenefit.asym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1,    output.dir=dir, file.name="all")
 fixSA.netbenefit.asym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 1.02, output.dir=dir, file.name="all")
 fixSA.netbenefit.asym(n=50, n.ben.tofix = 10^4, max.mutations = 10^7,  m=0.05, R = 0.98, output.dir=dir, file.name="all")
