### This file runs the simulations analyses of the manuscript.


################################################
## Wright - Fisher simulations
################################################

source("fixSA.netbenefit.R") # load the simulation function
source("num_pred&visualise.R") # load the numerical prediction and visualise

################################################
## Run simulations
################################################

dir <-"results/" #output directory
#simulations
set.seed(123)
for (R in seq(0.98,1.02,by=0.02)){
  for (n in c(50,25,10)){
    for (z in c(2,1,0.5)){
      print(paste(R,n,z))
      fixSA.netbenefit(n=n, n.ben.tofix = 10^3, max.mutations = 10^7,  m=0.05, z=z, R = R, output.dir=dir, file.name="all", visualize=TRUE)
    }  
  }
}

################################################
## Numerical predictions (n=50, z=1)
################################################

erf <- function(x){2*pnorm(x*sqrt(2)) - 1}

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

num_pred()

################################################
## Overlay numerical predictions with simulations (n=50, z=1)
################################################

# Main text figure (z = 1, n = 50)
dec <- read.csv(paste0(dir,"table_n50_m0.05_z1_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n50_m0.05_z1_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n50_m0.05_z1_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()

################################################
## Sup (n=50, z=2)
################################################

# Numerical predictions
n = 50 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 2 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

num_pred()

# Sup figures (z = 2)
dec <- read.csv(paste0(dir,"table_n50_m0.05_z2_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n50_m0.05_z2_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n50_m0.05_z2_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()

################################################
## Sup (n=50, z=0.5)
################################################

# Numerical predictions
n = 50 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 0.5 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

num_pred()

dec <- read.csv(paste0(dir,"table_n50_m0.05_z0.5_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n50_m0.05_z0.5_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n50_m0.05_z0.5_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()
################################################
## Sup (n=25, z=1)
################################################

# Numerical predictions
n = 25 #traits
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

num_pred()

dec <- read.csv(paste0(dir,"table_n25_m0.05_z1_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n25_m0.05_z1_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n25_m0.05_z1_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()

################################################
## Sup (n=25, z=2)
################################################

# Numerical predictions
n = 25 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 2 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

num_pred()
# Sup figures (z = 2)
dec <- read.csv(paste0(dir,"table_n25_m0.05_z2_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n25_m0.05_z2_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n25_m0.05_z2_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()


################################################
## Sup (n=25, z=0.5)
################################################

# Numerical predictions
n = 25 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 0.5 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

num_pred()

dec <- read.csv(paste0(dir,"table_n25_m0.05_z0.5_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n25_m0.05_z0.5_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n25_m0.05_z0.5_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()

################################################
## Sup (n=10, z=1)
################################################

# Numerical predictions
n = 10 #traits
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

num_pred()

dec <- read.csv(paste0(dir,"table_n10_m0.05_z1_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n10_m0.05_z1_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n10_m0.05_z1_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

plot(corr.sel.v2, Pr.SA.new, xlab = "Alignment of selection", ylab = "Proportion trade-off", type = "l", ylim = c(0, 1), lwd = 3, main="n=10, z=1")
lines(corr.sel.v2, Pr.SA.decline, col = "palegreen", lwd = 3)
points(x=dec$corr.sel,y=dec$Pr_SA,pch=19,col="palegreen")
lines(corr.sel.v2, Pr.SA.est, col = "green3", lwd = 3)
points(x=ct$corr.sel,y=ct$Pr_SA,pch=19,col="green3")
lines(corr.sel.v2, Pr.SA.expand, col="darkgreen", lwd = 3)
points(x=inc$corr.sel,y=inc$Pr_SA,pch=19,col="darkgreen")
points(x=ct$corr.sel,y=ct$Pr_SA_new,pch=19)


################################################
## Sup (n=10, z=2)
################################################

# Numerical predictions
n = 10 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 2 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

num_pred()

dec <- read.csv(paste0(dir,"table_n10_m0.05_z2_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n10_m0.05_z2_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n10_m0.05_z2_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()


################################################
## Sup (n=10, z=0.5)
################################################

# Numerical predictions
n = 10 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 0.5 #distance from optimum (assumed to be same in each sex)

corr.sel.v2 = 0:999/1000 #vector of cos(theta) values
Pr.SA.new = vector()
Pr.SA.est = vector()
Pr.SA.decline = vector()
Pr.SA.expand = vector()
R.decline = 0.98
R.expand = 1.02
s.min = max(0, (1 - R.decline)/R.decline)

num_pred()

dec <- read.csv(paste0(dir,"table_n10_m0.05_z0.5_R0.98.csv"))
dec$Pr_SA <- dec$SA.ben.fix/dec$netben.fix
inc <- read.csv(paste0(dir,"table_n10_m0.05_z0.5_R1.02.csv"))
inc$Pr_SA <- inc$SA.ben.fix/inc$netben.fix
ct <- read.csv(paste0(dir,"table_n10_m0.05_z0.5_R1.csv"))
ct$Pr_SA <- ct$SA.ben.fix/ct$netben.fix
ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

visualise()
