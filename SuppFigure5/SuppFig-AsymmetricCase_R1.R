# Make the Plot 
#
# Panel (A) - 2-d representation of the sampling scheme for asymmetric displacements
# Panel (B) - Plot overlaying the following:
#               * Numerical approxmations for symmetric displacements (z.1 = z.2, 0 < cos(theta_sel) < 1)
#               * Simulation results for symmetric displacements (z.1 = z.2, 0 < cos(theta_sel) < 1)
#               * Numerical approxmations for Asymmetric displacements (z.1 =/= z.2, 0.5 < cos(theta_sel) < 1)
#               * Simulation results for Assymmetric displacements (z.1 =/= z.2, 0.5 < cos(theta_sel) < 1)
#               * Analytic approximations for completely aligned but asymmetric displacements (z.1 =/= z.2, cos(theta_sel) = 1)


rm(list=ls())
# Source the dependencies
source('./suppFig-functions.R')

library(plotrix)


###########################
#save a pdf of the figure
pdf(file = "./SuppFig_Assym_Arc_n50.pdf",   # Directory + filename
    width = 11, # The width of the plot in inches
    height = 6) # The height of the plot in inches


# Colors
COL8  <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set plot layout
layout.mat  <-  matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
layout      <-  layout(layout.mat,respect=TRUE)

############################
# Panel (A): sampling scheme
z2  <-  1 # we are tracing the circle with r = 1 from O2

theta.sel  <-  seq(0, pi/3, length=10)
theta.1    <-  theta.sel
theta.2    <-  2*theta.sel
z1         <-  z1.len(theta=theta.sel)
xs  <-  cos(theta.2) + 1
ys  <-  -sin(theta.2)

par(omi=c(0.25, 0.5, 0.25, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(-2,2), 
                                             ylim = c(-2,2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        box()
        # Draw Circles
         draw.circle(0,0,1, lwd=1, border='blue', lty=2) # pos. sel.
         draw.circle(0,0,z1[5], lwd=0.5, border='blue', lty=2) # pos. sel.
         points(0,0, pch=21, col='blue', bg=adjustcolor('blue',alpha=0.6))
         proportionalLabel(0.46,  0.53,    expression(paste(italic(O[1]))), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
         draw.circle(1,0,1, lwd=1, border='red', lty=2) # pos. sel.
         points(1,0, pch=21, col='red', bg=adjustcolor('red',alpha=0.6))
         proportionalLabel(0.69,  0.53,    expression(paste(italic(O[2]))), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
         proportionalLabel(0.9,  0.29,    expression(paste(italic(A))), cex=0.75, adj=c(0.5, 0.5), xpd=NA)

        # Sampling points
        points(ys*1.05 ~ xs+1, pch=23, col=COL8[1], bg=adjustcolor(COL8[1], alpha=0.6))

        # z1 & z2 - draw geometry for double-angle
        draw.radial.line(0, 1, angle=0, col=adjustcolor(COL8[1], alpha=1))
        draw.radial.line(0, z1[5], angle=-theta.sel[5], col=adjustcolor(COL8[1], alpha=1))
        proportionalLabel(0.65,  0.4,    expression(paste(italic(z[1]))), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
        draw.radial.line(0, 1, center=c(1,0), angle=-theta.2[5], col=adjustcolor(COL8[1], alpha=1))
        proportionalLabel(0.825,  0.42,    expression(paste(italic(z[2]))), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
        draw.radial.line(0, 1, center=c(1,0), angle=0, col=adjustcolor(COL8[1], alpha=1), lty=2)

        # Draw in angles
        proportionalLabel(0.8075,  0.365,    expression(paste(theta[sel])), cex=0.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.58,  0.48,    expression(paste(theta[sel])), cex=0.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.785,  0.48,    expression(paste(2*italic(theta[sel]))), cex=0.5, adj=c(0.5, 0.5), xpd=NA)
        draw.arc(0,0, radius=0.5,angle1=0,angle2=-theta.sel[5])
        draw.arc(1,0, radius=0.42,angle1=0,angle2=-theta.2[5])
        draw.arc(x=xs[5],y=ys[5]*1.05, radius=0.45, angle1=(pi-theta.2[5]), angle2=(pi-theta.sel[5]))

        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5,  1.075,   expression(paste("Sampling Scheme in 2D")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.175,   expression(paste("Trait 1")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.175, 0.5,   expression(paste("Trait 2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

        proportionalLabel(0.08,  0.21,    expression(paste(italic(z[2])==1)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.125,  0.13,    expression(paste(italic(O[2])-italic(O[1])==1)), cex=0.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.21,  0.05,    expression(paste(italic(z[1])==2*sqrt(1 - sin^2*(theta[sel])))), cex=0.75, adj=c(0.5, 0.5), xpd=NA)


############################
# Panel (B): Numerical Results
dir = "./simResults/"

#inc  <-  read.csv(paste0(dir,"table_n50_m0.05_z1_R1.02.csv"))
#ct   <-  read.csv(paste0(dir,"table_n50_m0.05_z1_R1.csv"))
#dec  <-  read.csv(paste0(dir,"table_n50_m0.05_z1_R0.98.csv"))
inc  <-  read.csv(paste0(dir,"Symmetry_table_n50_m0.05_R1.02.csv"))
ct   <-  read.csv(paste0(dir,"Symmetry_table_n50_m0.05_R1.csv"))
dec  <-  read.csv(paste0(dir,"Symmetry_table_n50_m0.05_R0.98.csv"))

inc$Pr_SA    <- inc$SA.ben.fix/inc$netben.fix
ct$Pr_SA     <- ct$SA.ben.fix/ct$netben.fix
dec$Pr_SA    <- dec$SA.ben.fix/dec$netben.fix

ct$Pr_SA_new <- ct$netben.SA.mut/ct$netben.mut

# Numerical predictions 
#Parameters
n = 50 #traits
m = 0.05 #st. dev. mutation effects per trait
z = 1 #distance from optimum (assumed to be same in each sex)

corr.sel.v2   = 0:999/1000 #vector of cos(theta) values
Pr.SA.new     = vector()
Pr.SA.est     = vector()
Pr.SA.decline = vector()
Pr.SA.expand  = vector()
R.decline     = 0.98
R.expand      = 1.02
s.min         = max(0, (1 - R.decline)/R.decline)
num_pred() # calculate numerical predictions

# Read data for asymmetric case
tradeoffs_n50           <- read.csv(paste0(dir,"/assymetry_table_n50_m0.05_R1.csv"),h=T)
tradeoffs_n50$Pr_SA     <- tradeoffs_n50$SA.ben.fix/tradeoffs_n50$netben.fix
tradeoffs_n50$Pr_SA_new <- tradeoffs_n50$netben.SA.mut/tradeoffs_n50$netben.mut

# Expanding Populations, R = 1.02
tradeoffs_R1.02_n50           <- read.csv(paste0(dir,"/assymetry_table_n50_m0.05_R1.02.csv"),h=T)
tradeoffs_R1.02_n50$Pr_SA  <-  tradeoffs_R1.02_n50$SA.ben.fix/tradeoffs_R1.02_n50$netben.fix

# Declining Populations, R = 0.98
tradeoffs_R0.98_n50           <- read.csv(paste0(dir,"/assymetry_table_n50_m0.05_R0.98.csv"),h=T)
tradeoffs_R0.98_n50$Pr_SA  <-  tradeoffs_R0.98_n50$SA.ben.fix/tradeoffs_R0.98_n50$netben.fix


### Plot showing numerical results vs. simulation results
plot(NA, xlab = " ", ylab = " ", type = "l", ylim = c(0, 0.825),xlim=c(0,1), lwd = 3, main=" ",col="white")

# Plot Analytic approximations for Symmetric displacements
lines(corr.sel.v2, Pr.SA.new, col = "black", lwd = 3)
lines(corr.sel.v2, Pr.SA.expand, col="darkgreen", lwd = 3)
lines(corr.sel.v2, Pr.SA.est, col = "green3", lwd = 3)
lines(corr.sel.v2, Pr.SA.decline, col = "palegreen", lwd = 3)

# Plot main results for Symmetric displacements
points(x=ct$corr.sel,y=ct$Pr_SA_new, pch=21, col=adjustcolor("black",alpha=1),      bg=adjustcolor("black",alpha=0.8))
points(x=inc$corr.sel,y=inc$Pr_SA,   pch=21, col=adjustcolor("darkgreen",alpha=1),  bg=adjustcolor("darkgreen",alpha=0.8))
points(x=ct$corr.sel,y=ct$Pr_SA,     pch=21, col=adjustcolor("green",alpha=1),      bg=adjustcolor("green",alpha=0.8))
points(x=dec$corr.sel,y=dec$Pr_SA,   pch=21, col=adjustcolor("palegreen",alpha=1), bg=adjustcolor("palegreen",alpha=0.8))

# Plot Analytic approximations for Asymmetric displacements
num_pred_asym()
lines(corr.sel.v2[501:1000], Pr.SA.new[501:1000], col = "black", lwd = 3, lty=2)
lines(corr.sel.v2[501:1000], Pr.SA.expand[501:1000], col="darkgreen", lwd = 3, lty=2)
lines(corr.sel.v2[501:1000], Pr.SA.est[501:1000], col = "green3", lwd = 3, lty=2)
lines(corr.sel.v2[501:1000], Pr.SA.decline[501:1000], col = "palegreen", lwd = 3, lty=2)

# Approximations for completely aligned selection
x.B  <-  xBar(zBar=1.5)
x.M  <-  x.min(z.min=1)
pTradeoff.Aligned.new      <-  Pr.tradeoff.Aligned.New(xBar = x.B, x.min=x.M)
pTradeoff.Aligned.R1       <-  Pr.Tradeoff.Aligned.R1(xBar = x.B, x.min=x.M)
pTradeoff.Aligned.expand   <-  Pr.Tradeoff.Aligned.expand(R = 1.02, zBar = 1.5, xBar = x.B, x.min = x.M)
pTradeoff.Aligned.decline  <-  Pr.Tradeoff.Aligned.decline(R = 0.98, zBar = 1.5, xBar = x.B, x.min = x.M)
points(pTradeoff.Aligned.new ~ 1, pch=21, col=adjustcolor("black",alpha=1), bg=adjustcolor("black",alpha=1), cex=1.5)
points(pTradeoff.Aligned.expand ~ 1, pch=21, col=adjustcolor("darkgreen",alpha=1), bg=adjustcolor("darkgreen",alpha=1), cex=1.5)
points(pTradeoff.Aligned.R1 ~ 1, pch=21, col=adjustcolor("green",alpha=1), bg=adjustcolor("green",alpha=1), cex=1.5)
points(pTradeoff.Aligned.decline ~ 1, pch=21, col=adjustcolor("lightgreen",alpha=1), bg=adjustcolor("lightgreen",alpha=1), cex=1.5)
points(pTradeoff.Aligned.new ~ 1, pch=21, col=adjustcolor("white",alpha=1), bg=adjustcolor("white",alpha=1), cex=0.75)
points(pTradeoff.Aligned.expand ~ 1, pch=21, col=adjustcolor("white",alpha=1), bg=adjustcolor("white",alpha=1), cex=0.75)
points(pTradeoff.Aligned.R1 ~ 1, pch=21, col=adjustcolor("white",alpha=1), bg=adjustcolor("white",alpha=1), cex=0.75)
points(pTradeoff.Aligned.decline ~ 1, pch=21, col=adjustcolor("white",alpha=1), bg=adjustcolor("white",alpha=1), cex=0.75)

# Plot supp results for Asymmetric displacements
points(cos(tradeoffs_n50$theta)[-9],tradeoffs_n50$Pr_SA_new[-9],         pch=23, col=adjustcolor("black",alpha=1),     bg=adjustcolor("black",alpha=0.5))
points(cos(tradeoffs_R1.02_n50$theta)[-9],tradeoffs_R1.02_n50$Pr_SA[-9], pch=23, col=adjustcolor("darkgreen",alpha=1),  bg=adjustcolor("darkgreen",alpha=0.5))
points(cos(tradeoffs_n50$theta)[-9],tradeoffs_n50$Pr_SA[-9],             pch=23, col=adjustcolor("green",alpha=1), bg=adjustcolor("green",alpha=0.5))
points(cos(tradeoffs_R0.98_n50$theta)[-9],tradeoffs_R0.98_n50$Pr_SA[-9], pch=23, col=adjustcolor("lightgreen",alpha=1), bg=adjustcolor("lightgreen",alpha=0.5))

# Plot labels/annotations
proportionalLabel(0.5, 1.075, expression(paste("Symmetric vs. Asymmetric Displacements")),cex=1.4, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel(-0.175, 0.5, expression(paste("Proportion trade-off")),cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.175, expression(paste("Alignment of selection ", (cos(theta[sel])))),cex=1.2, adj=c(0.5, 0.5), xpd=NA)

#Legend
usr  <-  par('usr')
legend( x       =  usr[2]*0.8,
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
legend( x       =  usr[2]*0.7565,
        y       =  usr[4]*0.943,
        legend  =  c(" ", " ", " "),
        col     =   c('black'),
        pch     =  c(21, 23, 16),
        pt.bg   =  c(NA, NA, NA),
        pt.cex  =  c(1, 1, 1.5),
        cex     =  1,
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA)
legend( x       =  usr[2]*0.975,
        y       =  usr[4]*0.943,
        legend  =  c(expression(z[1]==z[2]),
                     expression(z[1]!=z[2]),
                     expression(cos(theta[sel]==1))),
        col     =   c('white'),
        pch     =  c(23, 16, 16),
        pt.bg   =  c(NA, NA, NA),
        pt.cex  =  c(0, 0, 0.75),
        cex     =  1,
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA)


# save the pdf
dev.off()