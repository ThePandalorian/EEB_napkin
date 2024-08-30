############################
# Compare the steady state densities of Gillespie, FPE, and Weak noise approx
# for the stochastic logistic equation.
# Written by Shikhara Bhat
# IISER Pune, India
# Date: Thu Jan 05 21:50:16 2023
###########################

library(GillespieSSA) #for exact stochastic simulations using Gillespie
library(deSolve) #for numerical integration of deterministic equations
library(sde) #for numerical integration of SDEs
library(ggplot2) #for plotting

set.seed(31415) #for reproducibility

#Model
b = 2 #birth rate
d = 1 #death rate
K = 10000 #carrying capacity

#for convenience, define difference and sum of birth and death rates
r=b-d
v=b+d

params_stoch <- c(b = b, d = d, K = K)       # Parameter list
tfinal <- 200                                # Final time
simName <- "Logistic growth (Gillespie)" 

#For the SSA
x0 <- c(N = K) #initial population. Start population at deterministic steady state (N=K)
nu <- matrix(c(+1, -1),ncol = 2) #transition matrix (specifies possible transitions)
a  <- c("b*N", "(d+(b-d)*N/K)*N") #update rules

#For the complete Fokker-Planck
SDE_drift <- expression(r*x*(1-x)) 
SDE_diffusion <- expression(sqrt((x*(v+r*x))/K))


#for the WNA. I have substituted alpha(t) = 1 (deterministic steady state)
WNA_drift <- expression(-r*x)
WNA_diffusion <- expression(sqrt(v+r))


steps = 100000 #step size for numerical integration

runs = 1000 #number of realizations to average over

sim_start_time <- Sys.time() #just to see how long code takes to run
SSA <- SDE <- WNA <- rep(0,runs)

for (i in 1:runs) {
  
  #Run the Gillespie algorithm (time intensive)
  stoch_out <- ssa(
    x0 = x0, #initial condition
    a = a,   #transition rates
    nu = nu, #possible transitions
    parms = params_stoch, #parameters
    tf = tfinal, #final time uptil which to run the model
    method =ssa.otl(), #simulation method. otl is optimized tau leap.
    simName = '', #unimportant
    verbose = FALSE, #dont print the output 
    consoleInterval = 1 #unimportant
  ) 
  
  SSA[i] <- as.data.frame(stoch_out[[1]])$N[tfinal]
  
  #The non-linear SDE and WNA are both in terms of pop density.
  #We need to multiply the output of both by K to convert back to pop size.
  
  #numerically integrate the non-linear SDE
  SDE[i] <- K*suppressMessages(sde.sim(
    X0=1, #initial condition
    T=tfinal, #time until which to integrate
    drift=SDE_drift, #drift term
    sigma=SDE_diffusion, #diffusion term
    N=steps))[steps+1] #partition size of the interval [0,T]). We only want the last entry, hence the [steps+1]
  
  #numerically integrate the WNA eqn
  WNA[i] <- K*(1 + (1/sqrt(K))*(suppressMessages(sde.sim( #WNA measures fluctuation y = sqrt(K)(x-alpha). We want x here.
    X0=0,
    T=tfinal,
    drift=WNA_drift,
    sigma=WNA_diffusion,
    N=steps))[steps+1]))
}

running_time = Sys.time() - sim_start_time
print(running_time) #see how long simulations took

################################################################################
#plotting

#collate the data
plot_data <- data.frame(sim=c(rep('Stochastic IbM (simulations)',length(SSA)),rep('Fokker-Planck equation (theory)',length(SDE)),rep('Weak Noise Approximation (theory)',length(WNA))),N=c(SSA,SDE,WNA))

#create the plot
p <- ggplot(data=plot_data,aes(x=N,group=sim,color=sim)) + geom_density(aes(linetype=sim,color=sim),adjust=1.7) + xlim(c(8500,11000)) + ylim(c(0,0.0035))
p <- p + scale_linetype_manual(name='',values=c("solid","dashed","dashed")) + scale_color_manual(name='',values=c('#000000','#ff0000','#0000ff'))

#aesthetics
p <- p + theme_light() + ylab('Steady state density') + xlab("Population size (N)")
p <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text = element_text(face = 'bold', color = 'black',size = 15))
p <- p + theme(axis.title = element_text(face = 'bold',color = 'black',size = 15))
p <- p + theme(legend.title = element_text(face = 'bold',color = 'black',size = 12))
p <- p + theme(legend.text = element_text(face = 'bold',color = 'black',size = 12))
p <- p + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',size=0.5))

#save to file
#ggsave(paste('./plots/1D_K_',K,'.svg',sep=''),p,width=12,height=8)


