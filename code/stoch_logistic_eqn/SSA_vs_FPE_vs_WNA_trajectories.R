############################
# Compare time series of Gillespie, Fokker-Planck, and weak noise approximation
# for the stochastic logistic equation.
# Written by Shikhara Bhat
# IISER Pune, India
# Date: Wed Jan 11 13:48:29 2023
###########################

library(GillespieSSA) #for stochastic simulations
library(deSolve) #for numerical integration of deterministic equations
library(sde) #for numerical integration of SDEs
library(ggplot2) #for plotting


set.seed(31415) #for reproducibility

#parameters for stochastic sims
b = 2 #birth rate
d = 1 #death rate
K = 10000 #carrying capacity


#for convenience, define difference and sum of birth and death rates
r=b-d
v=b+d

steps = 10000 #step size for numerical integration

params_stoch <- c(b = b, d = d, K = K)      # Parameter list
tfinal <- 50                                # Final time

x0 <- c(N = 100) #initial population
nu <- matrix(c(+1, -1),ncol = 2) #transition matrix (specifies possible transitions)
a  <- c("b*N", "(d+(b-d)*N/K)*N") #update rules

#Run the Gillespie simulation of the IbM
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

stoch_data <- as.data.frame(stoch_out[[1]]) #get output as dataframe
stoch_data <- cbind(stoch_data,state=rep('Exact stochastic IbM (SSA)',nrow(stoch_data)))

#############################################################################

logistic <- function(time,N,params){
  
  with(as.list(c(N, params)), {
    # "with" creates an "environment" within which R looks for values
    
    # our differential equation for growth
    dNdt <- r * N * (1 - (N/K))
    
    # making the function create output, which must be a list with
    # one element
    return(list( dNdt ))
  } )
}

params_deterministic <- c(r=r, K = K)
timelist <- seq(0,tfinal,by=tfinal/steps) #times for integration

#numerically integrate the logistic eqn for the deterministic trajectory
det_data <- as.data.frame(ode(x0,timelist,func=logistic,parms=params_deterministic))

#convert output to dataframe
names(det_data)[1] <- 't'
det_data <- cbind(det_data,state=rep('Deterministic trajectory',nrow(det_data)))

################################################################################
# Exact solution of the logistic eqn (for WNA)

alpha <- function(t){
  
  Yt <- K*x0/(x0+(K-x0)*exp(-r*t)) #exact solution (by integrating ODE)
  
  return (Yt)
}

################################################################################
#SDEs
dens0 <- x0/K

#For the complete Fokker-Planck
SDE_drift <- expression(r*x*(1-x)) 
SDE_diffusion <- expression(sqrt((x*(v+r*x))/K))


#for the WNA
WNA_drift <- expression(r*(1-2*(dens0/(dens0+(1-dens0)*exp(-r*t)))*x))
WNA_diffusion <- expression(sqrt((dens0/(dens0+(1-dens0)*exp(-r*t)))*(v+r*(dens0/(dens0+(1-dens0)*exp(-r*t))))))


#numerically integrate the non-linear SDE corresponding to the full Fokker-Planck eqn
SDE <- K*suppressMessages(sde.sim(
                                  X0=dens0, #initial condition
                                  T=tfinal, #time until which to integrate
                                  drift=SDE_drift, #drift term
                                  sigma=SDE_diffusion, #diffusion term
                                  N=steps)) #partition size of the interval [0,T]

#numerically integrate the linear SDE corresponding to the weak noise approximation
WNA_fluc <- (1/sqrt(K))*(suppressMessages(sde.sim(
                                  X0=0,
                                  T=tfinal,
                                  drift=WNA_drift,
                                  sigma=WNA_diffusion,
                                  N=steps)))

#WNA computes flucuations from the deterministic trajectory
#We want the complete dynamics, so add the fluctuations to the deterministic traj
#i.e. WNA measures the fluctuation y = sqrt(K)(x-alpha). We want x here. Transpose.
#we then get x = alpha + y/sqrt(K). Now convert from densities to pop numbers.
WNA <- rep(0,length(WNA_fluc))
for (i in 1:length(WNA)){
  WNA[i] <- alpha(timelist[i]) + K*as.vector(WNA_fluc)[i] #deterministic solution + scaled_fluctuations
}

SDE_data <- data.frame(cbind(t=timelist, N=as.vector(SDE),state=rep('Fokker-Planck equation (theory)',length(SDE))))
WNA_data <- data.frame(cbind(t=timelist, N=as.vector(WNA),state=rep('Weak Noise Approximation (theory)',length(SDE))))


################################################################################
#plotting

#collate the data for plotting
plot_data <- rbind(det_data,stoch_data,SDE_data,WNA_data)

#make sure everything is the right data type
plot_data$t <- as.numeric(plot_data$t)
plot_data$N <- as.numeric(plot_data$N)

#create the plot
p <- ggplot(data=plot_data,aes(x=t,y=N,group=state)) + geom_line(aes(linetype=state,color=state,alpha=state))
p <- p + scale_linetype_manual(name='',values=c("dashed", "solid","solid","solid")) + scale_color_manual(name='',values=c('#000000','#00ff00','#ff0000','#0000ff'))
p <- p + scale_alpha_manual(name='',values=c(1,0.6,0.6,0.6))

#aesthetics
p <- p + theme_light() + ylab('Population size') + xlab("Time")
p <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text = element_text(face = 'bold', color = 'black',size = 15))
p <- p + theme(axis.title = element_text(face = 'bold',color = 'black',size = 15))
p <- p + theme(legend.title = element_text(face = 'bold',color = 'black',size = 12))
p <- p + theme(legend.text = element_text(face = 'bold',color = 'black',size = 12))
p <- p + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',size=0.5))

#display the plot
p

#save to file
#ggsave(paste('./plots/timeseries_1D_K_',K,'.svg',sep=''),p,width=12,height=6,dpi=600)
