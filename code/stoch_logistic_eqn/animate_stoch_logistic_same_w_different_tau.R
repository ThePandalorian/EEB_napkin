############################
# Create an animation comparing stochastic logistic eqn for same fitness but
# different turnover rates
# Written by Shikhara Bhat
# IISER Pune, India
# Date: Wed Mar 28 13:44:35 2023
###########################

library(GillespieSSA) #for stochastic simulations
library(deSolve) #for numerical integration of deterministic equations
library(ggplot2) #for plotting
library(gganimate) #for animating ggplots

setwd('./animation/logistic_sim_same_w_different_tau/')

set.seed(2878) #for reproducibility

#parameters for stochastic sims
b1 = 4 #birth rate 1
d1 = 2 #death rate 1
b2 = 31 #birth rate 2
d2 = 29 #death rate 2
K1 = 150 #carrying capacity 1
K2 = 150 #carrying capacity 2
runs = 10 #number of realizations

#for convenience, define difference and sum of birth and death rates
r1=b1-d1
v1=b1+d1

r2=b2-d2
v2=b2+d2

steps = 10000 #step size for numerical integration

params_stoch_1 <- c(b = b1, d = d1, K = K1)      # Parameter list
params_stoch_2 <- c(b = b2, d = d2, K = K2)      # Parameter list
tfinal <- 30                                   # Final time

x0 <- c(N = 50) #initial population
nu <- matrix(c(+1, -1),ncol = 2) #transition matrix (specifies possible transitions)
a  <- c("b*N", "(d+(b-d)*N/K)*N") #update rules

#Run the Gillespie simulation of the IbM

stoch_out_1 <- ssa(
  x0 = x0, #initial condition
  a = a,   #transition rates
  nu = nu, #possible transitions
  parms = params_stoch_1, #parameters
  tf = tfinal, #final time uptil which to run the model
  method =ssa.otl(), #simulation method. otl is optimized tau leap.
  simName = '', #unimportant
  verbose = FALSE, #dont print the output 
  consoleInterval = 1 #unimportant
) 

stoch_out_2 <- ssa(
  x0 = x0, #initial condition
  a = a,   #transition rates
  nu = nu, #possible transitions
  parms = params_stoch_2, #parameters
  tf = tfinal, #final time uptil which to run the model
  method =ssa.otl(), #simulation method. otl is optimized tau leap.
  simName = '', #unimportant
  verbose = FALSE, #dont print the output 
  consoleInterval = 1 #unimportant
) 


stoch_data_1 <- as.data.frame(stoch_out_1[[1]]) #get output as dataframe
stoch_data_1 <- cbind(stoch_data_1,state=rep('Stochastic dynamics (model 1)',nrow(stoch_data_1)))
stoch_data_2 <- as.data.frame(stoch_out_2[[1]]) #get output as dataframe
stoch_data_2 <- cbind(stoch_data_2,state=rep('Stochastic dynamics (model 2)',nrow(stoch_data_2)))

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

params_deterministic_1 <- c(r=r1, K = K1)
params_deterministic_2 <- c(r=r2, K = K2)
timelist <- seq(0,tfinal,by=tfinal/steps) #times for integration

#numerically integrate the logistic eqn for the deterministic trajectory
det_data_1 <- as.data.frame(ode(x0,timelist,func=logistic,parms=params_deterministic_1))
det_data_2 <- as.data.frame(ode(x0,timelist,func=logistic,parms=params_deterministic_2))

#convert output to dataframe
names(det_data_1)[1] <- 't'
det_data_1 <- cbind(det_data_1,state=rep(('Deterministic limit (model 1)'),nrow(det_data_1)))
names(det_data_2)[1] <- 't'
det_data_2 <- cbind(det_data_2,state=rep(('Deterministic limit (model 2)'),nrow(det_data_2)))

################################################################################
#plotting

#collate the data for plotting
plot_data <- rbind(det_data_1,det_data_2,stoch_data_1,stoch_data_2)

#make sure everything is the right data type
plot_data$t <- as.numeric(plot_data$t)
plot_data$N <- as.numeric(plot_data$N)

#create the plot
p <- ggplot(data=plot_data,aes(x=t,y=N,group=state)) + geom_line(aes(linetype=state,color=state,alpha=state))
p <- p + scale_linetype_manual(name='',values=c("dashed","dashed","solid", "solid")) + scale_color_manual(name='',values=c('#000000','#000000','#ff0000','#0000ff'))
p <- p + scale_alpha_manual(name='',values=c(1,1,0.2,0.2))

#aesthetics
p <- p + theme_light() + ylab('Population size') + xlab("Time")
p <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text = element_text(face = 'bold', color = 'black',size = 15))
p <- p + theme(axis.title = element_text(face = 'bold',color = 'black',size = 15))
#p <- p + theme(legend.title = element_text(face = 'bold',color = 'black',size = 12))
#p <- p + theme(legend.text = element_text(face = 'bold',color = 'black',size = 12))
p <- p + theme(legend.position = 'none') #remove legend
p <- p + theme(aspect.ratio = 1) #set aspect ratio
p <- p + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',linewidth=0.5))

for (i in 1:(runs-1)){
  
  stoch_out_1 <- ssa(
    x0 = x0, #initial condition
    a = a,   #transition rates
    nu = nu, #possible transitions
    parms = params_stoch_1, #parameters
    tf = tfinal, #final time uptil which to run the model
    method =ssa.otl(), #simulation method. otl is optimized tau leap.
    simName = '', #unimportant
    verbose = FALSE, #dont print the output 
    consoleInterval = 1 #unimportant
  ) 
  
  stoch_out_2 <- ssa(
    x0 = x0, #initial condition
    a = a,   #transition rates
    nu = nu, #possible transitions
    parms = params_stoch_2, #parameters
    tf = tfinal, #final time uptil which to run the model
    method =ssa.otl(), #simulation method. otl is optimized tau leap.
    simName = '', #unimportant
    verbose = FALSE, #dont print the output 
    consoleInterval = 1 #unimportant
  ) 
  
  
  stoch_data_1 <- as.data.frame(stoch_out_1[[1]]) #get output as dataframe
  stoch_data_1 <- cbind(stoch_data_1,state=rep('Stochastic dynamics (model 1)',nrow(stoch_data_1)))
  stoch_data_2 <- as.data.frame(stoch_out_2[[1]]) #get output as dataframe
  stoch_data_2 <- cbind(stoch_data_2,state=rep('Stochastic dynamics (model 2)',nrow(stoch_data_2)))
  
  #collate the data for plotting
  plot_data <- rbind(stoch_data_1,stoch_data_2)
  
  #make sure everything is the right data type
  plot_data$t <- as.numeric(plot_data$t)
  plot_data$N <- as.numeric(plot_data$N)
  
  #create the plot
  p <- p + geom_line(data=plot_data,aes(x=t,y=N,group=state,linetype=state,color=state,alpha=state))
 # p <- p + scale_linetype_manual(aes(x=t,y=N,group=state),name='',values=c("solid", "solid")) + scale_color_manual(aes(x=t,y=N,group=state),name='',values=c('#0000ff','#ff0000'))
 # p <- p + scale_alpha_manual(aes(x=t,y=N,group=state),name='',values=c(0.2,0.2))
}

#for saving animations - species 1
setwd('./slower_species')

p <- p + scale_alpha_manual(name='',values=c(1,0,0.2,0))
#animate the plot
anim <- p + transition_reveal(t)
animate(anim,nframes=200)

#for saving animations - species 2
setwd('./faster_species')

p <- p + scale_alpha_manual(name='',values=c(0,1,0,0.2))
#animate the plot
anim <- p + transition_reveal(t)
animate(anim,nframes=200) 

#for saving animations - both species overlaid
setwd('./combined')

p <- p + scale_alpha_manual(name='',values=c(0,1,0.3,0.3))
#animate the plot
anim <- p + transition_reveal(t)
animate(anim,nframes=200) 

