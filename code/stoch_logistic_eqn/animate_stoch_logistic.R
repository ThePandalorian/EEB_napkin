############################
# Create an animation comparing stochastic logistic eqn with deterministic limit
# Written by Shikhara Bhat
# IISER Pune, India
# Date: Wed Apr 05 16:06:13 2023
###########################

library(GillespieSSA2) #for stochastic simulations
library(deSolve) #for numerical integration of deterministic equations
library(ggplot2) #for plotting
library(gganimate) #for animating ggplots

set.seed(1414) #for reproducibility

#parameters for stochastic sims
b = 2 #birth rate
d = 1 #death rate
K = 100 #carrying capacity
N0 = 30


#for convenience, define difference and sum of birth and death rates
r=b-d
v=b+d

steps = 10000 #step size for numerical integration

model_params <- c(b = b, d = d, K = K)      # Parameter list
final_time <- 50                                   # Final time

initial_state <- c(N = N0) #initial population

#update rules
reactions <- list(
  
  #               rate                    effect
  
  
  #birth
  reaction(      "b*N"      ,          c(N = +1)),
  
  #death
  reaction("(d+(b-d)*N/K)*N",          c(N = -1))
)

#Run the Gillespie simulation of the IbM

stoch_data <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = model_params,
  final_time = final_time,
  method = ssa_exact(),
  census_interval = 0, #how often to store data. 0 means store everything.
  sim_name = 'Stoch Logistic'
)


stoch_data <- as.data.frame(cbind(stoch_data[['time']], stoch_data[['state']], rep("Stochastic simulation", length(stoch_data[['time']]))))
colnames(stoch_data) <- c('time','N','state')
#############################################################################

#The deterministic logistic equation
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
timelist <- seq(0,final_time,by=final_time/steps) #times for integration

#numerically integrate the logistic eqn for the deterministic trajectory
det_data <- as.data.frame(ode(N0,timelist,func=logistic,parms=params_deterministic))
#convert output to dataframe
det_data <- cbind(det_data,state=rep(('Deterministic logistic equation'),nrow(det_data)))
colnames(det_data) <- c('time','N','state')

################################################################################
#plotting

#collate the data for plotting
plot_data <- rbind(det_data,stoch_data)

#make sure everything is the right data type
plot_data$time <- as.numeric(plot_data$time)
plot_data$N <- as.numeric(plot_data$N)

#create the plot
p <- ggplot(data=plot_data,aes(x=time,y=N,group=state)) + geom_line(aes(linetype=state,color=state,alpha=state))
p <- p + scale_linetype_manual(name='',values=c("dashed","solid")) + scale_color_manual(name='',values=c('#000000','#ff0000'))
p <- p + scale_alpha_manual(name='',values=c(1,0.6))

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


#for saving animations
#this will raise an error if the required directory doesn't exist
setwd(paste('./animation/K_',K,sep=''))

#animate the plot
anim <- p + transition_reveal(time)
animate(anim,nframes=200)


