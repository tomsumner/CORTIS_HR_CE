# This script calculates the RR needed to match the TB incidence in the CORTIS cohort
# analysis is done for whole cohort, those with no IPT, and those with no IPT stratified by ART at baseline
#
# The model is run for n_samp sets of parameters sampled from distributions in "Par_dist_gen.R"
# The likelihood is calculated for each set of parameters
# n_resamp parameter sets are resampled using likelihood as weighting

setwd("~/GitHub/COR_repeat_screen")

# Load libraries
library(ggplot2)
library(reshape2)
library(deSolve)

# Load the functions used
source("Generate_transition_matrix.R") # function to generate transition matrix
source("Generate_ICs.R")               # function to generate initial conditions
source("Par_dist_gen.R")       

# Set the number of samples to do
n_samp <- 100000 # Number of parameter samples to generate
n_resamp <- 1000 # Number of resamples to draw

# Observed incidence in CORTIS HR (median and 95% quantiles)
obs_inc <- rbind(c(0.3,0.9,1.6),  # all
                 c(0.2,1.3,2.4),  # no IPT
                 c(0.0,0.3,1.0),  # no IPT, on ART at enrollment
                 c(0.1,4.0,7.0))  # no IPT, not on ART at enrollment
# Fit normal distributions to these quantiles
temp1 <- get.norm.par(p=c(0.025,0.5,0.975),q=exp(obs_inc[1,]),show.output=FALSE,plot=FALSE) 
temp2 <- get.norm.par(p=c(0.025,0.5,0.975),q=exp(obs_inc[2,]),show.output=FALSE,plot=FALSE) 
temp3 <- get.norm.par(p=c(0.025,0.5,0.975),q=exp(obs_inc[3,]),show.output=FALSE,plot=FALSE) 
temp4 <- get.norm.par(p=c(0.025,0.5,0.975),q=exp(obs_inc[4,]),show.output=FALSE,plot=FALSE) 
inc_pars <- rbind(temp1,temp2,temp3,temp4)

# vectors to store parameters
spec_L_IGRA <- rep(0,n_samp)   # specificity of IGRA for infection
sens_L_IGRA <- spec_L_IGRA     # sensitivity of IGRA for infection
prev_L <- spec_L_IGRA          # prevalence of infection
ARI <- spec_L_IGRA             # Annual risk of infection    
RR <- spec_L_IGRA              # RR for TB in HIV+
p <- spec_L_IGRA               # protectiv effect of prior infection

# Set up some parameters that stay constant                                     
cum_TB <- c(0,0.05341,0.07431,0.08268,0.08624,0.08794) # cumulative developing TB by year post infection
temp <- diff(cum_TB)/4                                 # calculate proportion per year
tds <- rep(temp,each=4)                                # distribute across 3 monthly periods (assuming uniform)
pop_size <- 10000                                      # cohort size
LtoL20s <- (1/((1/0.000594)-5)/4)                      # rate of progression from L to L20
t_run <- 5                                             # in 3 monthly cycles - look at 1 year for incidence
a_age <- 35                                            # average age of cohort
prev_I <- 0.456                                        # prevalence of IGRA+ from trial 
E <-0                                                  # efficacy of treatment - not used here
ER <- 0

# Array to store incidence from model
Incidence <- array(0,c(t_run,n_samp))

# Run the model - for each run, sample the IGRA parameters (these determine the ARI)
for (nn in 1:n_samp){
  
  RR[nn] <- runif(1,1,50)   # randomly sample the RR for HIV+ TB
  
  spec_L_IGRA[nn] = rbeta(1,fit_spec_IGRA_LTBI[1],fit_spec_IGRA_LTBI[2])  # sample IGRA characteristics
  sens_L_IGRA[nn] = rbeta(1,fit_sens_IGRA_LTBI[1],fit_sens_IGRA_LTBI[2])   

  prev_L[nn] <- min(1,(prev_I - 1 + spec_L_IGRA[nn])/(sens_L_IGRA[nn] - 1 + spec_L_IGRA[nn]))  # calculate infetion prevalence
  ARI[nn] <- 100*-log(1-prev_L[nn])/a_age           # calculate ARI
  
  # sample the protection from prior infection
  p[nn] <- rbeta(1,fit_p[1],fit_p[2])   
  
  # 
  A <- ARI[nn]/(100*4)      # risk of infection (adjusted to 3-monthly)
  R <- A*p[nn]              # risk of reinfection 
  td <- RR[nn]*tds          # proportion entering each infection state
  LtoL20 <- RR[nn]*LtoL20s  # rate of transtion from remote infection
  Ld=1-sum(td)              # proportion entering remote state after infection 
  
  initial_cond <- rep(0,22)
  # Calculate initial conditons for the model
  # Assumes ARI has been constant for last 5 years
  # ******* note doesn't account for change in L and S over previous 5 years ********#
  #S, L1, ..., L20, L
  initial_cond[1] <- 1-prev_L[nn]
  initial_cond[21] <- A*((1-prev_L[nn]) + p*prev_L[nn])*td[20] + LtoL20*prev_L[nn] 
  for (ic in 20:2){
    initial_cond[ic] <- initial_cond[ic+1]+A*((1-prev_L[nn]) + p*prev_L[nn])*td[ic-1]
  }
  initial_cond[22] <- 1-sum(initial_cond[1:21])
  initial_cond <- initial_cond*pop_size
    
  # these parameters determine transitions on to PT for repeat screening - set to zero here to fit baseline incidence
  covP_S = 0  # in susceptible population
  covP_L = 0  # in latent population
  covP_I = 0  # in progressors
  # generate the initial conditions for the model
  init <- IC_gen(initial_cond,0,0,0,E)
  # create matrices to store model states
  run_base <- mat.or.vec(length(init),t_run)

  # Then iterate through time
  for(ii in 1:t_run){
      
    # set transition matrices to no PT version
    mat_base <- mat_gen(0,0,0,E,ER,A,R,length(init))

    if (ii==1){
        run_base[,ii] <- init%*%mat_base
    }
    if (ii > 1){
        run_base[,ii] <- run_base[,ii-1]%*%mat_base
    }
  }
  # Combine and rearrange results
  run_base <- as.data.frame(cbind(t(run_base),seq(1,t_run)))
  colnames(run_base) <- c(names(init),"time")
  outputs <- run_base
  outputs_m <- melt(outputs,id.vars=c("time"))
    
  # Calculate incidence 
  for (ii in 1:t_run){
    temp <- as.numeric(as.character(outputs_m[outputs_m$time==ii&outputs_m$variable=="had_TB",3]))
    Incidence[ii,nn] <- exp(100*as.numeric(as.character(temp))/10000)
  }
    
}

# calculate likelihood for each run 
L <- mat.or.vec(4,n_samp)
L_w <- L
ssd <- L
ssd_fit <- rep(0,4)
for (i in 1:n_samp){
  # All
  L[1,i] <- ((2*pi*(inc_pars[1,2]^2))^(-1/2))*exp(-(1/(2*(inc_pars[1,2]^2)))*
                                                    ((log(Incidence[5,i])-log(inc_pars[1,1]))^2))
  # No IPT
  L[2,i] <- ((2*pi*(inc_pars[2,2]^2))^(-1/2))*exp(-(1/(2*(inc_pars[2,2]^2)))*
                                                    ((log(Incidence[5,i])-log(inc_pars[2,1]))^2))
  # No IPT, ART
  L[3,i] <- ((2*pi*(inc_pars[3,2]^2))^(-1/2))*exp(-(1/(2*(inc_pars[3,2]^2)))*
                                                    ((log(Incidence[5,i])-log(inc_pars[3,1]))^2))
  # No IPT, no ART
  L[4,i] <- ((2*pi*(inc_pars[4,2]^2))^(-1/2))*exp(-(1/(2*(inc_pars[4,2]^2)))*
                                                    ((log(Incidence[5,i])-log(inc_pars[4,1]))^2))
  
}
# remove na and normalise
L[is.nan(L)] <- 0
L_w[1,] <- L[1,]/sum(L[1,])
L_w[2,] <- L[2,]/sum(L[2,])
L_w[3,] <- L[3,]/sum(L[3,])
L_w[4,] <- L[4,]/sum(L[4,])

# resample based on replacement
i_samp <- mat.or.vec(4,n_resamp)
i_samp[1,] <- sample(seq(1,n_samp),n_resamp,replace=TRUE,prob=L_w[1,])
i_samp[2,] <- sample(seq(1,n_samp),n_resamp,replace=TRUE,prob=L_w[2,])
i_samp[3,] <- sample(seq(1,n_samp),n_resamp,replace=TRUE,prob=L_w[3,])
i_samp[4,] <- sample(seq(1,n_samp),n_resamp,replace=TRUE,prob=L_w[4,])

# plot model incidence and fitted incidence
inc_out_L <- rbind(quantile(log(Incidence[5,i_samp[1,]]),probs=c(0.025,0.5,0.975)),
                 quantile(log(Incidence[5,i_samp[2,]]),probs=c(0.025,0.5,0.975)),
                 quantile(log(Incidence[5,i_samp[3,]]),probs=c(0.025,0.5,0.975)),
                 quantile(log(Incidence[5,i_samp[4,]]),probs=c(0.025,0.5,0.975)))

dat <- as.data.frame(rbind(cbind(obs_inc,c(1,2,3,4),"Trial"),
                           cbind(inc_out_L,c(1,2,3,4),"Model")))
colnames(dat) <- c("min","med","max","S","type")
plot_inc <- ggplot(dat,aes(x=S,y=as.numeric(as.character(med)),
                           ymin=as.numeric(as.character(min)),
                           ymax=as.numeric(as.character(max)),
                           color=as.factor(type)))+
  geom_errorbar(width=0.2,position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5))+
  expand_limits(y=0)

# get ranges for fitted things
RR_out <- rbind(quantile(RR[i_samp[1,]],probs=c(0.025,0.5,0.975)),
                quantile(RR[i_samp[2,]],probs=c(0.025,0.5,0.975)),
                quantile(RR[i_samp[3,]],probs=c(0.025,0.5,0.975)),
                quantile(RR[i_samp[4,]],probs=c(0.025,0.5,0.975)))

p_out <- rbind(quantile(p[i_samp[1,]],probs=c(0.025,0.5,0.975)),
               quantile(p[i_samp[2,]],probs=c(0.025,0.5,0.975)),
               quantile(p[i_samp[3,]],probs=c(0.025,0.5,0.975)),
               quantile(p[i_samp[4,]],probs=c(0.025,0.5,0.975)))

ARI_out <- rbind(quantile(ARI[i_samp[1,]],probs=c(0.025,0.5,0.975)),
               quantile(ARI[i_samp[2,]],probs=c(0.025,0.5,0.975)),
               quantile(ARI[i_samp[3,]],probs=c(0.025,0.5,0.975)),
               quantile(ARI[i_samp[4,]],probs=c(0.025,0.5,0.975)))

prev_L_out <- rbind(quantile(prev_L[i_samp[1,]],probs=c(0.025,0.5,0.975)),
                    quantile(prev_L[i_samp[2,]],probs=c(0.025,0.5,0.975)),
                    quantile(prev_L[i_samp[3,]],probs=c(0.025,0.5,0.975)),
                    quantile(prev_L[i_samp[4,]],probs=c(0.025,0.5,0.975)))

spec_L_IGRA_out <- rbind(quantile(spec_L_IGRA[i_samp[1,]],probs=c(0.025,0.5,0.975)),
                    quantile(spec_L_IGRA[i_samp[2,]],probs=c(0.025,0.5,0.975)),
                    quantile(spec_L_IGRA[i_samp[3,]],probs=c(0.025,0.5,0.975)),
                    quantile(spec_L_IGRA[i_samp[4,]],probs=c(0.025,0.5,0.975)))

sens_L_IGRA_out <- rbind(quantile(sens_L_IGRA[i_samp[1,]],probs=c(0.025,0.5,0.975)),
                         quantile(sens_L_IGRA[i_samp[2,]],probs=c(0.025,0.5,0.975)),
                         quantile(sens_L_IGRA[i_samp[3,]],probs=c(0.025,0.5,0.975)),
                         quantile(sens_L_IGRA[i_samp[4,]],probs=c(0.025,0.5,0.975)))