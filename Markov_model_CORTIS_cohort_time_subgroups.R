# This script implements a markov model for looking at different preventive therapy strategies
# It is set up to track time till TB to allow us to model CORTIS like tests
#
# This version simulates the ART naive, no IPT population (subset 4) assuming they start ART at t =0 and become like on ART population (subset 3) by t = 5

# set the working directory
setwd("~/GitHub/COR_repeat_screen")

# Load libraries
library(ggplot2)
library(reshape2)
library(deSolve)
library(epiR)

# load the fitted parameter sets already generated - need to have run "CORTIS_inc_check.R" first
load("C:/Users/eidetsum/Desktop/COR_repeat_data_CORTIS_fit.RData")

# Load the functions used
source("Generate_transition_matrix.R") # function to generate transition matrix
source("Generate_ICs.R")               # function to generate initial conditions
source("Par_dist_gen.R")               # function to fit distributions to parameter ranges - sens and spec

# Set the subset to use for the first year (R_sub1) and subsequent years (R_sub2)
# Alters RR and sens and spec of RISK11 at t = 5
# Only set different values to simulate ART naive cohort starting ART (R_sub1 = 4, R_sub2 = 3)
# 1 = all, 2 = no IPT, 3 = ART, no IPT, 4 = no ART, no IPT
R_sub1 <- 1
R_sub2 <- 1

# results in paper use n_resamp = 1000
n_resamp <- 1000

# efficacy of PT
E_base_v <- c(1)            # efficacy of PT - relative results don't depend on this 
E_RR_v <- seq(0,1,0.01)     # relative efficacy of PT in RISK11+ - primary analysis uses 1

# time points to store outputs
t_store <- c(5,20,40) # 1, 5, 10 years

# Set the colour palette - accessible
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# number of strategies (none, all, RISK11, opt, minimum)
n_strat <- 5    

# Set up some parameters that stay constant                                     
cum_TB <- c(0,0.05341,0.07431,0.08268,0.08624,0.08794) # cumulative developing TB by year post infection
temp <- diff(cum_TB)/4                                 # calculate proportion per year
tds <- rep(temp,each=4)                                # distribute across 3 monthly periods (assuming uniform)
pop_size <- 10000                                      # cohort size
LtoL20s <- (1/((1/0.000594)-5)/4)                      # rate of progression from L to L20
y_Markov <- 10                                         # run Markov model for 10 years
t_Markov <- y_Markov*4                                 # in 3 monthly cycles
a_age <- 35                                            # average age of cohort
prev_I <- 0.456                                        # prevalence of IGRA+ from trial 

# vectors for sampled parameters
spec_L_IGRA_r <- rep(0,n_resamp)
sens_L_IGRA_r <- spec_L_IGRA_r
sens_opt <- spec_L_IGRA_r   
spec_opt <- spec_L_IGRA_r 
sens_min <- spec_L_IGRA_r   
spec_min <- spec_L_IGRA_r 
prev_L <- spec_L_IGRA_r
ARI <- spec_L_IGRA_r
p_r <- spec_L_IGRA_r

# These can change over time
spec_RISK11 <- array(0,c(t_Markov,n_resamp))
sens_RISK11 <- spec_RISK11
RR_r <- spec_RISK11

# These are vectors of when to repeat screening
screen_at <- cbind(c(1,rep(t_Markov+9999,19)),                           # once
                   c(1,13,25,37,rep(t_Markov+9999,16)),                  # every 3 yrs
                   c(1,9,17,25,33,rep(t_Markov+9999,15)),                # every 2 yrs
                   c(1,5,9,13,17,21,25,29,33,37,rep(t_Markov+9999,10)))  # every yr

N_runs <- dim(screen_at)[2]*length(E_RR_v)*length(E_base_v)
screen <- rep(0,N_runs)
E_base <- screen
E_RR <- screen

# Create arrays for storing outputs
Cases <- array(0,c(N_runs,length(t_store),n_strat,n_resamp))
Cases_avert <- Cases
Diff_Cases_avert <- Cases
Ratio_Cases_avert <- Cases
HR <- Cases
Cum_screen <- Cases
Cum_PT <- Cases
NNT <- Cases
g_cost <- Cases
g_per <- Cases

# Run the model - for each run, sample the test parameters (these also affect the ARI)
for (nn in 1:n_resamp){
  
  # Get parameters from fitting to incidence 
  spec_L_IGRA_r[nn] = spec_L_IGRA[i_samp[R_sub1,nn]]
  sens_L_IGRA_r[nn] = sens_L_IGRA[i_samp[R_sub1,nn]]
  
  # try having a declining RR - drop by 50% over whole period
  #RR_1 <- RR[i_samp[R_sub1,nn]]
  #RR_d <- 1
  #tt <- approx(c(1,t_Markov),c(RR_1,RR_d*RR_1),seq(1,t_Markov))
  #RR_r[,nn] = tt$y
  RR_r[1:4,nn] = RR[i_samp[R_sub1,nn]]
  RR_r[5:t_Markov,nn] = RR[i_samp[R_sub2,nn]]

  # Get test parameters
  # RISK11
  spec_RISK11[1:4,nn] = rbeta(1,fit_spec[R_sub1,1],fit_spec[R_sub1,2])  
  sens_RISK11[1:4,nn] = rbeta(1,fit_sens[R_sub1,1],fit_sens[R_sub1,2]) 
  spec_RISK11[5:t_Markov,nn] = rbeta(1,fit_spec[R_sub2,1],fit_spec[R_sub2,2])  
  sens_RISK11[5:t_Markov,nn] = rbeta(1,fit_sens[R_sub2,1],fit_sens[R_sub2,2]) 
  # Minimum and optimum TPP
  sens_opt[nn] = 0.9     
  spec_opt[nn] = 0.9     
  sens_min[nn] = 0.75     
  spec_min[nn] = 0.75
  
  # Get p from distribution
  p_r[nn] <- p[i_samp[R_sub1,nn]]
    
  # Calculate prevalence of infection and ARI
  prev_L[nn] <- min(1,(prev_I - 1 + spec_L_IGRA_r[nn])/(sens_L_IGRA_r[nn] - 1 + spec_L_IGRA_r[nn]))
  ARI[nn] <- 100*-log(1-prev_L[nn])/a_age

  # index for model run
  kk <- 1
  initial_cond <- rep(0,22)

  for (freq_screen in 1:dim(screen_at)[2]){
  
    for (rr in 1:length(E_RR_v)){
      
      for (ee in 1:length(E_base_v)){
      
        screen[kk] <- freq_screen
        E_base[kk] <- E_base_v[ee]
        E_RR[kk] <- E_RR_v[rr]
        E <- E_base[kk]
        ER <- E_RR[kk]*E_base[kk]

        t_screen <- screen_at[,freq_screen]
        A <- ARI[nn]/(100*4)    # risk of infection
        R <- A*p_r[nn]                # risk of reinfection
        td <- RR_r[1,nn]*tds
        LtoL20 <- RR_r[1,nn]*LtoL20s
        Ld=1-sum(td)                                           # proportion entering latent state after infection 
        # Calculate initial conditons for the model
        # Assumes ARI has been constant for last 5 years
        # ******* but doesn't account for change in L and S over previous 5 years ********#
        #S, L1, ..., L20, L
        initial_cond[1] <- 1-prev_L[nn]
        initial_cond[21] <- A*((1-prev_L[nn]) + p*prev_L[nn])*td[20] + LtoL20*prev_L[nn] 
        for (ic in 20:2){
          initial_cond[ic] <- initial_cond[ic+1]+A*((1-prev_L[nn]) + p*prev_L[nn])*td[ic-1]
        }
        initial_cond[22] <- 1-sum(initial_cond[1:21])
        initial_cond <- initial_cond*pop_size

        # these parameters determine transitions on to PT for repeat screening - set to zero initially
        # will be changed over time in model run and depend on sens and spec above 
        covP_S = 0  # in susceptible population
        covP_L = 0  # in latent population
        covP_I = 0  # in progressors
        # generate the initial conditions for the Markov model
        init <- IC_gen(initial_cond,0,0,0,E)
        # create matrices to store Markov model states
        run_base <- mat.or.vec(length(init),t_Markov)
        run_all <- run_base
        run_RISK11 <- run_base
        run_opt <- run_base
        run_min <- run_base
        # numbers screened
        S_base <- rep(0,t_Markov)
        S_all <- S_base
        S_RISK11 <- S_base
        S_opt <- S_base
        S_min <- S_base
        # and vectors for the numbers given PT
        PT_base <- rep(0,t_Markov)
        PT_all <- PT_base
        PT_RISK11 <- PT_base
        PT_opt <- PT_base
        PT_min <- PT_base
     
        # Then iterate through time
        for(ii in 1:t_Markov){

          td <- RR_r[ii,nn]*tds
          LtoL20 <- RR_r[ii,nn]*LtoL20s
          Ld=1-sum(td)  
      
          # set transition matrices to no PT version
          mat_base <- mat_gen(0,0,0,E,ER,A,R,length(init))
          mat_all <- mat_base
          mat_RISK11 <- mat_base
          mat_opt <- mat_base
          mat_min <- mat_base
  
          # if time step where PT is given 
          if (ii %in% t_screen){ 
            #calculate new matrices
            mat_all <- mat_gen(1,1,1,E,ER,A,R,length(init))
            mat_RISK11 <- mat_gen((1-spec_RISK11[ii,nn]),(1-spec_RISK11[ii,nn]),(sens_RISK11[ii,nn]),E,ER,A,R,length(init))
            mat_opt <- mat_gen((1-spec_opt[nn]),(1-spec_opt[nn]),sens_opt[nn],E,ER,A,R,length(init))
            mat_min <- mat_gen((1-spec_min[nn]),(1-spec_min[nn]),sens_min[nn],E,ER,A,R,length(init))
            # Calculate number screened and number given PT - if in first time step use init, otherwise use current state matrix
            if (ii==1){
              S_RISK11[ii] <- init[1] + sum(init[4:19]) + sum(init[20:24])
              S_opt[ii] <- init[1] + sum(init[4:19]) + sum(init[20:24])
              S_min[ii] <- init[1] + sum(init[4:19]) + sum(init[20:24])
              PT_all[ii] <- (init[1] + sum(init[4:19]) + sum(init[20:24]))
              PT_RISK11[ii] <- ((1-spec_RISK11[ii,nn])*(init[1] + sum(init[4:19])) + sens_RISK11[ii,nn]*sum(init[20:24]))
              PT_opt[ii] <- ((1-spec_opt[nn])*(init[1] + sum(init[4:19])) + sens_opt[nn]*sum(init[20:24]))
              PT_min[ii] <- ((1-spec_min[nn])*(init[1] + sum(init[4:19])) + sens_min[nn]*sum(init[20:24]))
            }
      
          if (ii > 1){
            S_RISK11[ii] <- run_RISK11[1,ii-1] + sum(run_RISK11[4:19,ii-1]) + sum(run_RISK11[20:24,ii-1])
            S_opt[ii] <- run_opt[1,ii-1] + sum(run_opt[4:19,ii-1]) + sum(run_opt[20:24,ii-1])
            S_min[ii] <- run_min[1,ii-1] + sum(run_min[4:19,ii-1]) + sum(run_min[20:24,ii-1])
            PT_all[ii] <- (run_all[1,ii-1] + sum(run_all[4:19,ii-1]) + sum(run_all[20:24,ii-1]))
            PT_RISK11[ii] <- ((1-spec_RISK11[ii,nn])*(run_RISK11[1,ii-1] + sum(run_RISK11[4:19,ii-1])) + sens_RISK11[ii,nn]*sum(run_RISK11[20:24,ii-1]))
            PT_opt[ii] <- ((1-spec_opt[nn])*(run_opt[1,ii-1] + sum(run_opt[4:19,ii-1])) + sens_opt[nn]*sum(run_opt[20:24,ii-1]))
            PT_min[ii] <- ((1-spec_min[nn])*(run_min[1,ii-1] + sum(run_min[4:19,ii-1])) + sens_min[nn]*sum(run_min[20:24,ii-1]))
          }
          
        }

        if (ii==1){
          run_base[,ii] <- init%*%mat_base
          run_all[,ii] <- init%*%mat_all
          run_RISK11[,ii] <- init%*%mat_RISK11
          run_opt[,ii] <- init%*%mat_opt
          run_min[,ii] <- init%*%mat_min
        }
        if (ii > 1){
          run_base[,ii] <- run_base[,ii-1]%*%mat_base
          run_all[,ii] <- run_all[,ii-1]%*%mat_all
          run_RISK11[,ii] <- run_RISK11[,ii-1]%*%mat_RISK11
          run_opt[,ii] <- run_opt[,ii-1]%*%mat_opt
          run_min[,ii] <- run_min[,ii-1]%*%mat_min
        }
      }

      # Combine and rearrange results
      run_base <- as.data.frame(cbind(t(run_base),S_base,PT_base,seq(1,t_Markov),"No PT"))
      colnames(run_base) <- c(names(init),"n_screen","n_PT","time","strategy")
      run_all <- as.data.frame(cbind(t(run_all),S_all,PT_all,seq(1,t_Markov),"All"))
      colnames(run_all) <- c(names(init),"n_screen","n_PT","time","strategy")
      run_RISK11 <- as.data.frame(cbind(t(run_RISK11),S_RISK11,PT_RISK11,seq(1,t_Markov),"RISK11+"))
      colnames(run_RISK11) <- c(names(init),"n_screen","n_PT","time","strategy")
      run_opt <- as.data.frame(cbind(t(run_opt),S_opt,PT_opt,seq(1,t_Markov),"opt+"))
      colnames(run_opt) <- c(names(init),"n_screen","n_PT","time","strategy")
      run_min <- as.data.frame(cbind(t(run_min),S_min,PT_min,seq(1,t_Markov),"min+"))
      colnames(run_min) <- c(names(init),"n_screen","n_PT","time","strategy")
    
      outputs <- rbind(run_base,run_all,run_RISK11,run_opt,run_min)
      outputs_m <- melt(outputs,id.vars=c("time","strategy"))

      # Calculate cases averted etc
      for (ii in 1:length(t_store)){
        
        tt <- t_store[ii]
        
        Cases[kk,ii,,nn] <- as.numeric(as.character(outputs_m[outputs_m$time==tt&outputs_m$variable=="had_TB",4]))
        Cases_avert[kk,ii,,nn] <- Cases[kk,ii,1,nn]-Cases[kk,ii,,nn] # compared to baseline
        Diff_Cases_avert[kk,ii,,nn] <- Cases[kk,ii,2,nn]-Cases[kk,ii,,nn] # absolute difference compared to treating all
        Ratio_Cases_avert[kk,ii,,nn] <- Cases_avert[kk,ii,,nn]/Cases_avert[kk,ii,2,nn] # ratio compared to treating all  
        HR[kk,ii,,nn] <- Cases[kk,ii,,nn]/Cases[kk,ii,1,nn]
        Cum_screen[kk,ii,,nn] <- c(sum(S_base[1:tt]),sum(S_all[1:tt]),
                                 sum(S_RISK11[1:tt]),sum(S_opt[1:tt]),sum(S_min[1:tt]))
        Cum_PT[kk,ii,,nn] <- c(sum(PT_base[1:tt]),sum(PT_all[1:tt]),
                             sum(PT_RISK11[1:tt]),sum(PT_opt[1:tt]),sum(PT_min[1:tt]))
        NNT[kk,ii,,nn] <- Cum_PT[kk,ii,,nn]/Cases_avert[kk,ii,,nn]
        g_per[kk,ii,,nn] <- (1/Cum_screen[kk,ii,,nn])*((Cum_PT[kk,ii,2,nn]*Cases_avert[kk,ii,,nn]/Cases_avert[kk,ii,2,nn])-Cum_PT[kk,ii,,nn])
        g_cost[kk,ii,,nn] <- (1/Cum_screen[kk,ii,,nn])*(Cum_PT[kk,ii,2,nn]-Cum_PT[kk,ii,,nn])
      }
        
      kk <- kk+1
      }
    }
  } 
}


# Combine results for plotting
Freq <- rep(screen,times=n_resamp*length(t_store)*n_strat)
Efficacy <- rep(E_RR,times=n_resamp*length(t_store)*n_strat)
timestep <- rep(rep(t_store,each=N_runs),times=n_resamp*n_strat)
strategy <- rep(rep(c("None","All","RISK11","Optimum","Minimum"),each=N_runs*length(t_store)),times=n_resamp)  
samp <- rep(seq(1,n_resamp),each=n_strat*length(t_store)*N_runs)
  
dim(Cases_avert) <- n_resamp*N_runs*length(t_store)*n_strat
dim(Diff_Cases_avert) <- n_resamp*N_runs*length(t_store)*n_strat
dim(Ratio_Cases_avert) <- n_resamp*N_runs*length(t_store)*n_strat
dim(Cases) <- n_resamp*N_runs*length(t_store)*n_strat
dim(HR) <- n_resamp*N_runs*length(t_store)*n_strat
dim(Cum_screen) <- n_resamp*N_runs*length(t_store)*n_strat
dim(Cum_PT) <- n_resamp*N_runs*length(t_store)*n_strat
dim(NNT) <- n_resamp*N_runs*length(t_store)*n_strat
dim(g_cost) <- n_resamp*N_runs*length(t_store)*n_strat
dim(g_per) <- n_resamp*N_runs*length(t_store)*n_strat

outputs <- as.data.frame(cbind(Cases,Cases_avert,Diff_Cases_avert,Ratio_Cases_avert,Cum_PT,Cum_screen,NNT,g_cost,g_per,HR,
                               Freq,Efficacy,timestep,strategy,samp))
outputs <- melt(outputs,id.vars=c("Freq","Efficacy","timestep","strategy","samp"))

# check incidence of TB per 100 ps yrs over 12 months - compare to CORTIS-HR
temp <- outputs[outputs$variable=="Cases"&
                outputs$strategy=="None"&
                outputs$efficacy==1&
                outputs$Freq==1&
                outputs$timestep==5,"value"]
Incidence <- quantile(100*as.numeric(as.character(temp))/10000,probs=c(0.5,0.025,0.975),na.rm=TRUE)

# Get median and 95% CI of results by freq, efficacy, strategy and time
med_outputs <- as.data.frame(as.table(with(outputs,tapply(as.numeric(as.character(value)), list("Freq"=Freq, 
                                                                                                "Efficacy"=Efficacy,
                                                                                                "Strategy"=strategy,
                                                                                                "Time"=timestep,
                                                                                                "Variable"=variable),
                                         function(x) quantile(x, probs=c(0.5),na.rm = TRUE)))))
lo_outputs <- as.data.frame(as.table(with(outputs,tapply(as.numeric(as.character(value)), list("Freq"=Freq, 
                                                                                                "Efficacy"=Efficacy,
                                                                                                "Strategy"=strategy,
                                                                                                "Time"=timestep,
                                                                                                "Variable"=variable),
                                                          function(x) quantile(x, probs=c(0.025),na.rm = TRUE)))))
hi_outputs <- as.data.frame(as.table(with(outputs,tapply(as.numeric(as.character(value)), list("Freq"=Freq, 
                                                                                               "Efficacy"=Efficacy,
                                                                                               "Strategy"=strategy,
                                                                                               "Time"=timestep,
                                                                                               "Variable"=variable),
                                                         function(x) quantile(x, probs=c(0.975),na.rm = TRUE)))))

outputs_range <- cbind(med_outputs,lo_outputs[,6],hi_outputs[,6])
colnames(outputs_range) <- c(colnames(outputs_range)[1:5],"med","min","max")

# PLOTS ##########################################################################################################

# Plot main results - 100% relative efficacy
#                     5y and 10 y horizon
#                     once, 3yr, 3yr, 1yr frequency
#                     whole cohort
#                     RISK11, min spec, opt spec
# subset data
temp  <- outputs_range[outputs_range$Variable%in%c("Ratio_Cases_avert","Cum_PT","Cum_screen","g_cost","g_per")&
                      outputs_range$Efficacy==1&
                      outputs_range$Strategy%in%c("RISK11","Optimum","Minimum")&
                      outputs_range$Time%in%c(20,40),]

# rename stuff
temp$Freq_f = factor(temp$Freq, levels=c('1','2','3','4'))
levels(temp$Freq_f) <- c("Once","3yrs","2yrs","1yr")

temp$Variable_f = factor(temp$Variable, levels=c('Cases','Cases_avert','Diff_Cases_avert','Ratio_Cases_avert','Cum_screen','Cum_PT','NNT','g_cost','g_per','HR'))
levels(temp$Variable_f) <- c('Cases',
                         'Cases averted',
                         'Additional Cases\naverted',
                         'Ratio Cases\naverted',
                         'Cumulative\nscreened',
                         'Cumulative\ntreated',
                         'NNT',
                         'Max test cost\n(total cost)',
                         'Max test cost\n(/case averted)',
                         'Hazard ratio')

temp$Time_f = factor(temp$Time, levels=c(5,20,40))
levels(temp$Time_f) <- c("1yr","5yrs","10yrs")

# Plot 
dummy <- data.frame(Variable_f="Ratio Cases\naverted",z=1)  # add line at 1 in Ratio cases panel 
plot_results_freq <- ggplot(temp,aes(x=Freq_f,
                                     y=as.numeric(as.character(med)),
                                     ymin=as.numeric(as.character(min)),
                                     ymax=as.numeric(as.character(max)),
                                     color=Strategy))+
  geom_errorbar(width=0.2,position=position_dodge(width=0.5),size=1)+
  geom_point(position=position_dodge(width=0.5),size=1)+
  geom_hline(data = dummy, aes(yintercept = z),linetype="dashed")+
  theme_bw()+
  ylab("")+
  xlab("Screening interval")+
  expand_limits(y=0)+
  scale_colour_manual(values=cbbPalette)+
  facet_grid(Variable_f~Time_f,scales="free",switch="y")+
  labs(colour = "Test")+ 
  theme(strip.background =element_rect(fill="white"))+
  theme(text = element_text(size=14))

# Plot RISK11 results over 5 and 10yr horizon, annual screening as function of relative efficacy ##################
# subset data
temp  <- outputs_range[outputs_range$Variable%in%c("Ratio_Cases_avert","g_per")&
                         outputs_range$Freq%in%c(4)&
                         outputs_range$Strategy%in%c("RISK11")&
                         outputs_range$Time%in%c(20,40),]

temp$Time_f = factor(temp$Time, levels=c(5,20,40))
levels(temp$Time_f) <- c("1yr","5yrs","10yrs")

temp$Variable_f = factor(temp$Variable, levels=c('Cases','Cases_avert','Diff_Cases_avert','Ratio_Cases_avert','Cum_screen','Cum_PT','NNT','g_cost','g_per','HR'))
levels(temp$Variable_f) <- c('Cases',
                             'Cases averted',
                             'Additional Cases averted',
                             'Ratio Cases averted',
                             'Cumulative screened',
                             'Cumulative treated',
                             'NNT',
                             'Max test cost (total cost)',
                             'Max test cost (/case averted)',
                             'Hazard ratio')

dummy <- data.frame(Variable_f=c("Ratio Cases averted","Max test cost (/case averted)"),z=c(1,0))  # add line at 1 in Ratio cases panel 

plot_results_eff <- ggplot(temp,aes(x=as.numeric(as.character(Efficacy)),
                                     y=as.numeric(as.character(med)),
                                     ymin=as.numeric(as.character(min)),
                                     ymax=as.numeric(as.character(max)),fill=Time_f))+
geom_ribbon(alpha=0.5)+
  geom_line(aes(y=as.numeric(as.character(med))))+
  geom_hline(data = dummy, aes(yintercept = z),linetype="dashed")+
  theme_bw()+
  ylab("")+
  xlab("Relative efficacy of PT in RISK11+")+
  expand_limits(y=0)+
  scale_fill_manual(values=cbbPalette)+
  facet_wrap(~Variable_f,scales="free_y")+
  theme(strip.background =element_rect(fill="white"))+
  labs(fill = "Time horizon")+
  theme(legend.position="top")

# check for which efficacy ratio cases >1 and cost >0
thresholds <- c(min(as.numeric(as.character(temp[which(temp[temp$Variable=="Ratio_Cases_avert"&temp$Time==20,"min"]>1),"Efficacy"]))),
           min(as.numeric(as.character(temp[which(temp[temp$Variable=="Ratio_Cases_avert"&temp$Time==40,"min"]>1),"Efficacy"]))),
           min(as.numeric(as.character(temp[which(temp[temp$Variable=="g_per"&temp$Time==20,"min"]>0),"Efficacy"]))),
           min(as.numeric(as.character(temp[which(temp[temp$Variable=="g_per"&temp$Time==40,"min"]>0),"Efficacy"]))))

## Now calculate PRCCs on outputs for RISK11 strategy
## Look at 10yr horizon, annual screening, 100% relative efficacy
temp1 <- outputs[outputs$strategy=="RISK11"&
                 outputs$Efficacy==1&
                 outputs$variable%in%c("Ratio_Cases_avert","g_per")&
                 outputs$timestep%in%c(20,40),]
PRCCs <- c()
for(ff in c(1,2,3,4)){
  for(vv in c(c("Ratio_Cases_avert","g_per"))){
    for (tt in c(20,40)){ 
      
      temp2 <- temp1[temp1$Freq==ff&temp1$variable==vv&temp1$timestep==tt,] 
      temp3 <- cbind(spec_L_IGRA_r,sens_L_IGRA_r,p_r,sens_RISK11[1,],spec_RISK11[1,],RR_r[1,],temp2$value)
      temp3 <- as.data.frame(temp3)
      PRCCs <- rbind(PRCCs,cbind(epi.prcc(temp3, sided.test = 2),
                                 c("Specificity IGRA","Sensitivity IGRA","p","Sensitivity RISK11","Specificity RISK11","RR_HIV"),
                                 ff,vv,tt))
      
    }
  }
}
colnames(PRCCs) <- c("gamma","teststat","df","p","parameter","Frequency","variable","Time")
PRCCs$Time_f = factor(PRCCs$Time, levels=c(20,40))
levels(PRCCs$Time_f) <- c("5yrs","10yrs")
PRCCs$Freq_f = factor(PRCCs$Freq, levels=c('1','2','3','4'))
levels(PRCCs$Freq_f) <- c("Once","3yrs","2yrs","1yr")
PRCCs$Variable_f = factor(PRCCs$variable, levels=c('Ratio_Cases_avert','g_per'))
levels(PRCCs$Variable_f) <- c('Ratio Cases averted','Max test cost (/case averted)')

# Plot PRCCs
plot_PRCCs <- ggplot(PRCCs,aes(parameter,as.numeric(as.character(gamma)),fill=Variable_f))+
  geom_bar(stat="identity",position=position_dodge(),color="black")+
  facet_grid(Time_f~Freq_f)+
  coord_cartesian(ylim = c(-1, 1))+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_fill_manual(values=cbbPalette)+
  theme_bw()+
  theme(legend.position="right")+
  theme(text = element_text(size=20))+
  ylab("PRCC")+xlab("Parameter")+
  labs(fill = "")+
  theme(legend.position="top")
