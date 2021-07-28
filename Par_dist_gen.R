# Calculates parameters for distributions of model input parameters

library(rriskDistributions)

# Sens and spec of IGRA for infection
spec_IGRA_LTBI_data <- c(0.94,0.96,0.98) # Pai et al, Ann Intern Med, 2008 (Adults, BCG vaccinated, all low incidence countries) 
sens_IGRA_LTBI_data <- c(0.41,0.61,0.75) # Cattamanchi et al, JAIDS, 2011 (HIV+ populations)

fit_spec_IGRA_LTBI <- get.beta.par(p=c(0.025,0.5,0.975),q=spec_IGRA_LTBI_data,show.output=FALSE,plot=FALSE)
fit_sens_IGRA_LTBI <- get.beta.par(p=c(0.025,0.5,0.975),q=sens_IGRA_LTBI_data,show.output=FALSE,plot=FALSE)

# Sens and spec of RISK11 for incident TB (these are main primary results from CORTIS-HR trial)
spec_RISK11_TB_data <- c(0.653,0.689,0.723) # Mendelsohn et al, Lancet GH, 2021, HIV+ adults, 2 sample positive incident TB over 15 mnths
sens_RISK11_TB_data <- c(0.435,0.886,0.987) # Mendelsohn et al, Lancet GH, 2021, HIV+ adults, 2 sample positive incident TB over 15 mnths
fit_spec_RISK11_TB <- get.beta.par(p=c(0.025,0.5,0.975),q=spec_RISK11_TB_data,show.output=FALSE,plot=FALSE)
fit_sens_RISK11_TB <- get.beta.par(p=c(0.025,0.5,0.975),q=sens_RISK11_TB_data,show.output=FALSE,plot=FALSE)

# Among those who got no IPT - sensitivity was 100% as no incident cases among RISK11- who got no IPT
spec_no_IPT <- c(0.586,0.641,0.692)
spec_ART_no_IPT <- c(0.640,0.702,0.757)
spec_no_ART_no_IPT <- c(0.349,0.459,0.574)
fit_spec_no_IPT <- get.beta.par(p=c(0.025,0.5,0.975),q=spec_no_IPT,show.output=FALSE,plot=FALSE)
fit_spec_ART_no_IPT <- get.beta.par(p=c(0.025,0.5,0.975),q=spec_ART_no_IPT,show.output=FALSE,plot=FALSE)
fit_spec_no_ART_no_IPT <- get.beta.par(p=c(0.025,0.5,0.975),q=spec_no_ART_no_IPT,show.output=FALSE,plot=FALSE)

fit_spec <- rbind(fit_spec_RISK11_TB,fit_spec_no_IPT,fit_spec_ART_no_IPT,fit_spec_no_ART_no_IPT)
fit_sens <- rbind(fit_sens_RISK11_TB,c(1,0),c(1,0),c(1,0))

# Protective effect of prior infection
p_data <- c(0.14,0.21,0.30) # Andrews et al
fit_p <- get.beta.par(p=c(0.025,0.5,0.975),q=p_data,show.output=FALSE,plot=FALSE)

# sens and spec of IGRA for incident TB from CORTIS-HR
spec_IGRA_TB_data <- c(0.525,0.562,0.599) # Mendelsohn et al, Lancet GH, 2021, HIV+ adults, 2 sample positive incident TB over 15 mnths
sens_IGRA_TB_data <- c(0.259,0.621,0.885) # Mendelsohn et al, Lancet GH, 2021, HIV+ adults, 2 sample positive incident TB over 15 mnths
fit_spec_IGRA_TB <- get.beta.par(p=c(0.025,0.5,0.975),q=spec_IGRA_TB_data,show.output=FALSE,plot=FALSE)
fit_sens_IGRA_TB <- get.beta.par(p=c(0.025,0.5,0.975),q=sens_IGRA_TB_data,show.output=FALSE,plot=FALSE)


