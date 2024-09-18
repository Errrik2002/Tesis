```{r, results='hide'}
install.packages("deSolve")
library(deSolve)

#Variables
#Parametros
trypanosoma <- function(t, state, parameter){
  with(as.list(c(state, parameter)), {
    
  })
}
install.packages("deSolve")
library(deSolve)


zika <- function(t, state, parameter){
  with(as.list(c(state, parameter)), {
    dS_h <- mu_h-mu_h*S_h -betah*S_h*I_v
    dI_h <- betah*S_h*I_v - mu_h*I_h - gamma*I_h
    dR_h <- gamma*I_h - mu_h*R_h
    dS_v <- mu_v - betav*S_v*I_v - mu_v*S_v
    dE_v <- betav*S_v*I_h - mu_v*E_v - gamma*E_v
    dI_v <- gamma*E_v - mu_v*I_v
    list(c(dS_h, dI_h, dR_h, dS_v, dE_v, dI_v))
  })
}

#####################


#Variables
#Parametros
trypanosoma <- function(t, state, parameter){
  with(as.list(c(state, parameter)), {
  d_TL <- -(alpha1*((TNF^h)/(n^h*(TL)*(TNF)+TNF^h))*((n^h*(M)*(IL10))/(n^h*(M)*(IL10)+IL10^h)))*TL*M-alpha2*TL*Cn-(mu1*((TNF^h)/(n^h*(TL)*(TNF)+TNF^h))*((IFN^h)/(n^h*(TL)*(IFN)+IFN^h))*((n^h*(TL)*(IL10))/(n^h*(TL)*(IL10)+IL10^h)))*TL
  })
}











