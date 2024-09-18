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

##################### ECUACIONES DE TESIS ###############


#Variables

#Parametros

trypanosoma <- function(t, state, parameter){
  with(as.list(c(state, parameter)), {
  d_TL <- -(alpha1*((TNF^h)/(n1^h1*(TL)*(TNF)+TNF^h1))*((n1^h1*(M)*(IL10))/(n1^h1*(M)*(IL10)+IL10^h1)))*TL*M-alpha2*TL*Cn-(mu1*((TNF^h1)/(n1^h1*(TL)*(TNF)+TNF^h1))*((IFN^h1)/(n1^h1*(TL)*(IFN)+IFN^h1))*((n1^h1*(TL)*(IL10))/(n1^h1*(TL)*(IL10)+IL10^h1)))*TL
  d_M <- (nu2*((TNF^h2)/(n2^h2*(TL)*(TNF)+TNF^h2))((n2^h2*(M)*(IL10))/(n2^h2*(M)*(IL10)+IL10^h2)))*(M-M0)-(alpha1*(((TNF^h2)/(n2^h2*(M)*(TNF)+TNF^h2))((n2^h2*(M)*(IL10))/(n2^h2*(M)*(IL10)+IL10^h2))))*TL*M-(mu2*(((IFN^h2)/(n2^h2*(M)*(IFN)+IFN^h2))((n2^h2*(M)*(IL10))/(n2^h2*(M)*(IL10)+IL10^h2))))*M
  d_Cn <- -alpha2*TL*Cn-(mu3*(((IFN^h3)/(n3^h3*(Cn)*(IFN)+IL10^h3))((n3^h3*(Cn)*(IL10))/(n3^h3*(Cn)*(IL10)+IL10^h3))))*Cn
  d_Ti <- 
  
  
  })
}



#Curvas solucion proxima semana deSolve
#Revisarlas con lupa
#escribir ecuaciones con descriptores












