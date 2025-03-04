library(deSolve)


#Parametros
alpha1 <- 5.8*10**-8 #1 o 3.14, ELISA D H
alpha2 <- 2.5*10**-9 #2
#mu1 <- 0.06 # REFERENCIA DE YANG 2015 (HAY OTRO VALOR)
mu1 <- 0.06 #3
nu1 <- 1 # REFERENCIA DE FREITAS, LA OTRA OPCION ES MUY DRASTICA  #4 
mu2 <- 5*10**-1 #5
mu3 <- 0 #MISMO VALOR QUE EL DE LOS MACROFAGOS, MU2 #6
mu4 <- 6.5*10**-4 #ES LA SUMA DE PARAMETRO DE MUERTE DEL PARASITO Y SU ELIMINACION POR CD8 #7
alpha3 <- 0 #replicacion Mi #8 PUEDE VARIAR SI ES M1 O M2
alpha4 <- 90 #Replicacion Ci #9
mu5 <- 5*10**-1 #10 #NO AFECTAN, PUEDE SER POR LAS CITOCINAS
mu6 <- 1*10**-6 #11
alpha5 <- 14.4 #12
alpha7 <- 0.45 #SECRECION DE IL10 #13
mu7 <- 24 #TASA DE degradacion DE TNF #en base a la literatura CADA 18 MIN #14
mu9 <- 19.2 #15 Degradacion de IL10
n1 <- 17.4 #TNF-IL10
n2 <- 560 #IL10-IL6
h1 <- 3 #TNF-IL10
h2 <- 3.68 #IL10-IL6
M0 <-209# NUMERO EQUIS
qTNF <- 0.14 #16
qIL10 <- 0.15 #17  
n3 <- 100 #MTNF
h3 <- 3.16 #MTNF
n4 <- 4.35 #MIL10
h4 <- 0.3#MIL10


    
#    return(list(c(d_TL,d_M,d_Cn,d_Mi,d_Ci,d_TNF,d_IL10)))


trypanosoma_29_NOV <- function(t,state,parameter){
  with(as.list(c(state,parameter)),  {
    d_TL <- (alpha1*((TNF**h3)/((n3**h3)+TNF**h3))*(n4**h4)/((n4**h4)+IL10**h4))*TL*M-alpha2*TL*Cn+(mu1*((TNF**h3)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*TL+(alpha3*((TNF**h3)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*Mi+alpha4*Ci
    d_M <- (nu1*((TNF**h1)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*(M0-M)+(alpha1*((TNF**h3)/((n3**h3)+TNF**h3))*(n4**h4)/((n4**h4)+IL10**h4))*TL*M+(mu2*((TNF**h3)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*M
    d_Cn <- -alpha2*TL*Cn+(mu3*((TNF**h3)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*Cn
    d_Mi <- (alpha1*((TNF**h3)/((n3**h3)+TNF**h3))*(n4**h4)/((n4**h4)+IL10**h4))*TL*M+(mu5*((TNF**h3)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*Mi
    d_Ci <- alpha2*TL*Cn+(mu6*((TNF**h3)/((n3**h3)+TNF**h3))*((n4**h4)/((n4**h4)+IL10**h4)))*Ci
    d_TNF <- (alpha5*((n1**h1)/((n1**h1)+IL10**h1)))*Mi-mu7*(TNF-qTNF)+(alpha5*((n1**h1)/((n1**h1)+IL10**h1)))*Ci
    d_IL10 <- alpha7*Mi - mu9*(IL10-qIL10)

    return(list(c(d_TL,d_M,d_Cn,d_Mi,d_Ci,d_TNF,d_IL10)))      
  })
}


parameter <- c(alpha1,alpha2, mu1, nu1, mu2, mu3, mu4, alpha3, alpha4, mu5, mu6, alpha5, alpha7, mu7, mu9, n1, n2, h1, h2, qIL10, qTNF,M0) 
state <- c(TL=50, M=209, Cn= 136, Mi= 0, Ci = 0, TNF = 0.14, IL10 =0.15)
times <- seq(0,3650, by=1)
cuatro_de_marzo_2025 <- ode(y= state, times = times, func= trypanosoma_29_NOV, parms = parameter)

plot(cuatro_de_marzo_2025,col='red')



   #############Guardar las imagenes ###########

png("images/Param_alpha5_200_25FEB.png")
plot(Param_alpha5_200,col='red')
dev.off()
   ##################################



matplot(out29_NOV[,1], out29_NOV[,2:7], type = "l", xlab = "Dias (10 years)", ylab="valores",
        main="Enfermedad de Chagas", lwd = 2)
legend("topright", c("d_TL", "D_M", "d_Cn", "d_Mi", "d_Ci","d_TNF","d_IL10" ), col = 1:7, lty=1:7, cex = 0.5)
##Se ve muy grande los parasitos, lo demas por ser valores peque no se ven o se ve plano literal


                                 #### ANALISIS DE SENSIBILIDAD PARAMETRICA ####
Param_Normal #Modelado con parametros normales 
saveRDS(Param_Normal, "objetosRDS/Parametrosnormales.rds")
saveRDS(Param_alpha1_3p14, "objetosRDS/Paramtro_alpha1_3p14.rds")
saveRDS(Param_mu1_1x10m1, "objetosRDS/parametro_mu1_1x10m1.rds")
saveRDS(Param_nu1_1080728, "objetosRDS/parametro_nu_1080728.rds")
saveRDS(Param_mu2_0p0019, "objetosRDS/parametro_mu2_0p0019.rds")
saveRDS(Param_alpha3_1p95, "objetosRDS/parametro_alpha3_1p95.rds")
saveRDS(Param_alpha3_20, "objetosRDS/parametro_alpha3_20.rds")
saveRDS(Param_alpha5_200,"objetosRDS/parametro_alpha5_200.rds" )
saveRDS(Param_alpha7_1p29, "objetosRDS/parametro_alpha7_1p29.rds")




out29_NOV[3650,]
#Cambiar los parametros y hacer las diferencias 
diferencias <- (out29_NOV[3650,])-(out29_NOV[2000,]) #esta son las diferencias


out29_NOV[3650,]

plot(x = out29_NOV[3650,1], y = out29_NOV[3650,2])

#HACER UNA BASE DE DATOS DE LOS VALORES DE LAS VARIABLES CON LOS CMBIOS DE LOS PARAMETROS 

plot(out29_NOV[,7])
  
  
  

  
  
  
