install.packages("earlywarnings")
install.packages("EWS")
install.packages("EWSmethods")


############## 30 de Mayo 2024 #############


#punto critico conveniente
#alterar sistema y no retorno, distinguir
#entre situaciones, condiciones cambio y estado del sistema. cambios pqeue;os modifica poquito

#modificar parametros y ver dinamica del estilo, grandes y ver efecto grandes

#condiciones y sistem del estado cambia aburptamente. TP1 tipping point, con un poco de cambio 
#pasa a otro punto de equilibrio

#se;ales de alerta

#oerturbacion de punto critico con varias e;ales para avisar de un punto critico

#punto de no retorno
#alta y baja resiliencia que tanto nos recuperamos y aplicaciones a nivel de psicologia
#

#DOS CLASES DE ALERTA,

#Resultado de alentamiento oo se;ales asociadas al cambio de estabilidad. de atractor a repulsor

#indicadores para cada una de estas

#alentamiento el tiempo de recuperacion  CRITICAL SLOWING DOWN

#Unos son consecuecnias y uno es que cambia el estado de ellos

#CTIRICAL SLOWING DOWN tasa de recuperacion muy baja. 

#cerca de punto critico tarda mas en recuperarse, en perturbaciones peque;as



#CARACTERISTICA ESTADISTICA
#COrrelacion, que tanto se parecen esas dos variables, correlacion cero, no se parecen en nada

#Autocorrelacion, correlacionar tiempo despues misma variables correlacionada tiempo despues

#si no cambia mucho no hay se;al. 

#correlacion si es -1 es negativo pq es totalmente opuesta


######### Cuando se correlaciona parece bolita, si no es hace como nube




##########################

#Incremento de la varianza, distribucion de datos, con el ancho de la distribucion

#LEJOS DEL PUNTO CRITICO, DISTRIBUCION MUY POCA
#CERCA DEL PUNTO CRITICO LA DISTRIBUCION SE ENSANCHA

#VALORES DE DIVERSIDAD


########################

#SE;ALES DE ALERTA TEMPRANA, si hay perturbacion se retarda mas en recuperarse.

#CAMBIO DEL ESTABILITY LANDSCAPE

#cambio de paisaje y una perturbacion nos lleve a otro

##########################
#FAMILIA INCREASING SKEWNESS


#incrementa simestria de la distribucion

############################
#CAMBIO EN EL LANDSCAPE. FLICKERING




library(earlywarnings)
library(EWS)
library(EWSmethods)

####

data("simTransComms")

View(bacalao)
data("CODrecovery")
View(bacalao1)

View(CODrecovery)


uniEWS(CODrecovery$scenario2, metrics = "ar1")












































