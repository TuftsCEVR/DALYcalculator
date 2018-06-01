#         DALY calculator, function definitions         #
#                     Joanna Emerson                    #
# Center for the Evaluation of Value and Risk in Health #
#                  Tufts Medical Center                 #
#                      May 1, 2018                      #





#                                 CALCULATE DALYS FOR AN INDIVIDUAL                            #
#                                                                                              #
# Calculates disability-adjusted life years for an individual                                  # 
# Returns a vector of Years of Life Life (discounted), Years Lived in Disease, and total DALYs #
#                                                                                              #
# PARAMETERS:                                                                                  #
#   K = Age weighting modulation factor (1=use age weighting, 0=no age weighting)              #
#   C = contant (default = 0.16243)                                                            #
#   r = discount rate (between 0-1)                                                            #
#   beta = parameter of age weighting function (between 0-1)                                   #
#   a_death = Age of premature death due to disease (in years)                                 #
#   a_disability = Age of disease onset (in years)                                             #
#   YLL_L = Life expectancy at age of death (in years)                                         #
#   D = Disability weight (between 0-1)                                                        #

f_DALY<-function(K, C = 0.16243, r, beta, a_death, a_disability, YLL_L, D){

  #calculate YLD_L
  YLD_L <- a_death - a_disability
  
  #only r = 0
  if(r==0 & K!=0) {
    YLL<<-((K*C*exp(-beta*a_death))/(beta^2))*(exp(-beta*YLL_L)*((-beta)*(YLL_L+a_death)-1)-(-beta*a_death-1))+((1-K)*YLL_L)
    YLL_discounted<-YLL
    
    YLD<<-D*(((K*C*exp(-beta*a_disability))/(beta^2))*(exp(-beta*YLD_L)*(-beta*(YLD_L+a_disability)-1)-(-beta*a_disability-1))+((1-K)*YLD_L))
    
    DALY_total<<-YLL_discounted+YLD
  } 
  
  #both r= 0 and k = 0
  else if(r==0 & K==0) {
    YLL<<- ((1-K)/0.00000001)*(1-exp(-0.00000001*YLL_L))
    
    YLD<<- D*(((1-K)/0.00000001)*(1-exp(-0.00000001*YLD_L)))
    
    s<-a_death-a_disability
    YLL_discounted<-YLL*exp(-(0.00000001*s))
    
    DALY_total<<-YLL_discounted+YLD
  }
  
  #only k = 0
  else if(r!=0 & K==0) {
    YLL<<-((1-K)/r)*(1-exp(-r*YLL_L))
    
    YLD<<- D*(((1-K)/r)*(1-exp(-r*YLD_L)))
    
    s<-a_death-a_disability
    YLL_discounted<-YLL*exp(-(r*s))
    
    DALY_total<<-YLL_discounted+YLD
  }
  
  #neither r = 0 nor k = 0
  else if (r!=0 & K!=0){
    YLL<<-((K*C*exp(r*a_death))/((r+beta)^2))*(exp(-(r+beta)*(YLL_L+a_death))*(-(r+beta)*(YLL_L+a_death)-1)-exp(-(r+beta)*a_death)*(-(r+beta)*a_death-1))+((1-K)/r)*(1-exp(-r*YLL_L))
    
    YLD<<-D*((K*C*exp(r*a_disability))/((r+beta)^2))*(exp(-(r+beta)*(YLD_L+a_disability))*(-(r+beta)*(YLD_L+a_disability)-1)-exp(-(r+beta)*a_disability)*(-(r+beta)*a_disability-1))+((1-K)/r)*(1-exp(-r*YLD_L))
    
    s<-a_death-a_disability
    YLL_discounted<-YLL*exp(-(r*s))
    
    DALY_total<<-YLL_discounted+YLD
  }
  #return vector of years of life lost, years lived in disease, and total DALYs
  Amount<-c(YLL_discounted, YLD, DALY_total)
  return(Amount)
}

#                                   CALCULATE DALYS FOR A POPULATION                           #
#                                                                                              #
# Calculates disability-adjusted life years for a population                                   # 
# Returns a vector of Years of Life Life (discounted), Years Lived in Disease, and total DALYs #
#                                                                                              #
# PARAMETERS:                                                                                  #
#   K = Age weighting modulation factor (1=use age weighting, 0=no age weighting)              #
#   C = contant (default = 0.16243)                                                            #
#   r = discount rate (between 0-1)                                                            #
#   beta = parameter of age weighting function (between 0-1)                                   #
#   a_death = Average age of premature death due to disease for population (in years)          #
#   a_disability = Average age of disease onset for population (in years)                      #
#   YLL_L = Life expectancy at age of death (in years)                                         #
#   D = Disability weight (between 0-1)                                                        #
#   incident_deaths = average number of deaths for the population                              #
#   incident_cases = average number of cases in the population                                 #

f_DALYpop<-function(K, C = 0.16243, r, beta, a_death, a_disability, YLL_L, D, incident_deaths, incident_cases){
  #calculate YLD_L
  YLD_L <- a_death - a_disability

  #only r = 0
  if(r==0 & K!=0) {
    YLL<-((K*C*exp(-beta*a_death))/(beta^2))*(exp(-beta*YLL_L)*((-beta)*(YLL_L+a_death)-1)-(-beta*a_death-1))+((1-K)*YLL_L)
    popYLL_discounted<<-YLL*incident_deaths
    
    YLD<-D*(((K*C*exp(-beta*a_disability))/(beta^2))*(exp(-beta*YLD_L)*(-beta*(YLD_L+a_disability)-1)-(-beta*a_disability-1))+((1-K)*YLD_L))
    popYLD<<-YLD*incident_cases
    
    popDALY_total<<-popYLL_discounted+popYLD
  } 
  
  #both r = 0 and k = 0
  else if(r==0 & K==0) {
    YLL<- ((1-K)/0.00000001)*(1-exp(-0.00000001*YLL_L))
    s<-a_death-a_disability
    YLL_discounted<-YLL*exp(-(0.00000001*s))
    popYLL_discounted<<-YLL_discounted*incident_deaths
    
    YLD<<- D*(((1-K)/0.00000001)*(1-exp(-0.00000001*YLD_L)))
    popYLD<<-YLD*incident_cases
    
    popDALY_total<<-popYLL_discounted+popYLD
  }
  
  #only k = 0
  else if(r!=0 & K==0) {
    YLL<<-((1-K)/r)*(1-exp(-r*YLL_L))
    s<-a_death-a_disability
    YLL_discounted<-YLL*exp(-(r*s))
    popYLL_discounted<<-YLL_discounted*incident_deaths
    
    YLD<<- D*(((1-K)/r)*(1-exp(-r*YLD_L)))
    popYLD<<-YLD*incident_cases
    
    popDALY_total<<-popYLL_discounted+popYLD
  }
  
  #neither r = 0 nor k = 0
  else if (r!=0 & K!=0){
    YLL<-((K*C*exp(r*a_death))/((r+beta)^2))*(exp(-(r+beta)*(YLL_L+a_death))*(-(r+beta)*(YLL_L+a_death)-1)-exp(-(r+beta)*a_death)*(-(r+beta)*a_death-1))+((1-K)/r)*(1-exp(-r*YLL_L))
    s<-a_death-a_disability
    YLL_discounted<-YLL*exp(-(r*s))
    popYLL_discounted<<-YLL_discounted*incident_deaths
    
    YLD<-D*((K*C*exp(r*a_disability))/((r+beta)^2))*(exp(-(r+beta)*(YLD_L+a_disability))*(-(r+beta)*(YLD_L+a_disability)-1)-exp(-(r+beta)*a_disability)*(-(r+beta)*a_disability-1))+((1-K)/r)*(1-exp(-r*YLD_L))
    popYLD<<-YLD*incident_cases
    
    popDALY_total<<-popYLL_discounted+popYLD
  }
  
  #return vector of years of life lost, years lived in disease, and total DALYs
  Amount<-c(popYLL_discounted, popYLD, popDALY_total)
  return(Amount)
}