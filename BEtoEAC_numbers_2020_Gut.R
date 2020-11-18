
##############################################################################
#####   Code to produce MSCE-EAC model Results for Curtius et al. Gut 2020  ##
##############################################################################

source('/BaseCaseII/legauss.R')
dyn.load('/BaseCaseII/legauss.so')
library('deSolve')

phidot = function(s, phi, parms){
  phidot=numeric(6)
  with(as.list(parms),{
  tau = t - s
  #RR=4
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))#/1.5
  nu = nu0*((1-GERDpr)+RR*GERDpr)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## Armitage-Doll transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}

phidot_GERD = function(s, phi, parms){
  phidot=numeric(6)
  with(as.list(parms),{
  tau = t - s
  #RR=4
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))
  #nu = nu0*((1-GERDpr)+RR*GERDpr)
  nu = nu0*(RR*1)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## A-D transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}


h_EAC <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,-ode$phi10[2]/ode$phi9[2])
   }
  return(besurv.ode)
}  



s_EAC <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi9[2])
   }
  return(besurv.ode)
}  


s_EAC_GERD <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi9[2])
   }
  return(besurv.ode)
}  

h_EAC_GERD <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,-ode$phi10[2]/ode$phi9[2])
   }
  return(besurv.ode)
}  


BEdensity <- function(screen_age,gerd1,gerd2,gerd3,nu0_temp,RR_temp){
  ages=1:screen_age
  nu <- rep(0,length(ages))
  nu_cum <- rep(0,length(ages))
  GERDpr1 <- rep(0,length(ages))
  for ( i in 1:length(ages)){
    GERDpr1[i]=1-exp(-gerd1*min(gerd3,ages[i])-gerd2*max(0,ages[i]-gerd3))
    nu[i] <- nu0_temp*((1-GERDpr1[i])+RR_temp*GERDpr1[i])
    GLquad1<-legauss(0 ,ages[i],20)
    wg <- GLquad1$weights
    xg <- GLquad1$mesh
    int1 <- rep(0,20)
    for ( j in 1:20){
      GERDpr=1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,xg[j]-gerd3))
      nu_temp <-nu0_temp*((1-GERDpr)+RR_temp*GERDpr)
      int1[j] <- nu_temp
    }
    nu_cum[i]<-1-exp(-sum(wg*int1))
  }
  nu_dens <- nu*(1-nu_cum)
  return(list(nu=nu,nu_cum=nu_cum,nu_dens=nu_dens,GERDpr=GERDpr1))
}

BEdensity_GERD <- function(screen_age,gerd1,gerd2,gerd3,nu0_temp,RR_temp){
  ages=1:screen_age
  nu <- rep(0,length(ages))
  nu_cum <- rep(0,length(ages))
  GERDpr1 <- rep(0,length(ages))
  for ( i in 1:length(ages)){
    GERDpr1[i]=1-exp(-gerd1*min(gerd3,ages[i])-gerd2*max(0,ages[i]-gerd3))
    nu[i] <- nu0_temp*(RR_temp*1)
    GLquad1<-legauss(0 ,ages[i],20)
    wg <- GLquad1$weights
    xg <- GLquad1$mesh
    int1 <- rep(0,20)
    for ( j in 1:20){
      GERDpr=1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,xg[j]-gerd3))
      nu_temp <-nu0_temp*(RR_temp*1)
      int1[j] <- nu_temp
    }
    nu_cum[i]<-1-exp(-sum(wg*int1))
  }
  nu_dens <- nu*(1-nu_cum)
  return(list(nu=nu,nu_cum=nu_cum,nu_dens=nu_dens,GERDpr=GERDpr1))
}


## US CENSUS DATA FOR 2010 ####

m_num_2010 = c(10393977,11209085,10933274,9523648,8077500,5852547,4243972,3182388,2294374,1273867,424387,82263,9162)
f_num_2010 = c(10496987,11499506,11364851,10141157,8740424,6582716,5034194,4135407,3448953,2346592,1023979,288981,44202)
totalpop_2010 = sum(m_num_2010,f_num_2010)
print(totalpop_2010)
## 142,648,393


## NUMBER OF BE CASES
b_c1 = 1950
birth_cohort = b_c1
RR=5
age_groups = seq(42.5,102.5,by=5)
source('am_parameter_list_2020.R')

nu_cum_m = BEdensity(age_groups[length(age_groups)],gerd1,gerd2,gerd3,nu0,RR)$nu_cum
source('/Users/curtiu01/Dropbox/BaseCaseII/af_parameter_list.R')

nu_cum_f = BEdensity(age_groups[length(age_groups)],gerd1,gerd2,gerd3,nu0,RR)$nu_cum
## Choose BE prevalence (unconditional on EAC) for 5 year age groups 40+ to 100+
age_group_ind = seq(42,102,by=5)

BECases_2010 = sum(nu_cum_m[age_group_ind]*m_num_2010)+sum(nu_cum_f[age_group_ind]*f_num_2010)

## Total BE cases in 40+ population (prev taken at 42 - 102 by 5 year age groups)
print(BECases_2010)
#  2,229,757

## Percentage of adult population with BE in 2010
print(BECases_2010/sum(m_num_2010,f_num_2010))
#  0.01563114

## total US population from 40 to 100+ in 2010
print(sum(m_num_2010,f_num_2010))
#  142,648,393

#########################################################################################################
###############################     FOR NUMBERS PAPER RESULTS         ###################################
################### SEER Explorer population incidence rates for expected cases in 2010    ##############
#########################################################################################################

## EAC cases expected ages 40-90
expcases_2010 = (0.63*sum(f_num_2010[1:2]) + 5.7*sum(m_num_2010[1:2])+1.06*sum(f_num_2010[3:5]) + 9.57*sum(m_num_2010[3:5])
  +2.92*sum(f_num_2010[6:7]) + 22.92*sum(m_num_2010[6:7])+4.26*sum(f_num_2010[8:10]) + 28.17*sum(m_num_2010[8:10]))/100000
print(expcases_2010)
## 9399.355


####################################################
##  PARAMETER FUNCTION GIVEN SEX AND BIRTH COHORT ##
####################################################
param_list <- function(X, nu, alphaP, alphaM, rho, gP,gM,mu0,mu1,mu2){
    betaM <- alphaM-gM-rho
    mu2eff <- mu2*(1-betaM/alphaM)        ## effective rate to pre-clinical cancer 
    betaeff<- alphaP-gP-mu2eff
    betaP<- alphaP-gP-mu2
    peff <- (1/2)*(-alphaP+betaeff + mu2eff - sqrt((alphaP+betaeff+mu2eff)^2-4*alphaP*betaeff))
    qeff <- (1/2)*(-alphaP+betaeff + mu2eff + sqrt((alphaP+betaeff+mu2eff)^2-4*alphaP*betaeff))
    pP <- (1/2)*(-alphaP+betaP + mu2 - sqrt((alphaP+betaP+mu2)^2-4*alphaP*betaP))
    qP <- (1/2)*(-alphaP+betaP + mu2 + sqrt((alphaP+betaP+mu2)^2-4*alphaP*betaP))
    pM <- (1/2)*(-alphaM+betaM + rho - sqrt((alphaM+betaM+rho)^2-4*alphaM*betaM))
    qM <- (1/2)*(-alphaM+betaM + rho + sqrt((alphaM+betaM+rho)^2-4*alphaM*betaM))
    zeta=-qM/pM 
    tlag <- log((zeta/(1+zeta)))/pM             ## exact T_2
    return(list(X=X,nu=nu,mu0=mu0,mu1=mu1,mu2=mu2, rho=rho,alphaP=alphaP,alphaM=alphaM, betaP=betaP,betaM=betaM, gP=gP,gM=gM, tlag=tlag, mu2eff=mu2eff,betaeff=betaeff,pP=pP, qP=qP,pM=pM,qM=qM,peff=peff,qeff=qeff))
  }


get_params = function(gender,bc){
  if (gender=="M"){
    am_g0<- 0.99061E-01        
    am_g1<- 0.50880E+00       
    am_g2<- 0.53768E-01    
    am_trefc<- 0.19125E+04 
    am_gC0<- 0.75000E+00     

    gerd1<- 0.00061422    
    gerd2<- 0.0070447     
    gerd3<- 26.002   
    nu0 <-amnu0<- 0.36494E-03 

    #ammu0<- ammu1 <-  0.79942E-03
    ## ADJUSTED VALUES FOR BASE CASE 2 need a factor change for BE segment length distribution, mean length of beta(.5,4)*15+1
    crypts_per_avBE<- (1+15*(.5/4.5))*10*75/15*1000
    ammu0<- 0.79942E-03*250000/crypts_per_avBE
    ammu1 <-  0.79942E-03
    ammu2<-    0.45439E-04 

    kstem<- 4           ## stem cells/crypt
    #X<- 250000*kstem       ## total stem cells in a 5 cm BE segment assuming 250,000 crypts/5cm
    X <- crypts_per_avBE*kstem

    amrho<- 1.0000E-09        ## detection rate of clinical cancers
    amalphaP0<- 10        ## premalignant cell division rate
    amalphaM0 <- 150      ## pre-clinical malignant cell division rate
    birth_cohort<- bc

    am_gP<- am_g0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
    am_gM<- am_gC0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
    #am_gM=am_gM*1.25
    amalphaP<- amalphaP0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
    amalphaM<- amalphaM0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
    allmales<- param_list(X, amnu0, amalphaP, amalphaM, amrho, am_gP,am_gM, ammu0, ammu1, ammu2)
    #amparams<- c(allmales$mu0,allmales$mu1,allmales$betaP,allmales$mu2,allmales$alphaP,allmales$betaM,allmales$rho,allmales$alphaM,allmales$pM,allmales$qM,allmales$tlag)
    maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
    return(maleparms)
  }
  if (gender=="F"){
    af_g0<- 0.12323E+00   
    af_g1<- 0.63999E+00   
    af_g2<- 0.29781E-01    
    af_trefc<- 0.19453E+04   
    af_gC0<- 0.75000E+00     

    gerd1<- 0.18012E-03     
    gerd2<- 0.71326E-02     
    gerd3<- 0.25808E+02  
    nu0<-afnu0<- 0.74828E-04 

    ## ADJUSTED VALUES FOR BASE CASE 2 need a factor change for BE segment length distribution, mean length of beta(.5,4)
    crypts_per_avBE<- (1+15*(.5/4.5))*10*75/15*1000
    afmu0<- 0.70504E-03*250000/crypts_per_avBE
    afmu1 <-  0.70504E-03
    afmu2<-   0.68903E-04 

    kstem<- 4           ## stem cells/crypt
    #X<- 250000*kstem       ## total stem cells in a 5 cm BE segment assuming 250,000 crypts/5cm
    X<- crypts_per_avBE*kstem       

    afrho<- 1.0000E-09       ## detection rate of clinical cancers
    afalphaP0<- 10        
    afalphaM0 <- 150      ## pre-clinical malignant cell division rate

    birth_cohort<-bc

    af_gP<- af_g0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
    af_gM<- af_gC0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
    afalphaP<- afalphaP0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
    afalphaM<- afalphaM0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
    allfemales<- param_list(X, afnu0, afalphaP, afalphaM, afrho, af_gP,af_gM, afmu0, afmu1, afmu2)
    #afparams<- c(allfemales$mu0,allfemales$mu1,allfemales$betaP,allfemales$mu2,allfemales$alphaP,allfemales$betaM,allfemales$rho,allfemales$alphaM,allfemales$pM,allfemales$qM,allfemales$tlag)
    femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)
    return(femaleparms)
  }
}


age_groups = seq(42.5,87.5,by=5)
Lambda_m = rep(0,length(age_groups))
source('am_parameter_list_2020.R')
for (i in 1:length(Lambda_m)){
  parms_temp = get_params("M",(2010-age_groups[i]))
  Lambda_m[i] =  h_EAC(parms_temp,age_groups[i])*m_num_2010[i]
  print(h_EAC(parms_temp,age_groups[i])*100000)
}

Lambda_f = rep(0,length(age_groups))
source('af_parameter_list_2020.R')
for (i in 1:length(Lambda_f)){
  parms_temp = get_params("F",(2010-age_groups[i]))
  Lambda_f[i] =  h_EAC(parms_temp,age_groups[i])*f_num_2010[i]
  print(h_EAC(parms_temp,age_groups[i])*100000)
}
totalcases_2010=sum(Lambda_f,Lambda_m)
print(totalcases_2010)
## [1] 9966.073




#########################################################################################################
### MAIN RESULTS FOR corresponding BE prevalence-  VALIDATION WITH CORI FOR MALES/FEMALES             ###
#########################################################################################################

#####################################################################################
#################   FIGURE 3 for differing RR values            #####################
#####################################################################################

w=seq(0,20, .25)
RR=5
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
  strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
  strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50=strat1_w_GERD_noMortality
strat1_w_noMortality50=strat1_w_noMortality

strat1_w_f_noMortality50=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50 = strat1_w_GERD_f_noMortality



w=seq(0,20, .25)
RR=6
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}



for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}



for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
    strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50_6=strat1_w_GERD_noMortality
strat1_w_noMortality50_6=strat1_w_noMortality

strat1_w_f_noMortality50_6=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50_6 = strat1_w_GERD_f_noMortality



w=seq(0,20, .25)
RR=2
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}



for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
    strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50_2=strat1_w_GERD_noMortality
strat1_w_noMortality50_2=strat1_w_noMortality

strat1_w_f_noMortality50_2=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50_2 = strat1_w_GERD_f_noMortality


##  FIGURE 3 A-B
dev.new()
par(mfrow=c(1,2),xaxs="i", yaxs="i", las = 1, mar=c(5,5,2,3), bty='l')

plot(t_s[seq(16,66,1)],strat1_w_f_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="orange",lty=1,pch=18,cex=1,cex.axis=2, cex.lab=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="green3",lty=1,pch=17,cex=1)
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_f_noMortality50_2[seq(16,66),5]*100,strat1_w_f_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "orange", alpha.f = 0.1),border=adjustcolor( "orange", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_noMortality50_2[seq(16,66),5]*100,strat1_w_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "green3", alpha.f = 0.1),border=adjustcolor( "green3", alpha.f = 0.1))

legend('topleft',c('Model: Males', 'CORI Males: general', 'Model: Females','CORI Females: general' ), col=c("green3", "green3", "orange","orange"),cex=1.5, lwd = 3, bty='n',lty=c(1,3,1,3))


plot(t_s[seq(16,66,1)],strat1_w_GERD_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="blue",lty=1,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_GERD_f_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="red",lty=1,pch=16,cex=1)


polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_noMortality50_2[seq(16,66),5]*100,strat1_w_GERD_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "blue", alpha.f = 0.1),border=adjustcolor( "blue", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_f_noMortality50_2[seq(16,66),5]*100,strat1_w_GERD_f_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "red", alpha.f = 0.1),border=adjustcolor( "red", alpha.f = 0.1))

legend('topleft',c('Model: GERD Males', 'CORI Males: screening', 'Model: GERD Females','CORI Females: screening' ), col=c("blue", "blue", "red","red"),cex=1, lwd = 3, bty='n',lty=c(1,3,1,3))



##  CHANGING GERD PREV, include in Supplementary material , Figure S1
w=seq(0,20, .25)
RR=5
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
GERD_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$GERDpr[t_s]
#nu_t = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu[t_s]

source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)
GERD_cumf = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$GERDpr[t_s]
#nu_tf = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]

GERD_cum_low = GERD_cum*1.5
GERD_cumf_low = GERD_cumf*1.5
GERD_cum_high = GERD_cum*.5
GERD_cumf_high = GERD_cumf*.5



## CHANGING GERD PREV INPUT DIRECTLY
BEdensity <- function(screen_age,gerd1,gerd2,gerd3,nu0_temp,RR_temp){
    ages=1:screen_age
    nu <- rep(0,length(ages))
    nu_cum <- rep(0,length(ages))
    GERDpr1 <- rep(0,length(ages))
    for ( i in 1:length(ages)){
        GERDpr1[i]=(1-exp(-gerd1*min(gerd3,ages[i])-gerd2*max(0,ages[i]-gerd3)))*1.5
        nu[i] <- nu0_temp*((1-GERDpr1[i])+RR_temp*GERDpr1[i])
        GLquad1<-legauss(0 ,ages[i],20)
        wg <- GLquad1$weights
        xg <- GLquad1$mesh
        int1 <- rep(0,20)
        for ( j in 1:20){
            GERDpr=(1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,xg[j]-gerd3)))*1.5
            nu_temp <-nu0_temp*((1-GERDpr)+RR_temp*GERDpr)
            int1[j] <- nu_temp
        }
        nu_cum[i]<-1-exp(-sum(wg*int1))
    }
    nu_dens <- nu*(1-nu_cum)
    return(list(nu=nu,nu_cum=nu_cum,nu_dens=nu_dens,GERDpr=GERDpr1))
}


phidot = function(s, phi, parms){
    phidot=numeric(6)
    with(as.list(parms),{
        tau = t - s
        #RR=4
        GERDpr=(1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3)))*1.5
        nu = nu0*((1-GERDpr)+RR*GERDpr)
        phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
        phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
        phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
        phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
        phidot[5]=mu1*phi[5]*(phi[3]-1)
        phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
        phidot[7]=mu0*phi[7]*(phi[5]-1)
        phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
        ## Armitage-Doll transition
        phidot[9]=nu*(phi[7]-phi[9])
        phidot[10]=nu*(phi[8]-phi[10])
        list(c(phidot))
    })
}


w=seq(0,20, .25)
RR=5
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))

strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}

for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
}


strat1_w_noMortality50_high=strat1_w_noMortality
strat1_w_f_noMortality50_high=strat1_w_f_noMortality

BEdensity <- function(screen_age,gerd1,gerd2,gerd3,nu0_temp,RR_temp){
    ages=1:screen_age
    nu <- rep(0,length(ages))
    nu_cum <- rep(0,length(ages))
    GERDpr1 <- rep(0,length(ages))
    for ( i in 1:length(ages)){
        GERDpr1[i]=(1-exp(-gerd1*min(gerd3,ages[i])-gerd2*max(0,ages[i]-gerd3)))*.5
        nu[i] <- nu0_temp*((1-GERDpr1[i])+RR_temp*GERDpr1[i])
        GLquad1<-legauss(0 ,ages[i],20)
        wg <- GLquad1$weights
        xg <- GLquad1$mesh
        int1 <- rep(0,20)
        for ( j in 1:20){
            GERDpr=(1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,xg[j]-gerd3)))*.5
            nu_temp <-nu0_temp*((1-GERDpr)+RR_temp*GERDpr)
            int1[j] <- nu_temp
        }
        nu_cum[i]<-1-exp(-sum(wg*int1))
    }
    nu_dens <- nu*(1-nu_cum)
    return(list(nu=nu,nu_cum=nu_cum,nu_dens=nu_dens,GERDpr=GERDpr1))
}


phidot = function(s, phi, parms){
    phidot=numeric(6)
    with(as.list(parms),{
        tau = t - s
        #RR=4
        GERDpr=(1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3)))*.5
        nu = nu0*((1-GERDpr)+RR*GERDpr)
        phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
        phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
        phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
        phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
        phidot[5]=mu1*phi[5]*(phi[3]-1)
        phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
        phidot[7]=mu0*phi[7]*(phi[5]-1)
        phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
        ## Armitage-Doll transition
        phidot[9]=nu*(phi[7]-phi[9])
        phidot[10]=nu*(phi[8]-phi[10])
        list(c(phidot))
    })
}


w=seq(0,20, .25)
RR=5
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))

strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}

for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
}


strat1_w_noMortality50_low=strat1_w_noMortality
strat1_w_f_noMortality50_low=strat1_w_f_noMortality

## FOR GERD PREV, Figure S1
dev.new()
par(mfrow=c(2,2),xaxs="i", yaxs="i", las = 1, mar=c(5,5,2,3), bty='l')
plot(t_s[seq(11,66,1)],GERD_cum[seq(11,66,1)]*100,type='l',lwd=3,col="blue",lty=1,pch=17,cex=1.5,cex.axis=1.5, cex.lab=1.5, ylim=c(0,50),xlim=c(20,80),ylab= "GERD prevalence(%)", xlab="Age")
lines(t_s[seq(11,66,1)],GERD_cum_high[seq(11,66,1)]*100,type='l',lwd=3,col="goldenrod",lty=1,pch=17,cex=1)
lines(t_s[seq(11,66,1)],GERD_cum_low[seq(11,66,1)]*100,type='l',lwd=3,col="goldenrod",lty=4,pch=17,cex=1)
polygon(c(t_s[seq(11,66,1)],t_s[seq(66,11)]),c(GERD_cum_low[seq(11,66)]*100,GERD_cum_high[seq(66,11)]*100),col= adjustcolor( "goldenrod", alpha.f = 0.1),border=adjustcolor( "goldenrod", alpha.f = 0.1))
legend('topleft',c('GERD Males: 50% increase', 'GERD Males: baseline prev', 'GERD Males: 50% decrease' ), col=c("goldenrod", "blue","goldenrod"),cex=1, lwd = 3, bty='n',lty=c(4,1,1))

plot(t_s[seq(11,66,1)],GERD_cumf[seq(11,66,1)]*100,type='l',lwd=3,col="red",lty=1,pch=18,cex=1.5,cex.axis=1.5, cex.lab=1.5, ylim=c(0,50),xlim=c(20,80),ylab= "GERD prevalence(%)", xlab="Age")
lines(t_s[seq(11,66,1)],GERD_cumf_high[seq(11,66,1)]*100,type='l',lwd=3,col="purple",lty=1,pch=18,cex=1,cex.axis=1.5, cex.lab=1.5, ylim=c(0,50),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(11,66,1)],GERD_cumf_low[seq(11,66,1)]*100,type='l',lwd=3,col="purple",lty=4,pch=18,cex=1,cex.axis=2, cex.lab=2, ylim=c(0,50),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
polygon(c(t_s[seq(11,66,1)],t_s[seq(66,11)]),c(GERD_cumf_low[seq(11,66)]*100,GERD_cumf_high[seq(66,11)]*100),col= adjustcolor( "purple", alpha.f = 0.1),border=adjustcolor( "purple", alpha.f = 0.1))
legend('topleft',c('GERD Females: 50% increase', 'GERD Femles: baseline prev', 'GERD Females: 50% decrease' ), col=c("purple", "red","purple"),cex=1, lwd = 3, bty='n',lty=c(4,1,1))


plot(t_s[seq(16,66,1)],strat1_w_noMortality50_low[seq(16,66,1),5]*100,type='l',lwd=3,col="goldenrod",lty=1,pch=17,cex=1.5,cex.axis=1.5, cex.lab=1.5, ylim=c(0,6),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_noMortality50_high[seq(16,66,1),5]*100,type='l',lwd=3,col="goldenrod",lty=4,pch=17,cex=1.5,cex.axis=1.5, cex.lab=1.5, ylim=c(0,6),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")

polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_noMortality50_low[seq(16,66),5]*100,strat1_w_noMortality50_high[seq(66,16),5]*100),col= adjustcolor( "goldenrod", alpha.f = 0.1),border=adjustcolor( "goldenrod", alpha.f = 0.1))


plot(t_s[seq(16,66,1)],strat1_w_f_noMortality50_low[seq(16,66,1),5]*100,type='l',lwd=3,col="purple",lty=1,pch=17,cex=1.5,cex.axis=1.5, cex.lab=1.5, ylim=c(0,6),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_f_noMortality50_high[seq(16,66,1),5]*100,type='l',lwd=3,col="purple",lty=4,pch=17,cex=1.5,cex.axis=1.5, cex.lab=1.5, ylim=c(0,6),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")

polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_f_noMortality50_low[seq(16,66),5]*100,strat1_w_f_noMortality50_high[seq(66,16),5]*100),col= adjustcolor( "purple", alpha.f = 0.1),border=adjustcolor( "purple", alpha.f = 0.1))




## How many EACs does the model predict if BE prev in GERD was lower: Figure S2
## Let's first try for RR=3

w=seq(0,20, .25)
RR=3
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}



for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
    strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50_3=strat1_w_GERD_noMortality
strat1_w_noMortality50_3=strat1_w_noMortality

strat1_w_f_noMortality50_3=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50_3 = strat1_w_GERD_f_noMortality

## Now RR=4
w=seq(0,20, .25)
RR=4
birth_cohort = 1950
source('am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
    strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50_4=strat1_w_GERD_noMortality
strat1_w_noMortality50_4=strat1_w_noMortality

strat1_w_f_noMortality50_4=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50_4 = strat1_w_GERD_f_noMortality


## NEW FIGURE 3 A-B for second revision, Figure S2
dev.new()
par(mfrow=c(1,2),xaxs="i", yaxs="i", las = 1, mar=c(5,5,2,3), bty='l')

plot(t_s[seq(16,66,1)],strat1_w_f_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="orange",lty=1,pch=18,cex=1,cex.axis=2, cex.lab=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="green3",lty=1,pch=17,cex=1)
lines(t_s[seq(16,66,1)],strat1_w_noMortality50_3[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=3,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_f_noMortality50_3[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=3,pch=16,cex=1)
lines(t_s[seq(16,66,1)],strat1_w_noMortality50_4[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=2,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_f_noMortality50_4[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=2,pch=16,cex=1)


polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_f_noMortality50_2[seq(16,66),5]*100,strat1_w_f_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "orange", alpha.f = 0.1),border=adjustcolor( "orange", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_noMortality50_2[seq(16,66),5]*100,strat1_w_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "green3", alpha.f = 0.1),border=adjustcolor( "green3", alpha.f = 0.1))



plot(t_s[seq(16,66,1)],strat1_w_GERD_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="blue",lty=1,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_GERD_f_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="red",lty=1,pch=16,cex=1)
lines(t_s[seq(16,66,1)],strat1_w_GERD_noMortality50_3[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=3,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_GERD_f_noMortality50_3[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=3,pch=16,cex=1)
lines(t_s[seq(16,66,1)],strat1_w_GERD_noMortality50_4[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=2,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE Prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_GERD_f_noMortality50_4[seq(16,66,1),5]*100,type='l',lwd=3,col="darkgrey",lty=2,pch=16,cex=1)


polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_noMortality50_2[seq(16,66),5]*100,strat1_w_GERD_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "blue", alpha.f = 0.1),border=adjustcolor( "blue", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_f_noMortality50_2[seq(16,66),5]*100,strat1_w_GERD_f_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "red", alpha.f = 0.1),border=adjustcolor( "red", alpha.f = 0.1))




## Compute number of EACs from much lower RR =3,4
phidot = function(s, phi, parms){
  phidot=numeric(6)
  with(as.list(parms),{
  tau = t - s
  #RR=4
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))#/1.5
  nu = nu0*((1-GERDpr)+RR*GERDpr)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## Armitage-Doll transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}

RR=3
age_groups = seq(42.5,87.5,by=5)
Lambda_m = rep(0,length(age_groups))
source('am_parameter_list_2020.R')
for (i in 1:length(Lambda_m)){
  parms_temp = get_params("M",(2010-age_groups[i]))
  Lambda_m[i] =  h_EAC(parms_temp,age_groups[i])*m_num_2010[i]
  print(h_EAC(parms_temp,age_groups[i])*100000)
}

Lambda_f = rep(0,length(age_groups))
source('af_parameter_list_2020.R')
for (i in 1:length(Lambda_f)){
  parms_temp = get_params("F",(2010-age_groups[i]))
  Lambda_f[i] =  h_EAC(parms_temp,age_groups[i])*f_num_2010[i]
  print(h_EAC(parms_temp,age_groups[i])*100000)
}
totalcases_2010=sum(Lambda_f,Lambda_m)
print(totalcases_2010)
## [1] 9842.859
RR=4
age_groups = seq(42.5,87.5,by=5)
Lambda_m = rep(0,length(age_groups))
source('am_parameter_list_2020.R')
for (i in 1:length(Lambda_m)){
  parms_temp = get_params("M",(2010-age_groups[i]))
  Lambda_m[i] =  h_EAC(parms_temp,age_groups[i])*m_num_2010[i]
  print(h_EAC(parms_temp,age_groups[i])*100000)
}

Lambda_f = rep(0,length(age_groups))
source('af_parameter_list_2020.R')
for (i in 1:length(Lambda_f)){
  parms_temp = get_params("F",(2010-age_groups[i]))
  Lambda_f[i] =  h_EAC(parms_temp,age_groups[i])*f_num_2010[i]
  print(h_EAC(parms_temp,age_groups[i])*100000)
}
totalcases_2010=sum(Lambda_f,Lambda_m)
print(totalcases_2010)
## [1] 9904.47

ind55= which(t_s[seq(16,66,1)]==55)
ind75 = which(t_s[seq(16,66,1)]==75)

print(strat1_w_GERD_noMortality50_3[seq(16,66,1),5][ind55]*100)
print(strat1_w_GERD_noMortality50_4[seq(16,66,1),5][ind55]*100)

print(strat1_w_GERD_noMortality50_3[seq(16,66,1),5][ind75]*100)
print(strat1_w_GERD_noMortality50_4[seq(16,66,1),5][ind75]*100)


