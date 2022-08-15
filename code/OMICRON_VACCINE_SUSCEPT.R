########################################################################################################
################################ - OMICRON VACCINE SUSCEPT - #########################################
########################################################################################################
library(optiSolve)
library(quadprog)
library(polynom)
library(logitnorm)
library(stats)
library(dplyr)
library(matrixStats)
library(lubridate)
library(data.table)


rm(list=ls())
gc()

skage.groups_new = as.vector(read.csv('korea_population.csv',header=T)$x)
skage = read.csv('skage.csv',header=T)
skage = skage[1:101,5]
skage = as.numeric(sapply(skage, function(x) gsub(',','',x)))
skage_0_2 = sum(skage[1:3])
skage_3_4 = sum(skage[4:5])
skage_5_6 = sum(skage[6:7])
skage_7_9 = sum(skage[8:10])
skage_10_12 = sum(skage[11:13])

skage_10_11 = sum(skage[11:12])
skage_12_14 = sum(skage[13:15])
skage_13_14 = sum(skage[14:15])
skage_15 = skage[16]
skage_16_18 = sum(skage[17:19])

skage_15_17 = sum(skage[16:18])
skage_18_19 = sum(skage[19:20])

skage_19 = skage[20]
skage_18_19 = sum(skage[19:20])
skage_20_29 = sum(skage[21:30])

vaccine_omi = read.csv('vaccinated_omi_eng.csv',fileEncoding = 'euc-kr')
vaccine_omi = vaccine_omi[,1:57]

vaccine_res = copy(vaccine_omi)

########### Under 17 #################
num_17_1 = 384430
num_17_2 = 348070

num_16_1 = 320217
num_16_2 = 278506

num_15_1 = 315174
num_15_2 = 195138

num_14_1 = 258822
num_14_2 = 169633

num_13_1 = 247346
num_13_2 = 141436

num_12_1 = 161669
num_12_2 = 79587

vaccine_res[c(1,4,7,9),4:ncol(vaccine_omi)] = round(vaccine_omi[c(1,4,7,9),4:ncol(vaccine_omi)]*(num_14_1+num_13_1+num_12_1)/(num_17_1+num_16_1+num_15_1+num_14_1+num_13_1+num_12_1))
vaccine_res[c(1+11,4+11,7+11,9+11),4:ncol(vaccine_omi)] = round(vaccine_omi[c(1,4,7,9),4:ncol(vaccine_omi)]*(num_17_1+num_16_1+num_15_1)/(num_17_1+num_16_1+num_15_1+num_14_1+num_13_1+num_12_1))
vaccine_res[c(1+11,4+11,7+11,9+11),4:ncol(vaccine_omi)] = vaccine_res[c(1+11,4+11,7+11,9+11),4:ncol(vaccine_omi)] + vaccine_omi[c(1+11,4+11,7+11,9+11),4:ncol(vaccine_omi)] 
vaccine_res[c(2,5,10),4:ncol(vaccine_omi)] = round(vaccine_omi[c(2,5,10),4:ncol(vaccine_omi)]*(num_14_2+num_13_2+num_12_2)/(num_17_2+num_16_2+num_15_2+num_14_2+num_13_2+num_12_2))
vaccine_res[c(2+11,5+11,10+11),4:ncol(vaccine_omi)] = round(vaccine_omi[c(2,5,10),4:ncol(vaccine_omi)] *(num_17_2+num_16_2+num_15_2)/(num_17_2+num_16_2+num_15_2+num_14_2+num_13_2+num_12_2))
vaccine_res[c(2+11,5+11,10+11),4:ncol(vaccine_omi)] = vaccine_res[c(2+11,5+11,10+11),4:ncol(vaccine_omi)] + vaccine_omi[c(2+11,5+11,10+11),4:ncol(vaccine_omi)]

############## 70-80 ##########################################

vaccine_res[100:110,4:ncol(vaccine_res)] = vaccine_res[100:110,4:ncol(vaccine_res)] + vaccine_res[111:121,4:ncol(vaccine_res)]
vaccine_res = vaccine_res[1:110,]
vaccine_res[100:110,1] = '75+'
############## 60-70 #########################################

AZ_1 = vaccine_res %>% filter(X.1=='AZ' & X.2 =='Dose 1')

AZ_2 = vaccine_res %>% filter(X.1=='AZ' & (X.2 =='Dose 2' |X.2 == 'Dose 2 ') )

AZ_3 = vaccine_res %>% filter(X.1=='AZ' & (X.2 =='Dose 3' |X.2 == 'Dose 3 ') )

PF_1 = vaccine_res %>% filter(X.1=='P' & X.2 =='Dose 1')

PF_2 = vaccine_res %>% filter(X.1=='P' & (X.2 =='Dose 2'))

PF_3 = vaccine_res %>% filter(X.1=='P' & (X.2 =='Dose 3'))

M_1 = vaccine_res %>% filter(X.1=='M' & X.2 =='Dose 1')

M_2 = vaccine_res %>% filter(X.1=='M' & (X.2 =='Dose 2'))

M_3 = vaccine_res %>% filter(X.1=='M' & (X.2 =='Dose 3'))


JJ_1 = vaccine_res %>% filter(X.1=='JJ' & (X.2 =='Dose 1'))
JJ_3 = vaccine_res %>% filter(X.1=='JJ' & (X.2 =='Dose 3'))


AZ_1_res = matrix(0,nrow=500,ncol=16)
AZ_2_res = matrix(0,nrow=500,ncol=16)
AZ_3_res = matrix(0,nrow=500,ncol=16)
PF_1_res = matrix(0,nrow=500,ncol=16)
PF_2_res = matrix(0,nrow=500,ncol=16)
PF_3_res = matrix(0,nrow=500,ncol=16)
M_1_res = matrix(0,nrow=500,ncol=16)
M_2_res = matrix(0,nrow=500,ncol=16)
M_3_res = matrix(0,nrow=500,ncol=16)
JJ_1_res = matrix(0,nrow=500,ncol=16)
JJ_3_res = matrix(0,nrow=500,ncol=16)
start_day = as.Date('2021-02-26')-as.Date('2021-01-01') #56


k = start_day

for( j in 4:ncol(AZ_1)){
  if(j==4) {
    AZ_1_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(AZ_1[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    AZ_2_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(AZ_2[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    AZ_3_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(AZ_3[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    PF_1_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(PF_1[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    PF_2_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(PF_2[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    PF_3_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(PF_3[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    M_1_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(M_1[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    M_2_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(M_2[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    M_3_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(M_3[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    JJ_1_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(JJ_1[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    JJ_3_res[(k:(k+8)),c(3,4,13,14,15,16)] = matrix(round(JJ_3[c(1,2,7,8,9,10),5]/9),nrow=9,ncol=6,byrow=T)
    for(i in 3:6){
      AZ_1_res[(k:(k+8)),(2*i-1)] = round((AZ_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      AZ_1_res[(k:(k+8)),(2*i)] = round((AZ_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      AZ_2_res[(k:(k+8)),(2*i-1)] = round((AZ_2[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      AZ_2_res[(k:(k+8)),(2*i)] = round((AZ_2[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      AZ_3_res[(k:(k+8)),(2*i-1)] = round((AZ_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      AZ_3_res[(k:(k+8)),(2*i)] = round((AZ_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      
      PF_1_res[(k:(k+8)),(2*i-1)] = round((PF_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      PF_1_res[(k:(k+8)),(2*i)] = round((PF_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      PF_2_res[(k:(k+8)),(2*i-1)] = round((PF_2[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      PF_2_res[(k:(k+8)),(2*i)] = round((PF_2[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      PF_3_res[(k:(k+8)),(2*i-1)] = round((PF_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      PF_3_res[(k:(k+8)),(2*i)] = round((PF_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      
      
      M_1_res[(k:(k+8)),(2*i-1)] = round((M_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      M_1_res[(k:(k+8)),(2*i)] = round((M_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      M_2_res[(k:(k+8)),(2*i-1)] = round((M_2[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      M_2_res[(k:(k+8)),(2*i)] = round((M_2[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      M_3_res[(k:(k+8)),(2*i-1)] = round((M_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      M_3_res[(k:(k+8)),(2*i)] = round((M_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      
      
      JJ_1_res[(k:(k+8)),(2*i-1)] = round((JJ_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      JJ_1_res[(k:(k+8)),(2*i)] = round((JJ_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      JJ_3_res[(k:(k+8)),(2*i-1)] = round((JJ_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
      JJ_3_res[(k:(k+8)),(2*i)] = round((JJ_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/9)
    }
    k = k+9
  }
  else{
    AZ_1_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(AZ_1[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    AZ_2_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(AZ_2[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    AZ_3_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(AZ_3[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    PF_1_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(PF_1[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    PF_2_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(PF_2[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    PF_3_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(PF_3[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    M_1_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(M_1[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    M_2_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(M_2[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    M_3_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(M_3[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    JJ_1_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(JJ_1[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    JJ_3_res[(k:(k+6)),c(3,4,13,14,15,16)] = matrix(round(JJ_3[c(1,2,7,8,9,10),j]/7),nrow=7,ncol=6,byrow=T)
    for(i in 3:6){
      AZ_1_res[(k:(k+6)),(2*i-1)] = round((AZ_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      AZ_1_res[(k:(k+6)),(2*i)] = round((AZ_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      AZ_2_res[(k:(k+6)),(2*i-1)] = round((AZ_2[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      AZ_2_res[(k:(k+6)),(2*i)] = round((AZ_2[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      AZ_3_res[(k:(k+6)),(2*i-1)] = round((AZ_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      AZ_3_res[(k:(k+6)),(2*i)] = round((AZ_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      
      PF_1_res[(k:(k+6)),(2*i-1)] = round((PF_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      PF_1_res[(k:(k+6)),(2*i)] = round((PF_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      PF_2_res[(k:(k+6)),(2*i-1)] = round((PF_2[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      PF_2_res[(k:(k+6)),(2*i)] = round((PF_2[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      PF_3_res[(k:(k+6)),(2*i-1)] = round((PF_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      PF_3_res[(k:(k+6)),(2*i)] = round((PF_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      
      
      M_1_res[(k:(k+6)),(2*i-1)] = round((M_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      M_1_res[(k:(k+6)),(2*i)] = round((M_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      M_2_res[(k:(k+6)),(2*i-1)] = round((M_2[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      M_2_res[(k:(k+6)),(2*i)] = round((M_2[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      M_3_res[(k:(k+6)),(2*i-1)] = round((M_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      M_3_res[(k:(k+6)),(2*i)] = round((M_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      
      JJ_1_res[(k:(k+6)),(2*i-1)] = round((JJ_1[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      JJ_1_res[(k:(k+6)),(2*i)] = round((JJ_1[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      JJ_3_res[(k:(k+6)),(2*i-1)] = round((JJ_3[i,j]*skage.groups_new[2*i-1]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      JJ_3_res[(k:(k+6)),(2*i)] = round((JJ_3[i,j]*skage.groups_new[2*i]/(skage.groups_new[2*i-1]+skage.groups_new[2*i]))/7)
      
    }
    k = k+7    
  }
  
}


AZ_1_4_eff = 0.177
AZ_2_3_eff = 0.489
AZ_2_7_eff = 0.337
AZ_2_12_eff = 0.286
AZ_2_17_eff = 0.178
AZ_2_22_eff = 0.04
AZ_2_27_eff = 0.0
AZ_3_AZ_1_eff = 0.577
AZ_3_AZ_3_eff = 0.556
AZ_3_AZ_7_eff = 0.467
AZ_3_PF_1_eff = 0.588
AZ_3_PF_3_eff = 0.624
AZ_3_PF_7_eff = 0.529
AZ_3_PF_12_eff = 0.396
AZ_3_M_1_eff = 0.68
AZ_3_M_3_eff = 0.701
AZ_3_M_7_eff = 0.609


AZ_1_eff = function(t){
  if(t<=14){
    a=0
  }else{
    a = max(((AZ_1_4_eff-(AZ_1_4_eff+0.1))/(43-18))*(t-18)+0.277,0)
  }
  return(a)
}

AZ_2_eff = function(t){
  if(t<=14){
    a=0
  }
  else if(t>=15 & t<50){
    a = max(((AZ_2_7_eff-AZ_2_3_eff)/(50-22))*(t-50)+AZ_2_7_eff,0)
  }
  else if(t>=50 & t<85){
    a = max(((AZ_2_12_eff-AZ_2_7_eff)/(85-50))*(t-85)+AZ_2_12_eff,0)
  }
  else if(t>=85 & t<120){
    a = max(((AZ_2_17_eff-AZ_2_12_eff)/(120-85))*(t-120)+AZ_2_17_eff,0)
  }
  else if(t>=120 & t<155){
    a = max(((AZ_2_22_eff-AZ_2_17_eff)/(155-120))*(t-155)+AZ_2_22_eff,0)
  }
  else{
    a = max(((AZ_2_27_eff-AZ_2_22_eff)/(190-155))*(t-190)+AZ_2_27_eff,0)
  }
  return(a)
}

PF_1_3_eff = 0.428 #0.451
PF_1_4_eff = 0.315 #0.331
PF_2_3_eff = 0.655 #0.670
PF_2_7_eff = 0.487 #0.502
PF_2_12_eff = 0.301 #0.315
PF_2_17_eff = 0.154 #0.166
PF_2_22_eff = 0.115 #0.129
PF_2_27_eff = 0.088 #0.105
PF_3_1_eff = 0.669  #0.676
PF_3_3_eff = 0.672 #0.678
PF_3_7_eff = 0.55 #0.558
PF_3_12_eff = 0.457 #0.467

# PF_1_3_eff = 0.403 #0.451
# PF_1_4_eff = 0.299 #0.331
# PF_2_3_eff = 0.639 #0.670
# PF_2_7_eff = 0.471 #0.502
# PF_2_12_eff = 0.287 #0.315
# PF_2_17_eff = 0.142 #0.166
# PF_2_22_eff = 0.101 #0.129
# PF_2_27_eff = 0.070 #0.105
# PF_3_1_eff = 0.661  #0.676
# PF_3_3_eff = 0.665 #0.678
# PF_3_7_eff = 0.542 #0.558
# PF_3_12_eff = 0.447 #0.467

PF_1_eff = function(t){
  if(t<=14){
    a=0
  }else{
    a = max(((PF_1_4_eff-PF_1_3_eff)/(43-18))*(t-43)+PF_1_4_eff,0)
  }
  return(a)
}

PF_2_eff = function(t){
  if(t<=14){
    a=0
  }
  else if(t>=15 & t<50){
    a = max(((PF_2_7_eff-PF_2_3_eff)/(50-22))*(t-50)+PF_2_7_eff,0)
  }
  else if(t>=50 & t<85){
    a = max(((PF_2_12_eff-PF_2_7_eff)/(85-50))*(t-85)+PF_2_12_eff,0)
  }
  else if(t>=85 & t<120){
    a = max(((PF_2_17_eff-PF_2_12_eff)/(120-85))*(t-120)+PF_2_17_eff,0)
  }
  else if(t>=120 & t<155){
    a = max(((PF_2_22_eff-PF_2_17_eff)/(155-120))*(t-155)+PF_2_22_eff,0)
  }
  else{
    a = max(((PF_2_27_eff-PF_2_22_eff)/(190-155))*(t-190)+PF_2_27_eff,0)
  }
  return(a)
}



M_1_3_eff = 0.479 #0.523
M_1_4_eff = 0.319 #0.361

M_2_3_eff = 0.751 #0.787
M_2_7_eff = 0.528 #0.571
M_2_12_eff = 0.356 #0.384
M_2_17_eff = 0.253 #0.274
M_2_22_eff = 0.15 #0.182
M_2_27_eff = 0.149 #0.247

M_3_1_eff = 0.681 #0.705
M_3_3_eff = 0.663 #0.688

PF_3_1_eff = 0.669  #0.676
PF_3_3_eff = 0.672 #0.678
PF_3_7_eff = 0.55 #0.558
PF_3_12_eff = 0.457 #0.467


M_PF_3_1_eff = PF_3_1_eff*(2/3) + M_3_1_eff*(1/3)
M_PF_3_3_eff = PF_3_3_eff*(2/3) + M_3_3_eff*(1/3)
M_PF_3_7_eff = PF_3_7_eff
M_PF_3_12_eff = PF_3_12_eff


M_1_eff = function(t){
  if(t<=14){
    a=0
  }else{
    a = max(((M_1_4_eff-M_1_3_eff)/(43-18))*(t-43)+M_1_4_eff,0)
  }
  return(a)
}

M_2_eff = function(t){
  if(t<=14){
    a=0
  }
  else if(t>=15 & t<50){
    a = max(((M_2_7_eff-M_2_3_eff)/(50-22))*(t-50)+M_2_7_eff,0)
  }
  else if(t>=50 & t<85){
    a = max(((M_2_12_eff-M_2_7_eff)/(85-50))*(t-85)+M_2_12_eff,0)
  }
  else if(t>=85 & t<120){
    a = max(((M_2_17_eff-M_2_12_eff)/(120-85))*(t-120)+M_2_17_eff,0)
  }
  else if(t>=120 & t<155){
    a = max(((M_2_22_eff-M_2_17_eff)/(155-120))*(t-155)+M_2_22_eff,0)
  }
  else{
    a = max(((M_2_27_eff-M_2_22_eff)/(190-155))*(t-190)+M_2_27_eff,0)
  }
  return(a)
}

M_PF_3_eff = function(t){
  if(t<=7){
    a=0
  }
  else if(t>=8 & t<22){
    a = max(((M_PF_3_3_eff - M_PF_3_1_eff)/(22-11))*(t-22)+M_PF_3_3_eff,0)
  }
  else if(t>=22 & t<50){
    a = max(((M_PF_3_7_eff-M_PF_3_3_eff)/(50-22))*(t-50)+M_PF_3_7_eff,0)
  }
  else{
    a = max(((M_PF_3_12_eff-M_PF_3_7_eff)/(85-50))*(t-85)+M_PF_3_12_eff,0)
  }
  return(a)
}




AZ_1_index = rep(start_day,16)
AZ_2_index = rep(start_day,16)
PF_1_index = rep(start_day,16)
PF_2_index = rep(start_day,16)
M_1_index = rep(start_day,16)
M_2_index = rep(start_day,16)
JJ_1_index = rep(start_day,16)


AZ_1_res_temp = AZ_1_res
AZ_2_res_temp = AZ_2_res
AZ_3_res_temp = AZ_3_res

PF_1_res_temp = PF_1_res
PF_2_res_temp = PF_2_res
PF_3_res_temp = PF_3_res

M_1_res_temp = M_1_res
M_2_res_temp = M_2_res

JJ_1_res_temp = JJ_1_res
JJ_3_res_temp = JJ_3_res

vac_3_res_temp = AZ_3_res + PF_3_res + M_3_res + JJ_3_res
vac_3_res = vac_3_res_temp

vaccine_sus = matrix(0, nrow=500,ncol=16)

for(t in start_day:500){
  for(j in 3:16){
    while((sum(AZ_1_res_temp[1:AZ_1_index[j],j])<=AZ_2_res[t-15,j]) & AZ_2_res[t-15,j]!=0 & AZ_1_index[j]<=(t-15-15)){
      AZ_1_index[j] = AZ_1_index[j] + 1
    }
    AZ_1_res_temp[AZ_1_index[j],j] = sum(AZ_1_res_temp[1:AZ_1_index[j],j]) - AZ_2_res[t-15,j]
    AZ_1_res_temp[1:(AZ_1_index[j]-1),j]=0
    
    while((sum(M_1_res_temp[1:M_1_index[j],j])<=M_2_res[t-15,j]) & M_2_res[t-15,j]!=0  & M_1_index[j]<=(t-15-15) ){
      M_1_index[j] = M_1_index[j] + 1
    }
    M_1_res_temp[M_1_index[j],j] = sum(M_1_res_temp[1:M_1_index[j],j]) - M_2_res[t-15,j]
    M_1_res_temp[1:(M_1_index[j]-1),j]=0
    
    
    if(t<start_day+30){
      while((sum(PF_1_res_temp[1:PF_1_index[j],j])<=PF_2_res[t-15,j]) & PF_2_res[t-15,j]!=0 ){
        PF_1_index[j] = PF_1_index[j] + 1
      }
    }
    else{
      while((sum(PF_1_res_temp[1:PF_1_index[j],j])<=PF_2_res[t-15,j]) & PF_2_res[t-15,j]!=0  & PF_1_index[j]<=(t-15-15) ){
        PF_1_index[j] = PF_1_index[j] + 1
      }
    }
    
    if(sum(PF_1_res_temp[1:PF_1_index[j],j]) >= PF_2_res[t-15,j]){
      PF_1_res_temp[PF_1_index[j],j] = sum(PF_1_res_temp[1:PF_1_index[j],j]) - PF_2_res[t-15,j]
      PF_1_res_temp[1:(PF_1_index[j]-1),j]=0
    }
    else{
      PF_1_res_temp[PF_1_index[j],j] = 0
      PF_1_res_temp[1:(PF_1_index[j]-1),j]=0
      while((sum(AZ_1_res_temp[1:AZ_1_index[j],j])<= - sum(PF_1_res_temp[1:PF_1_index[j],j]) + PF_2_res[t-15,j] ) & AZ_1_index[j]<=(t-15-15) ){
        AZ_1_index[j] = AZ_1_index[j] + 1
      }
      temp_vac = sum(AZ_1_res_temp[1:AZ_1_index[j],j]) + sum(PF_1_res_temp[1:PF_1_index[j],j]) - PF_2_res[t-15,j]
      if(temp_vac>=0){
        AZ_1_res_temp[AZ_1_index[j],j] = sum(AZ_1_res_temp[1:AZ_1_index[j],j]) + sum(PF_1_res_temp[1:PF_1_index[j],j]) - PF_2_res[t-15,j]
        AZ_1_res_temp[1:(AZ_1_index[j]-1),j]=0
      }
      else{
        AZ_1_res_temp[AZ_1_index[j],j] = 0
        AZ_1_res_temp[1:(AZ_1_index[j]-1),j]=0
        
        while((sum(M_1_res_temp[1:M_1_index[j],j])<= -temp_vac) & M_1_index[j]<=(t-15-15) ){
          M_1_index[j] = M_1_index[j] + 1
        }
        if(sum(M_1_res_temp[1:M_1_index[j],j]) + temp_vac >=0 ){
          M_1_res_temp[M_1_index[j],j] = sum(M_1_res_temp[1:M_1_index[j],j]) + temp_vac
          M_1_res_temp[1:(M_1_index[j]-1),j]=0
        }else{
          M_1_res_temp[M_1_index[j],j] = 0
          M_1_res_temp[1:(M_1_index[j]-1),j]=0
        }
      }
    }
    
    ###################################Dose 3
    temp_vac = vac_3_res[t-8,j]
    while(temp_vac > 0){
      temp_AZ_vac =  AZ_2_res_temp[AZ_2_index[j],j] - temp_vac
      temp_vac = -temp_AZ_vac
      if(temp_AZ_vac>=0){
        AZ_2_res_temp[AZ_2_index[j],j] = temp_AZ_vac
      }else{
        AZ_2_res_temp[AZ_2_index[j],j] = 0
        temp_JJ_vac =  JJ_1_res_temp[JJ_1_index[j],j] + temp_AZ_vac
        temp_vac = - temp_JJ_vac
        if(temp_JJ_vac>=0){
          JJ_1_res_temp[JJ_1_index[j],j] = temp_JJ_vac
        }else{
          JJ_1_res_temp[JJ_1_index[j],j] = 0
          
          temp_M_vac =  M_2_res_temp[M_2_index[j],j] + temp_JJ_vac
          temp_vac = - temp_M_vac
          if(temp_M_vac>=0){
            M_2_res_temp[M_2_index[j],j] = temp_M_vac
          }else{
            M_2_res_temp[M_2_index[j],j] = 0
            
            temp_PF_vac =  PF_2_res_temp[PF_2_index[j],j] + temp_M_vac
            temp_vac = - temp_PF_vac
            if(temp_PF_vac>=0){
              PF_2_res_temp[PF_2_index[j],j] = temp_PF_vac
            }else{
              PF_2_res_temp[PF_2_index[j],j] = 0
              AZ_2_index[j] = AZ_2_index[j] + 1
              JJ_1_index[j] = JJ_1_index[j] + 1
              M_2_index[j] = M_2_index[j] + 1
              PF_2_index[j] = PF_2_index[j] + 1
            }
          }
        }
      }
    }
  }
  for(s in as.numeric(start_day):t){
    vaccine_sus[t,] = vaccine_sus[t,] + AZ_1_res_temp[s,]*AZ_1_eff(t-s) + AZ_2_res_temp[s,]*AZ_2_eff(t-s) + PF_1_res_temp[s,]*PF_1_eff(t-s) + PF_2_res_temp[s,]*PF_2_eff(t-s)+ M_1_res_temp[s,]*M_1_eff(t-s) + M_2_res_temp[s,]*M_2_eff(t-s) + JJ_1_res_temp[s,]*AZ_1_eff(t-s) + vac_3_res_temp[s,]*M_PF_3_eff(t-s)
  }
  
}

