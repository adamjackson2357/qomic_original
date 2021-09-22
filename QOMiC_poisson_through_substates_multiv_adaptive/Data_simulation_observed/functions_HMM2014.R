# (W_t)_t denotes the observed chain
# (Z_t)_t denotes the hidden chain


## P(Z_t=S)
Compute_H_S=function(P_S_I, years_vector){
  H_S=rep(NA, length(years_vector))
  H_S[1]=1
  for (k in 2:length(years_vector)){
    H_S[k]=(1-P_S_I)^(k-1)
  }
  return(H_S)
}


## P(Z_t=I_1)
Compute_H_1=function(P_S_I, years_vector){
  H_1=rep(NA, length(years_vector))
  H_1[1]=0
  for (k in 2:length(years_vector)){
    H_1[k]=P_S_I*(1-P_S_I)^(k-2)
  }
  return(H_1)
}


## P(Z_t=I_2)
Compute_H_2=function(P_S_I, years_vector){
  H_2=rep(NA, length(years_vector))
  H_2[1:2]=0
  for (k in 3:length(years_vector)){
    H_2[k]=0.5*P_S_I*(1-P_S_I)^(k-3)
  }
  return(H_2)
}


## P(W_t=SI)
Compute_H_SI=function(P_S_I, years_vector){
  H_S=Compute_H_S(P_S_I, years_vector)
  H_1=Compute_H_1(P_S_I, years_vector)
  H_2=Compute_H_2(P_S_I, years_vector)
  
  H_SI=rep(NA, length(years_vector))
  for (k in 1:length(years_vector)){
    H_SI=H_S+H_1+H_2
  }
  return(H_SI)
}


## P(W_t=SI inter W_t+1=R)
Compute_P_SI_inter_R=function(P_S_I, years_vector) {
  P_SI_inter_R=rep(NA, length(years_vector))
  ti=years_vector[1]
  for (k in 1:length(years_vector)){
    year=years_vector[k]
    if (year<ti+1) {
      P_SI_inter_R[k]=0
    } else if (year==ti+1) {
      P_SI_inter_R[k]=(1-P_S_I)^(year-1-ti) * P_S_I/2
    } else if (year-2>=ti) {
      P_SI_inter_R[k]=(1-P_S_I)^(year-1-ti) * P_S_I/2 + (1-P_S_I)^(year-2-ti) * P_S_I/2
    }
  }
  return(P_SI_inter_R)
}


## P(W_t+1=R | W_t=SI)
Compute_P_SI_R=function(P_S_I, years_vector){
  P_SI_inter_R=Compute_P_SI_inter_R(P_S_I, years_vector)
  H_SI=Compute_H_SI(P_S_I, years_vector)
  
  P_SI_R=P_SI_inter_R/H_SI
  return(P_SI_R)
}


ComputeLikelihood=function(transitions_SI_SI, transitions_SI_R, mu, years_vector){
  P_S_I=exp(mu)/(1+exp(mu))
  P_SI_R=Compute_P_SI_R(P_S_I, years_vector)
  start_year=years_vector[1]
  P_SI_R=c(rep(0,start_year), P_SI_R)
  tot_log_likelihood=0
  for (i in 1:nrow(transitions_SI_SI)){
    log_lik=sum(log((1-P_SI_R)^transitions_SI_SI[i,]))+sum(log(P_SI_R^transitions_SI_R[i,]))
    tot_log_likelihood=tot_log_likelihood+log_lik
  }
  return(tot_log_likelihood)
}


ComputeLikelihoodFormer=function(transitions_SI_SI, transitions_SI_R, mu, years_vector){
  P_S_I=exp(mu)/(1+exp(mu))
  P_SI_R=Compute_P_SI_inter_R(P_S_I, years_vector) # typo in previous QOMiC version
  start_year=years_vector[1]
  P_SI_R=c(rep(0,start_year), P_SI_R)
  tot_log_likelihood=0
  for (i in 1:nrow(transitions_SI_SI)){
    log_lik=sum(log((1-P_SI_R)^transitions_SI_SI[i,]))+sum(log(P_SI_R^transitions_SI_R[i,]))
    tot_log_likelihood=tot_log_likelihood+log_lik
  }
  return(tot_log_likelihood)
}


SimulateDODiagMu=function(start_year, end_year, true_mu){
  ## Computing true probability of SI-R transition
  years_vector=seq(start_year, end_year)
  
  P_S_I=exp(true_mu)/(1+exp(true_mu))
  
  P_diag=Compute_P_SI_inter_R(P_S_I=P_S_I, years_vector=years_vector)
  # P_diag=P_diag[-1]
  P_diag=c(P_diag, 1-sum(P_diag))
  names(P_diag)=c(years_vector, 1900-init_date)
  
  ## Simulation of dates of diagnosis
  date_diag=sample(as.numeric(names(P_diag))+init_date, size=1, replace=TRUE, prob=P_diag)
  return(date_diag)
}


SimulateTrajectoriesMu=function(start_year, end_year, init_year, true_mu){
  ## Computing true probability of SI-R transition
  years_vector=seq(start_year, end_year)
  
  P_S_I=exp(true_mu)/(1+exp(true_mu))
  P_Ii_Ij=matrix(c(0,0.5,0.5,0,
                   0,0,1,0,
                   0,0,1,0,
                   0,0,0,1), ncol=4, byrow=TRUE)
  
  W_t=rep("S",end_year-start_year+2)
  names(W_t)=start_year:(end_year+1) # last time point to show transitions happening in last year (2010)
  for (k in 2:(length(start_year:end_year)+1)){ # should be generalised for no mortality and R absorbing state
    # print(k)
    rand=runif(1)
    if (W_t[k-1]=="S"){
      # print(W_t[k-1])
      if (rand<P_S_I){
        # print("S-I transition")
        W_t[k]="I_1"
      } 
    }
    if (W_t[k-1]=="I_1"){
      if (rand<P_Ii_Ij[1,1]){
        # print("I_1-I_1 transition")
        W_t[k]="I_1"
      }
      if ((P_Ii_Ij[1,1]<rand)&(rand<(P_Ii_Ij[1,1]+P_Ii_Ij[1,2]))){
        # print("I_1-I_2 transition")
        W_t[k]="I_2"
      } 
      if (rand>(P_Ii_Ij[1,1]+P_Ii_Ij[1,2])){
        # print("I_1-R transition")
        W_t[k]="R"
      }
    }
    if (W_t[k-1]=="I_2"){
      if (rand<P_Ii_Ij[2,2]){
        # print("I_2-I_2 transition")
        W_t[k]="I_2"
      } 
      if ((P_Ii_Ij[2,2]<rand)&(rand<(P_Ii_Ij[2,2]+P_Ii_Ij[2,3]))){
        # print("I_2-R transition")
        W_t[k]="R"
      } 
    }
    if (W_t[k-1]=="R"){
      W_t[k]="R"
    }
  }
  # W[i,]=W_t
  # if ("R"%in%W_t){
  #   date_diag[i]=names(W_t)[min(which(W_t=="R"))-1]
  # }
  fill=rep("S",start_year-init_year)
  names(fill)=init_year:(start_year-1)
  W_t=c(fill, W_t)
  
  return(W_t)
}



