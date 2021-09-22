rm(list=ls())

library(R.utils)
library(data.table)


## Parameters
args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
execpath=as.character(args$execpath)
dates_filepath=as.character(args$dates)
details_filepath=as.character(args$specs)
trajectoriesFile=as.character(args$trajectories)
nameRun=as.character(args$output)
filepath=as.character(args$filepath)
logpath=as.character(args$logpath)
seed=as.numeric(args$seed)

set.seed(seed)
source(paste0(execpath, "functions_HMM2014.R"))
sink(paste0(logpath, "/simul_summary_",nameRun,".txt"))

## Loading original data
dates=read.table(dates_filepath, skip=2)
dates_years=apply(apply(dates, 2, substring, first=1, last=4), 2, as.numeric)
details=read.table(details_filepath, skip=2)
details1=details # removing country and center to keep in format of HMM 2014 code
details=details[,-c(4,5)] # removing country and center to keep in format of HMM 2014 code


## Dates and years
init_date=min(dates_years[,1])
end_date=max(dates_years[,3])
end_year=end_date-init_date

start_dates=dates_years[,2] # date of recruitment
start_years=start_dates-init_date
init_year=0

if (is.null(trajectoriesFile)){
  # For simulations using R language (only mu, no deaths)
  true_mu=-3
  W=matrix(NA, ncol=end_year+2, nrow=length(start_years)) # + 1 for year 0 and + 1 for year end_year+1 (states not transitions)
  for (i in 1:length(start_years)){
    start_year=start_years[i]
    W[i,]=SimulateTrajectoriesMu(start_year=start_year, end_year=end_year, init_year=init_year, true_mu=true_mu)
  }
} else {
  # Simulations using simulate_data in C++
  trajectories=as.matrix(data.frame(fread(trajectoriesFile)))
  trajectories=ifelse(trajectories==0, yes="S", no=trajectories)
  trajectories=ifelse(trajectories==1, yes="I_1", no=trajectories)
  trajectories=ifelse(trajectories==2, yes="I_2", no=trajectories)
  trajectories=ifelse(trajectories==3, yes="R", no=trajectories)
  trajectories=ifelse(trajectories==4, yes="M", no=trajectories)
  W=trajectories
}
colnames(W)=paste0("Y",init_year:(end_year+1))

print(getwd())
print(trajectoriesFile)
print(dim(W))
print("States at end of follow-up:")
table(W[,ncol(W)])


## Compute DODiag
dodiag_years=apply(W, 1, function(W_t){
  if ("R"%in%W_t){
      gsub("Y","",names(W_t)[min(which(W_t=="R"))-1]) # Transition the year before entering the state
  } else {
    1900-init_date
  }})
dodiags=as.numeric(dodiag_years)+init_date
dodiags_dates=as.numeric(paste0(dodiags, "1215")) # arbitrary month and day


## Create data in same format
simul_dates=dates
simul_dates[,4]=dodiags_dates
simul_details=details
simul_details[,2]=ifelse(dodiags==1900, yes=0, no=1)
simul_details1=details1
simul_details1[,2]=ifelse(dodiags==1900, yes=0, no=1)


## Remove deaths
simul_dates[,3]=19000101
simul_details[,4]=1
simul_details1[,6]=1


## Compute DOD
dod_years=apply(W, 1, function(W_t){
  if ("M"%in%W_t){
    gsub("Y","",names(W_t)[min(which(W_t=="M"))-1]) # Transition the year before entering the state
  } else {
    1900-init_date
  }})
dod=as.numeric(dod_years)+init_date
dod_dates=as.numeric(paste0(dod, "1215")) # arbitrary month and day
simul_details1[,6]=ifelse(dod_years>0, yes=2, no=1)
simul_dates[,3]=dod_dates

cat("\n")
print(simul_dates[13,3])
cat("\n")


## Save simulated data
write.table(matrix(c(nrow(simul_dates), ncol(simul_dates)), ncol=1), paste0(filepath, "dates_index.txt"), 
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(simul_dates, paste0(filepath, "pre_dates_simul_", seed, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
system(paste0("cat ",filepath, "dates_index.txt ", filepath, "pre_dates_simul_", seed, ".txt > ", filepath, "dates_simul_", nameRun, "_", seed, ".txt"))
system(paste0("rm ", filepath, "pre_dates_simul_", seed, ".txt"))
system(paste0("rm ", filepath, "dates_index.txt"))

# write.table(matrix(c(nrow(simul_details), ncol(simul_details)), ncol=1), paste0(filepath, "details_index.txt"), 
#             quote=FALSE, row.names=FALSE, col.names=FALSE)
# write.table(simul_details, paste0("", filepath, "pre_details_simul_", seed, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
# system(paste0("cat ", filepath, "details_index.txt ", filepath, "pre_details_simul_", seed, ".txt > ", filepath, "details_simul_", nameRun, "_", seed, ".txt"))
# system(paste0("rm ", filepath, "pre_details_simul_", seed, ".txt"))
# system(paste0("rm ", filepath, "details_index.txt"))

write.table(matrix(c(nrow(simul_details1), ncol(simul_details1)), ncol=1), paste0(filepath, "details_index1.txt"), 
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(simul_details1, paste0("", filepath, "pre_details1_simul_", seed, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
system(paste0("cat ", filepath, "details_index1.txt ", filepath, "pre_details1_simul_", seed, ".txt > ", filepath, "details1_simul_", nameRun, "_", seed, ".txt"))
system(paste0("rm ", filepath, "pre_details1_simul_", seed, ".txt"))
system(paste0("rm ", filepath, "details_index1.txt"))










