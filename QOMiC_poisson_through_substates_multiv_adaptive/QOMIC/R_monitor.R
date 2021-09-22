rm(list=ls())

library(R.utils)
library(grDevices)
library(colorspace)
library(MASS)
library(scales)
library(pheatmap)
library(RColorBrewer)


##### Arguments

args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
execpath=as.character(args$execpath)
nameFolder=as.character(args$namefolder)
nameRun=as.character(args$name)
NIter=as.character(args$iter)
BI=as.character(args$burn_in)
seeds=as.character(args$seed)

execpath = "/Users/adamjackson/Documents/summer_project/general_model/QOMiC_poisson_through_substates_multiv_adaptive/QOMIC/"
filepath = "/Users/adamjackson/Documents/summer_project/general_model/"
nameFolder=paste0(filepath, "Results_simulation_2_1")
nameRun="2_1"
NIter=as.character(format(100000, scientific = F))
BI=as.character(format(10000, scientific = F))
seeds=1

source(paste0(execpath, "/monitoring_functions.R"))

print(paste("Folder:", nameFolder))
print(paste("Name of the run:", nameRun))
print(paste("Number of iterations:", NIter))
print(paste("Burn in:", BI))
print(paste("seed:", seeds))

legCol <- c('darkgreen', 'orange', 'darkred', 'skyblue')

# nameFolder="Results/C_output/Run_adaptive_MWG/Run_poisson/Results_simulation_mu_minus2_gamma_3_1/"
# nameRun="simul_mu_minus2_gamma_3_expo"
# NIter="10000"
# BI=1000
# Nseed=1


##### Loading QOMiC outputs

mymodels=NULL
for (k in seeds){
  # for (k in 1){
  print(k)
  tmp <- read.table(paste(nameFolder,"/",nameRun,'_',k,'_History_iter_', NIter, '_k_2.txt',sep=''),header=T,stringsAsFactors=F)
  assign(paste0("Res",k), tmp)
  mymodels=c(mymodels, paste0("Res",k))
}


##### Path to outputs

figureFolder=paste0(nameFolder,"/Figures/")
dir.create(paste0(figureFolder), showWarnings=FALSE)
print(paste("Output path:", figureFolder))


##### Parameters and likelihood over MCMC iterations

{pdf(paste(figureFolder, '/Monitor_params_',nameRun,'_',paste(seeds, collapse='_'),'_',NIter,'iter.pdf',sep=''),height=15,width=20)
  par(mfrow=c(3,3), mar=c(5,5,2,2))
  layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
  
  MonitorParams(models=mymodels, param="Mu", name="mu", init=1, legend=TRUE)
    for (k in 1:4){
      MonitorParams(models=mymodels, param=paste0("Lambda",k), name=paste0("lambda[",k,"]"), init=1, legend=FALSE)
    }
    MonitorParams(models=mymodels, param="Gamma", name=expression(gamma), init=1, legend=FALSE, log_scale = FALSE)
    MonitorParams(models=mymodels, param="L", name="L", init=1, legend=FALSE, sigma=FALSE)
  dev.off()
}


{pdf(paste(figureFolder, '/Monitor_params_',nameRun,'_',paste(seeds, collapse='_'),'_afterBI_',NIter,'iter.pdf',sep=''),height=15,width=20)
  par(mfrow=c(3,3), mar=c(5,5,2,2))
  layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
  
  MonitorParams(models=mymodels, param="Mu", name="mu", init=BI, legend=TRUE)
    for (k in 1:4){
      MonitorParams(models=mymodels, param=paste0("Lambda",k), name=paste0("lambda[",k,"]"), init=BI, legend=FALSE)
    }
    MonitorParams(models=mymodels, param="Gamma", name=expression(gamma), init=BI, legend=FALSE, log_scale = FALSE)
    MonitorParams(models=mymodels, param="L", name="L", init=BI, legend=FALSE, sigma=FALSE)
  dev.off()
}


##### Posterior distribution (starting from BI)

{pdf(paste(figureFolder, '/Posterior_',nameRun,'_',paste(seeds, collapse='_'),'_',NIter,'iter.pdf',sep=''),height=15,width=20)
  par(mfrow=c(3,3), mar=c(5,5,2,2))
  layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
  
  GetPosterior(models=mymodels, param="Mu", name=expression(mu), init=BI, legend=TRUE)
  for (k in 1:4){
    GetPosterior(models=mymodels, param=paste0("Lambda",k), name=substitute(lambda[num], list(num=k)), init=BI, legend=TRUE)
  }
  GetPosterior(models=mymodels, param="Gamma", name=expression(gamma), init=BI, legend=TRUE)
  GetPosterior(models=mymodels, param="L", name="L", init=BI, legend=TRUE)
  dev.off()
}


##### Likelihood as a function of parameters 

{png(paste(figureFolder, '/Monitor_Likelihood_',nameRun,'_',paste(seeds, collapse='_'),'_',NIter,'iter.png',sep=''),
     height=10,width=20, unit='in', res=300)
  par(mfrow=c(2,3), mar=c(5,5,2,2))
  
  MonitorLikelihood(models=mymodels, param="Mu", name=expression(mu), init=BI, legend=TRUE)
  for (k in 1:4){
    MonitorLikelihood(models=mymodels, param=paste0("Lambda",k), name=substitute(lambda[num], list(num=k)), init=BI, legend=TRUE)
  }
  MonitorLikelihood(models=mymodels, param="Gamma", name=expression(gamma), init=BI, legend=TRUE)
  dev.off()
}


##### Correlation between parameters

Res1=eval(parse(text=mymodels[1]))
Res1$Gamma_log=log(Res1$Gamma)

#params=matrix(c("Mu","Lambda1","Mu","Gamma","Lambda1","Gamma"),ncol=2,byrow=TRUE)
#params_names=matrix(c("mu","lambda[1]","mu","gamma","lambda[1]","gamma"),ncol=2,byrow=TRUE)

params=matrix(c("Mu","Lambda1","Mu","Gamma_log","Lambda1","Gamma_log"),ncol=2,byrow=TRUE)
params_names=matrix(c("mu","lambda[1]","mu","log(gamma)","lambda[1]","log(gamma)"),ncol=2,byrow=TRUE)

loglik_cuts=as.numeric(as.character(cut(Res1$L[BI:nrow(Res1)], breaks = 100, labels = 1:100)))
loglik_col=colorRampPalette(c("gold","darkred"))(100)

{png(paste(figureFolder, '/Scatter_plot_params_',nameRun,'_',paste(seeds, collapse='_'),'_',NIter,'iter.png',sep=''),
     height=7,width=20, unit='in', res=300)
  par(mfrow=c(1,3), mar=c(5,5,1,1))
  for (k in 1:nrow(params)){
    mybandwidth=0.3
    myngrid=200
    mylwd=1
    x=eval(parse(text=paste0("Res1$",params[k,1],"[",BI,":",nrow(Res1),"]")))
    y=eval(parse(text=paste0("Res1$",params[k,2],"[",BI,":",nrow(Res1),"]")))
    z <- kde2d(x,y,
               n=myngrid)
    contour(z, las=1, main="", cex.main=2,
            xlab=eval(parse(text=paste0("expression(",params_names[k,1],")"))),
            ylab=eval(parse(text=paste0("expression(",params_names[k,2],")"))),
            cex.axis=1.5, cex.lab=1.5,
            lwd=mylwd, cex.lab=1.5, col=darken(legCol[1],amount=0.5),
            pch=19, cex=0.5, las=1, lty=1, drawlabels=FALSE)
    abline(h=axTicks(2),lty=3,col="grey",lwd=1.5)
    abline(v=axTicks(1),lty=3,col="grey",lwd=1.5)
    points(x,y,
           #col=alpha(lighten(legCol[1], amount=0.0001),0.4),
           col=alpha(loglik_col[loglik_cuts],0.4),
           pch=16, cex=0.5)
    contour(z, las=1, main="", cex.main=2,
            lwd=mylwd, cex.lab=1.5, col=darken(legCol[1],amount=0.5),
            pch=19, cex=0.5, las=1, lty=1, drawlabels=FALSE, add=TRUE)
    #mycor=formatC(cor(x,y, method="spearman"), format = "f", digits=2)
    mycor=formatC(cor(x,y), format = "f", digits=2)
    legend("topleft", legend=eval(parse(text=paste0("expression(rho*'=",mycor,"')"))),
           bty="n", cex=2)
  }
  dev.off()}


##### Correlation of the proposal

mycov=matrix(as.numeric(Res1[nrow(Res1),grep("sigma", colnames(Res1))]),
ncol=sqrt(length(grep("sigma", colnames(Res1)))), byrow=TRUE)

pheatmap(mycov, cluster_rows=FALSE, cluster_cols=FALSE, border=NA, display_numbers=TRUE, fontsize=20, number_format = "%.2e",
filename=paste(figureFolder, '/Heatmap_proposal_cov_',nameRun,'_',paste(seeds, collapse='_'),'_',NIter,'iter.png',sep=''))

pheatmap(cov2cor(mycov), cluster_rows=FALSE, cluster_cols=FALSE, border=NA, display_numbers=TRUE, fontsize=20,
filename=paste(figureFolder, '/Heatmap_proposal_cor_',nameRun,'_',paste(seeds, collapse='_'),'_',NIter,'iter.png',sep=''))

