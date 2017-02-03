# RedOakABC

### 1/ PROGRAMS
 
- To compile mscalc, the calculator of summary statistics (Ross-Ibarra et al. 2008; 2009; Roux et al. 2011):
 
 gcc *.c -lm -o mscalc
 
- To compile msnsam, the coalescent simulator (Ross-Ibarra et al. 2008):
 
 ./clms
 
- priorgen4_recentbottle.py, a prior generator.

To have more details: priorgen4_recentbottle.py -h

This generator is an extended version of the priorgen.py prior generator of Leroy et al. submitted (6 models: SI, AM, IM, SC, PAM, PSC; homo/hetero Ne; homo/hetero Nm; with or without recent bottlenecks*) (the Leroy et al. version of the software is also available: https://github.com/ThibaultLeroyFr/WhiteOaksABC)

Note that this version assumes recent bottlenecks. More precisely, bottlenecks are assumed to occur between the most recent event and present (i.e. posterior to the last secondary contact or ancient migration) .  
 
### 2/ DATASETS
 
All oak datasets used for our ABC analyses are available. For each pair of clusters, a bpfile, a spinput.txt and a file containing the summary statistics for the real dataset are available. To perform the multilocus coalescent simulations, the bpfile & spinput.txt files are required (i.e. need to be in your current directory). 
  
The target file contains the 19 summary statistics calculated on each real dataset.

### 3/ EXAMPLE: Multilocus coalescent simulations:
 - Introduction

25000 multilocus simulations assuming an AM scenario between cluster1 & cluster2 [i.e. number of SNPs (=68) x number of simulations (=25000) = 1,360,000]

Number of SNPs = 2nd line of the spinput.txt file
 - Bash script (note that here all programs are assumed to be in your bin directory):

> mknod myfifo p
 
> priorgen4_recentbottle.py bpfile=bpfile n1=0 n1=100 n2=0 n2=100 nA=0 nA=100 tau=0 tau=100 bottleneck=N taubottle=0 taubottle=10 alpha1=1 alpha1=1 alpha2=1 alpha2=1 M1=0 M1=100 M2=0 M2=100 shape1=0 shape1=100 shape2=0 shape2=500 model=AM nreps=20000 Nvariation=hetero Mvariation=hetero symMig=asym parameters=priorfile | msnsam tbs 1360000 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs >myfifo &
 
> mscalc < myfifo
 
### 4/ EXAMPLE: ABC Model Choice (Rscript):

here for each model, we assume : 200 directories containing an ABCstat.txt file containing 20,000 lines of summary statistics =4,000,000 simulations/model

> library(nnet)

> ### Import cv4ABC (Csillery et al. 2012)
> setwd("/home/tleroy/work2/ABC/")
> source("cv4abc.R")

> ### Real dataset: import summary statistics
> setwd("YOURPATH")
> target=read.table("sumstat_realdata_introduite1vs3.txt")
> ss=c(2:20)
> target=as.numeric(target)
> Print(target)

> ### Import each line of summary statistics corresponding to your simulations

> M1heteroNheteroMnobottle=M2heteroNheteroMnobottle=M3heteroNheteroMnobottle=M4heteroNheteroMnobottle=M5heteroNheteroMnobottle=M6heteroNheteroMnobottle=NULL

> for(i in 1:200){
  M1heteroNheteroMnobottle=rbind(M1heteroNheteroMnobottle, matrix(scan(paste("cluster1_cluster2_strictisolation_nobottleneck_Nhetero_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])
  M2heteroNheteroMnobottle=rbind(M2heteroNheteroMnobottle, matrix(scan(paste("cluster1_cluster2_ancientmig_nobottleneck_Nhetero_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])
  M3heteroNheteroMnobottle=rbind(M3heteroNheteroMnobottle, matrix(scan(paste("cluster1_cluster2_island_nobottleneck_Nhetero_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])
  M4heteroNheteroMnobottle=rbind(M4heteroNheteroMnobottle, matrix(scan(paste("cluster1_cluster2_2ndarycont_nobottleneck_Nhetero_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])
  M5heteroNheteroMnobottle=rbind(M5heteroNheteroMnobottle, matrix(scan(paste("cluster1_cluster2_multiancientmig_nobottleneck_Nhetero_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])
  M6heteroNheteroMnobottle=rbind(M6heteroNheteroMnobottle, matrix(scan(paste("cluster1_cluster2_multi2ndarycont_nobottleneck_Nhetero_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])
}
print(nrow(M1heteroNheteroMnobottle))
print(nrow(M2heteroNheteroMnobottle))
print(nrow(M3heteroNheteroMnobottle))
print(nrow(M4heteroNheteroMnobottle))
print(nrow(M5heteroNheteroMnobottle))
print(nrow(M6heteroNheteroMnobottle))

> ### Replace all "NaN" [!!! The number of lines of summary statistics need to be the same for the 6 models !!!]
> for(i in 1:ncol(M1heteroNheteroMnobottle)){
  M1heteroNheteroMnobottle[which(M1heteroNheteroMnobottle[,i]=="NaN"),i]=mean(M1heteroNheteroMnobottle[,i], na.rm=T)
  M2heteroNheteroMnobottle[which(M2heteroNheteroMnobottle[,i]=="NaN"),i]=mean(M2heteroNheteroMnobottle[,i], na.rm=T)
  M3heteroNheteroMnobottle[which(M3heteroNheteroMnobottle[,i]=="NaN"),i]=mean(M3heteroNheteroMnobottle[,i], na.rm=T)
  M4heteroNheteroMnobottle[which(M4heteroNheteroMnobottle[,i]=="NaN"),i]=mean(M4heteroNheteroMnobottle[,i], na.rm=T)
  M5heteroNheteroMnobottle[which(M5heteroNheteroMnobottle[,i]=="NaN"),i]=mean(M5heteroNheteroMnobottle[,i], na.rm=T)
  M6heteroNheteroMnobottle[which(M6heteroNheteroMnobottle[,i]=="NaN"),i]=mean(M6heteroNheteroMnobottle[,i], na.rm=T)
}

> ### To perform several ABC analyses, duplicate your real dataset  (here 20x)*

> obs=matrix(rep(target[ss],20), byrow=T, nrow=20)

> print(obs)

> ### Generate a long vector of numbers from the simulations (1=model1, 2=model2..., 4=model4), used as the dependent variable for the regression

> x=c(rep(1:6, each=nrow(M1heteroNheteroMnobottle)))

> print(summary(x))

> ### Perform your ABC analysis [note that tol is the most important paramters = required proportion of points nearest the target values (here 30,000/24,000,000 = 0,00125 best simulations)

> res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(M1heteroNheteroMnobottle,M2heteroNheteroMnobottle,M3heteroNheteroMnobottle,M4heteroNheteroMnobottle,M5heteroNheteroMnobottle,M6heteroNheteroMnobottle), tol=30000/(6*nrow(M1heteroNheteroMnobottle)), noweight=F, rejmethod=F, nb.nnet=20, size.nnet=8, output="OBS_introduite_SI-AM-CM-SC-PAM-PSC_nobottleneck_heteroM_10ABC_tol30000_190616")

### 5/ EXAMPLE: generate posterior distribution under the best model (Rscript):

> ### Import cv4estimations.R ((Csillery et al. 2012))
>  source("cv4estimations.R")

> ### Real dataset: import summary statistics (real data set)

> setwd("YOURPATH/")

> ss=c(2:20)

> ncv=10 #number of validations to perform
 
> tmp=as.numeric(read.table("sumstat_realdata_introduite1vs3.txt")[ss])
 
> ### Simulations: Import summary statistics 

> [edit the number of directory (for i in 1:X..) to compile all your priorfile & ABCstat.txt files; please check the number of columns in your files!!! ]

> target=NULL

> tmp.stat=NULL

> tmp.parameters=NULL
 
> for(i in 1:X){ 
   tmp.stat=matrix(scan("MYSUMSTATSFILE"), byrow=T, ncol=20)[,ss] 
   tmp.parameters=matrix(scan("MYPRIORFILE"), byrow=T, ncol=13) # ncol à modifier
   if(nrow(tmp.stat)==nrow(tmp.parameters)){
     prior.tmp=cbind(tmp.stat, tmp.parameters)
     prior=rbind(prior, prior.tmp)
   } 
  
  > ### Generate posteriors 
  
  > prior=na.omit(prior)
  
  > abc_nnet_multivar(target=target, x=prior[, -(1:length(ss))], sumstat=prior[, 1:length(ss)], tol=10000/(nrow(prior)), rejmethod=F, noweight=F, transf=rep("logit", ncol(prior)-length(ss)), bb=rbind(apply(prior[,-(1:length(ss))], MARGIN=2, FUN="min"), apply(prior[,-(1:length(ss))], MARGIN=2, FUN="max")), nb.nnet=25, size.nnet=10, trace=T, output="Posterior_10000_SC_MHETERO_v2_")
