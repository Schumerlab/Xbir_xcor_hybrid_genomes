###########################
###new ABC - reviewer 1
##########################

#install.packages("KScorrect")
library("KScorrect")

sim<-1:1e6
popsize<-round(runif(2e6,0,8000))
gen<-round(runif(2e6,10,300))
admixprop<-runif(2e6,0.6,1)
mig1<-rlunif(2e6, 1e-6, 0.025, base = exp(10))
mig2<-rlunif(2e6, 1e-6, 0.025, base = exp(10))

params<-cbind(sim,popsize,gen,admixprop,mig1,mig2)

write.table(params,file="ABC_1M_initprop_0.5-1_sim_params_mig_revisions_v2.csv",sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
