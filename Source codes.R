##Codes for running lasso-cox regression
library(glmnet)
lasso <- read.table("training features.csv", sep = ",", head = T,row.names=1)

fitx <- as.matrix(lasso[,1:90])
fity <- as.matrix(lasso[,91:92])

fit <- glmnet(fitx,fity,family = "cox")
plot(fit,xvar = "lambda",label = T)
set.seed(1111)
cvfit <- cv.glmnet(fitx,fity,family="cox",type.measure="deviance",nfolds=10)
plot(cvfit)

coef(cvfit,s=c(cvfit$lambda.min))
cvfit$lambda.min
#---------------------------------------------------------------------------------------------------------------------------------------------------


##Codes for cut-point of PSGC determination
cutpoint <- read.table("cut-point.csv", sep = ",", head = T,row.names=1)
library(survminer)
cut_point <- surv_cutpoint(data = cutpoint, time = "time", event = "status", variables = "sig",
                           minprop = 0.1, progressbar = TRUE)
plot(cut_point, palette = c("#CD3333", "#1874CD"), legend = "")
cut_point
#---------------------------------------------------------------------------------------------------------------------------------------------------

##Codes for generating survival curves
training_data <- read.table("training data.csv", sep = ",", head = T,row.names=1)
validation_data <- read.table("validation data.csv", sep = ",", head = T,row.names=1)
library(survminer)
library(survival)
Signature_OS_training <- survfit(Surv(OS,OSstatus)~Sub,data=training_data)
ggsurvplot(Signature_OS_training,training_data,size=1,
           palette = "nejm",pval = TRUE,risk.table = TRUE,risk.table.col = "strata",
           break.x.by=10,break.y.by=0.2,xlim=c(0,72),axes.offset=T)

Signature_DFS_training <- survfit(Surv(DFS,DFSstatus)~Sub,data=training_data)
ggsurvplot(Signature_DFS_training,training_cohort,size=1,
           palette = "nejm",pval = TRUE,risk.table = TRUE,risk.table.col = "strata",
           break.x.by=10,break.y.by=0.2,xlim=c(0,72),axes.offset=T)

Signature_OS_validation <- survfit(Surv(OS,OSstatus)~Sub,data=validation_data)
ggsurvplot(Signature_OS_validation,validation_data,size=1,
           palette = "nejm",pval = TRUE,risk.table = TRUE,risk.table.col = "strata",
           break.x.by=10,break.y.by=0.2,xlim=c(0,72),axes.offset=T)

Signature_DFS_validation <- survfit(Surv(DFS,DFSstatus)~Sub,data=validation_data)
ggsurvplot(Signature_DFS_validation,validation_data,size=1,
           palette = "nejm",pval = TRUE,risk.table = TRUE,risk.table.col = "strata",
           break.x.by=10,break.y.by=0.2,xlim=c(0,72),axes.offset=T)
#-------------------------------------------------------------------------------------------------------------------------------------------------

##Codes for development and validation of nomogram
training_data <- read.table("training data.csv", sep = ",", head = T,row.names=1)
validation_data <- read.table("validation data.csv", sep = ",", head = T,row.names=1)
library(rms)
library(survival)
dd_training <- datadist(training_data)
options(datadist="dd_training")

#Codes for constructing and plotting nomogram
f_OS <- cph(Surv(OS,OSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data, x = T, y = T, surv = T, time.inc = 60)
surv_OS <- Survival(f_OS)

nom_OS <- nomogram(f_OS, fun = list(function(x) surv_OS(24, x),
                                    function(x) surv_OS(36, x),
                                    function(x) surv_OS(60, x)),
                   fun.at = c(seq(.1,.9,by = .1),.95),
                   funlabel = c("2-year probability of OS","3-year probability of OS","5-year probability of OS"),lp = 0)
plot(nom_OS)

f_DFS <- cph(Surv(DFS,DFSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data, x = T, y = T, surv = T, time.inc = 60)
surv_DFS <- Survival(f_DFS)

nom_DFS <- nomogram(f_DFS, fun = list(function(x) surv_DFS(24, x),
                                      function(x) surv_DFS(36, x),
                                      function(x) surv_DFS(60, x)),
                    fun.at = c(.05,seq(.1,.9,by = .1),.95),
                    funlabel = c("2-year probability of DFS","3-year probability of DFS","5-year probability of DFS"),lp = 0)
plot(nom_DFS)

#Codes for PH assumption test
library("CoxPhLb")
fit_OS <- coxphlb(Surv(OS,OSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data,
                   method = "Bootstrap")
ptest_OS <- coxphlb.phtest(fit_OS, data = training_data, spec.p = NULL)
coxphlb.phtest.plot(ptest_OS)

fit_DFS <- coxphlb(Surv(DFS,DFSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data,
                  method = "Bootstrap")
ptest_DFS <- coxphlb.phtest(fit_DFS, data = training_data, spec.p = NULL)
coxphlb.phtest.plot(ptest_DFS)

#Codes for calculating C-index in training cohort
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS), data = training_data)
rcorrcens(Surv(DFS,DFSstatus) ~ predict(f_DFS), data = training_data)

#Codes for calculating C-index in validation cohort
f_OS_val <- cph(Surv(OS,OSstatus)~predict(f_OS, newdata = validation_data), x = T, y = T, surv = T,time.inc = 60, data = validation_data)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS, newdata = validation_data), data = validation_data)

f_DFS_val <- cph(Surv(DFS,DFSstatus)~predict(f_DFS, newdata = validation_data), x = T, y = T, surv = T,time.inc = 60, data = validation_data)
rcorrcens(Surv(DFS,OSstatus) ~ predict(f_DFS, newdata = validation_data), data = validation_data)

#Codes for calibration curve in training cohort
f_OS_36 <- cph(Surv(OS,OSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data, x = T, y = T, surv = T, time.inc = 36)
f_OS_24 <- cph(Surv(OS,OSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data, x = T, y = T, surv = T, time.inc = 24)

cal_OS_training_24 <-calibrate(f_OS_24, cmethod = "KM", method = "boot", u = 24, m = 70, B = 1000)
cal_OS_training_36 <-calibrate(f_OS_36, cmethod = "KM", method = "boot", u = 36, m = 70, B = 1000)
cal_OS_training <-calibrate(f_OS, cmethod = "KM", method = "boot", u = 60, m = 70, B = 1000)
plot(cal_OS_training,xlim = c(0,1),ylim = c(0,1), col =  c("#CD2626"), errbar.col =  c("#CD2626"))
plot(cal_OS_training_36,col =  c("#458B00"), errbar.col =  c("#458B00") , add = TRUE)
plot(cal_OS_training_24,col =  c("#1874CD"), errbar.col =  c("#1874CD") , add = TRUE)

#Codes for calibration curve in validation cohort
f_DFS_36 <- cph(Surv(DFS,DFSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data, x = T, y = T, surv = T, time.inc = 36)
f_DFS_24 <- cph(Surv(DFS,DFSstatus) ~ Tstage+Nstage+Mstage+Signature, data = training_data, x = T, y = T, surv = T, time.inc = 24)

cal_DFS_training_24 <-calibrate(f_DFS_24, cmethod = "KM", method = "boot", u = 24, m = 70, B = 1000)
cal_DFS_training_36 <-calibrate(f_DFS_36, cmethod = "KM", method = "boot", u = 36, m = 70, B = 1000)
cal_DFS_training <-calibrate(f_DFS, cmethod = "KM", method = "boot", u = 60, m = 70, B = 1000)

plot(cal_DFS_training,xlim = c(0,1),ylim = c(0,1), col =  c("#CD2626"), errbar.col =  c("#CD2626"))
plot(cal_DFS_training_36,col =  c("#458B00"), errbar.col =  c("#458B00") , add = TRUE)
plot(cal_DFS_training_24,col =  c("#1874CD"), errbar.col =  c("#1874CD") , add = TRUE)

#Codes for time-independent ROC curve and AUROC comparison in training cohort
library(riskRegression)
library(survival)
Srv_OS <- Surv(training_data$OS, training_data$OSstatus)
coxmod_nomogram_OS <- coxph(Srv_OS ~ Tstage+Nstage+Mstage+Signature, data=training_data,x=TRUE)
coxmod_TNM_OS <- coxph(Srv_OS ~ Tstage+Nstage+Mstage, data=training_data,x=TRUE)
coxmod_Sig_OS <- coxph(Srv_OS ~ Signature, data=training_data,x=TRUE)
coxmod_Tstage_OS <- coxph(Srv_OS ~ Tstage, data=training_data,x=TRUE)
coxmod_Nstage_OS <- coxph(Srv_OS ~ Nstage, data=training_data,x=TRUE)
coxmod_Mstage_OS <- coxph(Srv_OS ~ Mstage, data=training_data,x=TRUE)

ROC_training_OS <- Score(list("Model1"=coxmod_nomogram_OS,"Model2"=coxmod_TNM_OS, "Model3"=coxmod_Sig_OS,
                              "Model4"=coxmod_Tstage_OS,"Model5"=coxmod_Nstage_OS,"Model6"=coxmod_Mstage_OS),
                         formula=Hist(OS,OSstatus)~1,data=training_data,
                         times=60,plots="ROC",summary="risk",contrasts=TRUE)
ROC_training_OS
plotROC(ROC_training_OS)

Srv_DFS <- Surv(training_data$DFS, training_data$DFSstatus)
coxmod_nomogram_DFS <- coxph(Srv_DFS ~ Tstage+Nstage+Mstage+Signature, data=training_data,x=TRUE)
coxmod_TNM_DFS <- coxph(Srv_DFS ~ Tstage+Nstage+Mstage, data=training_data,x=TRUE)
coxmod_Sig_DFS <- coxph(Srv_DFS ~ Signature, data=training_data,x=TRUE)
coxmod_Tstage_DFS <- coxph(Srv_DFS ~ Tstage, data=training_data,x=TRUE)
coxmod_Nstage_DFS <- coxph(Srv_DFS ~ Nstage, data=training_data,x=TRUE)
coxmod_Mstage_DFS <- coxph(Srv_DFS ~ Mstage, data=training_data,x=TRUE)

ROC_training_DFS <- Score(list("Model1"=coxmod_nomogram_DFS,"Model2"=coxmod_TNM_DFS, "Model3"=coxmod_Sig_DFS,
                               "Model4"=coxmod_Tstage_DFS,"Model5"=coxmod_Nstage_DFS,"Model6"=coxmod_Mstage_DFS),
                        formula=Hist(DFS,DFSstatus)~1,data=training_data,
                        times=60,plots="ROC",summary="risk",contrasts=TRUE)
ROC_training_DFS
plotROC(ROC_training_DFS)

#Codes for time-independent ROC curve and AUROC comparison in validation cohort
ROC_validation_OS <- Score(list("Model1"=coxmod_nomogram_OS,"Model2"=coxmod_TNM_OS, "Model3"=coxmod_Sig_OS,
                              "Model4"=coxmod_Tstage_OS,"Model5"=coxmod_Nstage_OS,"Model6"=coxmod_Mstage_OS),
                         formula=Hist(OS,OSstatus)~1,data=validation_data,
                         times=60,plots="ROC",summary="risk",contrasts=TRUE)
ROC_validation_OS
plotROC(ROC_validation_OS)

ROC_validation_DFS <- Score(list("Model1"=coxmod_nomogram_DFS,"Model2"=coxmod_TNM_DFS, "Model3"=coxmod_Sig_DFS,
                               "Model4"=coxmod_Tstage_DFS,"Model5"=coxmod_Nstage_DFS,"Model6"=coxmod_Mstage_DFS),
                          formula=Hist(DFS,DFSstatus)~1,data=validation_data,
                          times=60,plots="ROC",summary="risk",contrasts=TRUE)
ROC_validation_DFS
plotROC(ROC_validation_DFS)

#Codes for C-index comparison in training cohort
library(compareC)
f_OS_TNM <- cph(Surv(OS,OSstatus) ~ Tstage+Nstage+Mstage, data = training_data, x = T, y = T, surv = T, time.inc = 60)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS_TNM), data = training_data)
f_DFS_TNM <- cph(Surv(DFS,DFSstatus) ~ Tstage+Nstage+Mstage, data = training_data, x = T, y = T, surv = T, time.inc = 60)
rcorrcens(Surv(DFS,DFSstatus) ~ predict(f_DFS_TNM), data = training_data)

compareC(timeX = training_data$OS, statusX = training_data$OSstatus,
         scoreY = predict(f_OS), scoreZ = predict(f_OS_TNM))
compareC(timeX = training_data$DFS, statusX = training_data$DFSstatus,
         scoreY = predict(f_DFS), scoreZ = predict(f_DFS_TNM))

#Codes for C-index comparison in validation cohort
f_OS_TNM_val <- cph(Surv(OS,OSstatus)~predict(f_OS_TNM, newdata = validation_data), x = T, y = T, surv = T,time.inc = 60, data = validation_data)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS_Nstage_val, newdata = validation_data), data = validation_data)

f_DFS_TNM_val <- cph(Surv(DFS,DFSstatus)~predict(f_DFS_TNM, newdata = validation_data), x = T, y = T, surv = T,time.inc = 60, data = validation_data)
rcorrcens(Surv(DFS,DFSstatus) ~ predict(f_DFS_TNM_val, newdata = validation_data), data = validation_data)

compareC(timeX = training_data$OS, statusX = training_data$OSstatus,
         scoreY = predict(f_OS_val), scoreZ = predict(f_OS_TNM_val))
compareC(timeX = training_data$DFS, statusX = training_data$DFSstatus,
         scoreY = predict(f_DFS_val), scoreZ = predict(f_DFS_TNM_val))

#Codes for decision curve analysis in training cohort
source("stdca.R")
library(survival)
Srv_OS <- Surv(training_data$OS, training_data$OSstatus)
coxmod_nomogram_OS <- coxph(Srv_OS ~ Tstage+Nstage+Mstage+Signature, data=training_data,x=TRUE)
coxmod_TNM_OS <- coxph(Srv_OS ~ Tstage+Nstage+Mstage, data=training_data,x=TRUE)
coxmod_Sig_OS <- coxph(Srv_OS ~ Signature, data=training_data,x=TRUE)

training_data$nomogram_OS <- c(1- (summary(survfit(coxmod_nomogram_OS,newdata=training_data), times=60)$surv))
training_data$TNM_OS <- c(1- (summary(survfit(coxmod_TNM_OS,newdata=training_data), times=60)$surv))
training_data$Sig_OS <- c(1- (summary(survfit(coxmod_Sig_OS,newdata=training_data), times=60)$surv))

dcr_training_nomogram_OS <- stdca(data=training_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                                  predictors="nomogram_OS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_training_TNM_OS <- stdca(data=training_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                             predictors="TNM_OS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_training_Sig_OS <- stdca(data=training_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                             predictors="Sig_OS", cmprsk=FALSE, smooth=TRUE, xstop=0.75)
plot(dcr_training_nomogram_OS$net.benefit$none,xlim=c(0, 90), ylim=c(-0.05, 0.43), 
     xlab="Threshold probability (%)", ylab="Net benefit",type="l", lwd=1,col="black") 
lines(dcr_training_nomogram_OS$net.benefit$all, type="l", lwd=1,col="#BF3EFF") 
lines(dcr_training_nomogram_OS$net.benefit$nomogram_OS, type="l",lwd=1,col="#1874CD")
lines(dcr_training_TNM_OS$net.benefit$TNM_OS, type="l",lwd=1,col="#FF8C00")
lines(dcr_training_Sig_OS$net.benefit$Sig_OS, type="l",lwd=1,col="#CD3333")

Srv_DFS <- Surv(training_data$DFS, training_data$DFSstatus)
coxmod_nomogram_DFS <- coxph(Srv_DFS ~ Tstage+Nstage+Mstage+Signature, data=training_data,x=TRUE)
coxmod_TNM_DFS <- coxph(Srv_DFS ~ Tstage+Nstage+Mstage, data=training_data,x=TRUE)
coxmod_Sig_DFS <- coxph(Srv_DFS ~ Signature, data=training_data,x=TRUE)

training_data$nomogram_DFS <- c(1- (summary(survfit(coxmod_nomogram_DFS,newdata=training_data), times=60)$surv))
training_data$TNM_DFS <- c(1- (summary(survfit(coxmod_TNM_DFS,newdata=training_data), times=60)$surv))
training_data$Sig_DFS <- c(1- (summary(survfit(coxmod_Sig_DFS,newdata=training_data), times=60)$surv))

dcr_training_nomogram_DFS <- stdca(data=training_data, outcome="DFSstatus", ttoutcome="DFS", timepoint=60, 
                                  predictors="nomogram_DFS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_training_TNM_DFS <- stdca(data=training_data, outcome="DFSstatus", ttoutcome="DFS", timepoint=60, 
                             predictors="TNM_DFS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_training_Sig_DFS <- stdca(data=training_data, outcome="DFSstatus", ttoutcome="DFS", timepoint=60, 
                             predictors="Sig_DFS", cmprsk=FALSE, smooth=TRUE, xstop=0.75)
plot(dcr_training_nomogram_DFS$net.benefit$none,xlim=c(0, 90), ylim=c(-0.05, 0.45), 
     xlab="Threshold probability (%)", ylab="Net benefit",type="l", lwd=1,col="black") 
lines(dcr_training_nomogram_DFS$net.benefit$all, type="l", lwd=1,col="#BF3EFF") 
lines(dcr_training_nomogram_DFS$net.benefit$nomogram_DFS, type="l",lwd=1,col="#1874CD")
lines(dcr_training_TNM_DFS$net.benefit$TNM_DFS, type="l",lwd=1,col="#FF8C00")
lines(dcr_training_Sig_DFS$net.benefit$Sig_DFS, type="l",lwd=1,col="#CD3333")

#Codes for decision curve analysis in validation cohort
validation_data$nomogram_OS <- c(1- (summary(survfit(coxmod_nomogram_OS,newdata=validation_data), times=60)$surv))
validation_data$TNM_OS <- c(1- (summary(survfit(coxmod_TNM_OS,newdata=validation_data), times=60)$surv))
validation_data$Sig_OS <- c(1- (summary(survfit(coxmod_Sig_OS,newdata=validation_data), times=60)$surv))

dcr_validation_nomogram_OS <- stdca(data=validation_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                                    predictors="nomogram_OS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_validation_TNM_OS <- stdca(data=validation_data, outcome="OSSstatus", ttoutcome="OS", timepoint=60, 
                               predictors="TNM_OS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_validation_Sig_OS <- stdca(data=validation_data, outcome="OSSstatus", ttoutcome="OS", timepoint=60, 
                               predictors="Sig_OS", cmprsk=FALSE, smooth=TRUE, xstop=0.8)
plot(dcr_validation_nomogram_OS$net.benefit$none,xlim=c(0, 90), ylim=c(-0.05, 0.52), 
     xlab="Threshold probability (%)", ylab="Net benefit",type="l", lwd=1,col="black") 
lines(dcr_validation_nomogram_OS$net.benefit$all, type="l", lwd=1,col="#BF3EFF") 
lines(dcr_validation_nomogram_OS$net.benefit$nomogram_OS, type="l",lwd=1,col="#1874CD")
lines(dcr_validation_TNM_OS$net.benefit$TNM_OS, type="l",lwd=1,col="#FF8C00")
lines(dcr_validation_Sig_OS$net.benefit$Sig_OS, type="l",lwd=1,col="#CD3333")

validation_data$nomogram_DFS <- c(1- (summary(survfit(coxmod_nomogram_DFS,newdata=validation_data), times=60)$surv))
validation_data$TNM_DFS <- c(1- (summary(survfit(coxmod_TNM_DFS,newdata=validation_data), times=60)$surv))
validation_data$Sig_DFS <- c(1- (summary(survfit(coxmod_Sig_DFS,newdata=validation_data), times=60)$surv))

dcr_validation_nomogram_DFS <- stdca(data=validation_data, outcome="DFSstatus", ttoutcome="DFS", timepoint=60, 
                                    predictors="nomogram_DFS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_validation_TNM_DFS <- stdca(data=validation_data, outcome="DFSstatus", ttoutcome="DFS", timepoint=60, 
                               predictors="TNM_DFS", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_validation_Sig_DFS <- stdca(data=validation_data, outcome="DFSstatus", ttoutcome="DFS", timepoint=60, 
                               predictors="Sig_DFS", cmprsk=FALSE, smooth=TRUE, xstop=0.8)
plot(dcr_validation_nomogram_DFS$net.benefit$none,xlim=c(0, 90), ylim=c(-0.05, 0.54), 
     xlab="Threshold probability (%)", ylab="Net benefit",type="l", lwd=1,col="black") 
lines(dcr_validation_nomogram_DFS$net.benefit$all, type="l", lwd=1,col="#BF3EFF") 
lines(dcr_validation_nomogram_DFS$net.benefit$nomogram_DFS, type="l",lwd=1,col="#1874CD")
lines(dcr_validation_TNM_DFS$net.benefit$TNM_DFS, type="l",lwd=1,col="#FF8C00")
lines(dcr_validation_Sig_DFS$net.benefit$Sig_DFS, type="l",lwd=1,col="#CD3333")




