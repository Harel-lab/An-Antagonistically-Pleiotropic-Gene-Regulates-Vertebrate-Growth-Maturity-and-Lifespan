###Code for BaSTA Analysis and Weibull regression analyses 

##load in required libraries
library(devtools)
#install_github("fercol/BaSTA/pkg")
library(BaSTA) #we use the multibasta function to run multiple BaSTA models on the same data set, perform model comparison and visualise the results of the multiple runs
library(snowfall) ##this allows us to run multiple BaSTA simulations in parallel
library(dplyr)
library(tidyverse)
library(ggthemes)
library(ggplot2)
library(patchwork)
library(survival)
library(survminer)
library(Rcpp)
library(ggsurvfit)

##set working directory 
setwd("")

##load in survival data for log rank and coxph analysis and get in correct format
NewestSurv2<-read.table("VGLL3-New_edit_cox.txt", header=T)
head(NewestSurv2)
str(NewestSurv2)
NewestSurv2$ID<-as.factor(NewestSurv2$ID)
NewestSurv2$Death<-as.numeric(NewestSurv2$Death)
NewestSurv2$Treatment<-as.factor(NewestSurv2$Treatment)
NewestSurv2$Gen.background<-as.factor(NewestSurv2$Gen.background)
NewestSurv2$Sex<-as.factor(NewestSurv2$Sex)
str(NewestSurv2)

##run log rank test for pairwise comparisons of treatments
##first subset data by sex
Females<-subset(NewestSurv2,Sex=="female")
Males<-subset(NewestSurv2,Sex=="male")

##compare WT and mutant for females then for males
logrankFemales<-survdiff(Surv(Death,Status)~Gen.background,data=Females)
logrankFemales 
logrankMales<-survdiff(Surv(Death,Status)~Gen.background,data=Males)
logrankMales

##Now run parametric survival analyses with the Weibull regression model
Weibullmodel<-survreg(Surv(Death,Status)~Sex*Gen.background,data=NewestSurv2,dist='weibull')
summary(Weibullmodel)

## Install SurvRegCensCov package for more in depth Weibull analysis
install.packages("SurvRegCensCov")
library(SurvRegCensCov)
ConvertWeibull(Weibullmodel,conf.level = 0.95) ##to get more relevant parameters from the Weibull model


##We can alternatively fit the Weibull regression model like this (combining survreg and ConvertWeibull)
WeibullmodelB<-WeibullReg(Surv(Death,Status)~Sex*Gen.background,data=NewestSurv2,conf.level=0.95)
WeibullmodelB ##this gives all outputs from before

##Diagnostic plots to test adequacy of Weibull model- stratified Kaplan Meier curves - want linear and parallel lines for an adequate model
WeibullDiag(Surv(Death,Status)~Sex,data=NewestSurv2) ##generally parallel and linear
WeibullDiag(Surv(Death,Status)~Gen.background,data=NewestSurv2) ##generally parallel and linear


##the interaction is non-significant so can remove from the model and run a separate model on each sex separately
Femalesurv<-subset(NewestSurv2,Sex=="female")
Malesurv<-subset(NewestSurv2,Sex=="male")
Weibullfemales<-WeibullReg(Surv(Death,Status)~Gen.background,data=Femalesurv,conf.level=0.95)
Weibullfemales
WeibullDiag(Surv(Death,Status)~Gen.background,data=Femalesurv)
Weibullmales<-WeibullReg(Surv(Death,Status)~Gen.background,data=Malesurv,conf.level=0.95)
Weibullmales
WeibullDiag(Surv(Death,Status)~Gen.background,data=Malesurv)


##Now subset by genotype
WTsurv<-subset(NewestSurv2,Gen.background=="WT")
mutantsurv<-subset(NewestSurv2,Gen.background="mutant")
WeibullWT<-WeibullReg(Surv(Death,Status)~Sex,data=WTsurv,conf.level=0.95)
WeibullWT
Weibullmutant<-WeibullReg(Surv(Death,Status)~Sex,data=mutantsurv,conf.level=0.95)
Weibullmutant

###Alternatively, model Weibull regression via eha package (Event History Analysis)
install.packages("eha")
library(eha)
Weibullaltmodel<-weibreg(Surv(Death,Status)~Sex*Gen.background, data 
= NewestSurv2) 
Weibullaltmodel

##Run separate models for each sex, as the interaction is non-significant, using the strata command
Weibullaltmodelfemale<-weibreg(Surv(Death,Status)~Gen.background,data=Femalesurv) 
Weibullaltmodelfemale
Weibullaltmodelmale<-weibreg(Surv(Death,Status)~Gen.background,data=Malesurv) 
Weibullaltmodelmale
fit<-weibreg(Surv(Death,Status)~strata(Gen.background), shape=~strata(Gen.background),data=Femalesurv)
summary(fit)
fit2<-weibreg(Surv(Death,Status)~strata(Gen.background), shape=~strata(Gen.background),data=Malesurv)
summary(fit2)
fit3<-weibreg(Surv(Death,Status)~strata(Treatment), shape=~strata(Treatment),data=NewestSurv2)
summary(fit3)
plot(fit3, fn="haz")
ggplot(fit3)
plot(fit3)
fit11<-weibreg(Surv(Death,Status)~Gen.background, shape=~Gen.background,data=Femalesurv)
summary(fit11)
fit12<-weibreg(Surv(Death,Status)~Gen.background, shape=~Gen.background,data=Malesurv)
summary(fit12)
waldtest(fit)

##plot
fita<-phreg(Surv(Death,Status)~strata(Gen.background), shape=~strata(Gen.background),dist="weibull", data=Femalesurv)
fita
plot(fita,fn = "haz")
plot(fita, fn="haz",main="Female- Weibull hazard function",xlab="days",col=c("purple"),lwd=2)
legend("topleft", legend=c("vgll3mutant","WT"),col=c("purple","black"),lty=c(1,2),lwd=c(2,1))
fitaaa<-phreg(Surv(Death,Status)~strata(Treatment),shape=~strata(Treatment),dist="weibull",data=NewestSurv2)
fitaaa
plot(fitaaa, fn="haz",main="Weibull hazard function",xlab="days",ylab="hazard",col=c("purple","lightpink1","blue","red"),lwd=2,lty=1,printLegend = FALSE)
legend("topleft", legend=c("vgll3 female","vgll3 male", "WT female", "WT male"),col=c("purple","lightpink1","blue","red"),lty=1,lwd=2)
plot(fitaaa)
plot(fit2, fn="haz")
plot(fit3, fn="haz", new.data = strata(NewestSurv2$Treatment))
legend("topleft",legend = NewestSurv2$Treatment,lwd=3)
plot(fitaaa,main="Weibull hazard function",xlab="days",ylab="hazard",col=c("purple","lightpink1","blue","red"),lwd=2,lty=1,printLegend=FALSE)

##try to run a Wald test on the shape parameters between genotypes, separately for the female and male models
library(aod)
library(sandwich)
library(lmtest)
wald.test(b = coef(fit), Sigma = vcov(fit), Terms = c(2,4))



#####
#############
##now preparing the data for BaSTA analysis 

##load in survival data with the columns: ID, birth, death (will add treatment info afterwards)...we have omitted matricide and censored data
Surv3 <- read.table("Short_VGLL3-New_edit.txt", header=T)
head(Surv3)
str(Surv3)
summary(Surv3)


##turning it into a dataframe with the number of columns equal to the maximum age at death plus 3 for the ID, birth and death columns, and the number of rows matching that in the original dataframe
bastasurv3 <- as.data.frame(matrix(ncol = max(Surv3$Death) + 3, nrow = nrow(Surv3)))

##loop to generate 0s and 1s ...for each row, have a 0, repeat 1 for as many times as the value of death day minus 2 (to not include the first day or last day), then repeat 0 for the number of times between the maximum death day and the death day for that individual; then make the final dataframe is the original dataframe bound together with these columns of 0's and 1's 
for (i in 1:nrow(Surv3)){
  p<-c(0,rep(1, length = Surv3$Death[i] - 2),0,rep(0, (max(Surv3$Death)-Surv3$Death[i])))
  bastasurv3[i,]<-as.data.frame(rbind(c(Surv3[i,],p)))}

##create df
df3<-as.data.frame(lapply(bastasurv3, unlist))
colnames(df3) = c("id", "birth", "death",1:max(Surv3$Death)) ##give column names
str(df3)
newdata3<-DataCheck(df3, studyStart = 1, studyEnd = max(Surv3$Death)) ##no problems were detected with the data

##now export df into Excel so that we can add treatment info
library(writexl)
write_xlsx(df3, 'file-path-to-data')

##now load in the final dataframe
dfnew3<-read.table("df3data_full.txt", header=T)
head(dfnew3)
str(dfnew3)
newdata4<-DataCheck(dfnew3, studyStart = 1, studyEnd = max(dfnew3$death)) ##no problems were detected with the data

##Run multi version of BaSTA on the data, using the best parameters from last time's multimodel comparison
outnew3 <- multibasta(dfnew3, studyStart = 1, studyEnd = max(dfnew3$death), models = c("GO", "WE", "LO", "EX"), shapes = c("simple", "Makeham", "bathtub"), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 600000, burnin = 60001, thinning = 600)

##Print results
summary(outnew3, digits = 3) ##Weibull Makeham and Weibull Bathtub did not converge, so try re-running with new parameters

##load gt package for displaying tables
library(gt)
gt(outnew3$DICs) ##gives the DICs in a fancy table

outnew3

##Re-run multi version of BaSTA on the data, using updated parameters
outnew4 <- multibasta(dfnew3, studyStart = 1, studyEnd = max(dfnew3$death), models = c("GO", "WE", "LO", "EX"), shapes = c("simple", "Makeham", "bathtub"), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 1000000, burnin = 100001, thinning = 1000)

##Print results
summary(outnew4, digits = 3) ##Weibull Makeham and Weibull Bathtub did not converge, so try re-running with new parameters

##load gt package for displaying tables
library(gt)
gt(outnew4$DICs) ##gives the DICs in a fancy table

outnew4

##Final re-run multi version of BaSTA on the data, using updated parameters
outnew5 <- multibasta(dfnew3, studyStart = 1, studyEnd = max(dfnew3$death), models = c("GO", "WE", "LO", "EX"), shapes = c("simple", "Makeham", "bathtub"), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 1500000, burnin = 150001, thinning = 1500)

##Print results
summary(outnew5, digits = 3)

##load gt package for displaying tables
library(gt)
gt(outnew5$DICs) ##gives the DICs in a fancy table

outnew5

## Run single version of BaSTA with the best model/shape/covariate structure (this was WE simple in our case) , can use the same parameters as above since it converged and serial autocorrelation was <5%
outbest <- basta(dfnew, studyStart = 1, studyEnd = max(dfnew$death), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 150000, burnin = 15001, thinning = 150, shape = "simple", model = "WE") 
save(outbest, file="model1.Rdata") ##this saves the model into the current working directory (survivaldata folder in one drive)
##to load the model, use load() command
## Print results: 
summary(outbest, digits = 3) ##want all serial autocorrelation values below 5%, not yet the case

##try another run of the model
outbest2 <- basta(dfnew, studyStart = 1, studyEnd = max(dfnew$death), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 600000, burnin = 60001, thinning = 600, shape = "simple", model = "WE") 
save(outbest2, file="model1.Rdata")

summary(outbest2, digits = 3) ##all serial autocorrelation less than 5% now

##try one more run with more iterations
outbest5 <- basta(dfnew, studyStart = 1, studyEnd = max(dfnew$death), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 1000000, burnin = 100001, thinning = 1000, shape = "simple", model = "WE") ##ran fine and survival parameters converged appropriately
save(outbest5, file="model5.Rdata")## this is the best one to use

##Try running a single Weibull Makeham model
outbest3 <- basta(dfnew, studyStart = 1, studyEnd = max(dfnew$death), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 600000, burnin = 60001, thinning = 600, shape = "Makeham", model = "WE") ##this didn't converge for some parameters unfortunately

##try a different run
outbest4 <- basta(dfnew, studyStart = 1, studyEnd = max(dfnew$death), updateJumps = TRUE, parallel = TRUE, ncpus = 8, nsim = 4, niter = 1000000, burnin = 100001, thinning = 1000, shape = "Makeham", model = "WE")
##this also did not converge for some parameters


##Plot survival and mortality curves: 
#export as PDF
load("model5.Rdata")
plot(outbest5, plot.trace = FALSE, noCI = FALSE, fancy=T, col = c("red","pink2","blue","#984EA3"), cex=2,names.legend= c("WT male", "vgll3 mutant male", "WT female", "vgll3 mutant female") )
plot(outbest5, densities = TRUE)


##plot KLDC from best model add various filters depending on choice. 
KLDC <- read.table("",header=T) 
KLDC
head(KLDC)
str(KLDC)
KLDC$Treatment_Comparisons<-as.factor(KLDC$Treatment_Comparisons)
KLDC$Parameter<-as.factor(KLDC$Parameter)
str(KLDC)
KLDC$Treatment_Comparisons
summary(KLDC)

library(cowplot)
library(stringr)
pl1<- KLDC %>% ggplot(aes(x = Treatment_Comparisons, y = KLDC_values, colour = Parameter)) + 
  geom_point(size = 5) + 
  geom_hline(yintercept = 0.85, linetype=2) + 
  scale_x_discrete(limits=c("WT_female-WT_male","mutant_female-mutant_male","mutant_female-WT_female","mutant_male-WT_male"),labels=c("WT_female-WT_male"="WT_female vs WT_male","mutant_female-mutant_male"="mutant_female vs mutant_male","mutant_female-WT_female"="mutant_female vs WT_female","mutant_male-WT_male"="mutant_male vs WT_male"))+
  xlab("Treatment Comparisons")+
  ylab("KLDC values")+
  labs(colour="Parameter")+
  theme(axis.text = element_text(size=20),axis.title=element_text(size=22))+
  theme_cowplot()
pl1

#####
#####
##plot b1 estimates and CI
b1estimates <- read.table("",header=T) 
b1estimates
str(b1estimates)
b1estimates$Treatment<-as.factor(b1estimates$Treatment)
str(b1estimates)
summary(b1estimates)

geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))

library(cowplot)
library(stringr)
pl2<- ggplot(b1estimates, aes(x = Treatment, y = Estimate)) + 
  geom_point(size = 5) + 
  geom_errorbar(aes(ymin= Lower, ymax = Upper))+
  scale_x_discrete(limits=c("WT_male","mutant_male","WT_female","mutant_female"),labels=c("WT_male"="WT male","mutant_male"="vgll3 mutant male","WT_female"="WT female","mutant_female"="vgll3 mutant female"))+
  xlab("Treatment")+
  ylim=c(0,0.05)+
  ylab("b1 estimates (with 95% CI)")+
  labs(colour= c("red","pink2","blue","#984EA3"))+
  theme(axis.text = element_text(size=20),axis.title=element_text(size=22))
pl2



#####

