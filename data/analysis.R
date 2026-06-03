# PCORI
rm(list=ls())

# Function to calculate influence function (for betas) for orm:

getBetaInfluenceFunctions <- function(fit, getC1 = FALSE) {
  # Quick checks
  if(!(all(class(fit) %in% c("orm", "rms")))) stop("Fitted object must be of class 'orm' from the rms package")
  if(is.null(fit$mscore)) stop("orm must be fit with 'mscore = TRUE'")
  if(strsplit(as.character(packageVersion("rms")), "[.]")[[1]][1] < 8) stop("Requires an orm fitted object from rms version 8.0.0 or higher")
  
  Na <- fit$non.slopes # Number of intercepts
  Np <- length(fit$coefficients) - Na # Number of slope parameters
  
  a  <- Matrix::bandSparse(Na, k = c(0, 1), diagonals = fit$info.matrix$a, symmetric = TRUE)
  ab <- fit$info.matrix$ab                        
  b  <- fit$info.matrix$b
  B1 <- fit$mscore[, 1:Na]
  B2 <- fit$mscore[, (Na + 1):(Na + Np)]
  
  a_inv_ab <- Matrix::solve(a, ab)
  C2       <- solve(b - t(ab) %*% a_inv_ab, t(B2 - B1 %*% a_inv_ab))
  
  if(getC1) {
    C1 <- -t(Matrix::solve(a, t(B1) - ab %*% C2))
  } else {
    C1 <- NULL
  }
  return(list(C1 = C1, C2 = -t(C2)))
}

setwd("~/Dropbox/chun-stuff/eric-kawaguchi/orm-pcori")
alldata <- read.csv("synthetic-mother-child-data.csv")

library(survival)
library(survey)
library(rms)


################################################################
#### Now doing proposed analysis with synthetic data
################################################################

# maternal weight gain ~ maternal BMI + age + race + ethnicity + depression + insurance + smoking; 
## exclude those that are not singleton

# Y = Estimated maternal weight change per week
# X = Maternal BMI
#     Age
#     Race
#     Ethnicity
#     Depression
#     Insurance
#     Singleton (remove)
#     Smoking
#     Number prior births
#     EGA

# orm

### Limiting the dataset to only singletons
d<-alldata[alldata$singleton.ph1==1 & (is.na(alldata$singleton) | alldata$singleton==1),]


### Fitting the naive model from which we will extract naive influence functions
mod.orm.n <- orm(estWtChangePerWk.ph1 ~ est_bmi_preg_mother.ph1+
                   mat_age_delivery.ph1+mat4race.ph1+ hispanic.ph1 + 
                   tobacPregr.ph1+depressionr.ph1+insurancer.ph1,
                 data = d,
                 family = "logistic", mscore = TRUE)
naiveInfl.orm <- getBetaInfluenceFunctions(mod.orm.n)$C2
d$nInfl.bmi<-naiveInfl.orm[,1]
d$nInfl.age<-naiveInfl.orm[,2]
d$nInfl.race1<-naiveInfl.orm[,3]
d$nInfl.race2<-naiveInfl.orm[,4]
d$nInfl.race3<-naiveInfl.orm[,5]
d$nInfl.hisp<-naiveInfl.orm[,6]
d$nInfl.tobac<-naiveInfl.orm[,7]
d$nInfl.depr<-naiveInfl.orm[,8]
d$nInfl.insur<-naiveInfl.orm[,9]

### Creating the sampling design object and calibrating weights based on the naive influence function
designO <- twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~R==1,
                    data=d,method="approx")
Ninflcal.orm <-calibrate(designO,
                         formula=~nInfl.bmi+nInfl.age+nInfl.race1+nInfl.race2+nInfl.race3+
                           nInfl.hisp+nInfl.tobac+nInfl.depr+nInfl.insur+wave4Strata,
                         phase=2,calfun="raking") 

design <- Ninflcal.orm
if(any(weights(design) < 0)) stop("weights must be non-negative")
sdat <- model.frame(design)
sdat$.survey.prob.weights <- (1 / design$prob) / mean(1 / design$prob)

### Fitting CPM with calibrated weights (generalized raking)
mod.orm.GR <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurancer+
                   egaWk, 
              family = "logistic",
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights)	

orig.var.GR <- infoMxop(mod.orm.GR$info.matrix, i = "x")
inflFun.GR  <- getBetaInfluenceFunctions(mod.orm.GR, getC1 = FALSE)
twophase.var.GR <- twophasevar(inflFun.GR$C2, Ninflcal.orm)
sqrt(diag(orig.var.GR))
sqrt(diag(twophase.var.GR))


### Now IPW analysis for comparison
if(any(weights(designO) < 0)) stop("weights must be non-negative")
sdatO <- model.frame(designO)
sdatO$.survey.prob.weights <- (1 / designO$prob) / mean(1 / designO$prob)

### Fitting CPM with inverse probability weights
mod.orm.IPW <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurancer+
                   egaWk, 
              family = "logistic",
              data = sdatO,
              mscore = TRUE,
              weights = .survey.prob.weights)	

orig.var.IPW <- infoMxop(mod.orm.IPW$info.matrix, i = "x")
inflFun.IPW  <- getBetaInfluenceFunctions(mod.orm.IPW, getC1 = FALSE)
twophase.var.IPW <- twophasevar(inflFun.IPW$C2, designO)          
sqrt(diag(orig.var.IPW))
sqrt(diag(twophase.var.IPW))



### A bootstrap procedure to obtain confidence intervals. This takes a little while
###  to run, so it has been commented out.
# ## Bootstrapping by sampling with replacement from each stratum
# 
# set.seed(123)
# nboot<-1000
# coef.GR<-coef.IPW<-matrix(NA, nrow=nboot, ncol=10)
# strat.names<-unique(d$wave4Strata)
# 
# for (i in 1:nboot){
#   boot.ids1<-boot.ids0<-NULL
#   boot.ids1<-sample(which(d$R==1 & d$wave4Strata==strat.names[1]), sum(d$R==1 & d$wave4Strata==strat.names[1]), replace=TRUE)
#   boot.ids0<-sample(which(d$R==0 & d$wave4Strata==strat.names[1]), sum(d$R==0 & d$wave4Strata==strat.names[1]), replace=TRUE)
#   for (j in 2:length(strat.names)) {
#     boot.ids1<-c(boot.ids1,sample(which(d$R==1 & d$wave4Strata==strat.names[j]), sum(d$R==1 & d$wave4Strata==strat.names[j]), replace=TRUE))
#     boot.ids0<-c(boot.ids0,sample(which(d$R==0 & d$wave4Strata==strat.names[j]), sum(d$R==0 & d$wave4Strata==strat.names[j]), replace=TRUE))
#   }
#   d.boot<-d[c(boot.ids1,boot.ids0),]
#   
#   mod.orm.n <- orm(estWtChangePerWk.ph1 ~ est_bmi_preg_mother.ph1+
#                    mat_age_delivery.ph1+mat4race.ph1+ hispanic.ph1 + 
#                    tobacPregr.ph1+depressionr.ph1+insurancer.ph1,
#                  data = d.boot,
#                  family = "logistic", mscore = TRUE)
#   naiveInfl.orm <- getBetaInfluenceFunctions(mod.orm.n)$C2
#   d.boot$nInfl.bmi<-naiveInfl.orm[,1]
#   d.boot$nInfl.age<-naiveInfl.orm[,2]
#   d.boot$nInfl.race1<-naiveInfl.orm[,3]
#   d.boot$nInfl.race2<-naiveInfl.orm[,4]
#   d.boot$nInfl.race3<-naiveInfl.orm[,5]
#   d.boot$nInfl.hisp<-naiveInfl.orm[,6]
#   d.boot$nInfl.tobac<-naiveInfl.orm[,7]
#   d.boot$nInfl.depr<-naiveInfl.orm[,8]
#   d.boot$nInfl.insur<-naiveInfl.orm[,9]
# 
#             ## Errors if one of the strata are empty in the bootstrap sample. I skip these
#             ##  replications when it happens.
#   if (min(table(d.boot$wave4Strata,d.boot$R)) >0){
#   designO <- twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~R==1,
#                       data=d.boot,method="approx")
# 
#   Ninflcal.orm <-calibrate(designO,
#                            formula=~nInfl.bmi+nInfl.age+nInfl.race1+nInfl.race2+nInfl.race3+
#                              nInfl.hisp+nInfl.tobac+nInfl.depr+nInfl.insur+wave4Strata,
#                            phase=2,calfun="raking") 
# 
#   design <- Ninflcal.orm
#   if(any(weights(design) < 0)) stop("weights must be non-negative")
#   sdat <- model.frame(design)
#   sdat$.survey.prob.weights <- (1 / design$prob) / mean(1 / design$prob)
# 
#   mod.orm.GR <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
#                      mat_age_delivery+mat4race+ hispanic + 
#                      tobacPregr+depressionr+insurancer+
#                    egaWk, 
#               family = "logistic",
#               data = sdat,
#               weights = .survey.prob.weights)	
# 
#   if(any(weights(designO) < 0)) stop("weights must be non-negative")
#   sdatO <- model.frame(designO)
#   sdatO$.survey.prob.weights <- (1 / designO$prob) / mean(1 / designO$prob)
# 
#   mod.orm.IPW <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
#                    mat_age_delivery+mat4race+ hispanic + 
#                    tobacPregr+depressionr+insurancer+
#                    egaWk, 
#               family = "logistic",
#               data = sdatO,
#               weights = .survey.prob.weights)	
# 
#   coef.GR[i,]<-tail(mod.orm.GR$coefficients,10)
#   coef.IPW[i,]<-tail(mod.orm.IPW$coefficients,10)
#   } 
# }
# 
# sum(is.na(coef.GR[,1]))
# 
# apply(coef.GR,2,sd,na.rm=TRUE)
# apply(coef.IPW,2,sd,na.rm=TRUE)
# 
# ses.GR
# ses.IPW
# 
# ## The standard errors using the sampling 
# ##  with replacement from strata-by-phase2 status were fairly similar to those that only sampled with 
# ##  replacement by phase2 status, but there were no empty strata (which caused some problems when only 
# ##  sampling with replacement based on phase2 status).



################################################################
################################################################
### Now doing proposed analysis with the real data. The real data
###  cannot be provided, but an interested investigator could apply 
###  this code (likely with a few minor tweaks) to the synthetic 
###  data provided above.
################################################################
################################################################



rm(list=ls())
alldata<-readRDS(file = "analysisData20240523-bes.rds")

library(rms)
library(survey)

getBetaInfluenceFunctions <- function(fit, getC1 = FALSE) {
  # Quick checks
  if(!(all(class(fit) %in% c("orm", "rms")))) stop("Fitted object must be of class 'orm' from the rms package")
  if(is.null(fit$mscore)) stop("orm must be fit with 'mscore = TRUE'")
  if(strsplit(as.character(packageVersion("rms")), "[.]")[[1]][1] < 8) stop("Requires an orm fitted object from rms version 8.0.0 or higher")
  
  Na <- fit$non.slopes # Number of intercepts
  Np <- length(fit$coefficients) - Na # Number of slope parameters
  
  a  <- Matrix::bandSparse(Na, k = c(0, 1), diagonals = fit$info.matrix$a, symmetric = TRUE)
  ab <- fit$info.matrix$ab                        
  b  <- fit$info.matrix$b
  B1 <- fit$mscore[, 1:Na]
  B2 <- fit$mscore[, (Na + 1):(Na + Np)]
  
  a_inv_ab <- Matrix::solve(a, ab)
  C2       <- solve(b - t(ab) %*% a_inv_ab, t(B2 - B1 %*% a_inv_ab))
  
  if(getC1) {
    C1 <- -t(Matrix::solve(a, t(B1) - ab %*% C2))
  } else {
    C1 <- NULL
  }
  return(list(C1 = C1, C2 = -t(C2)))
}


d<-alldata[alldata$singleton.ph1==1 & (is.na(alldata$singleton) | alldata$singleton==1),]
d$R<-d$in.Ophase2
d$mat4race.ph1<-factor(d$mat4race.ph1, levels=c("Asian","Black","White","Other"))
d$mat4race<-factor(d$mat4race, levels=c("Asian","Black","White","Other"))
d$mat4race.ph1 = relevel(as.factor(d$mat4race.ph1),ref="White")
d$mat4race = relevel(as.factor(d$mat4race),ref="White")


d$insurance.priv.ph1 = with(d, ifelse(insurancer.ph1=="Private",1,0))
d$insurance.priv = with(d, ifelse(insurancer=="Private",1,0))



hist(d$estWtChangePerWk[d$R==1])
summary(d$estWtChangePerWk[d$R==1])
hist(d$estWtChangePerWk.ph1)
summary(d$estWtChangePerWk.ph1)

pdf("weight-change-boxplot.pdf", width=5)
boxplot(d$estWtChangePerWk.ph1, d$estWtChangePerWk[d$R==1],
        ylab="Weight Change during Pregnancy (kg/week)", xlab="Cohort")
mtext("Phase-1", side=1, at=1,line=1)
mtext("Phase-2", side=1, at=2,line=1)
dev.off()


junk1<-hist(d$estWtChangePerWk.ph1, nclass=100)
junk2<-hist(d$estWtChangePerWk,nclass=100)

pdf("weight-change-correlation.pdf",height=5,width=5)
mar=c(5,2,1,1)
min.change<-with(d,min(c(estWtChangePerWk.ph1,estWtChangePerWk),na.rm=TRUE))
max.change<-with(d,max(c(estWtChangePerWk.ph1,estWtChangePerWk),na.rm=TRUE))
plot(c(min.change,max.change),c(min.change,max.change),type="n",
     xlab="Weight Change (kg/wk) (Phase 1)",ylab="Weight Change (kg/wk) (Phase 2)")
with(d, points(estWtChangePerWk.ph1,estWtChangePerWk))
abline(0,1,col=gray(.8))
for(i in 1:length(junk1$density)){
  if (junk1$density[i]>0){
    lines(rep(junk1$mids[i],2), c(-0.8, -0.8+junk1$density[i]*.1) ,lwd=2)
  }
}
for(i in 1:length(junk2$density)){
  if (junk2$density[i]>0){
    lines(c(-0.8, -0.8+junk2$density[i]*.1), rep(junk2$mids[i],2) ,lwd=2)
  }
}
dev.off()

with(d, cor(estWtChangePerWk.ph1[R==1],estWtChangePerWk[R==1]))

round(quantile(d$estWtChangePerWk.ph1, c(.5,.25,.75)),2)
round(quantile(d$estWtChangePerWk[d$R==1], c(.5,.25,.75)),2)
mean(d$estWtChangePerWk.ph1[d$R==1]!=d$estWtChangePerWk[d$R==1])
round(summary(d$estWtChangePerWk.ph1[d$R==1]-d$estWtChangePerWk[d$R==1])[c(3,1,6)],2)

round(quantile(d$est_bmi_preg_mother.ph1, c(.5,.25,.75)),1)
round(quantile(d$est_bmi_preg_mother[d$R==1], c(.5,.25,.75)),1)
mean(d$est_bmi_preg_mother.ph1[d$R==1]!=d$est_bmi_preg_mother[d$R==1])
round(summary(d$est_bmi_preg_mother.ph1[d$R==1]-d$est_bmi_preg_mother[d$R==1])[c(3,1,6)],2)
with(d, cor(est_bmi_preg_mother.ph1[R==1],est_bmi_preg_mother[R==1]))
with(d, plot(est_bmi_preg_mother.ph1[R==1],est_bmi_preg_mother[R==1]))


table(d$mat4race.ph1)
round(table(d$mat4race.ph1)/sum(table(d$mat4race.ph1)),2)
table(d$mat4race)
round(table(d$mat4race)/sum(table(d$mat4race)),2)
round(mean(d$mat4race[d$R==1]!=d$mat4race.ph1[d$R==1]),2)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1=="Asian"]=="Asian"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1!="Asian"]!="Asian"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1=="Black"]=="Black"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1!="Black"]!="Black"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1=="White"]=="White"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1!="White"]!="White"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1=="Other"]=="Other"),3)
round(mean(d$mat4race[d$R==1 & d$mat4race.ph1!="Other"]!="Other"),3)

sum(d$hispanic.ph1)
round(mean(d$hispanic.ph1),2)
sum(d$hispanic[d$R==1])
round(mean(d$hispanic[d$R==1]),2)
round(mean(d$hispanic.ph1[d$R==1]!=d$hispanic[d$R==1]),2)
round(mean(d$hispanic[d$R==1 & d$hispanic.ph1==1]==1),3)
round(mean(d$hispanic[d$R==1 & d$hispanic.ph1!=1]!=1),3)

sum(d$tobacPregr.ph1)
round(mean(d$tobacPregr.ph1),2)
sum(d$tobacPregr[d$R==1])
round(mean(d$tobacPregr[d$R==1]),2)
round(mean(d$tobacPregr.ph1[d$R==1]!=d$tobacPregr[d$R==1]),2)
round(mean(d$tobacPregr[d$R==1 & d$tobacPregr.ph1==1]==1),3)
round(mean(d$tobacPregr[d$R==1 & d$tobacPregr.ph1!=1]!=1),3)

sum(d$depressionr.ph1)
round(mean(d$depressionr.ph1),2)
sum(d$depressionr[d$R==1])
round(mean(d$depressionr[d$R==1]),2)
round(mean(d$depressionr.ph1[d$R==1]!=d$depressionr[d$R==1]),2)
round(mean(d$depressionr[d$R==1 & d$depressionr.ph1==1]==1),3)
round(mean(d$depressionr[d$R==1 & d$depressionr.ph1!=1]!=1),3)

sum(d$insurance.priv.ph1)
round(mean(d$insurance.priv.ph1),2)
sum(d$insurance.priv[d$R==1])
round(mean(d$insurance.priv[d$R==1]),2)
round(mean(d$insurance.priv.ph1[d$R==1]!=d$insurance.priv[d$R==1]),2)
round(mean(d$insurance.priv[d$R==1 & d$insurance.priv.ph1==1]==1),3)
round(mean(d$insurance.priv[d$R==1 & d$insurance.priv.ph1!=1]!=1),3)

sum(d$insurancer.ph1=="Public/Other")
round(mean(d$insurancer.ph1=="Public/Other"),2)
sum(d$insurancer[d$R==1]=="Public/Other")
round(mean(d$insurancer[d$R==1]=="Public/Other"),2)
round(mean(d$insurancer.ph1[d$R==1]!=d$insurancer[d$R==1]),2)
round(mean(d$insurancer[d$R==1 & d$insurancer.ph1=="Public/Other"]=="Public/Other"),3)
round(mean(d$insurancer[d$R==1 & d$insurancer.ph1!="Public/Other"]!="Public/Other"),3)

round(quantile(d$egaWk[d$R==1], c(.5,.25,.75)),1)


dd  <- with(d, datadist(est_bmi_preg_mother, mat_age_delivery, mat4race, hispanic, tobacPregr,
                depressionr, insurance.priv, insurancer, egaWk,
                est_bmi_preg_mother.ph1, mat_age_delivery.ph1, mat4race.ph1, hispanic.ph1, tobacPregr.ph1,
                depressionr.ph1, insurance.priv.ph1, insurancer.ph1)) 
options(datadist='dd')

mod.orm.n <- orm(estWtChangePerWk.ph1 ~ est_bmi_preg_mother.ph1+
                   mat_age_delivery.ph1+mat4race.ph1+ hispanic.ph1 + 
                   tobacPregr.ph1+depressionr.ph1+insurance.priv.ph1,
                 data = d,
                 family = "logistic", mscore = TRUE)
naiveInfl.orm <- getBetaInfluenceFunctions(mod.orm.n)$C2
d$nInfl.bmi<-naiveInfl.orm[,1]
d$nInfl.age<-naiveInfl.orm[,2]
d$nInfl.race1<-naiveInfl.orm[,3]
d$nInfl.race2<-naiveInfl.orm[,4]
d$nInfl.race3<-naiveInfl.orm[,5]
d$nInfl.hisp<-naiveInfl.orm[,6]
d$nInfl.tobac<-naiveInfl.orm[,7]
d$nInfl.depr<-naiveInfl.orm[,8]
d$nInfl.insur<-naiveInfl.orm[,9]

## I am simplifying things here by only including those who were sampled as part of the obesity study.
##  Others were sampled and validated as part of the asthma study, and we also used their data in our 
##  original Biometrics paper. But those individuals had a separate sampling frame, so we had to do some
##  additional things to properly include those individuals. We can do it again in this study if needed,
##  but I currently prefer just using the phase 2 subsample that was sampled for the obesity study.

designO <- twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~R==1,
                    data=d,method="approx")

Ninflcal.orm <-calibrate(designO,
                         formula=~nInfl.bmi+nInfl.age+nInfl.race1+nInfl.race2+nInfl.race3+
                           nInfl.hisp+nInfl.tobac+nInfl.depr+nInfl.insur+wave4Strata,
                         phase=2,calfun="raking") 

design <- Ninflcal.orm
if(any(weights(design) < 0)) stop("weights must be non-negative")
sdat <- model.frame(design)
sdat$.survey.prob.weights <- (1 / design$prob) / mean(1 / design$prob)

mod.orm.GR <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurance.priv+
                   egaWk, 
              family = "logistic",
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights)	

# Calculate robust SE by mimicking what is done in svycoxph
orig.var.GR <- infoMxop(mod.orm.GR$info.matrix, i = "x")
inflFun.GR  <- getBetaInfluenceFunctions(mod.orm.GR, getC1 = FALSE)
twophase.var.GR <- twophasevar(inflFun.GR$C2, Ninflcal.orm)
sqrt(diag(orig.var.GR))
sqrt(diag(twophase.var.GR))


if(any(weights(designO) < 0)) stop("weights must be non-negative")
sdatO <- model.frame(designO)
sdatO$.survey.prob.weights <- (1 / designO$prob) / mean(1 / designO$prob)

mod.orm.IPW <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurance.priv+
                   egaWk, 
              family = "logistic",
              data = sdatO,
              mscore = TRUE,
              weights = .survey.prob.weights)	

# Calculate robust SE by mimicking what is done in svycoxph
orig.var.IPW <- infoMxop(mod.orm.IPW$info.matrix, i = "x")
inflFun.IPW  <- getBetaInfluenceFunctions(mod.orm.IPW, getC1 = FALSE)
twophase.var.IPW <- twophasevar(inflFun.IPW$C2, designO)          ### I think this is designO, but I'm just copying
sqrt(diag(orig.var.IPW))
sqrt(diag(twophase.var.IPW))


summary(diag(twophase.var.GR)/diag(twophase.var.IPW))

## Estimates

ests.orm.GR<-cbind(tail(mod.orm.GR$coefficients,10),
      tail(mod.orm.GR$coefficients,10)-1.96*sqrt(diag(twophase.var.GR)),
      tail(mod.orm.GR$coefficients,10)+1.96*sqrt(diag(twophase.var.GR)))
ests.GR<-ests.orm.GR
ests.GR["est_bmi_preg_mother",]<-5*ests.GR["est_bmi_preg_mother",]
ests.GR["mat_age_delivery",]<-10*ests.GR["mat_age_delivery",]

ests.orm.IPW<-cbind(tail(mod.orm.IPW$coefficients,10),
      tail(mod.orm.IPW$coefficients,10)-1.96*sqrt(diag(twophase.var.IPW)),
      tail(mod.orm.IPW$coefficients,10)+1.96*sqrt(diag(twophase.var.IPW)))
ests.IPW<-ests.orm.IPW
ests.IPW["est_bmi_preg_mother",]<-5*ests.IPW["est_bmi_preg_mother",]
ests.IPW["mat_age_delivery",]<-10*ests.IPW["mat_age_delivery",]

### I can't figure out how to easily extract the SEs from the orm fit!
ses.mod.orm.n<-c(0.0030,0.0031,0.0758,0.0685,0.0902,0.0510,0.0745,0.0645,0.0381)
ests.orm.n<-cbind(tail(mod.orm.n$coefficients,9),
      tail(mod.orm.n$coefficients,9)-1.96*ses.mod.orm.n,
      tail(mod.orm.n$coefficients,9)+1.96*ses.mod.orm.n)
ests.n<-ests.orm.n
ests.n["est_bmi_preg_mother.ph1",]<-5*ests.n["est_bmi_preg_mother.ph1",]
ests.n["mat_age_delivery.ph1",]<-10*ests.n["mat_age_delivery.ph1",]



## Figure
n.covs<-10
gap<-0.02
pdf("odds-ratios.pdf",height=5)
mar=c(5,2,1,1)
plot(c(-3.2,2),c(-0.02,0.92),type="n",xlab="",ylab="", axes=FALSE)
axis(1, at=log(c(0.2,0.5,1,2,5)), lab=c(0.2,0.5,1,2,5))
mtext("Odds of Greater Weight Gain", side=1, line=3, at=0)
abline(v=0,col=gray(.9))
points(ests.GR[,1],1-c(1:n.covs/n.covs), pch=3, cex=0.5)
for (i in 1:n.covs){
  lines(c(ests.GR[i,2],ests.GR[i,3]),1-rep(i/n.covs,each=2))
}

points(ests.IPW[,1],1-c(1:n.covs)/n.covs-gap, col=3, pch=3,cex=0.5)
for (i in 1:n.covs){
  lines(c(ests.IPW[i,2],ests.IPW[i,3]),1-rep(i/n.covs,each=2)-gap,col=3)
}

points(ests.n[,1],1-c(1:(n.covs-1))/n.covs-2*gap, col=2, pch=3,cex=0.5)
for (i in 1:(n.covs-1)){
  lines(c(ests.n[i,2],ests.n[i,3]),1-rep(i/n.covs,each=2)-2*gap, col=2)
}
text(rep(-3.5,n.covs),1-c(1:n.covs)/n.covs,
      c("BMI (per 5 kg/m2)","Age (per 10 years)","Asian Race (ref White)","Black Race","Other Race","Hispanic Ethnicity",
        "Tobacco","Depression","Private Insurance","Pregnancy length (weeks)"), pos=4)

legend(x="topright",legend=c("GR","IPW","naive"),lwd=c(1,1,1),col=c(1,3,2),bty="n")
dev.off()


###########
### Including splines


mod.orm.GR.ns <- orm(estWtChangePerWk ~ rcs(est_bmi_preg_mother,3)+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurance.priv+
                   rcs(egaWk,3), 
              family = "logistic",
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights)	

# Calculate robust SE by mimicking what is done in svycoxph
orig.var.GR.ns <- infoMxop(mod.orm.GR.ns$info.matrix, i = "x")
inflFun.GR.ns  <- getBetaInfluenceFunctions(mod.orm.GR.ns, getC1 = TRUE)
twophase.var.GR.ns<-twophasevar(cbind(inflFun.GR.ns$C1, inflFun.GR.ns$C2), Ninflcal.orm)
sqrt(diag(orig.var.GR.ns))
sqrt(diag(twophase.var.GR.ns))

coef.est.ns<-tail(mod.orm.GR.ns$coefficients,12)
lower.ns<-coef.est.ns-1.96*tail(sqrt(diag(twophase.var.GR.ns)),12)
upper.ns<-coef.est.ns+1.96*tail(sqrt(diag(twophase.var.GR.ns)),12)
wald.p.ns<-2*pnorm(-abs(coef.est.ns/tail(sqrt(diag(twophase.var.GR.ns)),12)))
cbind(coef.est.ns,lower.ns,upper.ns,wald.p.ns)
exp(cbind(coef.est.ns,lower.ns,upper.ns,wald.p.ns)["mat_age_delivery",]*10)

# Calculate median and 95% CIs as in orm
qu  <- Quantile(mod.orm.GR.ns)                       
med <- function(lp) qu(0.5, lp = lp)  
p  <- Predict(mod.orm.GR.ns, est_bmi_preg_mother = seq(14, 50, 1), fun = med); p

f<-mod.orm.GR.ns
nd <- as.data.frame(p[, c("est_bmi_preg_mother","mat_age_delivery","egaWk")])
B_bmi <- Hmisc::rcspline.eval(nd$est_bmi_preg_mother, 
                              knots = f$Design$parms$est_bmi_preg_mother,
                              inclx = TRUE)
B_egaWk <- Hmisc::rcspline.eval(nd$egaWk, 
                              knots = f$Design$parms$egaWk,
                              inclx = TRUE)

X = as.matrix(cbind(B_bmi, 27.37, 0, 0, 0, 0, 0, 0, 0, B_egaWk))
colnames(X) = c("est_bmi_preg_mother","est_bmi_preg_mother'",
                "mat_age_delivery","mat4race=Asian","mat4race=Black","mat4race=Other",
                "hispanic","tobacPregr","depressionr","insurance.priv",
                "egaWk","egaWk'")


ns   <- f$non.slopes                    # number of intercepts
ints <- f$coef[1:ns]                    # intercept vector (cutpoints)
bet  <- f$coef[-(1:ns)]                 # slopes
vals <- f$yunique; if(!length(vals)) vals <- names(f$freq)
vals <- as.numeric(vals)                  # ordered outcome values (numeric)
inverse <- eval(f$famfunctions[2])
cumprob <- eval(f$famfunctions[1])

all(names(bet) %in% colnames(X))  # Verify all coefficient names are correct and included in X
tau <- 0.5 # Quantile of interest

lp <- X %*% bet 
lb <- matrix(sapply(ints, `+`, lp), ncol = ns)
m_cdf <- cbind(1 - cumprob(lb), 1)

# Get interpolated fitted values
cp <- cbind(cumprob(lb), 0)

# build y-grid for interpolation, row by row
m_yvals <- matrix(NA, nrow(lb), ns + 2)
for (j in seq_len(nrow(lb))) {
  # normalized weights between first and last cutpoint
  ws <- c(0, (cp[j, -(ns+1)] - cp[j, 1]) / (cp[j, ns] - cp[j, 1]), 1)
  m_yvals[j, ] <- (1 - ws) * c(vals[1], vals) + ws * c(vals, vals[ns + 1])
}

# map tau onto the same probability scale and interpolate y
target <- cumprob(inverse(1 - tau))
z <- vapply(seq_len(nrow(lb)), function(i)
  approx(x = c(1, cp[i, ]), y = m_yvals[i, ], xout = target, rule = 2)$y,
  numeric(1)) 


p$yhat - z # Confirm predicted values are the same as in Predict function

# Compute standard errors
lb.se        <- matrix(NA, ncol = ns, nrow = nrow(X))
idx          <- which(names(c(ints, bet)) %in% colnames(X)); idx
dlb.dtheta   <- as.matrix(cbind(1, X))
info.inverse <- twophase.var.GR.ns #infoMxop(f$info.matrix, invert = TRUE) # Variance 

# Doing this the original way since info.inverse will be calculated via survey package
for(i in 1:ns){
  v.i <- info.inverse[c(i, idx), c(i, idx)]
  lb.se[, i] <- sqrt(diag(dlb.dtheta %*% v.i %*% t(dlb.dtheta)))
  # Compute (i, idx) portion of info inverse, multiplied by t(dlb.dtheta)
  #v.i        <- infoMxop(f$info.matrix, i=c(i, idx), B=t(dlb.dtheta))
  #lb.se[, i] <- sqrt(Matrix::diag(dlb.dtheta %*% v.i))
}
w <- qnorm((1 + 0.95) / 2);w
ci.ub <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] - w * lb.se[, i])}), ncol = ns)
ci.lb <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] + w * lb.se[, i])}), ncol = ns)

z.ub <- sapply(1:nrow(lb),
                 function(i) approx(c(1, 1 - ci.lb[i,], 0), m_yvals[i,],
                                    xout = cumprob(inverse(1 - tau)),
                                    rule = 2)$y)
z.lb <- sapply(1:nrow(lb),
                 function(i) approx(c(1, 1 - ci.ub[i,], 0), m_yvals[i,],
                                    xout = cumprob(inverse(1 - tau)),
                                    rule = 2)$y)



stuff<-hist(d$est_bmi_preg_mother.ph1, main = "", xlab = "", ylab = "", axes = FALSE,nclass=100)

pdf("median-weight-gain-by-BMI.pdf",height=5)
mar=c(5,2,1,1)
plot(c(min(p$est_bmi_preg_mother),max(p$est_bmi_preg_mother)),
     c(min(z.lb),max(z.ub)),type="n",xlab="Body Mass Index before Pregnancy (kg/m2)",ylab="Pregnancy Weight Gain (kg/wk)")
lines(p$est_bmi_preg_mother,p$yhat,col=1)
lines(p$est_bmi_preg_mother,z.lb,lty=2,col=1)
lines(p$est_bmi_preg_mother,z.ub,lty=2,col=1)

for(i in 1:length(stuff$density)){
  lines(rep(stuff$mids[i],2), c(min(z.lb), min(z.lb)+stuff$density[i]*2*(max(z.ub)-min(z.lb))) )
}
dev.off()




### Examining the relationship between EGA and Pregnancy Gain


# Calculate median and 95% CIs as in ?orm
qu  <- Quantile(mod.orm.GR.ns)                       
med <- function(lp) qu(0.5, lp = lp)  
p  <- Predict(mod.orm.GR.ns, egaWk = seq(23, 43, 1), fun = med); p

f<-mod.orm.GR.ns
nd <- as.data.frame(p[, c("est_bmi_preg_mother","mat_age_delivery","egaWk")])
B_bmi <- Hmisc::rcspline.eval(nd$est_bmi_preg_mother, 
                              knots = f$Design$parms$est_bmi_preg_mother,
                              inclx = TRUE)
B_egaWk <- Hmisc::rcspline.eval(nd$egaWk, 
                              knots = f$Design$parms$egaWk,
                              inclx = TRUE)

X = as.matrix(cbind(B_bmi, 27.37, 0, 0, 0, 0, 0, 0, 0, B_egaWk))
colnames(X) = c("est_bmi_preg_mother","est_bmi_preg_mother'",
                "mat_age_delivery","mat4race=Asian","mat4race=Black","mat4race=Other",
                "hispanic","tobacPregr","depressionr","insurance.priv",
                "egaWk","egaWk'")


ns   <- f$non.slopes                    # number of intercepts
ints <- f$coef[1:ns]                    # intercept vector (cutpoints)
bet  <- f$coef[-(1:ns)]                 # slopes
vals <- f$yunique; if(!length(vals)) vals <- names(f$freq)
vals <- as.numeric(vals)                  # ordered outcome values (numeric)
inverse <- eval(f$famfunctions[2])
cumprob <- eval(f$famfunctions[1])

all(names(bet) %in% colnames(X))  # Verify all coefficient names are correct and included in X
tau <- 0.5 # Quantile of interest

lp <- X %*% bet 
lb <- matrix(sapply(ints, `+`, lp), ncol = ns)
m_cdf <- cbind(1 - cumprob(lb), 1)

# Get interpolated fitted values
cp <- cbind(cumprob(lb), 0)

# build y-grid for interpolation, row by row
m_yvals <- matrix(NA, nrow(lb), ns + 2)
for (j in seq_len(nrow(lb))) {
  # normalized weights between first and last cutpoint
  ws <- c(0, (cp[j, -(ns+1)] - cp[j, 1]) / (cp[j, ns] - cp[j, 1]), 1)
  m_yvals[j, ] <- (1 - ws) * c(vals[1], vals) + ws * c(vals, vals[ns + 1])
}

# map tau onto the same probability scale and interpolate y
target <- cumprob(inverse(1 - tau))
z <- vapply(seq_len(nrow(lb)), function(i)
  approx(x = c(1, cp[i, ]), y = m_yvals[i, ], xout = target, rule = 2)$y,
  numeric(1)) 


p$yhat - z # Confirm predicted values are the same as in Predict function

# Compute standard errors
lb.se        <- matrix(NA, ncol = ns, nrow = nrow(X))
idx          <- which(names(c(ints, bet)) %in% colnames(X)); idx
dlb.dtheta   <- as.matrix(cbind(1, X))
info.inverse <- twophase.var.GR.ns #infoMxop(f$info.matrix, invert = TRUE) # Variance 

# Doing this the original way since info.inverse will be calculated via survey package
for(i in 1:ns){
  v.i <- info.inverse[c(i, idx), c(i, idx)]
  lb.se[, i] <- sqrt(diag(dlb.dtheta %*% v.i %*% t(dlb.dtheta)))
  # Compute (i, idx) portion of info inverse, multiplied by t(dlb.dtheta)
  #v.i        <- infoMxop(f$info.matrix, i=c(i, idx), B=t(dlb.dtheta))
  #lb.se[, i] <- sqrt(Matrix::diag(dlb.dtheta %*% v.i))
}
w <- qnorm((1 + 0.95) / 2);w
ci.ub <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] - w * lb.se[, i])}), ncol = ns)
ci.lb <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] + w * lb.se[, i])}), ncol = ns)

z.ub <- sapply(1:nrow(lb),
                 function(i) approx(c(1, 1 - ci.lb[i,], 0), m_yvals[i,],
                                    xout = cumprob(inverse(1 - tau)),
                                    rule = 2)$y)
z.lb <- sapply(1:nrow(lb),
                 function(i) approx(c(1, 1 - ci.ub[i,], 0), m_yvals[i,],
                                    xout = cumprob(inverse(1 - tau)),
                                    rule = 2)$y)


stuff<-hist(d$egaWk, main = "", xlab = "", ylab = "", axes = FALSE,nclass=100)

pdf("median-weight-gain-by-EGA.pdf",height=5)
mar=c(5,2,1,1)
plot(c(min(p$egaWk),max(p$egaWk)),
     c(min(z.lb),max(z.ub)),type="n",xlab="Length of Pregnancy (weeks)",ylab="Pregnancy Weight Gain (kg/wk)")
lines(p$egaWk,p$yhat,col=1)
lines(p$egaWk,z.lb,lty=2,col=1)
lines(p$egaWk,z.ub,lty=2,col=1)

for(i in 1:length(stuff$density)){
  lines(rep(stuff$mids[i],2), c(min(z.lb), min(z.lb)+stuff$density[i]*(max(z.ub)-min(z.lb))) )
}
dev.off()








## Now let's estimate P(weight change >0.4)



# Calculate Exceedence probabilities using ExProb

f<-mod.orm.GR.ns
nd <- data.frame(est_bmi_preg_mother=24, mat_age_delivery=27.37303, mat4race=c("White","Asian","Black","Other"), 
                 hispanic=0, tobacPregr=0, depressionr=0, insurance.priv=0, egaWk=39.28571)

d <- ExProb(f)
lp0 <- predict(f, newdata = nd)
w <- d(lp0); w
# w <- d(lp0, X = nd, conf.int = 0.95); w  ## Giving me an error, but I don't think I need it


# Manually doing this (w/ help from Hmisc)

B_bmi <- Hmisc::rcspline.eval(nd$est_bmi_preg_mother, 
                              knots = f$Design$parms$est_bmi_preg_mother,
                              inclx = TRUE)
B_egaWk <- Hmisc::rcspline.eval(nd$egaWk, 
                              knots = f$Design$parms$egaWk,
                              inclx = TRUE)
race.mat <- cbind(c(0,1,0,0), c(0,0,1,0), c(0,0,0,1))


X = as.matrix(cbind(B_bmi, 27.373, race.mat, 0, 0, 0, 0, B_egaWk))
colnames(X) = c("est_bmi_preg_mother","est_bmi_preg_mother'",
                "mat_age_delivery","mat4race=Asian","mat4race=Black","mat4race=Other",
                "hispanic","tobacPregr","depressionr","insurance.priv",
                "egaWk","egaWk'")

# Manual (most of this code is derived from ExProb)
# Get necessary components from fitted model "f"

ns   <- f$non.slopes                    # number of intercepts
ints <- f$coef[1:ns]                    # intercept vector (cutpoints)
bet  <- f$coef[-(1:ns)]                 # slopes
vals <- f$yunique; if(!length(vals)) vals <- names(f$freq)
vals <- as.numeric(vals)                  # ordered outcome values (numeric)
inverse <- eval(f$famfunctions[2])
cumprob <- eval(f$famfunctions[1])
yname <- f$yname
info <- f$info.matrix

all(names(bet) %in% colnames(X))  # Verify all coefficient names are correct and included in X
lp <- X %*% bet 
y <- NULL

# From here on out is from ExProb
prob <- cumprob(sapply(c(1e+30, ints), "+", lp)); prob


dim(prob) <- c(length(lp), length(vals))
if (!length(y)) {
  colnames(prob) <- paste("Prob(Y>=", vals, ")", sep = "")
  y <- vals
  result <- structure(list(y = vals, prob = prob, yname = yname), 
                      class = "ExProb")
}

conf.int = 0.95

index <- sapply(y, FUN = function(x) {
  if (x <= min(vals)) 
    result <- 1
  else if (x >= max(vals)) 
    result <- length(vals)
  else which(x <= vals)[1] - 1
})
idx <- which(names(c(ints, bet)) %in% colnames(X))
dlb.dtheta <- as.matrix(cbind(1, X))



### New code from Eric: 12/17/25
lb.se <- sapply(1:length(y), function(i) Matrix::diag(dlb.dtheta %*% twophase.var.GR.ns[c(index[i], idx), c(index[i], idx)] %*% t(dlb.dtheta)))

# lb.se <- sapply(1:length(y), function(i) Matrix::diag(dlb.dtheta %*% 
#                                                         infoMxop(info, i = c(index[i], idx), B = t(dlb.dtheta))))
# lb.se <- matrix(sqrt(lb.se), ncol = length(y))
m.alpha <- c(ints, bet)[index]
lb <- matrix(sapply(m.alpha, "+", lp), ncol = length(y))
ci.ub <- matrix(sapply(1:length(y), FUN = function(i) {
  cumprob(lb[, i] + qnorm((1 + conf.int)/2) * lb.se[, 
                                                    i])
}), ncol = length(y))
ci.lb <- matrix(sapply(1:length(y), FUN = function(i) {
  cumprob(lb[, i] - qnorm((1 + conf.int)/2) * lb.se[, 
                                                    i])
}), ncol = length(y))
ci.ub[, which(y <= min(vals))] <- ci.lb[, which(y <= 
                                                    min(vals))] <- 1
ci.ub[, which(y >= max(vals))] <- ci.lb[, which(y >= 
                                                    max(vals))] <- 0
if (length(y) > 1) 
  colnames(ci.ub) <- colnames(ci.lb) <- colnames(result$prob)
attr(result, "limits") <- list(lower = ci.lb, upper = ci.ub)

result


pdf("pr-unhealthy-wtgain-by-race.pdf",height=5)
mar=c(5,2,1,1)
plot(c(0,4),c(0.1,0.3), type="n", ylab="Probability of Unhealthy Weight Gain", xlab="Race", axes=FALSE)
axis(2)
axis(1, labels=c("White","Asian","Black","Other/Unknown"), at=c(0.5,1.5,2.5,3.5))
points(c(0.5,1.5,2.5,3.5),
       c(result$prob[1,min(which(result$y>=0.4))],
         result$prob[2,min(which(result$y>=0.4))],
         result$prob[3,min(which(result$y>=0.4))],
         result$prob[4,min(which(result$y>=0.4))]))
lines(c(0.5,0.5),c(ci.lb[1,min(which(result$y>=0.4))],ci.ub[1,min(which(result$y>=0.4))]))
lines(c(1.5,1.5),c(ci.lb[2,min(which(result$y>=0.4))],ci.ub[2,min(which(result$y>=0.4))]))
lines(c(2.5,2.5),c(ci.lb[3,min(which(result$y>=0.4))],ci.ub[3,min(which(result$y>=0.4))]))
lines(c(3.5,3.5),c(ci.lb[4,min(which(result$y>=0.4))],ci.ub[4,min(which(result$y>=0.4))]))
dev.off()





##### Fitting with probit link function
mod.orm.GR.ns.probit <- orm(estWtChangePerWk ~ rcs(est_bmi_preg_mother,3)+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurance.priv+
                   rcs(egaWk,3), 
              family = "probit",
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights)	

##### Fitting with loglog link function
mod.orm.GR.ns.loglog <- orm(estWtChangePerWk ~ rcs(est_bmi_preg_mother,3)+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurance.priv+
                   rcs(egaWk,3), 
              family = "loglog",
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights)	

##### Fitting with cloglog link function
mod.orm.GR.ns.cloglog <- orm(estWtChangePerWk ~ rcs(est_bmi_preg_mother,3)+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurance.priv+
                   rcs(egaWk,3), 
              family = "cloglog",
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights)	

logLik(mod.orm.GR.ns)
logLik(mod.orm.GR.ns.probit)
logLik(mod.orm.GR.ns.cloglog)
logLik(mod.orm.GR.ns.loglog)





## Bootstrapping

set.seed(123)
nboot<-1000
coef.GR<-coef.IPW<-matrix(NA, nrow=nboot, ncol=10)
strat.names<-unique(d$wave4Strata)

for (i in 1:nboot){
  boot.ids1<-boot.ids0<-NULL
  boot.ids1<-sample(which(d$R==1 & d$wave4Strata==strat.names[1]), sum(d$R==1 & d$wave4Strata==strat.names[1]), replace=TRUE)
  boot.ids0<-sample(which(d$R==0 & d$wave4Strata==strat.names[1]), sum(d$R==0 & d$wave4Strata==strat.names[1]), replace=TRUE)
  for (j in 2:length(strat.names)) {
    boot.ids1<-c(boot.ids1,sample(which(d$R==1 & d$wave4Strata==strat.names[j]), sum(d$R==1 & d$wave4Strata==strat.names[j]), replace=TRUE))
    if (sum(d$R==0 & d$wave4Strata==strat.names[j])==1) {
      boot.ids0<-c(boot.ids0, which(d$R==0 & d$wave4Strata==strat.names[j]))
      }
    else { 
      boot.ids0<-c(boot.ids0, sample(which(d$R==0 & d$wave4Strata==strat.names[j]), sum(d$R==0 & d$wave4Strata==strat.names[j]), replace=TRUE))
    }
  }
  d.boot<-d[c(boot.ids1,boot.ids0),]
  
  mod.orm.n <- orm(estWtChangePerWk.ph1 ~ est_bmi_preg_mother.ph1+
                   mat_age_delivery.ph1+mat4race.ph1+ hispanic.ph1 + 
                   tobacPregr.ph1+depressionr.ph1+insurancer.ph1,
                 data = d.boot,
                 family = "logistic", mscore = TRUE)
  naiveInfl.orm <- getBetaInfluenceFunctions(mod.orm.n)$C2
  d.boot$nInfl.bmi<-naiveInfl.orm[,1]
  d.boot$nInfl.age<-naiveInfl.orm[,2]
  d.boot$nInfl.race1<-naiveInfl.orm[,3]
  d.boot$nInfl.race2<-naiveInfl.orm[,4]
  d.boot$nInfl.race3<-naiveInfl.orm[,5]
  d.boot$nInfl.hisp<-naiveInfl.orm[,6]
  d.boot$nInfl.tobac<-naiveInfl.orm[,7]
  d.boot$nInfl.depr<-naiveInfl.orm[,8]
  d.boot$nInfl.insur<-naiveInfl.orm[,9]

            ## Errors if one of the strata are empty in the bootstrap sample. I skip these
            ##  replications when it happens.
  if (min(table(d.boot$wave4Strata,d.boot$R)) >0){
  designO <- twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~R==1,
                      data=d.boot,method="approx")

  Ninflcal.orm <-calibrate(designO,
                           formula=~nInfl.bmi+nInfl.age+nInfl.race1+nInfl.race2+nInfl.race3+
                             nInfl.hisp+nInfl.tobac+nInfl.depr+nInfl.insur+wave4Strata,
                           phase=2,calfun="raking") 

  design <- Ninflcal.orm
  if(any(weights(design) < 0)) stop("weights must be non-negative")
  sdat <- model.frame(design)
  sdat$.survey.prob.weights <- (1 / design$prob) / mean(1 / design$prob)

  mod.orm.GR <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
                     mat_age_delivery+mat4race+ hispanic + 
                     tobacPregr+depressionr+insurancer+
                   egaWk, 
              family = "logistic",
              data = sdat,
              weights = .survey.prob.weights)	

  if(any(weights(designO) < 0)) stop("weights must be non-negative")
  sdatO <- model.frame(designO)
  sdatO$.survey.prob.weights <- (1 / designO$prob) / mean(1 / designO$prob)

  mod.orm.IPW <- orm(estWtChangePerWk ~ est_bmi_preg_mother+
                   mat_age_delivery+mat4race+ hispanic + 
                   tobacPregr+depressionr+insurancer+
                   egaWk, 
              family = "logistic",
              data = sdatO,
              weights = .survey.prob.weights)	

  coef.GR[i,]<-tail(mod.orm.GR$coefficients,10)
  coef.IPW[i,]<-tail(mod.orm.IPW$coefficients,10)
  } 
}

sum(is.na(coef.GR[,1]))
ses.GR<-apply(coef.GR,2,sd,na.rm=TRUE)
ses.IPW<-apply(coef.IPW,2,sd,na.rm=TRUE)

ses.GR
ses.IPW

#> ses.GR
# [1] 0.008407604 0.009866176 0.273082385 0.361608537 0.280631567 0.212396270 0.272297270 0.317152253
# [9] 0.264299321 0.033011428
#> ses.IPW
# [1] 0.01864638 0.01663924 0.33479496 0.40703106 0.32639934 0.25506913 0.33387186 0.29940096 0.26818593
#[10] 0.03007202

lower.GR<-apply(coef.GR,2,quantile,.025,na.rm=TRUE)
upper.GR<-apply(coef.GR,2,quantile,.975,na.rm=TRUE)
lower.IPW<-apply(coef.IPW,2,quantile,.025,na.rm=TRUE)
upper.IPW<-apply(coef.IPW,2,quantile,.975,na.rm=TRUE)

cbind(tail(mod.orm.GR$coefficients,10),ses.GR,lower.GR,upper.GR)
# > cbind(tail(mod.orm.GR$coefficients,10),ses.GR,lower.GR,upper.GR)
#                                         ses.GR    lower.GR     upper.GR
# est_bmi_preg_mother     -0.0527771 0.008407604 -0.07106496 -0.037668773
# mat_age_delivery        -0.0403685 0.009866176 -0.04626717 -0.007294217
# mat4race=Black           0.2805301 0.273082385 -0.03119601  1.049900218
# mat4race=Other           0.5350032 0.361608537  0.21263904  1.579960510
# mat4race=White           0.5379938 0.280631567  0.43038699  1.550947837
# hispanic                 0.1165083 0.212396270 -0.86725441 -0.015302322
# tobacPregr              -0.4463744 0.272297270 -1.23965849 -0.153872621
# depressionr              0.6246101 0.317152253 -0.16831956  1.066538410
# insurancer=Public/Other -0.6853487 0.264299321 -0.61465042  0.395247135
# egaWk                    0.1831924 0.033011428  0.08622936  0.211088914


cbind(tail(mod.orm.IPW$coefficients,10),ses.IPW,lower.IPW,upper.IPW)
# > cbind(tail(mod.orm.IPW$coefficients,10),ses.IPW,lower.IPW,upper.IPW)
#                                        ses.IPW   lower.IPW    upper.IPW
# est_bmi_preg_mother     -0.07563880 0.01864638 -0.09396782 -0.021041166
# mat_age_delivery        -0.05083973 0.01663924 -0.06118130  0.003891894
# mat4race=Black           0.28487350 0.33479496 -0.06082040  1.237602386
# mat4race=Other          -0.34765224 0.40703106 -0.35051113  1.292547398
# mat4race=White           0.10325235 0.32639934  0.14722326  1.370312260
# hispanic                 0.58520627 0.25506913 -0.22496171  0.773622200
# tobacPregr               0.08548488 0.33387186 -1.02224309  0.190894949
# depressionr              0.46309766 0.29940096 -0.13680825  1.046053875
# insurancer=Public/Other -0.76246937 0.26818593 -0.69013914  0.341319360
# egaWk                    0.17872849 0.03007202  0.09671841  0.209963999





