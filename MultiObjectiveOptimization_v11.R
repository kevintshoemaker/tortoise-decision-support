#################################################
#Floreana Multi-objective optimization          #
#E. Hunter and K. Shoemaker, Summer 2019			  #
#################################################



#rm(list=ls())

ELIZABETHOFFICE = F
ELIZABETH = F
KEVINOFFICE = T
KEVIN = F


###########
# GLOBAL VARS
###########

objectives <- c("Abund","Gen_Div","Qvals")

###################################
##########  SET WORKING DIRECTORY

if(ELIZABETHOFFICE) rootDir <- "C:/Users/elizabethhunter/Dropbox/Galapagos/TortoiseRepatriation/Wolf/"
if(ELIZABETH) rootDir <- "~/Dropbox/Galapagos/TortoiseRepatriation/Wolf/"
if(KEVINOFFICE) rootDir <- "E:/Dropbox/Tortoise repatriation/Wolf/"
if(KEVIN) rootDir <- "C:/Users/Kevin/Dropbox/Tortoise repatriation/Wolf/"

ScriptDir <- paste(rootDir,"RCode",sep="")
FiguresDir <- paste(rootDir,"Figs",sep="")
NetLogoDir <- paste(rootDir,"NetLogo",sep="")
PresDir <- paste(rootDir,"Presentation/PresentationImages",sep="")

setwd(NetLogoDir)

getwd()

##################################
##########  LOAD PACKAGES

library(MASS)
library(optimization)
library(mgcv)
library(GenSA)

###################################
##########  READ IN DATA

dat <- read.csv("Wolf_netlogo_model17_MASTER_AllOptions.csv", header=T)
dat2 <- read.csv("Wolf_netlogo_model18_MASTER_AllOptions.csv", header=T)
dat <- rbind(dat, dat2)

rm(dat2)

#Change T/F variables to 0/1
dat$release.adult.hybrids.q <- as.factor(ifelse(dat$release.adult.hybrids==T, 1, 0))
dat$release.adult.hybrids.q <- ifelse(dat$release.adult.hybrids==T, 20, 0)

#Rescale large covariates (helps with optimization)
dat$release.end.s <- (dat$release.end-min(dat$release.end))/(max(dat$release.end)-min(dat$release.end))
dat$sex.ratio.s <- (dat$sex.ratio-min(dat$sex.ratio))/(max(dat$sex.ratio)-min(dat$sex.ratio))
dat$release.age.s <- (dat$release.age-min(dat$release.age))/(max(dat$release.age)-min(dat$release.age))
dat$Corral.capacity.s <- (dat$Corral.capacity-min(dat$Corral.capacity))/(max(dat$Corral.capacity)-min(dat$Corral.capacity))
dat$release.adult.hybrids.s <- (dat$release.adult.hybrids.q-min(dat$release.adult.hybrids.q))/(max(dat$release.adult.hybrids.q)-min(dat$release.adult.hybrids.q))

#Scaling each response variable to 0-1
dat$pop.norm <- (dat$pop.size.wild - min(dat$pop.size.wild)) / (max(dat$pop.size.wild) - min(dat$pop.size.wild))+0.5
dat$cost.norm <- (dat$cost - min(dat$cost)) / (max(dat$cost) - min(dat$cost))+0.5
dat$qval.norm <- (dat$final.mean.qval - min(dat$final.mean.qval)) / (max(dat$final.mean.qval) - min(dat$final.mean.qval))+0.5
dat$sar.norm <- (dat$final.sar - min(dat$final.sar)) / (max(dat$final.sar) - min(dat$final.sar))+0.5
dat$aar.norm <- (dat$final.aar - min(dat$final.aar)) / (max(dat$final.aar) - min(dat$final.aar))+0.5
dat$ho.norm <- (dat$final.ho - min(dat$final.ho)) / (max(dat$final.ho) - min(dat$final.ho))+0.5

nrow(dat)

summary(dat)


names(dat)

actions <- data.frame(
     label = c("Sex Ratio (female:male)", "Number of Breeders", "Translocation Age", "Program Duration (years)", "Direct Release of Adults"),
     unstd = c("sex.ratio","Corral.capacity","release.age","release.end","release.adult.hybrids.q"),
     std = c("sex.ratio.s","Corral.capacity.s","release.age.s","release.end.s","release.adult.hybrids.s"),
     stringsAsFactors = F
)

action_ranges <- lapply(actions$unstd,function(t) range(dat[[t]]))
names(action_ranges) <- actions$label

##########
# Exploratory visualizations
##########

#####
# predictor vars

graphics.off()

hist(dat$release.age) # release age
hist(dat$release.age.s)  # rescaled

hist(dat$release.end) # program duration
hist(dat$release.end.s)   # rescaled

hist(dat$sex.ratio) # sex ratio
hist(dat$sex.ratio.s)

hist(dat$release.adult.hybrids.q)    # release adult hybrids
hist(dat$release.adult.hybrids.s)     # refactored to 0/1

hist(dat$Corral.capacity)     # corral capacity
hist(dat$Corral.capacity.s)   # rescaled

# 5 real predictor vars: release.age.s, release.end.s, sex.ratio.s, release.adult.hybrids.q, Corral.capacity.s

######
# response vars

hist(dat$cost.norm)
hist(dat$sar.norm)
hist(dat$pop.norm)
hist(dat$qval.norm)

###########
# model population abundance objective
###########

best.mods <- list()

mod2.2 <- gam(pop.norm ~ release.adult.hybrids.s*Corral.capacity.s + Corral.capacity.s:release.end.s +
                Corral.capacity.s:release.age.s + Corral.capacity.s:sex.ratio.s + release.adult.hybrids.s:release.end.s +
                release.adult.hybrids.s:release.age.s + release.adult.hybrids.s:sex.ratio.s + 
              s(release.age.s,release.end.s,sex.ratio.s, k=25,bs="ts"),
              data=dat,method="ML",select=T,family="Gamma")  # minimize gcv score

summary(mod2.2)

mod2.3 <- gam(pop.norm ~ release.adult.hybrids.s*Corral.capacity.s + Corral.capacity.s:release.end.s +
                Corral.capacity.s:release.age.s +  Corral.capacity.s:sex.ratio.s + release.adult.hybrids.s:release.end.s +
                s(release.age.s,release.end.s,sex.ratio.s, k=30,bs="ts"),
              data=dat,method="ML",select=T,family=Gamma(link="inverse"))  # minimize gcv score   family = Gamma(link="inverse")

summary(mod2.3)  # best model

##########
# goodness of fit checks

graphics.off()
layout(matrix(1:4,nrow=2))
gam.check(mod2.3)
plot(fitted(mod2.3),residuals(mod2.3)) 

#########
# visualize model


setwd(FiguresDir)
filename=sprintf("gamPop_%s.pdf", Sys.Date())
pdf(filename, width=7, height=5.5)
layout(matrix(1:12,nrow=3,byrow=T))
par(mar=c(2,2,1,1))

vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","sex.ratio.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nDuration", ylab="\nSex Ratio")
vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","release.age.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nDuration", ylab="\nTranslocation Age")
vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","Corral.capacity.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nDuration", ylab="\nNumber Breeders")
vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","release.adult.hybrids.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nDuration", ylab="\nRelease Adults")

vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.age.s","sex.ratio.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nTranslocation Age", ylab="\nSex Ratio")
vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.age.s","Corral.capacity.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nTranslocation Age", ylab="\nNumber Breeders")
vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.age.s","release.adult.hybrids.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nTranslocation Age", ylab="\nRelease Adults")

vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","sex.ratio.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nNumber Breeders", ylab="\nSex Ratio")
vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.adult.hybrids.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nNumber Breeders", ylab="\nRelease Adults")

vis.gam(mod2.3,color="heat",plot.type = "persp",type="response",
        view=c("sex.ratio.s","release.adult.hybrids.s"),
        se=2,theta = 25, phi = 25, zlab="\nPopulation Size", xlab="\nSex Ratio", ylab="\nRelease Adults")

dev.off()

# for now, let's assume mod2.3 is the best model for pop size.

best.mods[[objectives[1]]] <- mod2.3


##########
# Model genetic diversity
##########

gmod2.4 <- gam(sar.norm ~ release.adult.hybrids.s*Corral.capacity.s + Corral.capacity.s:release.end.s +
                 Corral.capacity.s:release.age.s + Corral.capacity.s:sex.ratio.s + release.adult.hybrids.s:release.end.s +
                 release.adult.hybrids.s:release.age.s + release.adult.hybrids.s:sex.ratio.s + 
                 s(release.end.s,k=4) + s(sex.ratio.s,k=3) + s(release.age.s,k=3),
              data=dat,method="ML",select=T,family="Gamma")  # minimize gcv score


gmod2.5 <- gam(sar.norm ~ release.adult.hybrids.s*Corral.capacity.s + sex.ratio.s +
                 release.adult.hybrids.s:sex.ratio.s + 
                 s(release.end.s,k=4), #+ s(sex.ratio.s,k=3) + s(release.age.s,k=3),
               data=dat,method="ML",select=T)  # minimize gcv score  family="Gamma"

summary(gmod2.5)

plot(gmod2.5)

##########
# goodness-of-fit checks

graphics.off()
layout(matrix(1:4,nrow=2))
gam.check(gmod2.5)
plot(fitted(gmod2.5),residuals(gmod2.5))    # some heteroskedasticity 

#########
# visualize model


setwd(FiguresDir)
filename=sprintf("gamDiv_%s.pdf", Sys.Date())
pdf(filename, width=7, height=3.5)
layout(matrix(1:6,nrow=2,byrow=T))
par(mar=c(2,2,1,1))

vis.gam(gmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.end.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenetic Diversity", xlab="\nNumber Breeders", ylab="\nDuration")
vis.gam(gmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","sex.ratio.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenetic Diversity", xlab="\nNumber Breeders", ylab="\nSex Ratio")
vis.gam(gmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.adult.hybrids.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenetic Diversity", xlab="\nNumber Breeders", ylab="\nRelease Adults")
vis.gam(gmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("release.adult.hybrids.s","release.end.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenetic Diversity", xlab="\nRelease Adults", ylab="\nDuration")
vis.gam(gmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("release.adult.hybrids.s","sex.ratio.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenetic Diversity", xlab="\nRelease Adults", ylab="\nSex Ratio")
vis.gam(gmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","sex.ratio.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenetic Diversity", xlab="\nDuration", ylab="\nSex Ratio")


dev.off()

#AIC(gmod2.3,gmod2.4,gmod2.5)

best.mods[[objectives[2]]] <- gmod2.5     # best model explains 50% of the variance

##########
# Model floreana genome
##########

fmod2.4 <- gam(qval.norm ~ release.adult.hybrids.s*Corral.capacity.s + Corral.capacity.s:release.end.s +
                 Corral.capacity.s:release.age.s + Corral.capacity.s:sex.ratio.s + release.adult.hybrids.s:release.end.s +
                 release.adult.hybrids.s:release.age.s + release.adult.hybrids.s:sex.ratio.s + 
                 s(release.end.s,k=4) + s(sex.ratio.s,k=3) + s(release.age.s,k=3),
               data=dat,method="ML",select=T)  # minimize gcv score    family="Gamma"

summary(fmod2.4)

fmod2.5 <- gam(qval.norm ~ release.adult.hybrids.s+Corral.capacity.s+release.end.s + #sex.ratio.s + 
                 Corral.capacity.s:sex.ratio.s + release.adult.hybrids.s:release.end.s+ s(sex.ratio.s,k=3) ,
                # s(release.end.s,k=3) + s(sex.ratio.s,k=3), #+ s(release.age.s,k=3),
               data=dat,method="ML",select=T)  # minimize gcv score    family="Gamma"  family=Gamma(link="inverse")

summary(fmod2.5)

##########
# goodness of fit checks

graphics.off()
layout(matrix(1:4,nrow=2))
gam.check(fmod2.5)
plot(fitted(fmod2.5),residuals(fmod2.5))    # some heteroskedasticity 

#########
# visualize model


plot(fmod2.5)    # best model // not much of deviance really explained...

setwd(FiguresDir)
filename=sprintf("gamGen_%s.pdf", Sys.Date())
pdf(filename, width=7, height=3.5)
layout(matrix(1:6,nrow=2,byrow=T))
par(mar=c(2,2,1,1))

vis.gam(fmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.end.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenome", xlab="\nNumber Breeders", ylab="\nDuration")
vis.gam(fmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","sex.ratio.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenome", xlab="\nNumber Breeders", ylab="\nSex Ratio")
vis.gam(fmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.adult.hybrids.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenome", xlab="\nNumber Breeders", ylab="\nRelease Adults")
vis.gam(fmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("release.adult.hybrids.s","release.end.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenome", xlab="\nRelease Adults", ylab="\nDuration")
vis.gam(fmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("release.adult.hybrids.s","sex.ratio.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenome", xlab="\nRelease Adults", ylab="\nSex Ratio")
vis.gam(fmod2.5,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","sex.ratio.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nGenome", xlab="\nDuration", ylab="\nSex Ratio")

dev.off()
#AIC(fmod2.3,fmod2.4,fmod2.5)

best.mods[[objectives[3]]] <- fmod2.5     # best model explains 33% of the variance


#############
# Cost function
#############

cmod2.2 <- gam(cost.norm ~ release.adult.hybrids.s*Corral.capacity.s + Corral.capacity.s:release.end.s +
                Corral.capacity.s:release.age.s + release.adult.hybrids.s:release.end.s +
                release.adult.hybrids.s:release.age.s + 
                s(release.age.s,release.end.s,k=20,bs="ts"),
              data=dat,method="ML",select=T,family="Gamma")  # minimize gcv score
summary(cmod2.2)

cmod2.3 <- gam(cost.norm ~ release.adult.hybrids.s*Corral.capacity.s + Corral.capacity.s:release.end.s +
                 Corral.capacity.s:release.age.s + release.adult.hybrids.s:release.end.s +
                 release.adult.hybrids.s:release.age.s + 
                 s(release.age.s,release.end.s,k=20,bs="ts"),
               data=dat,method="ML",select=T,family=Gamma(link="inverse"))  # minimize gcv score family="Gamma"
summary(cmod2.3)   # best model

##########
# goodness of fit checks

graphics.off()
layout(matrix(1:4,nrow=2))
gam.check(cmod2.3)
plot(fitted(cmod2.3),residuals(cmod2.3))    # some heteroskedasticity 

#########
# visualize model

setwd(FiguresDir)
filename=sprintf("gamCost_%s.pdf", Sys.Date())
pdf(filename, width=7, height=3.5)
layout(matrix(1:6,nrow=2,byrow=T))
par(mar=c(2,2,1,1))

vis.gam(cmod2.3,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.end.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nCost", xlab="\nNumber Breeders", ylab="\nDuration")
vis.gam(cmod2.3,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.age.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nCost", xlab="\nNumber Breeders", ylab="\nTranslocation Age")
vis.gam(cmod2.3,color="heat",plot.type = "persp",type="response",
        view=c("Corral.capacity.s","release.adult.hybrids.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nCost", xlab="\nNumber Breeders", ylab="\nRelease Adults")
vis.gam(cmod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.adult.hybrids.s","release.end.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nCost", xlab="\nRelease Adults", ylab="\nDuration")
vis.gam(cmod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.adult.hybrids.s","release.age.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nCost", xlab="\nRelease Adults", ylab="\nTranslocation Age")
vis.gam(cmod2.3,color="heat",plot.type = "persp",type="response",
        view=c("release.end.s","release.age.s"),
        cond=list("release.adult.hybrids.q" = factor("1",levels = levels(dat$release.adult.hybrids.q)) ),
        se=2,theta = 25, phi = 25, zlab="\nCost", xlab="\nDuration", ylab="\nTranslocation Age")
dev.off()

#AIC(cmod2.3,cmod2.5,cmod2.4)    # model 2.3 is way better!

cost.mod <- cmod2.3

costfunc <- function(release.end.s=0.2,release.age.s=0.2,
                     release.adult.hybrids.s=0.2,
                     Corral.capacity.s=0.2){
    #cost <- 0.2+release.end.s*2+((9*release.age.s)*release.end.s)^1.1+as.numeric(release.adult.hybrids.q)*0.5
    temp <- data.frame(
      release.end.s=release.end.s,
      release.age.s=release.age.s,
      release.adult.hybrids.s=release.adult.hybrids.s,
      Corral.capacity.s=Corral.capacity.s
    )
    cost <- predict(cost.mod,newdata = temp,type="response")-0.5
    cost <- (cost * (max(dat$cost) - min(dat$cost)))+min(dat$cost)
    return(cost/1000000)
}

costfunc(0.1,0.1,0.1,0.1)    # in million $



#############
# Objective functions
#############

optvec <- c(
  sex.ratio.s = 0.5,
  Corral.capacity.s = 0.5,
  release.age.s = 0.5,
  release.end.s = 0.5,
  release.adult.hybrids.s = 0.5
)

optvec_l <- c(
  sex.ratio.s = 0,
  Corral.capacity.s = 0,
  release.age.s = -0.25,
  release.end.s = 0,
  release.adult.hybrids.s = 0
)

optvec_u <- c(
  sex.ratio.s = 1.5,
  Corral.capacity.s = 1.5,
  release.age.s = 1.5,
  release.end.s = 1.5,
  release.adult.hybrids.s = 1.5
)

control = list(t0 = 1000,
               nlimit = 100,
               r = 0.2,
               t_min = 0.1,
               dyn_rf = T,
               rf = 1
)

optvec

objfunc <- function(cost.threshold=1,objective=objectives[1]){ # optvec,
  
  # release.adult.hybrids.q<-factor(as.numeric(release.adult.hybrids),
  #                                 levels=c("0","1"))
  function(optvec){
    newdata <- data.frame(
      sex.ratio.s=optvec[1],
      Corral.capacity.s = optvec[2],
      release.age.s = optvec[3],
      release.end.s = optvec[4],
      release.adult.hybrids.s= optvec[5]
    )
    value1 <- max(0.01,predict(best.mods[[objective]],newdata=newdata,type="response")-0.5)
    value2 <- 1/value1
    cost <- costfunc(optvec[4],optvec[3],optvec[5],optvec[2])
    if(cost>cost.threshold) value2 = value2+20
   
    return(value2)
  }
} 

objfunc()(optvec)

multi_objfunc <- function(cost.threshold=3){ # optvec,

  function(optvec){
    newdata <- data.frame(
      sex.ratio.s=optvec[1],
      Corral.capacity.s = optvec[2],
      release.age.s = optvec[3],
      release.end.s = optvec[4],
      release.adult.hybrids.s= optvec[5]
    )
    values <- sapply(1:3,function(t) max(0.01,predict(best.mods[[t]],newdata=newdata,type="response")-0.5))
    values <- mean(values)   # could apply weighting scheme here
    value2 <- 1/values
    cost <- costfunc(optvec[4],optvec[3],optvec[5],optvec[2])
    if(cost>cost.threshold) value2 = value2+20
    
    return(value2)
  }
}

multi_objfunc(cost.threshold=3)(optvec)

# opt1 <- optim_sa(fun=objfunc(cost.threshold=3,objective=objectives[1]),
#                  start=optvec,trace=T,lower=optvec_l,upper=optvec_u)   

opt1 <- GenSA::GenSA(par=optvec,
                     fn=objfunc(cost.threshold=3,objective=objectives[1]),
                     lower = optvec_l,
                     upper = optvec_u,
                     control=list(
                       smooth=F,
                       max.time=40,
                       temperature=100000
                     )
)

opt2 <- GenSA::GenSA(par=optvec,
                     fn=multi_objfunc(cost.threshold=3),
                     lower = optvec_l,
                     upper = optvec_u,
                     control=list(
                       smooth=F,
                       max.time=40,
                       temperature=100000
                     )
)


multi_objfunc(cost.threshold=3)(opt2$par)


##############
# partial dependence plots for GAMS
##############

actions
objectives

# svg("Univariate_GAM_vis.svg",7,6)
# 
# layout(matrix(1:6,nrow=2,byrow=T))
# a=4
# for(a in 1:nrow(actions)){
#   plot(1,1,pch="",xlim=c(optvec_l[actions$std[a]],optvec_u[actions$std[a]]),
#        ylim=c(-0.1,2.1),xaxt="n",xlab=actions$label[a],ylab="Predicted success (standardized")  
#   xrange = seq(optvec_l[actions$std[a]],optvec_u[actions$std[a]],length=100)
#   axis(1,at=seq(optvec_l[actions$std[a]],optvec_u[actions$std[a]],length=5),
#        labels = seq(optvec_l[[actions$unstd[a]]],optvec_u[[actions$unstd[a]]],length=5))
#   otheractions <- actions$std[-a]
#   
#   nd <- data.frame(
#     temp <- xrange
#   )
#   
#   names(nd) <- actions$std[a]
#   for(a2 in otheractions){
#     nd[[a2]] <- 0.5
#   }
#   #nd 
#   o=1
#   cols <- rainbow(3)
#   for(o in 1:length(objectives)){
#     thismod <- best.mods[[objectives[o]]]
#     pred <- predict(thismod,newdata=nd,se.fit=T,type="response")
#     points(nd[[1]],pred$fit-0.5,type="l",col=cols[o],lwd=2)
#     
#   }
#   legend("topleft",col=cols,lwd=2,legend=objectives,bty="n")
# }
# 
# dev.off()


setwd(FiguresDir)
filename=sprintf("GAM_partials_%s.pdf", Sys.Date())
pdf(filename, width=7, height=5)
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.5, 0), las=1, tcl=-0.2, bty="l")
layout(matrix(1:6,nrow=2,byrow=T), widths=c(2.6, 2.2, 2.2))

col.vec <- c("firebrick2", "darkorange", "gold1")
line.types <- c(3,2,1)

#Program Duration
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.4, 0))
plot(1,1,pch="",xlim=c(optvec_l[actions$std[4]],optvec_u[actions$std[4]]),
     ylim=c(-0.1,1.2),xaxt="n",xlab=actions$label[4],ylab="Predicted success (standardized)")  
xrange = seq(optvec_l[actions$std[4]],optvec_u[actions$std[4]],length=100)
axis(1,at=seq(optvec_l[actions$std[4]],optvec_u[actions$std[4]],length=5),
     labels = seq(optvec_l[[actions$std[4]]]*(max(dat[[actions$unstd[4]]])-min(dat[[actions$unstd[4]]]))+min(dat[[actions$unstd[4]]]),optvec_u[[actions$std[4]]]*(max(dat[[actions$unstd[4]]])-min(dat[[actions$unstd[4]]]))+min(dat[[actions$unstd[4]]]),length=5))
otheractions <- actions$std[-4]
nd <- data.frame(
  temp <- xrange
)
names(nd) <- actions$std[4]
for(a2 in otheractions){
  nd[[a2]] <- 0.5
}
o=1
cols <- col.vec
for(o in 1:length(objectives)){
  thismod <- best.mods[[objectives[o]]]
  pred <- predict(thismod,newdata=nd,se.fit=T,type="response")
  points(nd[[1]],pred$fit-0.5,type="l",col=cols[o],lwd=2,lty=line.types[o])
}
legend("topright", "A", bty="n")

#Translocation Age
par(mar=c(3, 0, 1, 1), mgp=c(2, 0.5, 0))
plot(1,1,pch="",xlim=c(optvec_l[actions$std[3]],optvec_u[actions$std[3]]),
     ylim=c(-0.1,1.2),xaxt="n",xlab=actions$label[3],ylab="", yaxt="n")  
axis(1,at=seq(optvec_l[actions$std[3]],optvec_u[actions$std[3]],length=5),
     labels = seq(optvec_l[[actions$std[3]]]*(max(dat[[actions$unstd[3]]])-min(dat[[actions$unstd[3]]]))+min(dat[[actions$unstd[3]]]),optvec_u[[actions$std[3]]]*(max(dat[[actions$unstd[3]]])-min(dat[[actions$unstd[3]]]))+min(dat[[actions$unstd[3]]]),length=5))
xrange = seq(optvec_l[actions$std[3]],optvec_u[actions$std[3]],length=100)
otheractions <- actions$std[-3]
nd <- data.frame(
  temp <- xrange
)
names(nd) <- actions$std[3]
for(a2 in otheractions){
  nd[[a2]] <- 0.5
}
o=1
cols <- col.vec
for(o in 1:length(objectives)){
  thismod <- best.mods[[objectives[o]]]
  pred <- predict(thismod,newdata=nd,se.fit=T,type="response")
  points(nd[[1]],pred$fit-0.5,type="l",col=cols[o],lwd=2,lty=line.types[o])
}
legend("topright", "B", bty="n")

#Legend
plot.new()
legend("center",lty=c(3,2,1),lwd=c(2,2,2),col=c("firebrick2","darkorange","gold1"), legend=c("Population", "Diversity", "Genome"),bty="n", cex=1.5)

#Sex Ratio
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.4, 0))
plot(1,1,pch="",xlim=c(optvec_l[actions$std[1]],optvec_u[actions$std[1]]),
     ylim=c(-0.1,1.2),xaxt="n",xlab=actions$label[1],ylab="Predicted success (standardized)")  
xrange = seq(optvec_l[actions$std[1]],optvec_u[actions$std[1]],length=100)
axis(1,at=seq(optvec_l[actions$std[1]],optvec_u[actions$std[1]],length=4),
     labels = seq(min(dat[[actions$unstd[1]]]),optvec_u[[actions$std[1]]]*(max(dat[[actions$unstd[1]]])-min(dat[[actions$unstd[1]]]))+min(dat[[actions$unstd[1]]]),length=4))
otheractions <- actions$std[-1]
nd <- data.frame(
  temp <- xrange
)
names(nd) <- actions$std[1]
for(a2 in otheractions){
  nd[[a2]] <- 0.5
}
o=1
cols <- col.vec
for(o in 1:length(objectives)){
  thismod <- best.mods[[objectives[o]]]
  pred <- predict(thismod,newdata=nd,se.fit=T,type="response")
  points(nd[[1]],pred$fit-0.5,type="l",col=cols[o],lwd=2,lty=line.types[o])
}
legend("topright", "C", bty="n")

#Number of breeders
par(mar=c(3, 0, 1, 1), mgp=c(2, 0.5, 0))
plot(1,1,pch="",xlim=c(optvec_l[actions$std[2]],optvec_u[actions$std[2]]),
     ylim=c(-0.1,1.2),xaxt="n",xlab=actions$label[2],ylab="", yaxt="n")  
xrange = seq(optvec_l[actions$std[2]],optvec_u[actions$std[2]],length=100)
axis(1,at=seq(optvec_l[actions$std[2]],optvec_u[actions$std[2]],length=5),
     labels = seq(optvec_l[[actions$std[2]]]*(max(dat[[actions$unstd[2]]])-min(dat[[actions$unstd[2]]]))+min(dat[[actions$unstd[2]]]),optvec_u[[actions$std[2]]]*(max(dat[[actions$unstd[2]]])-min(dat[[actions$unstd[2]]]))+min(dat[[actions$unstd[2]]]),length=5))
otheractions <- actions$std[-2]
nd <- data.frame(
  temp <- xrange
)
names(nd) <- actions$std[2]
for(a2 in otheractions){
  nd[[a2]] <- 0.5
}
o=1
cols <- col.vec
for(o in 1:length(objectives)){
  thismod <- best.mods[[objectives[o]]]
  pred <- predict(thismod,newdata=nd,se.fit=T,type="response")
  points(nd[[1]],pred$fit-0.5,type="l",col=cols[o],lwd=2,lty=line.types[o])
}
legend("topright", "D", bty="n")

#Direct release
par(mar=c(3, 0, 1, 1), mgp=c(2, 0.5, 0))
plot(1,1,pch="",xlim=c(optvec_l[actions$std[5]],optvec_u[actions$std[5]]),
     ylim=c(-0.1,1.2),xaxt="n",xlab=actions$label[5],ylab="", yaxt="n")  
xrange = seq(optvec_l[actions$std[5]],optvec_u[actions$std[5]],length=100)
axis(1,at=seq(optvec_l[actions$std[5]],optvec_u[actions$std[5]],length=5),
     labels = seq(optvec_l[[actions$std[5]]]*(max(dat[[actions$unstd[5]]])-min(dat[[actions$unstd[5]]]))+min(dat[[actions$unstd[5]]]),optvec_u[[actions$std[5]]]*(max(dat[[actions$unstd[5]]])-min(dat[[actions$unstd[5]]]))+min(dat[[actions$unstd[5]]]),length=5))
otheractions <- actions$std[-5]
nd <- data.frame(
  temp <- xrange
)
names(nd) <- actions$std[5]
for(a2 in otheractions){
  nd[[a2]] <- 0.5
}
o=1
cols <- col.vec
for(o in 1:length(objectives)){
  thismod <- best.mods[[objectives[o]]]
  pred <- predict(thismod,newdata=nd,se.fit=T,type="response")
  points(nd[[1]],pred$fit-0.5,type="l",col=cols[o],lwd=2,lty=line.types[o])
}
legend("topright", "E", bty="n")

dev.off()

###########
# Loop through cost constraints
###########

bestvals <- list()
optims <- list()
#rahs <- list()

temp <- sapply(objectives,function(t) bestvals[[t]] <<- numeric(0) )
temp <- sapply(objectives,function(t) optims[[t]] <<- list() )
#temp <- sapply(objectives,function(t) rahs[[t]] <<- numeric(0) )


bestvals_mo <- numeric(0)     # add multi-objective solution
optims_mo <- list()
#rahs_mo <- numeric(0)

cts <- c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
ct=5
ctr=0
for(ct in cts){
  ctr <- ctr+1
  obj = objectives[1]
  for(obj in objectives){
    opt <- GenSA::GenSA(par=optvec,
                        fn=objfunc(cost.threshold=ct,objective=obj),
                        lower = optvec_l,
                        upper = optvec_u,
                        control=list(
                          smooth=F,
                          max.time=100,
                          temperature=100000
                        ) 
    )
    
    #best <- which.min(sapply(opts,function(t) t$function_value/100))
    #rahs[[obj]] <- c(rahs[[obj]],2-best)
    bestvals[[obj]] <- c(bestvals[[obj]],opt$value)
    optims[[obj]][[ctr]] <- opt
  }
  
  ########
  # multi-objective solution
  ########
  
  opt <- GenSA::GenSA(par=optvec,
                      fn=multi_objfunc(cost.threshold=ct),
                      lower = optvec_l,
                      upper = optvec_u,
                      control=list(
                        smooth=F,
                        max.time=100,
                        temperature=100000
                      ) 
  )   
  #best <- which.min(sapply(opts,function(t) t$function_value/100))
  #rahs_mo <- c(rahs_mo,2-best)
  bestvals_mo <- c(bestvals_mo,opt$value)
  optims_mo[[ctr]] <- opt
}

# plot(optims[[3]][[8]])
# plot(optims_mo[[3]])

temp <- sapply(objectives,function(t) names(optims[[t]]) <<- cts)
temp <- sapply(objectives,function(t) bestvals[[t]] <<- 1/bestvals[[t]])


#############
# Change useless params to lowest-cost possibility

term_in_model <- list()

obj=objectives[1]
for(obj in objectives){
  thisform <- as.character(best.mods[[obj]]$formula[3])
  term_in_model[[obj]] <- sapply(1:nrow(actions),function(t) grepl(pattern=actions$std[t],x=thisform) )
  names(term_in_model[[obj]]) <- actions$std
}
term_in_model

optims2 <- optims
obj=objectives[2]
for(obj in objectives){
  
  ct=1
  for(ct in 1:length(cts)){
    optims2[[obj]][[ct]]$par[!term_in_model[[obj]]] <- optvec_l[!term_in_model[[obj]]]
  }
}

optims2[[obj]][[ct]]$par

optims <- optims2




############
# Plot without error bars
############

graphics.off()

plot(bestvals[[objectives[1]]]~cts,xlab="Cost constraint (million USD)",
     ylab="Quasi-optimal achievable value",type="b",col="blue",lwd=2,ylim=c(0,1.5))
points(cts,bestvals[[objectives[2]]],col="red",type="b",lwd=2)
points(cts,bestvals[[objectives[3]]],col="green",type="b",lwd=2)

legend("topleft",pch=c(1,1,1),lty=c(1,1,1),lwd=c(2,2,2),col=c("blue","red","green"),legend=c("Abund","Ho","q-val"),bty="n")


##########
# Plot with error bars
##########

PI <- list()
temp <- sapply(1:3, function(t) PI[[objectives[t]]] <<- data.frame(
  lower = numeric(length(cts)),
  med = 0,
  mean = 0,
  upper = 0
) )

#ctr=0
ct=3
for(ct in 1:length(cts)){
  #ctr <- ctr+1
  obj = objectives[1]
  for(obj in objectives){
    vals <- optims[[obj]][[ct]]$par
    newdata <- data.frame(
      sex.ratio.s=vals[1],
      Corral.capacity.s = vals[2],
      release.age.s = vals[3],
      release.end.s = vals[4],
      release.adult.hybrids.s= vals[5]
    )
    value1 <- predict(best.mods[[obj]],newdata=newdata,type="response")-0.5
    
    ## extract parameter estimates and cov matrix...
    beta <- coef(best.mods[[obj]])
    Vb <- vcov(best.mods[[obj]])
    
    ## simulate replicate beta vectors from posterior...
    Cv <- chol(Vb)
    n.rep=100
    nb <- length(beta)
    
    br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta
    
    ## turn these into replicate linear predictors...
    Xp <- predict(best.mods[[obj]],newdata=newdata,type="lpmatrix")
    lp <- Xp%*%br
    
    if(best.mods[[obj]]$family$family=="Gamma"){
      fv <- 1/lp # exp(lp) ## ... finally, replicate expected value vectors
    }else{
      fv <- lp
    }
    
    ## now simulate from Gamma deviates with mean as in fv
    ## and estimated scale...
    if(best.mods[[obj]]$family$family=="Gamma"){
      yr <- matrix(rgamma(fv*0,shape=1/best.mods[[obj]]$scale,scale=fv*best.mods[[obj]]$scale),nrow(fv),ncol(fv))
    }else{
      yr <- matrix(rnorm(fv*0,mean=fv,sd=sqrt(best.mods[[obj]]$sig2)),nrow(fv),ncol(fv))
      
    }
    
    ## compute 95% prediction interval...
    temp <- apply(yr-0.5,1,quantile,prob=c(.025,0.5,0.975))
    PI[[obj]]$lower[ct] <- temp[1,1]
    PI[[obj]]$med[ct] <- temp[2,1]
    PI[[obj]]$mean[ct] <- value1
    PI[[obj]]$upper[ct] <- temp[3,1]
   
  }
  
}

setwd(FiguresDir)
filename=sprintf("QuasiOpt_Achievable_%s.pdf", Sys.Date())
pdf(filename, width=4.5, height=4)
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.5, 0), bty="l", las=1, tcl=-0.2)

ndx <- optvals[[objectives[1]]]<6

plot(bestvals[[objectives[1]]]~cts,xlab="Cost constraint (million USD)",
     ylab="Quasi-optimal achievable value",type="b",bg="firebrick2",lwd=1,ylim=c(0,2), pch=21)
x <- c(cts[ndx],rev(cts[ndx]))
x[c(1,length(x))] <- x[c(1,length(x))]-0.05
x[c(length(ndx),length(ndx)-1)] <- x[c(length(ndx),length(ndx)-1)]+0.05
polygon(x,
        c(PI[[objectives[1]]]$lower[ndx],rev(PI[[objectives[1]]]$upper[ndx])),border = F,col=rgb(238,44,44,75, maxColorValue = 255))
polygon(x,
        c(PI[[objectives[2]]]$lower[ndx],rev(PI[[objectives[2]]]$upper[ndx])),border = F,col=rgb(255,140,0,75, maxColorValue = 255))
polygon(x,
        c(PI[[objectives[3]]]$lower[ndx],rev(PI[[objectives[3]]]$upper[ndx])),border = F,col=rgb(255,215,0,75, maxColorValue = 255))


points(cts,bestvals[[objectives[1]]],bg="firebrick2",type="b",lwd=1, pch=21)
points(cts,bestvals[[objectives[2]]],bg="darkorange",type="b",lwd=1, pch=21)
points(cts,bestvals[[objectives[3]]],bg="gold1",type="b",lwd=1, pch=21)
points(0.5, 0.0472, pch=19)

legend("topleft",pch=c(21,21,21), pt.bg=c("firebrick2","darkorange","gold1"),legend=c("Population","Diversity","Genome"),bty="n")

dev.off()


###########
# add multi-objective solution

PI2 <- list()
temp <- sapply(1:3, function(t) PI2[[objectives[t]]] <<- data.frame(
  lower = numeric(length(cts)),
  med = 0,
  mean = 0,
  upper = 0
) )


#ctr=0
ct=3
for(ct in 1:length(cts)){
  #ctr <- ctr+1
  obj = objectives[1]
  for(obj in objectives){
    vals <- optims_mo[[ct]]$par
    newdata <- data.frame(
      sex.ratio.s=vals[1],
      Corral.capacity.s = vals[2],
      release.age.s = vals[3],
      release.end.s = vals[4],
      release.adult.hybrids.s= vals[5]
    )
    value1 <- predict(best.mods[[obj]],newdata=newdata,type="response")-0.5
    
    ## extract parameter estimates and cov matrix...
    beta <- coef(best.mods[[obj]])
    Vb <- vcov(best.mods[[obj]])
    
    ## simulate replicate beta vectors from posterior...
    Cv <- chol(Vb)
    n.rep=100
    nb <- length(beta)
    
    br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta
    
    ## turn these into replicate linear predictors...
    Xp <- predict(best.mods[[obj]],newdata=newdata,type="lpmatrix")
    lp <- Xp%*%br
    
    if(best.mods[[obj]]$family$family=="Gamma"){
      fv <- 1/lp # exp(lp) ## ... finally, replicate expected value vectors
    }else{
      fv <- lp
    }
    
    ## now simulate from Gamma deviates with mean as in fv
    ## and estimated scale...
    if(best.mods[[obj]]$family$family=="Gamma"){
      yr <- matrix(rgamma(fv*0,shape=1/best.mods[[obj]]$scale,scale=fv*best.mods[[obj]]$scale),nrow(fv),ncol(fv))
    }else{
      yr <- matrix(rnorm(fv*0,mean=fv,sd=sqrt(best.mods[[obj]]$sig2)),nrow(fv),ncol(fv))
      
    }
    
    ## compute 95% prediction interval...
    temp <- apply(yr-0.5,1,quantile,prob=c(.025,0.5,0.975))
    PI2[[obj]]$lower[ct] <- temp[1,1]
    PI2[[obj]]$med[ct] <- temp[2,1]
    PI2[[obj]]$mean[ct] <- value1
    PI2[[obj]]$upper[ct] <- temp[3,1]
    
  }
  
}


setwd(FiguresDir)
filename=sprintf("QuasiOpt_Achievable_3pane_%s.pdf", Sys.Date())
pdf(filename, width=7, height=3)
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.5, 0), bty="l", las=1, tcl=-0.2)

ndx <- optvals[[objectives[1]]]<6
x <- c(cts[ndx],rev(cts[ndx]))
x[c(1,length(x))] <- x[c(1,length(x))]-0.05
x[c(length(ndx),length(ndx)-1)] <- x[c(length(ndx),length(ndx)-1)]+0.05

layout(matrix(1:3,nrow=1,byrow=T), widths=c(2.6,2.2,2.2))
par(mar=c(3, 3, 1, 1))
plot(bestvals[[objectives[1]]]~cts,xlab="Cost constraint (million USD)",
     ylab="Quasi-optimal achievable value",type="b",col="firebrick2",lwd=1,ylim=c(0,1.7),main="Population Size")
polygon(x,
        c(PI2[[objectives[1]]]$lower[ndx]-0.02,rev(PI2[[objectives[1]]]$upper[ndx]-0.02)),border = F,col=rgb(230,230,230,150, maxColorValue = 255))
polygon(x,
        c(PI[[objectives[1]]]$lower[ndx],rev(PI[[objectives[1]]]$upper[ndx])),border = F, col=rgb(238,44,44,75, maxColorValue = 255))
points(cts,bestvals[[objectives[1]]],pch=21,bg="firebrick2",type="b",lwd=1)
points(cts[ndx],PI2[[objectives[1]]]$mean[ndx]-0.02,type="l",col="black",lwd=1,lty=2)

par(mar=c(3, 1, 1, 1))
plot(bestvals[[objectives[2]]]~cts,xlab="Cost constraint (million USD)",
     ylab="",type="b",col="darkorange",lwd=1,ylim=c(0,1.7),main="Overall Genetic Diversity", yaxt="n")
polygon(x,
        c(PI2[[objectives[2]]]$lower[ndx]-0.02,rev(PI2[[objectives[2]]]$upper[ndx]-0.02)),border = F,col=rgb(230,230,230,150, maxColorValue = 255))
polygon(x,
        c(PI[[objectives[2]]]$lower[ndx],rev(PI[[objectives[2]]]$upper[ndx])),border = F,col=rgb(255,140,0,75, maxColorValue = 255))
points(cts,bestvals[[objectives[2]]],pch=21,bg="darkorange",type="b",lwd=1)
points(cts[ndx],PI2[[objectives[2]]]$mean[ndx]-0.02,type="l",col="black",lwd=1,lty=2)

par(mar=c(3, 1, 1, 1))
plot(bestvals[[objectives[3]]]~cts,xlab="Cost constraint (million USD)",
     ylab="",type="b",col="gold1",lwd=1,ylim=c(0,1.7),main="Genome Representation", yaxt="n")
polygon(x,
        c(PI2[[objectives[3]]]$lower[ndx]-0.02,rev(PI2[[objectives[3]]]$upper[ndx]-0.02)),border = F,col=rgb(230,230,230,150, maxColorValue = 255))
polygon(x,
        c(PI[[objectives[3]]]$lower[ndx],rev(PI[[objectives[3]]]$upper[ndx])),border = F,col=rgb(255,215,0,75, maxColorValue = 255))
points(cts,bestvals[[objectives[3]]],pch=21,bg="gold1",type="b",lwd=1)
points(cts[ndx],PI2[[objectives[3]]]$mean[ndx]-0.02,type="l",col="black",lwd=1,lty=2)

dev.off()


###########
# Suboptimality matrix
###########

subopt_mat <- matrix(0,nrow=3,ncol=3)
colnames(subopt_mat) <- objectives     # model we optimize for
rownames(subopt_mat) <- objectives     # model we assess optimality for

# sub-optimality of (e.g.) best-abundance solution (column) for genetic diversity (row)

i=objectives[1];j=objectives[2]
for(i in objectives){      # loop through rows
  for(j in objectives){    # loop through columns  
    subopt_mat[i,j] <- mean(sapply(1:length(cts),function(t) 
        (1/(objfunc(cost.threshold=cts[t],objective=i)(optims[[i]][[t]]$par)) - 
        1/(objfunc(cost.threshold=cts[t],objective=i)(optims[[j]][[t]]$par)))/
          (1/(objfunc(cost.threshold=cts[t],objective=i)(optims[[i]][[t]]$par))) ))
    
  }
}

subopt_mat   # suboptimality matrix

## add a multi-objective solution

subopt_mat <- cbind(subopt_mat,c(0,0,0))
colnames(subopt_mat)[4] <- "Multi-obj"

for(i in objectives){
  subopt_mat[i,4] <- mean(sapply(1:length(cts),
              function(t) (1/(objfunc(cost.threshold=cts[t],objective=i)(optims[[i]][[t]]$par)) - 
                1/(objfunc( cost.threshold=cts[t],objective=i)(optims_mo[[t]]$par)))/
                (1/(objfunc(cost.threshold=cts[t],objective=i)(optims[[i]][[t]]$par))) ))
}

setwd(rootDir)
write.csv(subopt_mat,row.names=T,"subopt_mat.csv")

subopt_mat

library(raster)



# sprintf(c(10,1), fmt = '%#.1f')

###########
# Suboptimality matrix at 2 million dollars
###########

subopt_mat2 <- matrix(0,nrow=3,ncol=3)
colnames(subopt_mat2) <- objectives     # model we optimize for
rownames(subopt_mat2) <- objectives     # model we assess optimality for

# sub-optimality of (e.g.) best-abundance solution (column) for genetic diversity (row)

i=objectives[1];j=objectives[2]
for(i in objectives){      # loop through rows
  for(j in objectives){    # loop through columns  
    subopt_mat2[i,j] <- (1/(objfunc(cost.threshold=2,objective=i)(optims[[i]][[4]]$par)) - 
                                    1/(objfunc(cost.threshold=2,objective=i)(optims[[j]][[4]]$par)) )/
                                   (1/(objfunc(cost.threshold=2,objective=i)(optims[[i]][[4]]$par)))
    
  }
}

subopt_mat2   # suboptimality matrix


## add a multi-objective solution

subopt_mat2 <- cbind(subopt_mat2,c(0,0,0))
colnames(subopt_mat2)[4] <- "Multi-obj"

for(i in objectives){
  subopt_mat2[i,4] <- (1/(objfunc(cost.threshold=2,objective=i)(optims[[i]][[4]]$par)) - 
                         1/(objfunc(cost.threshold=2,objective=i)(optims_mo[[4]]$par)) )/
                         (1/(objfunc(cost.threshold=2,objective=i)(optims[[i]][[4]]$par)))    
}

setwd(rootDir)
write.csv(subopt_mat2,row.names=T,"subopt_mat2.csv")


setwd(FiguresDir)
filename=sprintf("Suboptimality_%s.pdf", Sys.Date())
pdf(filename, width=4.5, height=3)
par(mar=c(7, 2, 3, 1), mgp=c(2, 0.5, 0), bty="n", las=1, tcl=-0.2)
layout(matrix(1:2,nrow=1), widths=c(2.65,1.85))

par(mar=c(7, 6, 3, 1))
r <- raster(extent(0,4,0,3),res=1)
r <- setValues(r,subopt_mat)
image(r,col=grey(c(10:1)/10),main="Mean across all \ncost constraints",xlab="",
     ylab="",xaxt="n",yaxt="n",legend=F)
text(x=rep(c(0.5,1.5,2.5,3.5), times=3), y=rep(c(2.5,1.5,0.5), each=4), labels=round(getValues(r),digits=2)*100, col=ifelse(round(getValues(r),digits=2)>0.5,"white","black"))
axis(1,at=c(0.5,1.5,2.5,3.5),labels = c("Population", "Diversity", "Genome","Multi-obj."),las=2, line=0, lwd=0)
axis(2,at=c(2.5,1.5,0.5),labels = c("Population", "Diversity", "Genome"),las=2, lwd=0)
title(ylab="Loss of success for (%)", line=5)
title(xlab="Optimized for:", line=5)

par(mar=c(7, 2, 3, 1))
r <- raster(extent(0,4,0,3),res=1)
r <- setValues(r,subopt_mat2)
image(r,col=grey(c(10:1)/10),main="$2m cost \nconstraint",xlab="",
     ylab="",xaxt="n",yaxt="n",legend=F)
text(x=rep(c(0.5,1.5,2.5,3.5), times=3), y=rep(c(2.5,1.5,0.5), each=4), labels=round(getValues(r),digits=2)*100, col=ifelse(round(getValues(r),digits=2)>0.5,"white","black"))
axis(1,at=c(0.5,1.5,2.5,3.5),labels = c("Population", "Diversity", "Genome","Multi-obj."),las=2, line=0, lwd=0)
title(xlab="Optimized for:", line=5)

dev.off()

###############
# Plot out the optimal management actions for each cost constraint


optvals <- list()

best <- list()
temp <- sapply(actions$std,function(t) best[[t]] <<- list())

o=objectives[1]
for(o in objectives){
  best[[o]] <- list()
  a=5
  for(a in 1:length(actions$std)){
    best[[o]][[actions$label[a]]] <- sapply(1:length(cts),function(t) optims[[o]][[t]]$par[a] )
  }
  #best[[o]][[actions[5]]] <- rahs[[o]]
  optvals[[o]] <- sapply(1:length(cts),function(t) optims[[o]][[t]]$value )
}

best[["Multi-obj"]] <- list()
for(a in 1:length(actions$std)){
  best[["Multi-obj"]][[actions$label[a]]] <- sapply(1:length(cts),function(t) optims_mo[[t]]$par[a] )
}



# plot for each objective

# layout(matrix(1:3,nrow=1))
# o=1
# for(o in 1:length(objectives)){
#   ndx <- optvals[[objectives[o]]]<10
#   plot(jitter(best[[actions[1]]][[objectives[o]]][ndx],10)~jitter(cts[ndx]),xlab="Cost constraint (million USD)",
#        ylab="Quasi-optimal management action (std)",type="l",col="blue",lwd=2,
#        ylim=c(-0.5,3),xlim=c(1,4),main=objectives[o])
#   points(jitter(cts[ndx]),jitter(best[[actions[2]]][[objectives[o]]][ndx],10),col="red",type="l",lwd=2)
#   points(jitter(cts[ndx]),jitter(best[[actions[3]]][[objectives[o]]][ndx],10),col="green",type="l",lwd=2)
#   points(jitter(cts[ndx]),jitter(best[[actions[4]]][[objectives[o]]][ndx],10),col="black",type="l",lwd=2)
#   points(jitter(cts[ndx]),jitter(best[[actions[5]]][[objectives[o]]][ndx],1),col=gray(0.6),type="l",lwd=2)
#   legend("topleft",lty=c(1,1,1,1,1),lwd=c(1,1,1,1,1),col=c("blue","red","green","black",gray(0.6)),legend=actions,bty="n")
#   
# }


# plot for each objective 

setwd(FiguresDir)
filename=sprintf("QuasiOptActions_%s.pdf", Sys.Date())
pdf(filename, width=7, height=5)
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.5, 0), las=1, tcl=-0.2)
layout(matrix(1:6,nrow=2,byrow=T), widths=c(2.33, 2.33, 2.33))

#Program Duration
par(mar=c(3, 4, 1, 1), mgp=c(2.5, 0.4, 0))
ndx <- optvals[[objectives[1]]]<6
plot(jitter(best[[objectives[1]]][[4]][ndx])~jitter(cts[ndx]),xlab="", ylab="Quasi-optimal management action",type="l",col="firebrick2",lwd=2, ylim=c(0,1.7),xlim=c(1,4),main=actions$label[4], bty="l", xaxt="n", yaxt="n")
axis(2,at=seq(optvec_l[actions$std[4]],optvec_u[actions$std[4]],length=5),
     labels = seq(optvec_l[[actions$std[4]]]*(max(dat[[actions$unstd[4]]])-min(dat[[actions$unstd[4]]]))+min(dat[[actions$unstd[4]]]),optvec_u[[actions$std[4]]]*(max(dat[[actions$unstd[4]]])-min(dat[[actions$unstd[4]]]))+min(dat[[actions$unstd[4]]]),length=5))
points(jitter(cts[ndx]),jitter(best[[objectives[2]]][[4]][ndx])+0.02,col="darkorange",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[[objectives[3]]][[4]][ndx]),col="gold1",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[["Multi-obj"]][[4]][ndx])-0.01,col="black",type="l",lty=2,lwd=2)

#Translocation Age
par(mar=c(3, 4, 1, 1), mgp=c(2.5, 0.4, 0))
ndx <- optvals[[objectives[1]]]<6
plot(jitter(best[[objectives[1]]][[3]][ndx])~jitter(cts[ndx]),xlab="", ylab="",type="l",col="firebrick2",lwd=2, ylim=c(-0.25,1.7),xlim=c(1,4),main=actions$label[3], bty="l", yaxt="n", xaxt="n")
axis(2,at=seq(optvec_l[actions$std[3]],optvec_u[actions$std[3]],length=5),
     labels = seq(optvec_l[[actions$std[3]]]*(max(dat[[actions$unstd[3]]])-min(dat[[actions$unstd[3]]]))+min(dat[[actions$unstd[3]]]),optvec_u[[actions$std[3]]]*(max(dat[[actions$unstd[3]]])-min(dat[[actions$unstd[3]]]))+min(dat[[actions$unstd[3]]]),length=5))
points(jitter(cts[ndx]),jitter(best[[objectives[2]]][[3]][ndx])+0.02,col="darkorange",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[[objectives[3]]][[3]][ndx]),col="gold1",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[["Multi-obj"]][[3]][ndx])-0.01,col="black",type="l",lty=2,lwd=2)

plot.new()
legend("center",lty=c(1,1,1,2),lwd=c(2,2,2,2),col=c("firebrick2","darkorange","gold1","black"), legend=c("Population", "Diversity", "Genome","Multi-obj."),bty="n", cex=1.5)

#Sex Ratio
par(mar=c(3, 4, 1, 1), mgp=c(2.5, 0.4, 0))
ndx <- optvals[[objectives[1]]]<6
plot(jitter(best[[objectives[1]]][[1]][ndx])~jitter(cts[ndx]),xlab="", ylab="Quasi-optimal management action",type="l",col="firebrick2",lwd=2, ylim=c(-0.5,1.7),xlim=c(1,4),main=actions$label[1], bty="l", yaxt="n")
axis(2,at=seq(optvec_l[actions$std[1]],optvec_u[actions$std[1]],length=4),
     labels = seq(min(dat[[actions$unstd[1]]]),optvec_u[[actions$std[1]]]*(max(dat[[actions$unstd[1]]])-min(dat[[actions$unstd[1]]]))+min(dat[[actions$unstd[1]]]),length=4))
title(xlab="Cost constraint (million USD)", mgp=c(2, 0.4, 0))
points(jitter(cts[ndx]),jitter(best[[objectives[2]]][[1]][ndx])+0.02,col="darkorange",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[[objectives[3]]][[1]][ndx]),col="gold1",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[["Multi-obj"]][[1]][ndx])-0.01,col="black",type="l",lty=2,lwd=2)

#Number of breeders
par(mar=c(3, 4, 1, 1), mgp=c(2.5, 0.4, 0))
ndx <- optvals[[objectives[1]]]<6
plot(jitter(best[[objectives[1]]][[2]][ndx])~jitter(cts[ndx]),xlab="", ylab="",type="l",col="firebrick2",lwd=2, ylim=c(0,1.7),xlim=c(1,4),main=actions$label[2], bty="l", yaxt="n")
axis(2,at=seq(optvec_l[actions$std[2]],optvec_u[actions$std[2]],length=5),
     labels = seq(optvec_l[[actions$std[2]]]*(max(dat[[actions$unstd[2]]])-min(dat[[actions$unstd[2]]]))+min(dat[[actions$unstd[2]]]),optvec_u[[actions$std[2]]]*(max(dat[[actions$unstd[2]]])-min(dat[[actions$unstd[2]]]))+min(dat[[actions$unstd[2]]]),length=5))
title(xlab="Cost constraint (million USD)", mgp=c(2, 0.4, 0))
points(jitter(cts[ndx]),jitter(best[[objectives[2]]][[2]][ndx])+0.02,col="darkorange",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[[objectives[3]]][[2]][ndx]),col="gold1",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[["Multi-obj"]][[2]][ndx])-0.01,col="black",type="l",lty=2,lwd=2)

#Direct release
par(mar=c(3, 4, 1, 1), mgp=c(2.5, 0.4, 0))
ndx <- optvals[[objectives[1]]]<6
plot(jitter(best[[objectives[1]]][[5]][ndx])~jitter(cts[ndx]),xlab="", ylab="",type="l",col="firebrick2",lwd=2, ylim=c(0,1.7),xlim=c(1,4),main=actions$label[5], bty="l", yaxt="n")
axis(2,at=seq(optvec_l[actions$std[5]],optvec_u[actions$std[5]],length=5),
     labels = seq(optvec_l[[actions$std[5]]]*(max(dat[[actions$unstd[5]]])-min(dat[[actions$unstd[5]]]))+min(dat[[actions$unstd[5]]]),optvec_u[[actions$std[5]]]*(max(dat[[actions$unstd[5]]])-min(dat[[actions$unstd[5]]]))+min(dat[[actions$unstd[5]]]),length=5))
points(jitter(cts[ndx]),jitter(best[[objectives[2]]][[5]][ndx])+0.02,col="darkorange",type="l",lwd=2)
title(xlab="Cost constraint (million USD)", mgp=c(2, 0.4, 0))
points(jitter(cts[ndx]),jitter(best[[objectives[3]]][[5]][ndx]),col="gold1",type="l",lwd=2)
points(jitter(cts[ndx]),jitter(best[["Multi-obj"]][[5]][ndx])-0.01,col="black",type="l",lty=2,lwd=2)

dev.off()









################
# OLD CODE
################

mod1.1 <- gam(pop.norm ~ s(release.age,sex.ratio,k=4),data=dat)  # minimize gcv score
plot(mod1.1)
plot(mod1.1,residuals = T,pch=19)
vis.gam(mod1.1)
gam.check(mod1.1)
plot(fitted(mod1.1),residuals(mod1.1))

## check for unmodelled pattern in the residuals 
rsd <- residuals(mod1.1,type="deviance")
gam(rsd~s(release.age,sex.ratio,k=10)-1,data=dat,select=TRUE)



mod1.2 <- gam(pop.norm ~ te(release.age.s,release.end.s,sex.ratio, k=c(3,3,3)),data=dat)  # minimize gcv score
mod1.3 <- gam(pop.norm ~ release.adult.hybrids.q + te(release.age.s,release.end.s,sex.ratio,Corral.capacity.s, 
                                                      k=c(3,3,3,3)),data=dat,select=T,method="ML")  # minimize gcv score
mod1.4 <- gam(pop.norm ~ Corral.capacity.s + te(release.age.s,release.end.s,sex.ratio, k=c(4,4,3),bs="cs",by=release.adult.hybrids.q),
              data=dat,select=T,method="ML")  # minimize gcv score

plot(mod1.4,shade=T,scale=0)
plot(mod1.4,residuals = T,pch=19)
vis.gam(mod1.4,view=c("release.age.s","release.end.s"),se=1)
vis.gam(mod1.4,view=c("release.age.s","sex.ratio"),se=1)
gam.check(mod1.4)
plot(fitted(mod1.1),residuals(mod1.1))
summary(mod1.4)
anova(mod1.4)

AIC(mod1.4,mod1.1)


mod2.1 <- gam(pop.norm ~ s(release.age.s,k=4),data=dat)  # minimize gcv score
mod2.1 <- gam(pop.norm ~ s(release.end.s,k=4),data=dat)
mod2.1 <- gam(pop.norm ~ Corral.capacity.s,data=dat)

mod2.2 <- gam(pop.norm ~ s(release.age.s,release.end.s,sex.ratio,k=15),data=dat)  # minimize gcv score
#mod2.3 <- gam(pop.norm ~ release.adult.hybrids.q + s(release.age.s,release.end.s,sex.ratio,k=15,bs="cr"),data=dat)  # minimize gcv score


mod2.4 <- gam(pop.norm ~ Corral.capacity.s + s(release.age.s,release.end.s,sex.ratio, k=30,bs="ts",by=release.adult.hybrids.q),
              data=dat,method="ML",select=T)  # minimize gcv score

mod2.5 <- gam(pop.norm ~ release.adult.hybrids.q + s(release.age.s,k=4) + s(release.end.s,k=4) 
              + s(sex.ratio,k=3) + Corral.capacity.s,
              data=dat,method="ML",select=T)  # minimize gcv score

summary(mod1)
summary(mod2)
summary(mod1.2)
summary(mod2.2)
summary(mod1.3)
summary(mod2.3)

summary(mod2.4) 


AIC(mod2.3,mod2.4,mod2.5)

names(dat)

gmod2.3 <- gam(sar.norm ~ s(release.age.s,release.end.s,sex.ratio,Corral.capacity.s, k=30,bs="ts",by=release.adult.hybrids.q),
               data=dat,method="ML",select=T)
summary(gmod2.3)

gmod2.4 <- gam(sar.norm ~ Corral.capacity.s + s(release.age.s,release.end.s,sex.ratio, k=30,bs="ts",by=release.adult.hybrids.q),
               data=dat,method="ML",select=T)

names(dat)

fmod2.3 <- gam(qval.norm ~ s(release.age.s,release.end.s,sex.ratio,Corral.capacity.s, k=30,bs="ts",by=release.adult.hybrids.q),
               data=dat,method="ML",select=T)
summary(fmod2.3)

fmod2.4 <- gam(qval.norm ~ Corral.capacity.s + s(release.age.s,release.end.s,sex.ratio, k=30,bs="ts",by=release.adult.hybrids.q),
               data=dat,method="ML",select=T)

summary(fmod2.4)


cmod2.4 <- gam(cost.norm ~ Corral.capacity.s + s(release.age.s,release.end.s, k=20,bs="ts",by=release.adult.hybrids.q),
               data=dat,method="ML",select=T)
summary(cmod2.4)   

summary(cmod2.4)

vis.gam(cmod2.3,color="heat",plot.type = "persp",view=c("release.end.s","release.age.s"),se=1)
vis.gam(cmod2.3,color="heat",plot.type = "persp",view=c("release.end.s","Corral.capacity.s"))


cmod2.5 <- gam(cost.norm ~ release.adult.hybrids.q + s(release.age.s,k=4) + s(release.end.s,k=4) 
               + s(sex.ratio,k=3) + Corral.capacity.s,
               data=dat,method="ML",select=T)  # minimize gcv score
summary(cmod2.5)


#### old versions

multi_objfunc <- function(optvec,release.adult.hybrids=TRUE,cost.threshold=1){
  
  release.adult.hybrids.q<-factor(as.numeric(release.adult.hybrids),
                                  levels=levels(dat$release.adult.hybrids.q))
  newdata <- data.frame(
    release.end.s = optvec['release.end.s'],
    release.age.s = optvec['release.age.s'],
    release.adult.hybrids.q=release.adult.hybrids.q,
    sex.ratio.s=optvec['sex.ratio.s'],
    Corral.capacity.s = optvec['Corral.capacity.s']
  )
  values <- sapply(1:3,function(t) max(0.01,predict(best.mods[[t]],newdata=newdata)))
  values <- mean(values)   # could apply weighting scheme here
  value2 <- 1/values
  cost <- costfunc(optvec['release.end.s'],optvec['release.age.s'],release.adult.hybrids.q,optvec['Corral.capacity.s'])
  if(cost>cost.threshold) value2 = value2+20
  # if((optvec['release.end.s']<0.5)|(optvec['release.end.s']>1.5)) value2 = value2+2*abs(optvec['release.end.s']-0.5)
  # if((optvec['release.age.s']<0.5)|(optvec['release.age.s']>1.5)) value2 = value2+2*abs(optvec['release.age.s']-0.5)
  # if((optvec['Corral.capacity.s']<0.5)|(optvec['Corral.capacity.s']>1.5)) value2 = value2+5*abs(optvec['Corral.capacity.s']-0.5)
  # if((optvec['sex.ratio.s']<0.5)|(optvec['sex.ratio.s']>1.5)) value2 = value2+2*abs(optvec['sex.ratio.s']-0.5)
  return(value2)
} 

optim(par=init.parms,fn=objfunc,method="Nelder-Mead",cost.threshold=ct,
      release.adult.hybrids=FALSE,objective=obj,control=list(maxit=10000))


