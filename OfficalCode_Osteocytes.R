############################### ANALYSES FOR OSTEOCYTES LACUNAE PAPER #############################

######## I PERFORM THE FOLLOWING 5 ANALYSES ON THIS SCRIPT:
########### 1) 4-way ANOVA for each Lacunae Metric 
########### 2) Multilevel Linear Model for each Lacunae Metric (that was significant in anova)
########### 3) Semivariograms by mouse, loading, and strain 
########### 4) Semivariogram by strain 
########### 4) Semivariogram by strain and loading 




#loading necessary packages
library(car)
library("geoR")
library("gstat")
library("sp")
library("fields")
library(dplyr)
library(stringr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(base)
library(lme4)
library(lattice)
setwd("/Users/sandratredinnick/Downloads")


#LOADING IN DATA AND CONSTRUCTING DATAFRAME APPROPRIATELY
#read in data
lacundf= read.csv('FullDataset_OsteocyteCharacteristics.csv')
#make id column
df2 = lacundf %>% mutate(str_extract(Sample,"[0-9]+"))
lacundf$id = df2$`str_extract(Sample,"[0-9]+")`
df2 = lacundf %>% mutate(str_extract(Sample,"[0-9]+"))
lacundf$id = df2[,21]
#make location column
lacundf$location <- ifelse(grepl("PL", lacundf$Sample), "PL", "AM")
#make group column 
lacundf$group=ifelse(lacundf$id==252|lacundf$id==233|
                       lacundf$id==238|lacundf$id==250,0,1) #1 if flight
#make categorical predictor varible for linear models with 4 levels
###value 1 = AM and flight 
df111= lacundf%>% filter(group == 1 & location =="AM")
df111$Interaction = rep("AM_Flight",16703)
###value 2 = PL and flight 
df222= lacundf%>% filter(group == 1 & location =="PL")
df222$Interaction = rep("PL_Flight",15947)
###value 3 = AM and ground
df333= lacundf%>% filter(group == 0 & location =="AM")
df333$Interaction = rep("AM_Ground",20359)
###value 4 = PL and ground
df444= lacundf%>% filter(group == 0 & location =="PL")
df444$Interaction = rep("PL_Ground",15776)
lacundf <- rbind(df111, df222,df333,df444)
#code variables as factors 
lacundf$id=as.factor(lacundf$id)
lacundf$group=as.factor(lacundf$group)
lacundf$Interaction=as.factor(lacundf$Interaction)
lacundf$location=as.factor(lacundf$location)




################ 4-WAY ANOVA WITH BLOCKING VARIABLE ###############
#NOTE: transformations were performed to help normalize the residuals, allowing us to use anova. tranformatins were chosen
#based on which model most closely alligned with the assumptions of anova 

#the assumptions that need to be checked are:
##### 1) normality (dont use shapiro wilk because too sensitive to small deviations of normality)
##### 2) homogeneity (dont use levenes test because too sensitive to small deviations of homogeneity)
##### 3) outliers



####LOG VOLUME
rcbdVOL = lm(log(Volume..µm..) ~location*group + id, data = lacundf)
summary(rcbdVOL)
anova(rcbdVOL)
#check assumptions 
#1: homogeneity
plot(rcbdVOL, 1) 
#2: normality 
plot(rcbdVOL, 2)
ggqqplot(residuals(rcbdVOL))
#3: outliers 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logvol)
as.data.frame(out)# a few extreme but not bad
lacundf$logvol = log(lacundf$Volume..µm..) 


####RANK OBLATENESS
rcbdOB = lm(rank(Lacuna.Oblateness..Untested...NA.) ~location*group + id, data = lacundf)
summary(rcbdOB)
anova(rcbdOB)
#check assumptions 
#1: homogeneity
plot(rcbdOB, 1) 
#2: normality 
plot(rcbdOB, 2) #normality 
ggqqplot(residuals(rcbdOB))
#3: outliers 
lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(rankob) #no extreme 
lacundf$rankob= rank(lacundf$Lacuna.Oblateness..Untested...NA.)



####LOG SURFACE AREA 
rcbdSA = lm(log(Surface.Area..Lindblad.2005...µm..) ~location*group + id, data = lacundf)
summary(rcbdSA)
anova(rcbdSA)
#check assumptions 
#1: homogeneity
plot(rcbdSA, 1) 
#2: normality 
plot(rcbdSA, 2)
ggqqplot(residuals(rcbdSA))
#3: outliers 
lacundf$logSA = log(lacundf$Surface.Area..Lindblad.2005...µm..) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logSA)
as.data.frame(out)# one extreme 



####LOG ASPECT RATIO 
rcbdAR = lm(log(Aspect.Ratio..NA.) ~location*group + id, data = lacundf)
summary(rcbdAR)
anova(rcbdAR)
#check assumptions 
#1: homogeneity
plot(rcbdAR, 1) 
#2: normality 
plot(rcbdAR, 2)
ggqqplot(residuals(rcbdAR))
#3: outliers 
lacundf$logAR = log(lacundf$Aspect.Ratio..NA.) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logAR)
as.data.frame(out)# one extreme 




####RANK STRETCH 
rcbdST = lm(rank(Lacuna.Stretch..Untested...NA.) ~location*group + id, data = lacundf)
summary(rcbdST)
anova(rcbdST)
#check assumptions 
#1: homogeneity
plot(rcbdST, 1) 
#2: normality 
plot(rcbdST, 2)
ggqqplot(residuals(rcbdST))
#3: outliers 
lacundf$rankST = rank(lacundf$Lacuna.Stretch..Untested...NA.) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(rankST)
as.data.frame(out)# no extreme 




####LOG WIDTH 
rcbdWI = lm(log(Lacuna.Width..Untested...µm.) ~location*group + id, data = lacundf)
summary(rcbdWI)
anova(rcbdWI)
#check assumptions 
#1: homogeneity
plot(rcbdWI, 1) 
#2: normality 
plot(rcbdWI, 2)
ggqqplot(residuals(rcbdWI))
#3: outliers 
lacundf$logWI = log(lacundf$Lacuna.Width..Untested...µm.) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logWI)
as.data.frame(out) # few extreme  



####RANK FITTED ELLIPSOID VOLUME 
rcbdEVOL = lm(rank(Fitted.Ellipsoid.Volume..Untested...µm..) ~location*group + id, data = lacundf)
summary(rcbdEVOL)
anova(rcbdEVOL)
#check assumptions 
#1: homogeneity
plot(rcbdEVOL, 1) 
#2: normality 
plot(rcbdEVOL, 2)
ggqqplot(residuals(rcbdEVOL))
#3: outliers 
lacundf$rankEVOL = rank(lacundf$Fitted.Ellipsoid.Volume..Untested...µm..) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(rankEVOL)
as.data.frame(out)#no extreme 



####LOG LENGTH 
rcbdLE = lm(log(Lacuna.Length..Untested...µm.) ~location*group + id, data = lacundf)
summary(rcbdLE)
anova(rcbdLE)
#check assumptions 
#1: homogeneity
plot(rcbdLE, 1) 
#2: normality 
plot(rcbdLE, 2)
ggqqplot(residuals(rcbdLE))
#3: outliers 
lacundf$logLE = log(lacundf$Lacuna.Length..Untested...µm.) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logLE)
as.data.frame(out)# a few extreme



####LOG HEIGHT 
rcbdHE = lm(log(Lacuna.Height..Untested...µm.) ~location*group + id, data = lacundf)
summary(rcbdHE)
anova(rcbdHE)
#check assumptions 
#1: homogeneity
plot(rcbdHE, 1) 
#2: normality 
plot(rcbdHE, 2)
ggqqplot(residuals(rcbdHE))
#3: outliers 
lacundf$logHE = log(lacundf$Lacuna.Height..Untested...µm.) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logHE)
as.data.frame(out)# a few extreme 



####LOG FITTED ELLIPSOID SURFACE AREA 
rcbdESA = lm(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..) ~location*group + id, data = lacundf)
summary(rcbdESA)
anova(rcbdESA)
#check assumptions 
#1: homogeneity
plot(rcbdESA, 1) 
#2: normality 
plot(rcbdESA, 2)
ggqqplot(residuals(rcbdESA))
#3: outliers 
lacundf$logESA = log(lacundf$Fitted.Ellipsoid.Surface.Area..Untested...µm..) 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(logESA)
as.data.frame(out)# one extreme 



####SPHERICITY 
rcbdSPH = lm(Sphericity..Untested...NA. ~location*group + id, data = lacundf)
summary(rcbdSPH)
anova(rcbdSPH)
#check assumptions 
#1: homogeneity
plot(rcbdSPH, 1) 
#2: normality 
plot(rcbdSPH, 2)
ggqqplot(residuals(rcbdSPH))
#3: outliers 
out= lacundf %>%
  group_by(location,group,id) %>%
  identify_outliers(Sphericity..Untested...NA.)
as.data.frame(out)# a few extreme 


















################ MULTILEVEL LINEAR MODELING ###############


#####LOG VOLUME
###one: varying intercept 
lmevolume1G = lmer(log(Volume..µm..)~Interaction+(1|id),data=lacundf)
summary(lmevolume1G)
display(lmevolume1G)
confint(lmevolume1G,level=0.95)
###two: varying intercept and slope , DIVERGENT... 
lmevolume2G = lmer(log(Volume..µm..)~Interaction+(1+Interaction|id),data=lacundf)
summary(lmevolume2G)
display(lmevolume2G)
confint(lmevolume2G,level=0.95)
#check assumptions 
### 1: linearity
lin<-plot(resid(lmevolume1G),lacundf$logvol)
### 2: homogeneity
plot(lmevolume1G)
### 3: normality 
qqmath(lmevolume1G, id = 0.05)



#####RANK OBLATENESS
###one: varying intercept 
lmeRO = lmer(rank(Lacuna.Oblateness..Untested...NA.)~Interaction+(1|id),data=lacundf)
summary(lmeRO)
display(lmeRO)
confint(lmeRO,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeRO),lacundf$rankob)
### 2: homogeneity
plot(lmeRO)
### 3: normality 
qqmath(lmeRO, id = 0.05)




#####LOG SURFACE AREA 
###one: varying intercept 
lmeSA = lmer(log(Surface.Area..Lindblad.2005...µm..)~Interaction+(1|id),data=lacundf)
summary(lmeSA)
display(lmeSA)
confint(lmeSA,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeSA),lacundf$logSA)
### 2: homogeneity
plot(lmeSA)
### 3: normality 
qqmath(lmeSA, id = 0.05)



#####LOG ASPECT RATIO
###one: varying intercept 
lmeAR = lmer(log(Aspect.Ratio..NA.)~Interaction+(1|id),data=lacundf)
summary(lmeAR)
display(lmeAR)
confint(lmeAR,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeAR),lacundf$logAR)
### 2: homogeneity
plot(lmeAR)
### 3: normality 
qqmath(lmeAR, id = 0.05)





#####RANK STRETCH
###one: varying intercept 
lmeST = lmer(rank(Lacuna.Stretch..Untested...NA.)~Interaction+(1|id),data=lacundf)
summary(lmeST)
display(lmeST)
confint(lmeST,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeST),lacundf$rankST)
### 2: homogeneity
plot(lmeST)
### 3: normality 
qqmath(lmeST, id = 0.05)




#####LOG WIDTH
###one: varying intercept 
lmeWI = lmer(log(Lacuna.Width..Untested...µm.)~Interaction+(1|id),data=lacundf)
summary(lmeWI)
display(lmeWI)
confint(lmeWI,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeWI),lacundf$logWI)
### 2: homogeneity
plot(lmeWI)
### 3: normality 
qqmath(lmeWI, id = 0.05)




#####RANK FITTED ELLIPSOID VOLUME
###one: varying intercept 
lmeEV = lmer(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~Interaction+(1|id),data=lacundf)
summary(lmeEV)
display(lmeEV)
confint(lmeEV,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeEV),lacundf$logEVOL)
### 2: homogeneity
plot(lmeEV)
### 3: normality 
qqmath(lmeEV, id = 0.05)




#####LOG LENGTH
###one: varying intercept 
lmeLE = lmer(log(Lacuna.Length..Untested...µm.)~Interaction+(1|id),data=lacundf)
summary(lmeLE)
display(lmeLE)
confint(lmeLE,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeLE),lacundf$logLE)
### 2: homogeneity
plot(lmeLE)
### 3: normality 
qqmath(lmeLE, id = 0.05)





#####LOG HEIGHT
###one: varying intercept 
lmeHE = lmer(log(Lacuna.Height..Untested...µm.)~Interaction+(1|id),data=lacundf)
summary(lmeHE)
display(lmeHE)
confint(lmeHE,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeHE),lacundf$logHE)
### 2: homogeneity
plot(lmeHE)
### 3: normality 
qqmath(lmeHE, id = 0.05)





#####LOG FITTED ELLIPSOID SURFACE AREA 
###one: varying intercept 
lmeFSA = lmer(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~Interaction+(1|id),data=lacundf)
summary(lmeFSA)
display(lmeFSA)
confint(lmeFSA,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeFSA),lacundf$logFSA)
### 2: homogeneity
plot(lmeFSA)
### 3: normality 
qqmath(lmeFSA, id = 0.05)




#####SPHERICITY
###one: varying intercept 
lmeSP = lmer(Sphericity..Untested...NA.~Interaction+(1|id),data=lacundf)
summary(lmeSP)
display(lmeSP)
confint(lmeSP,level=0.95)
#checking assumptions 
### 1: linearity
lin<-plot(resid(lmeSP),lacundf$Sphericity..Untested...NA.)
### 2: homogeneity
plot(lmeSP)
### 3: normality 
qqmath(lmeSP, id = 0.05)































################ SEMIVARIOGRAM BY ID, LOADING, AND STRAIN ###############


######################MOUSE 233
lacundfGPL233= lacundf%>%filter(id=="233"&
                                  location =="PL"&
                                  group==0)
lacundfGAM233= lacundf%>%filter(id=="233"&
                                  location =="AM"&
                                  group==0)
coordinates(lacundfGPL233)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM233)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
quilt.plot(x=lacundfGPL233$Center.Of.Mass.X..µm., y=lacundfGPL233$Center.Of.Mass.Y..µm., 
           z= log(lacundfGPL233$Volume..µm..), main="Heat Map of the Log Volume of Osteocytes for Mouse 233 in the PL Region",xlab='Center of Mass X',ylab='Center of Mass Y')


####LOG VOLUME
#PL
v.cloudvol233PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL233, cloud = F)
plot(v.cloudvol233PL, main = "Variogram - Log Volume: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 20 first x axis point
#sill = .70 #first maximum 
#range = 205 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelmess = vgm(psill = .17, nugget=25 , range= 200, model= "Sph")
#now fit
fitmess = fit.variogram(object= v.cloudvol233PL, model =modelmess )
plot(v.cloudvol233PL, model=fitmess)
#AM
v.cloudvol233AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM233, cloud = F)
plot(v.cloudvol233AM, main = "Variogram - Log Volume: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelvol233AM = vgm(psill = .18, nugget=25 , range= 210, model= "Exp")
modelvol233AMsph = vgm(psill = .18, nugget=25 , range= 210, model= "Sph")
#now fit
fitvol233AM = fit.variogram(object= v.cloudvol233AM, model =modelvol233AMsph )
plot(v.cloudvol233AM, model=fitvol233AM)



####RANK OBLATENESS
#PL
v.cloud233PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLOB, main = "Variogram - Rank Oblateness: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB233PLsph = vgm(psill = 7e5, nugget=15 , range= 200, model= "Sph")
#now fit
fitOB233PL = fit.variogram(object= v.cloud233PLOB, model =modelOB233PLsph)
plot(v.cloud233PLOB, model=fitOB233PL)
#AM
v.cloud233AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMOB, main = "Variogram - Rank Oblateness: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB233AMsph = vgm(psill = 1201000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB233AM = fit.variogram(object= v.cloud233AMOB, model =modelOB233AMsph)
plot(v.cloud233AMOB, model=fitOB233AM)



####LOG SURFACE AREA
#PL
v.cloud233PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLSA, main = "Variogram - Log Surface Area: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA233PLsph = vgm(psill = .13, nugget=25 , range= 200, model= "Sph")
#now fit
fitSA233PL = fit.variogram(object= v.cloud233PLSA, model =modelSA233PLsph)
plot(v.cloud233PLSA, model=fitSA233PL)
#AM
v.cloud233AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMSA, main = "Variogram - Log Surface Area: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA233AMsph = vgm(psill = .13, nugget=25 , range= 205, model= "Sph")
#now fit
fitSA233AM = fit.variogram(object= v.cloud233AMSA, model =modelSA233AMsph)
plot(v.cloud233AMSA, model=fitSA233AM)



####LOG ASPECT RATIO
#PL
v.cloud233PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR233PLsph = vgm(psill = .20, nugget=15 , range= 200, model= "Sph")
#now fit
fitAR233PL = fit.variogram(object= v.cloud233PLAR, model =modelAR233PLsph)
plot(v.cloud233PLAR, model=fitAR233PL)
#AM
v.cloud233AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR233AMsph = vgm(psill = .24, nugget=15 , range= 200, model= "Sph")
#now fit
fitAR233AM = fit.variogram(object= v.cloud233AMAR, model =modelAR233AMsph)
plot(v.cloud233AMAR, model=fitAR233AM)




####RANK STRETCH
#PL
v.cloud233PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLST, main = "Variogram - Rank Stretch: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST233PLsph = vgm(psill = 7e5, nugget=15 , range= 200, model= "Sph")
#now fit
fitST233PL = fit.variogram(object= v.cloud233PLST, model =modelST233PLsph)
plot(v.cloud233PLST, model=fitST233PL)
#AM
v.cloud233AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMST, main = "Variogram - Rank Stretch: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST233AMsph = vgm(psill = 1205000, nugget=15 , range= 200, model= "Sph")
#now fit
fitST233AM = fit.variogram(object= v.cloud233AMST, model =modelST233AMsph)
plot(v.cloud233AMST, model=fitST233AM)


####LOG WIDTH
#PL
v.cloud233PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLWI, main = "Variogram - Log Width: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI233PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitWI233PL = fit.variogram(object= v.cloud233PLWI, model =modelWI233PLsph)
plot(v.cloud233PLWI, model=fitWI233PL)
#AM
v.cloud233AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMWI, main = "Variogram - Log Width: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI233AMsph = vgm(psill = .065, nugget=10 , range= 210, model= "Sph")
#now fit
fitWI233AM = fit.variogram(object= v.cloud233AMWI, model =modelWI233AMsph)
plot(v.cloud233AMWI, model=fitWI233AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud233PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV233PLsph = vgm(psill = 7e5, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV233PL = fit.variogram(object= v.cloud233PLFEV, model =modelFEV233PLsph)
plot(v.cloud233PLFEV, model=fitFEV233PL)
#AM
v.cloud233AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV233AMsph = vgm(psill = 12001000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV233AM = fit.variogram(object= v.cloud233AMFEV, model =modelFEV233AMsph)
plot(v.cloud233AMFEV, model=fitFEV233AM)



####LOG LENGTH
#PL
v.cloud233PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLLE, main = "Variogram - Log Length: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE233PLsph = vgm(psill = .04, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE233PL = fit.variogram(object= v.cloud233PLLE, model =modelLE233PLsph)
plot(v.cloud233PLLE, model=fitLE233PL)
#AM
v.cloud233AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMLE, main = "Variogram - Log Length: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE233AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE233AM = fit.variogram(object= v.cloud233AMLE, model =modelLE233AMsph)
plot(v.cloud233AMLE, model=fitLE233AM)


####LOG HEIGHT
#PL
v.cloud233PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLHE, main = "Variogram - Log Height: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE233PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE233PL = fit.variogram(object= v.cloud233PLHE, model =modelHE233PLsph)
plot(v.cloud233PLHE, model=fitHE233PL)
#AM
v.cloud233AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMHE, main = "Variogram - Log Height: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE233AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE233AM = fit.variogram(object= v.cloud233AMHE, model =modelHE233AMsph)
plot(v.cloud233AMHE, model=fitHE233AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud233PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA233PLsph = vgm(psill = .10, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA233PL = fit.variogram(object= v.cloud233PLFSA, model =modelFSA233PLsph)
plot(v.cloud233PLFSA, model=fitFSA233PL)
#AM
v.cloud233AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA233AMsph = vgm(psill = .12, nugget=10 , range= 210, model= "Sph")
#now fit
fitFSA233AM = fit.variogram(object= v.cloud233AMFSA, model =modelFSA233AMsph)
plot(v.cloud233AMFSA, model=fitFSA233AM)


####SPHERICITY
#PL
v.cloud233PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL233, cloud = F)
plot(v.cloud233PLFSP, main = "Variogram - Sphericity: PL Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH233PLsph = vgm(psill = .0030, nugget=10 , range= 200, model= "Sph")
#now fit
fitSPH233PL = fit.variogram(object= v.cloud233PLFSP, model =modelSPH233PLsph)
plot(v.cloud233PLSPH, model=fitSPH233PL)
#AM
v.cloud233AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM233, cloud = F)
plot(v.cloud233AMFSP, main = "Variogram - Sphericity: AM Mouse 233 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH233AMsph = vgm(psill = .0025, nugget=10 , range= 210, model= "Sph")
#now fit
fitSPH233AM = fit.variogram(object= v.cloud233AMFSP, model =modelSPH233AMsph)
plot(v.cloud233AMFSP, model=fitSPH233AM)
















######################MOUSE 238
lacundfGPL238= lacundf%>%filter(id=="238"&
                                  location =="PL"&
                                  group==0)
lacundfGAM238= lacundf%>%filter(id=="238"&
                                  location =="AM"&
                                  group==0)
coordinates(lacundfGPL238)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM238)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol238PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL238, cloud = F)
plot(v.cloudvol238PL, main = "Variogram - Log Volume: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 20 first x axis point
#sill = .25 #first maximum 
#range = 205 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelfit238PLvol = vgm(psill = .25, nugget=25 , range= 200, model= "Sph")
#now fit
fit238PLvol = fit.variogram(object= v.cloudvol238PL, model =modelfit238PLvol )
plot(v.cloudvol238PL, model=fit238PLvol)
#AM
v.cloudvol238AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM238, cloud = F)
plot(v.cloudvol238AM, main = "Variogram - Log Volume: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol238AMsph = vgm(psill = .25, nugget=15 , range= 205, model= "Sph")
#now fit
fitvol238AM = fit.variogram(object= v.cloudvol238AM, model =modelvol238AMsph )
plot(v.cloudvol238AM, model=fitvol238AM)



####RANK OBLATENESS
#PL
v.cloud238PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLOB, main = "Variogram - Rank Oblateness: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB238PLsph = vgm(psill = 1700000, nugget=15 , range= 200, model= "Sph")
#now fit
fitOB238PL = fit.variogram(object= v.cloud238PLOB, model =modelOB238PLsph)
plot(v.cloud238PLOB, model=fitOB238PL)
#AM
v.cloud238AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMOB, main = "Variogram - Rank Oblateness: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB238AMsph = vgm(psill = 1600000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB238AM = fit.variogram(object= v.cloud238AMOB, model =modelOB238AMsph)
plot(v.cloud238AMOB, model=fitOB238AM)


####LOG SURFACE AREA
#PL
v.cloud238PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLSA, main = "Variogram - Log Surface Area: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA238PLsph = vgm(psill = .13, nugget=15 , range= 205, model= "Sph")
#now fit
fitSA238PL = fit.variogram(object= v.cloud238PLSA, model =modelSA238PLsph)
plot(v.cloud238PLSA, model=fitSA238PL)
#AM
v.cloud238AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMSA, main = "Variogram - Log Surface Area: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA238AMsph = vgm(psill = .15, nugget=15 , range= 200, model= "Sph")
#now fit
fitSA238AM = fit.variogram(object= v.cloud238AMSA, model =modelSA238AMsph)
plot(v.cloud238AMSA, model=fitSA238AM)


####LOG ASPECT RATIO
#PL
v.cloud238PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR238PLsph = vgm(psill = .27, nugget=15 , range= 210, model= "Sph")
#now fit
fitAR238PL = fit.variogram(object= v.cloud238PLAR, model =modelAR238PLsph)
plot(v.cloud238PLAR, model=fitAR238PL)
#AM
v.cloud238AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR238AMsph = vgm(psill = .25, nugget=15 , range= 200, model= "Sph")
#now fit
fitAR238AM = fit.variogram(object= v.cloud238AMAR, model =modelAR238AMsph)
plot(v.cloud238AMAR, model=fitAR238AM)



####RANK STRETCH
#PL
v.cloud238PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLST, main = "Variogram - Rank Stretch: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST238PLsph = vgm(psill = 1700000, nugget=15 , range= 200, model= "Sph")
#now fit
fitST238PL = fit.variogram(object= v.cloud238PLST, model =modelST238PLsph)
plot(v.cloud238PLST, model=fitST238PL)
#AM
v.cloud238AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMST, main = "Variogram - Rank Stretch: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST238AMsph = vgm(psill = 1700000, nugget=15 , range= 200, model= "Sph")
#now fit
fitST238AM = fit.variogram(object= v.cloud238AMST, model =modelST238AMsph)
plot(v.cloud238AMST, model=fitST238AM)


####LOG WIDTH
#PL
v.cloud238PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLWI, main = "Variogram - Log Width: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI238PLsph = vgm(psill = .08, nugget=10 , range= 210, model= "Sph")
#now fit
fitWI238PL = fit.variogram(object= v.cloud238PLWI, model =modelWI238PLsph)
plot(v.cloud238PLWI, model=fitWI238PL)
#AM
v.cloud238AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMWI, main = "Variogram - Log Width: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI238AMsph = vgm(psill = .065, nugget=10 , range= 205, model= "Sph")
#now fit
fitWI238AM = fit.variogram(object= v.cloud238AMWI, model =modelWI238AMsph)
plot(v.cloud238AMWI, model=fitWI238AM)



####Rank Fitted Ellipsoid Volume
#PL
v.cloud238PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV238PLsph = vgm(psill = 16000000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV238PL = fit.variogram(object= v.cloud238PLFEV, model =modelFEV238PLsph)
plot(v.cloud238PLFEV, model=fitFEV238PL)
#AM
v.cloud238AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV238AMsph = vgm(psill = 17000000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV238AM = fit.variogram(object= v.cloud238AMFEV, model =modelFEV238AMsph)
plot(v.cloud238AMFEV, model=fitFEV238AM)



####LOG LENGTH
#PL
v.cloud238PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLLE, main = "Variogram - Log Length: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE238PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE238PL = fit.variogram(object= v.cloud238PLLE, model =modelLE238PLsph)
plot(v.cloud238PLLE, model=fitLE238PL)
#AM
v.cloud238AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMLE, main = "Variogram - Log Length: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE238AMsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE238AM = fit.variogram(object= v.cloud238AMLE, model =modelLE238AMsph)
plot(v.cloud238AMLE, model=fitLE238AM)


####LOG HEIGHT
#PL
v.cloud238PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLHE, main = "Variogram - Log Height: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE238PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE238PL = fit.variogram(object= v.cloud238PLHE, model =modelHE238PLsph)
plot(v.cloud238PLHE, model=fitHE238PL)
#AM
v.cloud238AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMHE, main = "Variogram - Log Height: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE238AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE238AM = fit.variogram(object= v.cloud238AMHE, model =modelHE238AMsph)
plot(v.cloud238AMHE, model=fitHE238AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud238PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA238PLsph = vgm(psill = .12, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA238PL = fit.variogram(object= v.cloud238PLFSA, model =modelFSA238PLsph)
plot(v.cloud238PLFSA, model=fitFSA238PL)
#AM
v.cloud238AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA238AMsph = vgm(psill = .13, nugget=10 , range= 210, model= "Sph")
#now fit
fitFSA238AM = fit.variogram(object= v.cloud238AMFSA, model =modelFSA238AMsph)
plot(v.cloud238AMFSA, model=fitFSA238AM)


####SPHERICITY
#PL
v.cloud238PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL238, cloud = F)
plot(v.cloud238PLFSP, main = "Variogram - Sphericity: PL Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH238PLsph = vgm(psill = .0025, nugget=10 , range= 200, model= "Sph")
#now fit
fitSPH238PL = fit.variogram(object= v.cloud238PLFSP, model =modelSPH238PLsph)
plot(v.cloud238PLFSP, model=fitSPH238PL)
#AM
v.cloud238AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM238, cloud = F)
plot(v.cloud238AMFSP, main = "Variogram - Sphericity: AM Mouse 238 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH238AMsph = vgm(psill = .0020, nugget=10 , range= 200, model= "Sph")
#now fit
fitSPH238AM = fit.variogram(object= v.cloud238AMFSP, model =modelSPH238AMsph)
plot(v.cloud238AMFSP, model=fitSPH238AM)













######################MOUSE 250
lacundfGPL250= lacundf%>%filter(id=="250"&
                                  location =="PL"&
                                  group==0)
lacundfGAM250= lacundf%>%filter(id=="250"&
                                  location =="AM"&
                                  group==0)
coordinates(lacundfGPL250)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM250)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol250PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL250, cloud = F)
plot(v.cloudvol250PL, main = "Variogram - Log Volume: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 20 first x axis point
#sill = .25 #first maximum 
#range = 205 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelfit250PLvol = vgm(psill = .25, nugget=25 , range= 200, model= "Sph")
#now fit
fit250PLvol = fit.variogram(object= v.cloudvol250PL, model =modelfit250PLvol )
plot(v.cloudvol250PL, model=fit250PLvol)
#AM
v.cloudvol250AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM250, cloud = F)
plot(v.cloudvol250AM, main = "Variogram - Log Volume: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol250AMsph = vgm(psill = .20, nugget=20 , range= 210, model= "Sph")
#now fit
fitvol250AM = fit.variogram(object= v.cloudvol250AM, model =modelvol250AMsph )
plot(v.cloudvol250AM, model=fitvol250AM)


####RANK OBLATENESS
#PL
v.cloud250PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLOB, main = "Variogram - Rank Oblateness: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB250PLsph = vgm(psill = 1000000, nugget=15 , range= 200, model= "Sph")
#now fit
fitOB250PL = fit.variogram(object= v.cloud250PLOB, model =modelOB250PLsph)
plot(v.cloud250PLOB, model=fitOB250PL)
#AM
v.cloud250AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMOB, main = "Variogram - Rank Oblateness: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB250AMsph = vgm(psill = 2700000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB250AM = fit.variogram(object= v.cloud250AMOB, model =modelOB250AMsph)
plot(v.cloud250AMOB, model=fitOB250AM)



####LOG SURFACE AREA
#PL
v.cloud250PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLSA, main = "Variogram - Log Surface Area: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA250PLsph = vgm(psill = .11, nugget=15 , range= 200, model= "Sph")
#now fit
fitSA250PL = fit.variogram(object= v.cloud250PLSA, model =modelSA250PLsph)
plot(v.cloud250PLSA, model=fitSA250PL)
#AM
v.cloud250AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMSA, main = "Variogram - Log Surface Area: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA250AMsph = vgm(psill = .13, nugget=15 , range= 210, model= "Sph")
#now fit
fitSA250AM = fit.variogram(object= v.cloud250AMSA, model =modelSA250AMsph)
plot(v.cloud250AMSA, model=fitSA250AM)




####LOG ASPECT RATIO
#PL
v.cloud250PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR250PLsph = vgm(psill = .25, nugget=15 , range= 200, model= "Sph")
#now fit
fitAR250PL = fit.variogram(object= v.cloud250PLAR, model =modelAR250PLsph)
plot(v.cloud250PLAR, model=fitAR250PL)
#AM
v.cloud250AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR250AMsph = vgm(psill = .22, nugget=15 , range= 200, model= "Sph")
#now fit
fitAR250AM = fit.variogram(object= v.cloud250AMAR, model =modelAR250AMsph)
plot(v.cloud250AMAR, model=fitAR250AM)


####RANK STRETCH
#PL
v.cloud250PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLST, main = "Variogram - Rank Stretch: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST250PLsph = vgm(psill = 1200000, nugget=15 , range= 200, model= "Sph")
#now fit
fitST250PL = fit.variogram(object= v.cloud250PLST, model =modelST250PLsph)
plot(v.cloud250PLST, model=fitST250PL)
#AM
v.cloud250AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMST, main = "Variogram - Rank Stretch: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST250AMsph = vgm(psill = 2600000, nugget=15 , range= 220, model= "Sph")
#now fit
fitST250AM = fit.variogram(object= v.cloud250AMST, model =modelST250AMsph)
plot(v.cloud250AMST, model=fitST250AM)


####LOG WIDTH
#PL
v.cloud250PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLWI, main = "Variogram - Log Width: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI250PLsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitWI250PL = fit.variogram(object= v.cloud250PLWI, model =modelWI250PLsph)
plot(v.cloud250PLWI, model=fitWI250PL)
#AM
v.cloud250AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMWI, main = "Variogram - Log Width: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI250AMsph = vgm(psill = .05, nugget=10 , range= 220, model= "Sph")
#now fit
fitWI250AM = fit.variogram(object= v.cloud250AMWI, model =modelWI250AMsph)
plot(v.cloud250AMWI, model=fitWI250AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud250PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV250PLsph = vgm(psill = 12000000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV250PL = fit.variogram(object= v.cloud250PLFEV, model =modelFEV250PLsph)
plot(v.cloud250PLFEV, model=fitFEV250PL)
#AM
v.cloud250AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV250AMsph = vgm(psill = 2600000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV250AM = fit.variogram(object= v.cloud250AMFEV, model =modelFEV250AMsph)
plot(v.cloud250AMFEV, model=fitFEV250AM)


####LOG LENGTH
#PL
v.cloud250PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLLE, main = "Variogram - Log Length: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE250PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE250PL = fit.variogram(object= v.cloud250PLLE, model =modelLE250PLsph)
plot(v.cloud250PLLE, model=fitLE250PL)
#AM
v.cloud250AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMLE, main = "Variogram - Log Length: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE250AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE250AM = fit.variogram(object= v.cloud250AMLE, model =modelLE250AMsph)
plot(v.cloud250AMLE, model=fitLE250AM)


####LOG HEIGHT
#PL
v.cloud250PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLHE, main = "Variogram - Log Height: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE250PLsph = vgm(psill = .052, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE250PL = fit.variogram(object= v.cloud250PLHE, model =modelHE250PLsph)
plot(v.cloud250PLHE, model=fitHE250PL)
#AM
v.cloud250AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMHE, main = "Variogram - Log Height: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE250AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE250AM = fit.variogram(object= v.cloud250AMHE, model =modelHE250AMsph)
plot(v.cloud250AMHE, model=fitHE250AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud250PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA250PLsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA250PL = fit.variogram(object= v.cloud250PLFSA, model =modelFSA250PLsph)
plot(v.cloud250PLFSA, model=fitFSA250PL)
#AM
v.cloud250AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA250AMsph = vgm(psill = .11, nugget=10 , range= 210, model= "Sph")
#now fit
fitFSA250AM = fit.variogram(object= v.cloud250AMFSA, model =modelFSA250AMsph)
plot(v.cloud250AMFSA, model=fitFSA250AM)


####SPHERICITY
#PL
v.cloud250PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL250, cloud = F)
plot(v.cloud250PLFSP, main = "Variogram - Sphericity: PL Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH250PLsph = vgm(psill = .0015, nugget=10 , range= 200, model= "Sph")
#now fit
fitSPH250PL = fit.variogram(object= v.cloud250PLFSP, model =modelSPH250PLsph)
plot(v.cloud250PLFSP, model=fitSPH250PL)
#AM
v.cloud250AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM250, cloud = F)
plot(v.cloud250AMFSP, main = "Variogram - Sphericity: AM Mouse 250 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH250AMsph = vgm(psill = .0020, nugget=10 , range= 210, model= "Sph")
#now fit
fitSPH250AM = fit.variogram(object= v.cloud250AMFSP, model =modelSPH250AMsph)
plot(v.cloud250AMFSP, model=fitSPH250AM)












######################MOUSE 252
lacundfGPL252= lacundf%>%filter(id=="252"&
                                  location =="PL"&
                                  group==0)
lacundfGAM252= lacundf%>%filter(id=="252"&
                                  location =="AM"&
                                  group==0)
coordinates(lacundfGPL252)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM252)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol252PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL252, cloud = F)
plot(v.cloudvol252PL, main = "Variogram - Log Volume: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 20 first x axis point
#sill = .25 #first maximum 
#range = 205 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelfit252PLvol = vgm(psill = .20, nugget=25 , range= 200, model= "Sph")
#now fit
fit252PLvol = fit.variogram(object= v.cloudvol252PL, model =modelfit252PLvol )
plot(v.cloudvol252PL, model=fit252PLvol)
#AM
v.cloudvol252AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM252, cloud = F)
plot(v.cloudvol252AM, main = "Variogram - Log Volume: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol252AMsph = vgm(psill = .30, nugget=15 , range= 210, model= "Sph")
#now fit
fitvol252AM = fit.variogram(object= v.cloudvol252AM, model =modelvol252AMsph )
plot(v.cloudvol252AM, model=fitvol252AM)



####RANK OBLATENESS
#PL
v.cloud252PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLOB, main = "Variogram - Rank Oblateness: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB252PLsph = vgm(psill = 1500000, nugget=10 , range= 200, model= "Sph")
#now fit
fitOB252PL = fit.variogram(object= v.cloud252PLOB, model =modelOB252PLsph)
plot(v.cloud252PLOB, model=fitOB252PL)
#AM
v.cloud252AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMOB, main = "Variogram - Rank Oblateness: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB252AMsph = vgm(psill = 3000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB252AM = fit.variogram(object= v.cloud252AMOB, model =modelOB252AMsph)
plot(v.cloud252AMOB, model=fitOB252AM)


####LOG SURFACE AREA
#PL
v.cloud252PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLSA, main = "Variogram - Log Surface Area: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA252PLsph = vgm(psill = .11, nugget=15 , range= 200, model= "Sph")
#now fit
fitSA252PL = fit.variogram(object= v.cloud252PLSA, model =modelSA252PLsph)
plot(v.cloud252PLSA, model=fitSA252PL)
#AM
v.cloud252AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMSA, main = "Variogram - Log Surface Area: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA252AMsph = vgm(psill = .19, nugget=15 , range= 210, model= "Sph")
#now fit
fitSA252AM = fit.variogram(object= v.cloud252AMSA, model =modelSA252AMsph)
plot(v.cloud252AMSA, model=fitSA252AM)


####LOG ASPECT RATIO
#PL
v.cloud252PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR252PLsph = vgm(psill = .23, nugget=15 , range= 200, model= "Sph")
#now fit
fitAR252PL = fit.variogram(object= v.cloud252PLAR, model =modelAR252PLsph)
plot(v.cloud252PLAR, model=fitAR252PL)
#AM
v.cloud252AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR252AMsph = vgm(psill = .30, nugget=15 , range= 210, model= "Sph")
#now fit
fitAR252AM = fit.variogram(object= v.cloud252AMAR, model =modelAR252AMsph)
plot(v.cloud252AMAR, model=fitAR252AM)





####RANK STRETCH
#PL
v.cloud252PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLST, main = "Variogram - Rank Stretch: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST252PLsph = vgm(psill = 1500000, nugget=15 , range= 200, model= "Sph")
#now fit
fitST252PL = fit.variogram(object= v.cloud252PLST, model =modelST252PLsph)
plot(v.cloud252PLST, model=fitST252PL)
#AM
v.cloud252AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMST, main = "Variogram - Rank Stretch: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST252AMsph = vgm(psill = 3000000, nugget=15 , range= 220, model= "Sph")
#now fit
fitST252AM = fit.variogram(object= v.cloud252AMST, model =modelST252AMsph)
plot(v.cloud252AMST, model=fitST252AM)


####LOG WIDTH
#PL
v.cloud252PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLWI, main = "Variogram - Log Width: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI252PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitWI252PL = fit.variogram(object= v.cloud252PLWI, model =modelWI252PLsph)
plot(v.cloud252PLWI, model=fitWI252PL)
#AM
v.cloud252AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMWI, main = "Variogram - Log Width: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI252AMsph = vgm(psill = .07, nugget=10 , range= 210, model= "Sph")
#now fit
fitWI252AM = fit.variogram(object= v.cloud252AMWI, model =modelWI252AMsph)
plot(v.cloud252AMWI, model=fitWI252AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud252PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV252PLsph = vgm(psill = 1300000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV252PL = fit.variogram(object= v.cloud252PLFEV, model =modelFEV252PLsph)
plot(v.cloud252PLFEV, model=fitFEV252PL)
#AM
v.cloud252AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV252AMsph = vgm(psill = 3000000, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV252AM = fit.variogram(object= v.cloud252AMFEV, model =modelFEV252AMsph)
plot(v.cloud252AMFEV, model=fitFEV252AM)



####LOG LENGTH
#PL
v.cloud252PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLLE, main = "Variogram - Log Length: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE252PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE252PL = fit.variogram(object= v.cloud252PLLE, model =modelLE252PLsph)
plot(v.cloud252PLLE, model=fitLE252PL)
#AM
v.cloud252AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMLE, main = "Variogram - Log Length: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE252AMsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE252AM = fit.variogram(object= v.cloud252AMLE, model =modelLE252AMsph)
plot(v.cloud252AMLE, model=fitLE252AM)


####LOG HEIGHT
#PL
v.cloud252PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLHE, main = "Variogram - Log Height: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE252PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE252PL = fit.variogram(object= v.cloud252PLHE, model =modelHE252PLsph)
plot(v.cloud252PLHE, model=fitHE252PL)
#AM
v.cloud252AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMHE, main = "Variogram - Log Height: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE252AMsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE252AM = fit.variogram(object= v.cloud252AMHE, model =modelHE252AMsph)
plot(v.cloud252AMHE, model=fitHE252AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud252PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA252PLsph = vgm(psill = .10, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA252PL = fit.variogram(object= v.cloud252PLFSA, model =modelFSA252PLsph)
plot(v.cloud252PLFSA, model=fitFSA252PL)
#AM
v.cloud252AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA252AMsph = vgm(psill = .17, nugget=10 , range= 210, model= "Sph")
#now fit
fitFSA252AM = fit.variogram(object= v.cloud252AMFSA, model =modelFSA252AMsph)
plot(v.cloud252AMFSA, model=fitFSA252AM)


####SPHERICITY
#PL
v.cloud252PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL252, cloud = F)
plot(v.cloud252PLFSP, main = "Variogram - Sphericity: PL Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH252PLsph = vgm(psill = .00175, nugget=10 , range= 200, model= "Sph")
#now fit
fitSPH252PL = fit.variogram(object= v.cloud252PLFSP, model =modelSPH252PLsph)
plot(v.cloud252PLFSP, model=fitSPH252PL)
#AM
v.cloud252AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM252, cloud = F)
plot(v.cloud252AMFSP, main = "Variogram - Sphericity: AM Mouse 252 (GROUND)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSPH252AMsph = vgm(psill = .00230, nugget=10 , range= 210, model= "Sph")
#now fit
fitSPH252AM = fit.variogram(object= v.cloud252AMFSP, model =modelSPH252AMsph)
plot(v.cloud252AMFSP, model=fitSPH252AM)





#MAKING AVERAGE DF TO GET ONE SEMIVARIOGRAM BY STRAIN AND GROUP

###LOG VOLUME

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfvolGPL = rbind(fit252PLvol, fit238PLvol, fitmess, 
                    fit250PLvol)
avgdfvolGG<- avgdfvolGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                "ang1",
            "ang2",
            "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
mess1 = rbind(v.cloudvol252PL,v.cloudvol233PL, 
              v.cloudvol238PL, v.cloudvol250PL)
#create fixed variogram fit based on averaged parameters 
vgmm = vgm(psill = .09683775, nugget=.11224837 , range= 196.9409, model= "Sph")
#plot
plot(mess1, vgmm,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Volume with Spherical Fit for Ground Control Group in PL Region")



#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfvolGAM = rbind(fitvol252AM, fitvol250AM, fitvol238AM, 
                    fitvol233AM)
avgdfvolGGAM<- avgdfvolGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatvolGAM = rbind(v.cloudvol252AM,v.cloudvol233AM, 
              v.cloudvol238AM, v.cloudvol250AM)
#create fixed variogram fit based on averaged parameters 
vgmvolGAM = vgm(psill = .08719647, nugget=.14959328 , range= 173.476, model= "Sph")
#plot
plot(semidatvolGAM, vgmvolGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Volume with Spherical Fit for Ground Control Group in AM Region")





###RANK OBLATENESS

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfOBGPL = rbind(fitOB252PL, fitOB238PL, fitOB233PL, 
                   fitOB250PL)
avgdfOBGG<- avgdfOBGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatOBGPL = rbind(v.cloud252PLOB,v.cloud233PLOB, 
              v.cloud238PLOB, v.cloud250PLOB)
#create fixed variogram fit based on averaged parameters 
vgmOBGPL = vgm(psill = 134221.2, nugget=1291879.6 , range= 1108.06, model= "Sph")
#plot
plot(semidatOBGPL, vgmOBGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Oblateness with Spherical Fit for Ground Control Group in PL Region")



#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfOBGAM = rbind(fitOB252AM, fitOB238AM, fitOB233AM, 
                 fitOB250AM)
avgdfOBGGAM<- avgdfOBGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatOBGAM = rbind(v.cloud252AMOB,v.cloud233AMOB, 
                     v.cloud238AMOB, v.cloud250AMOB)
#create fixed variogram fit based on averaged parameters 
vgmOBGAM = vgm(psill = 232855.3, nugget=1982501.8 , range= 131.8117, model= "Sph")
#plot
plot(semidatOBGAM, vgmOBGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Oblateness with Spherical Fit for Ground Control Group in AM Region")






###LOG SURFACE AREA 

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfSAGPL = rbind(fitSA252PL, fitSA238PL, fitSA233PL, 
                   fitSA250PL)
avgdfSAGG<- avgdfSAGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSAGPL = rbind(v.cloud252PLSA,v.cloud233PLSA, 
                     v.cloud238PLSA, v.cloud250PLSA)
#create fixed variogram fit based on averaged parameters 
vgmSAGPL = vgm(psill = 0.08170282, nugget=.07481537 , range= 811.074, model= "Sph")
#plot
plot(semidatSAGPL, vgmSAGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Surface Area with Spherical Fit for Ground Group in PL Region")





#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfSAGAM = rbind(fitSA252AM, fitSA238AM, fitSA233AM, 
                   fitSA250AM)
avgdfSAGGAM<- avgdfSAGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSAGAM = rbind(v.cloud252AMSA,v.cloud233AMSA, 
                     v.cloud238AMSA, v.cloud250AMSA)
#create fixed variogram fit based on averaged parameters 
vgmSAGAM = vgm(psill = 0.04777692, nugget=0.09614153 , range= 179.9277, model= "Sph")
#plot
plot(semidatSAGAM, vgmSAGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Surface Area with Spherical Fit for Ground Group in AM Region")




###LOG ASPECT RATIO 

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfARGPL = rbind(fitAR252PL, fitAR238PL, fitAR233PL, 
                   fitAR250PL)
avgdfARGG<- avgdfARGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatARGPL = rbind(v.cloud252PLAR,v.cloud233PLAR, 
                     v.cloud238PLAR, v.cloud250PLAR)
#create fixed variogram fit based on averaged parameters 
vgmARGPL = vgm(psill = 0.02723733, nugget=.21822528 , range= 381.5971, model= "Sph")
#plot
plot(semidatARGPL, vgmARGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Aspect Ratio with Spherical Fit for Ground Control Group in PL Region")




#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfARGAM = rbind(fitAR252AM, fitAR238AM, fitAR233AM, 
                   fitAR250AM)
avgdfARGGAM= avgdfARGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatARGAM = rbind(v.cloud252AMAR,v.cloud233AMAR, 
                     v.cloud238AMAR, v.cloud250AMAR)
#create fixed variogram fit based on averaged parameters 
vgmARGAM = vgm(psill = 0.02256842, nugget=.22817826 , range= 189.7254, model= "Sph")
#plot
plot(semidatARGAM, vgmARGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Aspect Ratio with Spherical Fit for Ground Group in AM Region")







###RANK STRETCH 

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfSTGPL = rbind(fitST252PL, fitST238PL, fitST233PL, 
                   fitST250PL)
avgdfSTGG<- avgdfSTGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSTGPL = rbind(v.cloud252PLST,v.cloud233PLST, 
                     v.cloud238PLST, v.cloud250PLST)
#create fixed variogram fit based on averaged parameters 
vgmSTGPL = vgm(psill = 351230.3, nugget=950305.5 , range= 122.2938, model= "Sph")
#plot
plot(semidatSTGPL, vgmSTGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Stretch with Spherical Fit for Ground Control Group in PL Region")




#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfSTGAM = rbind(fitST252AM, fitST238AM, fitST233AM, 
                   fitST250AM)
avgdfSTGGAM<- avgdfSTGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSTGAM = rbind(v.cloud252AMST,v.cloud233AMST, 
                     v.cloud238AMST, v.cloud250AMST)
#create fixed variogram fit based on averaged parameters 
vgmSTGAM = vgm(psill = 332363.3, nugget=1872011.7 , range= 135.4579, model= "Sph")
#plot
plot(semidatSTGAM, vgmSTGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Stretch with Spherical Fit for Ground Control Group in AM Region")






###LOG WIDTH 

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfWIGPL = rbind(fitWI252PL, fitWI238PL, fitWI233PL, 
                   fitWI250PL)
avgdfWIGG<- avgdfWIGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatWIGPL = rbind(v.cloud252PLWI,v.cloud233PLWI, 
                     v.cloud238PLWI, v.cloud250PLWI)
#create fixed variogram fit based on averaged parameters 
vgmWIGPL = vgm(psill = 0.02738793, nugget=0.03335808 , range= 153.7728, model= "Sph")
#plot
plot(semidatWIGPL, vgmWIGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Width with Spherical Fit for Ground Control Group in PL Region")



#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfWIGAM = rbind(fitWI252AM, fitWI238AM, fitWI233AM, 
                   fitWI250AM)
avgdfWIGGAM<- avgdfWIGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatWIGAM = rbind(v.cloud252AMWI,v.cloud233AMWI, 
                     v.cloud238AMWI, v.cloud250AMWI)
#create fixed variogram fit based on averaged parameters 
vgmWIGAM = vgm(psill = 0.02125010, nugget=0.03745973 , range= 167.3083, model= "Sph")
#plot
plot(semidatWIGAM, vgmWIGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Width with Spherical Fit for Ground Control Group in AM Region")






###RANK FITTED ELLIPSOID VLUME  

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfFEVGPL = rbind(fitFEV252PL, fitFEV238PL, fitFEV233PL, 
                   fitFEV250PL)
avgdfFEVGG<- avgdfFEVGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFEVGPL = rbind(v.cloud252PLFEV,v.cloud233PLFEV, 
                     v.cloud238PLFEV, v.cloud250PLFEV)
#create fixed variogram fit based on averaged parameters 
vgmFEVGPL = vgm(psill = 658450.0, nugget=620764.2 , range= 146.7021, model= "Sph")
#plot
plot(semidatFEVGPL, vgmFEVGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Fitted Ellipsoid Volume with Spherical Fit for Ground Control Group in PL Region")



#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfFEVGAM = rbind(fitFEV252AM, fitFEV238AM, fitFEV233AM, 
                    fitFEV250AM)
avgdfFEVGGAM<- avgdfFEVGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFEVGAM = rbind(v.cloud252AMFEV,v.cloud233AMFEV, 
                      v.cloud238AMFEV, v.cloud250AMFEV)
#create fixed variogram fit based on averaged parameters 
vgmFEVGAM = vgm(psill = 984977.4, nugget=1186761.9 , range= 140.6721, model= "Sph")
#plot
plot(semidatFEVGAM, vgmFEVGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Fitted Ellipsoid Volume with Spherical Fit for Ground Control Group in AM Region")






###LOG LENGTH

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfLEGPL = rbind(fitLE252PL, fitLE238PL, fitLE233PL, 
                    fitLE250PL)
avgdfLEGG<- avgdfLEGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatLEGPL = rbind(v.cloud252PLLE,v.cloud233PLLE, 
                      v.cloud238PLLE, v.cloud250PLLE)
#create fixed variogram fit based on averaged parameters 
vgmLEGPL = vgm(psill = 0.01418780, nugget=0.05254127 , range= 917.2437, model= "Sph")
#plot
plot(semidatLEGPL, vgmLEGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Length with Spherical Fit for Ground Control Group in PL Region")



#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfLEGAM = rbind(fitLE252AM, fitLE238AM, fitLE233AM, 
                   fitLE250AM)
avgdfLEGGAM<- avgdfLEGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatLEGAM = rbind(v.cloud252AMLE,v.cloud233AMLE, 
                     v.cloud238AMLE, v.cloud250AMLE)
#create fixed variogram fit based on averaged parameters 
vgmLEGAM = vgm(psill = 0.01189084, nugget=0.05554982 , range= 170.7729, model= "Sph")
#plot
plot(semidatLEGAM, vgmLEGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Length with Spherical Fit for Ground Control Group in AM Region")





###LOG HEIGHT

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfHEGPL = rbind(fitHE252PL, fitHE238PL, fitHE233PL, 
                   fitHE250PL)
avgdfHEGG<- avgdfHEGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatHEGPL = rbind(v.cloud252PLHE,v.cloud233PLHE, 
                     v.cloud238PLHE, v.cloud250PLHE)
#create fixed variogram fit based on averaged parameters 
vgmHEGPL = vgm(psill = 0.01241881, nugget=0.04766225 , range=499.1472, model= "Sph")
#plot
plot(semidatHEGPL, vgmHEGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Height with Spherical Fit for Ground Control Group in PL Region")


#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfHEGAM = rbind(fitHE252AM, fitHE238AM, fitHE233AM, 
                   fitHE250AM)
avgdfHEGGAM<- avgdfHEGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatHEGAM = rbind(v.cloud252AMHE,v.cloud233AMHE, 
                     v.cloud238AMHE, v.cloud250AMHE)
#create fixed variogram fit based on averaged parameters 
vgmHEGAM = vgm(psill = 0.01184656, nugget=0.05735947 , range=198.6737, model= "Sph")
#plot
plot(semidatHEGAM, vgmHEGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Height with Spherical Fit for Ground Control Group in AM Region")






###LOG Fitted Ellipsoid Surface Area

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfFSAGPL = rbind(fitFSA252PL, fitFSA238PL, fitFSA233PL, 
                   fitFSA250PL)
avgdfFSAGG<- avgdfFSAGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFSAGPL = rbind(v.cloud252PLFSA,v.cloud233PLFSA, 
                     v.cloud238PLFSA, v.cloud250PLFSA)
#create fixed variogram fit based on averaged parameters 
vgmFSAGPL = vgm(psill = 0.04007449, nugget=0.08235792 , range=496.4118, model= "Sph")
#plot
plot(semidatFSAGPL, vgmFSAGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Fitted Ellipsoid Surface Area with Spherical Fit for Ground Control Group in PL Region")



#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfFSAGAM= rbind(fitFSA252AM, fitFSA238AM, fitFSA233AM, 
                    fitFSA250AM)
avgdfFSAGGAM<- avgdfFSAGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFSAGAM = rbind(v.cloud252AMFSA,v.cloud233AMFSA, 
                      v.cloud238AMFSA, v.cloud250AMFSA)
#create fixed variogram fit based on averaged parameters 
vgmFSAGAM = vgm(psill = 0.03842507, nugget=0.09793628 , range=174.6563, model= "Sph")
#plot
plot(semidatFSAGAM, vgmFSAGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Fitted Ellipsoid Surface Area with Spherical Fit for Ground Control Group in AM Region")






###SPHERICITY

#GROUND CONTROL, PL REGION
#get avgeraged parameters
avgdfSPGPL = rbind(fitSPH252PL, fitSPH238PL, 
                    fitSPH250PL)
avgdfSPGG<- avgdfSPGPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSPGPL = rbind(v.cloud252PLFSP, v.cloud233PLFSP,
                     v.cloud250PLFSP, v.cloud238PLFSP)
#create fixed variogram fit based on averaged parameters 
vgmSPGPL = vgm(psill = 0.0002676828, nugget=0.0015660378 , range=193.8342, model= "Sph")
#plot
plot(semidatSPGPL, vgmSPGPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Sphericity with Spherical Fit for Ground Control Group in PL Region")


#GROUND CONTROL, AM REGION
#get avgeraged parameters
avgdfSPGAM = rbind(fitSPH252AM, fitSPH238AM, fitSPH233AM,
                   fitSPH250AM)
avgdfSPGGAM<- avgdfSPGAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSPGAM = rbind(v.cloud252AMFSP, v.cloud233AMFSP,
                     v.cloud250AMFSP, v.cloud238AMFSP)
#create fixed variogram fit based on averaged parameters 
vgmSPGAM = vgm(psill = 0.01191093, nugget=0.00191471 , range=5677.471, model= "Sph")
#plot
plot(semidatSPGAM, vgmSPGAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Sphericity with Spherical Fit for Ground Control Group in AM Region")





















######################MOUSE 202
lacundfGPL202= lacundf%>%filter(id=="202"&
                                  location =="PL"&
                                  group==1)
lacundfGAM202= lacundf%>%filter(id=="202"&
                                  location =="AM"&
                                  group==1)
coordinates(lacundfGPL202)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM202)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")


####LOG VOLUME
#PL
v.cloudvol202PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL202, cloud = F)
plot(v.cloudvol202PL, main = "Variogram - Log Volume: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol202PLsph = vgm(psill = .20, nugget=15 , range= 200, model= "Sph")
#now fit
fitvol202PL = fit.variogram(object= v.cloudvol202PL, model =modelvol202PLsph )
plot(v.cloudvol202PL, model=fitvol202PL)
#AM
v.cloudvol202AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM202, cloud = F)
plot(v.cloudvol202AM, main = "Variogram - Log Volume: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol202AMsph = vgm(psill = .20, nugget=25 , range= 210, model= "Sph")
#now fit
fitvol202AM = fit.variogram(object= v.cloudvol202AM, model =modelvol202AMsph )
plot(v.cloudvol202AM, model=fitvol202AM)



####RANK OBLATENESS
#PL
v.cloud202PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLOB, main = "Variogram - Rank Oblateness: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB202PLsph = vgm(psill = 900000, nugget=10 , range= 200, model= "Sph")
#now fit
fitOB202PL = fit.variogram(object= v.cloud202PLOB, model =modelOB202PLsph)
plot(v.cloud202PLOB, model=fitOB202PL)
#AM
v.cloud202AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMOB, main = "Variogram - Rank Oblateness: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB202AMsph = vgm(psill = 1700000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB202AM = fit.variogram(object= v.cloud202AMOB, model =modelOB202AMsph)
plot(v.cloud202AMOB, model=fitOB202AM)



####LOG SURFACE AREA
#PL
v.cloud202PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLSA, main = "Variogram - Log Surface Area: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA202PLsph = vgm(psill = .10, nugget=10 , range= 200, model= "Sph")
#now fit
fitSA202PL = fit.variogram(object= v.cloud202PLSA, model =modelSA202PLsph)
plot(v.cloud202PLSA, model=fitSA202PL)
#AM
v.cloud202AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMSA, main = "Variogram - Log Surface Area: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA202AMsph = vgm(psill = .12, nugget=10 , range= 210, model= "Sph")
#now fit
fitSA202AM = fit.variogram(object= v.cloud202AMSA, model =modelSA202AMsph)
plot(v.cloud202AMSA, model=fitSA202AM)




####LOG ASPECT RATIO
#PL
v.cloud202PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR202PLsph = vgm(psill = .24, nugget=5 , range= 200, model= "Sph")
#now fit
fitAR202PL = fit.variogram(object= v.cloud202PLAR, model =modelAR202PLsph)
plot(v.cloud202PLAR, model=fitAR202PL)
#AM
v.cloud202AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR202AMsph = vgm(psill = .25, nugget=5 , range= 210, model= "Sph")
#now fit
fitAR202AM = fit.variogram(object= v.cloud202AMAR, model =modelAR202AMsph)
plot(v.cloud202AMAR, model=fitAR202AM)



####RANK STRETCH
#PL
v.cloud202PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLST, main = "Variogram - Rank Stretch: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST202PLsph = vgm(psill = 9e5, nugget=10 , range= 200, model= "Sph")
#now fit
fitST202PL = fit.variogram(object= v.cloud202PLST, model =modelST202PLsph)
plot(v.cloud202PLST, model=fitST202PL)
#AM
v.cloud202AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMST, main = "Variogram - Rank Stretch: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST202AMsph = vgm(psill = 1600000, nugget=10 , range= 220, model= "Sph")
#now fit
fitST202AM = fit.variogram(object= v.cloud202AMST, model =modelST202AMsph)
plot(v.cloud202AMST, model=fitST202AM)




####LOG WIDTH
#PL
v.cloud202PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLWI, main = "Variogram - Log Width: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI202PLsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitWI202PL = fit.variogram(object= v.cloud202PLWI, model =modelWI202PLsph)
plot(v.cloud202PLWI, model=fitWI202PL)
#AM
v.cloud202AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMWI, main = "Variogram - Log Width: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI202AMsph = vgm(psill = .07, nugget=10 , range= 220, model= "Sph")
#now fit
fitWI202AM = fit.variogram(object= v.cloud202AMWI, model =modelWI202AMsph)
plot(v.cloud202AMWI, model=fitWI202AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud202PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV202PLsph = vgm(psill = 9e5, nugget=10 , range= 200, model= "Sph")
#now fit
fitFEV202PL = fit.variogram(object= v.cloud202PLFEV, model =modelFEV202PLsph)
plot(v.cloud202PLFEV, model=fitFEV202PL)
#AM
v.cloud202AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV202AMsph = vgm(psill = 15000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV202AM = fit.variogram(object= v.cloud202AMFEV, model =modelFEV202AMsph)
plot(v.cloud202AMFEV, model=fitFEV202AM)


####LOG LENGTH
#PL
v.cloud202PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLLE, main = "Variogram - Log Length: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE202PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE202PL = fit.variogram(object= v.cloud202PLLE, model =modelLE202PLsph)
plot(v.cloud202PLLE, model=fitLE202PL)
#AM
v.cloud202AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMLE, main = "Variogram - Log Length: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE202AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE202AM = fit.variogram(object= v.cloud202AMLE, model =modelLE202AMsph)
plot(v.cloud202AMLE, model=fitLE202AM)


####LOG HEIGHT
#PL
v.cloud202PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLHE, main = "Variogram - Log Height: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE202PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE202PL = fit.variogram(object= v.cloud202PLHE, model =modelHE202PLsph)
plot(v.cloud202PLHE, model=fitHE202PL)
#AM
v.cloud202AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMHE, main = "Variogram - Log Height: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE202AMsph = vgm(psill = .06, nugget=10 , range= 210, model= "Sph")
#now fit
fitHE202AM = fit.variogram(object= v.cloud202AMHE, model =modelHE202AMsph)
plot(v.cloud202AMHE, model=fitHE202AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud202PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA202PLsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA202PL = fit.variogram(object= v.cloud202PLFSA, model =modelFSA202PLsph)
plot(v.cloud202PLFSA, model=fitFSA202PL)
#AM
v.cloud202AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA202AMsph = vgm(psill = .11, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA202AM = fit.variogram(object= v.cloud202AMFSA, model =modelFSA202AMsph)
plot(v.cloud202AMFSA, model=fitFSA202AM)


####SPHERICITY
#PL
v.cloud202PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL202, cloud = F)
plot(v.cloud202PLFSP, main = "Variogram - Sphericity: PL Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP202PLsph = vgm(psill = .0015, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP202PL = fit.variogram(object= v.cloud202PLFSP, model =modelSP202PLsph)
plot(v.cloud202PLFSP, model=fitSP202PL)
#AM
v.cloud202AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM202, cloud = F)
plot(v.cloud202AMFSP, main = "Variogram - Sphericity: AM Mouse 202 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP202AMsph = vgm(psill = .0022, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP202AM = fit.variogram(object= v.cloud202AMFSP, model =modelSP202AMsph)
plot(v.cloud202AMFSP, model=fitSP202AM)












######################MOUSE 208
lacundfGPL208= lacundf%>%filter(id=="208"&
                                  location =="PL"&
                                  group==1)
lacundfGAM208= lacundf%>%filter(id=="208"&
                                  location =="AM"&
                                  group==1)
coordinates(lacundfGPL208)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM208)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol208PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL208, cloud = F)
plot(v.cloudvol208PL, main = "Variogram - Log Volume: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol208PLsph = vgm(psill = .20, nugget=15 , range= 200, model= "Sph")
#now fit
fitvol208PL = fit.variogram(object= v.cloudvol208PL, model =modelvol208PLsph )
plot(v.cloudvol208PL, model=fitvol208PL)
#AM
v.cloudvol208AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM208, cloud = F)
plot(v.cloudvol208AM, main = "Variogram - Log Volume: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol208AMsph = vgm(psill = .20, nugget=15 , range= 205, model= "Sph")
#now fit
fitvol208AM = fit.variogram(object= v.cloudvol208AM, model =modelvol208AMsph )
plot(v.cloudvol208AM, model=fitvol208AM)



####RANK OBLATENESS
#PL
v.cloud208PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLOB, main = "Variogram - Rank Oblateness: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB208PLsph = vgm(psill = 2000000, nugget=10 , range= 200, model= "Sph")
#now fit
fitOB208PL = fit.variogram(object= v.cloud208PLOB, model =modelOB208PLsph)
plot(v.cloud208PLOB, model=fitOB208PL)
#AM
v.cloud208AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMOB, main = "Variogram - Rank Oblateness: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB208AMsph = vgm(psill = 1800000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB208AM = fit.variogram(object= v.cloud208AMOB, model =modelOB208AMsph)
plot(v.cloud208AMOB, model=fitOB208AM)


####LOG SURFACE AREA
#PL
v.cloud208PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLSA, main = "Variogram - Log Surface Area: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA208PLsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitSA208PL = fit.variogram(object= v.cloud208PLSA, model =modelSA208PLsph)
plot(v.cloud208PLSA, model=fitSA208PL)
#AM
v.cloud208AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMSA, main = "Variogram - Log Surface Area: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA208AMsph = vgm(psill = .12, nugget=10 , range= 205, model= "Sph")
#now fit
fitSA208AM = fit.variogram(object= v.cloud208AMSA, model =modelSA208AMsph)
plot(v.cloud208AMSA, model=fitSA208AM)



####LOG ASPECT RATIO
#PL
v.cloud208PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR208PLsph = vgm(psill = .24, nugget=5 , range= 200, model= "Sph")
#now fit
fitAR208PL = fit.variogram(object= v.cloud208PLAR, model =modelAR208PLsph)
plot(v.cloud208PLAR, model=fitAR208PL)
#AM
v.cloud208AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR208AMsph = vgm(psill = .24, nugget=5 , range= 210, model= "Sph")
#now fit
fitAR208AM = fit.variogram(object= v.cloud208AMAR, model =modelAR208AMsph)
plot(v.cloud208AMAR, model=fitAR208AM)



####RANK STRETCH
#PL
v.cloud208PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLST, main = "Variogram - Rank Stretch: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST208PLsph = vgm(psill = 1700000, nugget=10 , range= 205, model= "Sph")
#now fit
fitST208PL = fit.variogram(object= v.cloud208PLST, model =modelST208PLsph)
plot(v.cloud208PLST, model=fitST208PL)
#AM
v.cloud208AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMST, main = "Variogram - Rank Stretch: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST208AMsph = vgm(psill = 13000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitST208AM = fit.variogram(object= v.cloud208AMST, model =modelST208AMsph)
plot(v.cloud208AMST, model=fitST208AM)


####LOG WIDTH
#PL
v.cloud208PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLWI, main = "Variogram - Log Width: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI208PLsph = vgm(psill = .07, nugget=10 , range= 205, model= "Sph")
#now fit
fitWI208PL = fit.variogram(object= v.cloud208PLWI, model =modelWI208PLsph)
plot(v.cloud208PLWI, model=fitWI208PL)
#AM
v.cloud208AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMWI, main = "Variogram - Log Width: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI208AMsph = vgm(psill = .07, nugget=10 , range= 205, model= "Sph")
#now fit
fitWI208AM = fit.variogram(object= v.cloud208AMWI, model =modelWI208AMsph)
plot(v.cloud208AMWI, model=fitWI208AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud208PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV208PLsph = vgm(psill = 17000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV208PL = fit.variogram(object= v.cloud208PLFEV, model =modelFEV208PLsph)
plot(v.cloud208PLFEV, model=fitFEV208PL)
#AM
v.cloud208AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV208AMsph = vgm(psill = 13000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV208AM = fit.variogram(object= v.cloud208AMFEV, model =modelFEV208AMsph)
plot(v.cloud208AMFEV, model=fitFEV208AM)


####LOG LENGTH
#PL
v.cloud208PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLLE, main = "Variogram - Log Length: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE208PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE208PL = fit.variogram(object= v.cloud208PLLE, model =modelLE208PLsph)
plot(v.cloud208PLLE, model=fitLE208PL)
#AM
v.cloud208AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMLE, main = "Variogram - Log Length: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE208AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE208AM = fit.variogram(object= v.cloud208AMLE, model =modelLE208AMsph)
plot(v.cloud208AMLE, model=fitLE208AM)


####LOG HEIGHT
#PL
v.cloud208PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLHE, main = "Variogram - Log Height: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE208PLsph = vgm(psill = .045, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE208PL = fit.variogram(object= v.cloud208PLHE, model =modelHE208PLsph)
plot(v.cloud208PLHE, model=fitHE208PL)
#AM
v.cloud208AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMHE, main = "Variogram - Log Height: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE208AMsph = vgm(psill = .06, nugget=10 , range= 210, model= "Sph")
#now fit
fitHE208AM = fit.variogram(object= v.cloud208AMHE, model =modelHE208AMsph)
plot(v.cloud208AMHE, model=fitHE208AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud208PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA208PLsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA208PL = fit.variogram(object= v.cloud208PLFSA, model =modelFSA208PLsph)
plot(v.cloud208PLFSA, model=fitFSA208PL)
#AM
v.cloud208AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA208AMsph = vgm(psill = .11, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA208AM = fit.variogram(object= v.cloud208AMFSA, model =modelFSA208AMsph)
plot(v.cloud208AMFSA, model=fitFSA208AM)


####SPHERICITY
#PL
v.cloud208PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL208, cloud = F)
plot(v.cloud208PLFSP, main = "Variogram - Sphericity: PL Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP208PLsph = vgm(psill = .0020, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP208PL = fit.variogram(object= v.cloud208PLFSP, model =modelSP208PLsph)
plot(v.cloud208PLFSP, model=fitSP208PL)
#AM
v.cloud208AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM208, cloud = F)
plot(v.cloud208AMFSP, main = "Variogram - Sphericity: AM Mouse 208 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP208AMsph = vgm(psill = .0022, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP208AM = fit.variogram(object= v.cloud208AMFSP, model =modelSP208AMsph)
plot(v.cloud208AMFSP, model=fitSP208AM)











######################MOUSE 215
lacundfGPL215= lacundf%>%filter(id=="215"&
                                  location =="PL"&
                                  group==1)
lacundfGAM215= lacundf%>%filter(id=="215"&
                                  location =="AM"&
                                  group==1)
coordinates(lacundfGPL215)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM215)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol215PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL215, cloud = F)
plot(v.cloudvol215PL, main = "Variogram - Log Volume: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol215PLsph = vgm(psill = .18, nugget=15 , range= 205, model= "Sph")
#now fit
fitvol215PL = fit.variogram(object= v.cloudvol215PL, model =modelvol215PLsph )
plot(v.cloudvol215PL, model=fitvol215PL)
#AM
v.cloudvol215AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM215, cloud = F)
plot(v.cloudvol215AM, main = "Variogram - Log Volume: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol215AMsph = vgm(psill = .25, nugget=15 , range= 205, model= "Sph")
#now fit
fitvol215AM = fit.variogram(object= v.cloudvol215AM, model =modelvol215AMsph )
plot(v.cloudvol215AM, model=fitvol215AM)



####RANK OBLATENESS
#PL
v.cloud215PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLOB, main = "Variogram - Rank Oblateness: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB215PLsph = vgm(psill = 900000, nugget=10 , range= 200, model= "Sph")
#now fit
fitOB215PL = fit.variogram(object= v.cloud215PLOB, model =modelOB215PLsph)
plot(v.cloud215PLOB, model=fitOB215PL)
#AM
v.cloud215AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMOB, main = "Variogram - Rank Oblateness: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB215AMsph = vgm(psill = 1900000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB215AM = fit.variogram(object= v.cloud215AMOB, model =modelOB215AMsph)
plot(v.cloud215AMOB, model=fitOB215AM)


####LOG SURFACE AREA
#PL
v.cloud215PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLSA, main = "Variogram - Log Surface Area: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA215PLsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitSA215PL = fit.variogram(object= v.cloud215PLSA, model =modelSA215PLsph)
plot(v.cloud215PLSA, model=fitSA215PL)
#AM
v.cloud215AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMSA, main = "Variogram - Log Surface Area: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA215AMsph = vgm(psill = .14, nugget=10 , range= 210, model= "Sph")
#now fit
fitSA215AM = fit.variogram(object= v.cloud215AMSA, model =modelSA215AMsph)
plot(v.cloud215AMSA, model=fitSA215AM)


####LOG ASPECT RATIO
#PL
v.cloud215PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR215PLsph = vgm(psill = .25, nugget=5 , range= 200, model= "Sph")
#now fit
fitAR215PL = fit.variogram(object= v.cloud215PLAR, model =modelAR215PLsph)
plot(v.cloud215PLAR, model=fitAR215PL)
#AM
v.cloud215AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR215AMsph = vgm(psill = .27, nugget=5 , range= 210, model= "Sph")
#now fit
fitAR215AM = fit.variogram(object= v.cloud215AMAR, model =modelAR215AMsph)
plot(v.cloud215AMAR, model=fitAR215AM)



####RANK STRETCH
#PL
v.cloud215PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLST, main = "Variogram - Rank Stretch: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST215PLsph = vgm(psill = 1000000, nugget=10 , range= 205, model= "Sph")
#now fit
fitST215PL = fit.variogram(object= v.cloud215PLST, model =modelST215PLsph)
plot(v.cloud215PLST, model=fitST215PL)
#AM
v.cloud215AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMST, main = "Variogram - Rank Stretch: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST215AMsph = vgm(psill = 16000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitST215AM = fit.variogram(object= v.cloud215AMST, model =modelST215AMsph)
plot(v.cloud215AMST, model=fitST215AM)



####LOG WIDTH
#PL
v.cloud215PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLWI, main = "Variogram - Log Width: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI215PLsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitWI215PL = fit.variogram(object= v.cloud215PLWI, model =modelWI215PLsph)
plot(v.cloud215PLWI, model=fitWI215PL)
#AM
v.cloud215AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMWI, main = "Variogram - Log Width: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI215AMsph = vgm(psill = .07, nugget=10 , range= 220, model= "Sph")
#now fit
fitWI215AM = fit.variogram(object= v.cloud215AMWI, model =modelWI215AMsph)
plot(v.cloud215AMWI, model=fitWI215AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud215PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV215PLsph = vgm(psill = 10000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV215PL = fit.variogram(object= v.cloud215PLFEV, model =modelFEV215PLsph)
plot(v.cloud215PLFEV, model=fitFEV215PL)
#AM
v.cloud215AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV215AMsph = vgm(psill = 16000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV215AM = fit.variogram(object= v.cloud215AMFEV, model =modelFEV215AMsph)
plot(v.cloud215AMFEV, model=fitFEV215AM)


####LOG LENGTH
#PL
v.cloud215PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLLE, main = "Variogram - Log Length: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE215PLsph = vgm(psill = .045, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE215PL = fit.variogram(object= v.cloud215PLLE, model =modelLE215PLsph)
plot(v.cloud215PLLE, model=fitLE215PL)
#AM
v.cloud215AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMLE, main = "Variogram - Log Length: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE215AMsph = vgm(psill = .06, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE215AM = fit.variogram(object= v.cloud215AMLE, model =modelLE215AMsph)
plot(v.cloud215AMLE, model=fitLE215AM)


####LOG HEIGHT
#PL
v.cloud215PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLHE, main = "Variogram - Log Height: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE215PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE215PL = fit.variogram(object= v.cloud215PLHE, model =modelHE215PLsph)
plot(v.cloud215PLHE, model=fitHE215PL)
#AM
v.cloud215AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMHE, main = "Variogram - Log Height: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE215AMsph = vgm(psill = .06, nugget=10 , range= 210, model= "Sph")
#now fit
fitHE215AM = fit.variogram(object= v.cloud215AMHE, model =modelHE215AMsph)
plot(v.cloud215AMHE, model=fitHE215AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud215PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA215PLsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA215PL = fit.variogram(object= v.cloud215PLFSA, model =modelFSA215PLsph)
plot(v.cloud215PLFSA, model=fitFSA215PL)
#AM
v.cloud215AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA215AMsph = vgm(psill = .12, nugget=10 , range= 210, model= "Sph")
#now fit
fitFSA215AM = fit.variogram(object= v.cloud215AMFSA, model =modelFSA215AMsph)
plot(v.cloud215AMFSA, model=fitFSA215AM)


####SPHERICITY
#PL
v.cloud215PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL215, cloud = F)
plot(v.cloud215PLFSP, main = "Variogram - Sphericity: PL Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP215PLsph = vgm(psill = .0015, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP215PL = fit.variogram(object= v.cloud215PLFSP, model =modelSP215PLsph)
plot(v.cloud215PLFSP, model=fitSP215PL)
#AM
v.cloud215AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM215, cloud = F)
plot(v.cloud215AMFSP, main = "Variogram - Sphericity: AM Mouse 215 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP215AMsph = vgm(psill = .0017, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP215AM = fit.variogram(object= v.cloud215AMFSP, model =modelSP215AMsph)
plot(v.cloud215AMFSP, model=fitSP215AM)









######################MOUSE 222
lacundfGPL222= lacundf%>%filter(id=="222"&
                                  location =="PL"&
                                  group==1)
lacundfGAM222= lacundf%>%filter(id=="222"&
                                  location =="AM"&
                                  group==1)
coordinates(lacundfGPL222)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfGAM222)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol222PL<-variogram(log(Volume..µm..)~1, data = lacundfGPL222, cloud = F)
plot(v.cloudvol222PL, main = "Variogram - Log Volume: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol222PLsph = vgm(psill = .14, nugget=15 , range= 200, model= "Sph")
#now fit
fitvol222PL = fit.variogram(object= v.cloudvol222PL, model =modelvol222PLsph )
plot(v.cloudvol222PL, model=fitvol222PL)
#AM
v.cloudvol222AM<-variogram(log(Volume..µm..)~1, data = lacundfGAM222, cloud = F)
plot(v.cloudvol222AM, main = "Variogram - Log Volume: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol222AMsph = vgm(psill = .15, nugget=15 , range= 210, model= "Sph")
#now fit
fitvol222AM = fit.variogram(object= v.cloudvol222AM, model =modelvol222AMsph )
plot(v.cloudvol222AM, model=fitvol222AM)



####RANK OBLATENESS
#PL
v.cloud222PLOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLOB, main = "Variogram - Rank Oblateness: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB222PLsph = vgm(psill = 1500000, nugget=10 , range= 200, model= "Sph")
#now fit
fitOB222PL = fit.variogram(object= v.cloud222PLOB, model =modelOB222PLsph)
plot(v.cloud222PLOB, model=fitOB222PL)
#AM
v.cloud222AMOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMOB, main = "Variogram - Rank Oblateness: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelOB222AMsph = vgm(psill = 1000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitOB222AM = fit.variogram(object= v.cloud222AMOB, model =modelOB222AMsph)
plot(v.cloud222AMOB, model=fitOB222AM)


####LOG SURFACE AREA
#PL
v.cloud222PLSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLSA, main = "Variogram - Log Surface Area: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA222PLsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitSA222PL = fit.variogram(object= v.cloud222PLSA, model =modelSA222PLsph)
plot(v.cloud222PLSA, model=fitSA222PL)
#AM
v.cloud222AMSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMSA, main = "Variogram - Log Surface Area: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSA222AMsph = vgm(psill = .09, nugget=10 , range= 210, model= "Sph")
#now fit
fitSA222AM = fit.variogram(object= v.cloud222AMSA, model =modelSA222AMsph)
plot(v.cloud222AMSA, model=fitSA222AM)



####LOG ASPECT RATIO
#PL
v.cloud222PLAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLAR, main = "Variogram - Log Aspect Ratio: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR222PLsph = vgm(psill = .26, nugget=5 , range= 200, model= "Sph")
#now fit
fitAR222PL = fit.variogram(object= v.cloud222PLAR, model =modelAR222PLsph)
plot(v.cloud222PLAR, model=fitAR222PL)
#AM
v.cloud222AMAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMAR, main = "Variogram - Log Aspect Ratio: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelAR222AMsph = vgm(psill = .25, nugget=5 , range= 210, model= "Sph")
#now fit
fitAR222AM = fit.variogram(object= v.cloud222AMAR, model =modelAR222AMsph)
plot(v.cloud222AMAR, model=fitAR222AM)


####RANK STRETCH
#PL
v.cloud222PLST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLST, main = "Variogram - Rank Stretch: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST222PLsph = vgm(psill = 1505000, nugget=10 , range= 200, model= "Sph")
#now fit
fitST222PL = fit.variogram(object= v.cloud222PLST, model =modelST222PLsph)
plot(v.cloud222PLST, model=fitST222PL)
#AM
v.cloud222AMST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMST, main = "Variogram - Rank Stretch: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelST222AMsph = vgm(psill = 12000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitST222AM = fit.variogram(object= v.cloud222AMST, model =modelST222AMsph)
plot(v.cloud222AMST, model=fitST222AM)



####LOG WIDTH
#PL
v.cloud222PLWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLWI, main = "Variogram - Log Width: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI222PLsph = vgm(psill = .065, nugget=10 , range= 200, model= "Sph")
#now fit
fitWI222PL = fit.variogram(object= v.cloud222PLWI, model =modelWI222PLsph)
plot(v.cloud222PLWI, model=fitWI222PL)
#AM
v.cloud222AMWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMWI, main = "Variogram - Log Width: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelWI222AMsph = vgm(psill = .055, nugget=10 , range= 220, model= "Sph")
#now fit
fitWI222AM = fit.variogram(object= v.cloud222AMWI, model =modelWI222AMsph)
plot(v.cloud222AMWI, model=fitWI222AM)


####Rank Fitted Ellipsoid Volume
#PL
v.cloud222PLFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV222PLsph = vgm(psill = 15000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV222PL = fit.variogram(object= v.cloud222PLFEV, model =modelFEV222PLsph)
plot(v.cloud222PLFEV, model=fitFEV222PL)
#AM
v.cloud222AMFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFEV222AMsph = vgm(psill = 10000000, nugget=10 , range= 210, model= "Sph")
#now fit
fitFEV222AM = fit.variogram(object= v.cloud222AMFEV, model =modelFEV222AMsph)
plot(v.cloud222AMFEV, model=fitFEV222AM)



####LOG LENGTH
#PL
v.cloud222PLLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLLE, main = "Variogram - Log Length: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE222PLsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE222PL = fit.variogram(object= v.cloud222PLLE, model =modelLE222PLsph)
plot(v.cloud222PLLE, model=fitLE222PL)
#AM
v.cloud222AMLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMLE, main = "Variogram - Log Length: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelLE222AMsph = vgm(psill = .05, nugget=10 , range= 200, model= "Sph")
#now fit
fitLE222AM = fit.variogram(object= v.cloud222AMLE, model =modelLE222AMsph)
plot(v.cloud222AMLE, model=fitLE222AM)


####LOG HEIGHT
#PL
v.cloud222PLHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLHE, main = "Variogram - Log Height: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE222PLsph = vgm(psill = .045, nugget=10 , range= 200, model= "Sph")
#now fit
fitHE222PL = fit.variogram(object= v.cloud222PLHE, model =modelHE222PLsph)
plot(v.cloud222PLHE, model=fitHE222PL)
#AM
v.cloud222AMHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMHE, main = "Variogram - Log Height: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelHE222AMsph = vgm(psill = .05, nugget=10 , range= 210, model= "Sph")
#now fit
fitHE222AM = fit.variogram(object= v.cloud222AMHE, model =modelHE222AMsph)
plot(v.cloud222AMHE, model=fitHE222AM)


####LOG FITTED ELLIPSOID SURFACE AREA
#PL
v.cloud222PLFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA222PLsph = vgm(psill = .09, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA222PL = fit.variogram(object= v.cloud222PLFSA, model =modelFSA222PLsph)
plot(v.cloud222PLFSA, model=fitFSA222PL)
#AM
v.cloud222AMFSA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMFSA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelFSA222AMsph = vgm(psill = .08, nugget=10 , range= 200, model= "Sph")
#now fit
fitFSA222AM = fit.variogram(object= v.cloud222AMFSA, model =modelFSA222AMsph)
plot(v.cloud222AMFSA, model=fitFSA222AM)


####SPHERICITY
#PL
v.cloud222PLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGPL222, cloud = F)
plot(v.cloud222PLFSP, main = "Variogram - Sphericity: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP222PLsph = vgm(psill = .0015, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP222PL = fit.variogram(object= v.cloud222PLFSP, model =modelSP222PLsph)
plot(v.cloud222PLFSP, model=fitSP222PL)
#AM
v.cloud222AMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfGAM222, cloud = F)
plot(v.cloud222AMFSP, main = "Variogram - Sphericity: AM Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= 25 first x axis point
#sill = .20 #first maximum 
#range = 210 #distance at which sill occurs (x value)
##get those values from empirical semi ^^
modelSP222AMsph = vgm(psill = .0020, nugget=10 , range= 200, model= "Sph")
#now fit
fitSP222AM = fit.variogram(object= v.cloud222AMFSP, model =modelSP222AMsph)
plot(v.cloud222AMFSP, model=fitSP222AM)





#MAKING AVERAGE DF TO GET ONE SEMIVARIOGRAM BY STRAIN AND GROUP

###LOG VOLUME

#FLIGHT, AM REGION
#get avgeraged parameters
avgdfvolFAM = rbind(fitvol202AM, fitvol208AM, fitvol215AM, 
                    fitvol222AM)
avgdfvolFGAM<- avgdfvolFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatvolFAM = rbind(v.cloudvol202AM,v.cloudvol208AM, 
                      v.cloudvol215AM, v.cloudvol222AM)
#create fixed variogram fit based on averaged parameters 
vgmvolFAM = vgm(psill = .07696824, nugget=.12431317 , range= 187.8835, model= "Sph")
#plot
plot(semidatvolFAM, vgmvolFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Volume with Spherical Fit for Flight Group in AM Region")



#FLIGHT, PL REGION
#get avgeraged parameters
avgdfvolFPL = rbind(fitvol202PL, fitvol208PL, fitvol215PL, 
                    fitvol222PL)
avgdfvolFGPL<- avgdfvolFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatvolFPL = rbind(v.cloudvol202PL,v.cloudvol208PL, 
                      v.cloudvol215PL, v.cloudvol222PL)
#create fixed variogram fit based on averaged parameters 
vgmvolFPL = vgm(psill = .07652503, nugget=.08506363, range= 183.1253, model= "Sph")
#plot
plot(semidatvolFPL, vgmvolFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Volume with Spherical Fit for Flight Group in PL Region")


#FLIGHT, AM REGION
#get avgeraged parameters
avgdfvolFAM = rbind(fitvol202AM, fitvol208AM, fitvol215AM, 
                    fitvol222AM)
avgdfvolFGAM<- avgdfvolFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatvolFAM = rbind(v.cloudvol202AM,v.cloudvol208AM, 
                      v.cloudvol215AM, v.cloudvol222AM)
#create fixed variogram fit based on averaged parameters 
vgmvolFAM = vgm(psill = .07696824, nugget=.12431317, range= 187.8835, model= "Sph")
#plot
plot(semidatvolFAM, vgmvolFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram with Spherical Fit for Flight Group in AM Region")





###RANK OBLATENESS

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfOBFPL = rbind(fitOB202PL, fitOB208PL, fitOB215PL, 
                    fitOB222PL)

avgdfOBFGPL<- avgdfOBFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatOBFPL = rbind(v.cloud202PLOB,v.cloud208PLOB, 
                      v.cloud215PLOB, v.cloud222PLOB)
#create fixed variogram fit based on averaged parameters 
vgmOBFPL = vgm(psill = 34414.58, nugget=1317450.39 , range= 178.1449, model= "Sph")
#plot
plot(semidatOBFPL, vgmOBFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Oblateness with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfOBFAM = rbind(fitOB202AM, fitOB208AM, fitOB215AM, 
                   fitOB222AM)

avgdfOBFGAM<- avgdfOBFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatOBFAM = rbind(v.cloud202AMOB,v.cloud208AMOB, 
                     v.cloud215AMOB, v.cloud222AMOB)
#create fixed variogram fit based on averaged parameters 
vgmOBFAM = vgm(psill = 111211.3, nugget=1338859.4 , range= 112.1159, model= "Sph")
#plot
plot(semidatOBFAM, vgmOBFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Oblateness with Spherical Fit for Flight Group in AM Region")







###LOG SURFACE AREA

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfSAFPL = rbind(fitSA202PL, fitSA208PL, fitSA215PL, 
                   fitSA222PL)

avgdfSAFGPL<- avgdfSAFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSAFPL = rbind(v.cloud202PLSA,v.cloud208PLSA, 
                     v.cloud215PLSA, v.cloud222PLSA)
#create fixed variogram fit based on averaged parameters 
vgmSAFPL = vgm(psill = .12467535, nugget=.059396856 , range= 1648.19, model= "Sph")
#plot
plot(semidatSAFPL, vgmSAFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Surface Area with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfSAFAM = rbind(fitSA202AM, fitSA208AM, fitSA215AM, 
                   fitSA222AM)

avgdfSAFGAM<- avgdfSAFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSAFAM = rbind(v.cloud202AMSA,v.cloud208AMSA, 
                     v.cloud215AMSA, v.cloud222AMSA)
#create fixed variogram fit based on averaged parameters 
vgmSAFAM = vgm(psill = .03982272, nugget=.07781163 , range= 677.1706, model= "Sph")
#plot
plot(semidatSAFAM, vgmSAFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Surface Area with Spherical Fit for Flight Group in AM Region")








###LOG ASPECT RATIO

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfARFPL = rbind(fitAR202PL, fitAR208PL, fitAR215PL, 
                   fitAR222PL)

avgdfARFGPL<- avgdfARFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatARFPL = rbind(v.cloud202PLAR,v.cloud208PLAR, 
                     v.cloud215PLAR, v.cloud222PLAR)
#create fixed variogram fit based on averaged parameters 
vgmARFPL = vgm(psill = .04260005, nugget=.19582008 , range= 120.6893, model= "Sph")
#plot
plot(semidatARFPL, vgmARFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Aspect Ratio with Spherical Fit for Flight Group in PL Region")




#FLIGHT, AM REGION
#get avgeraged parameters
avgdfARFAM = rbind(fitAR202AM, fitAR208AM, fitAR215AM, 
                   fitAR222AM)

avgdfARFGAM<- avgdfARFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatARFAM = rbind(v.cloud202AMAR,v.cloud208AMAR, 
                     v.cloud215AMAR, v.cloud222AMAR)
#create fixed variogram fit based on averaged parameters 
vgmARFAM = vgm(psill = .02637192, nugget=.22509407 , range= 133.2889, model= "Sph")
#plot
plot(semidatARFAM, vgmARFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Aspect Ratio with Spherical Fit for Flight Group in AM Region")






###RANK STRETCH

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfSTFPL = rbind(fitST202PL, fitST208PL, fitST215PL, 
                   fitST222PL)

avgdfSTFGPL<- avgdfSTFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSTFPL = rbind(v.cloud202PLST,v.cloud208PLST, 
                     v.cloud215PLST, v.cloud222PLST)
#create fixed variogram fit based on averaged parameters 
vgmSTFPL = vgm(psill = 585654.5, nugget=764749.8 , range= 119.172, model= "Sph")
#plot
plot(semidatSTFPL, vgmSTFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Stretch with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfSTFAM = rbind(fitST202AM, fitST208AM, fitST215AM, 
                   fitST222AM)

avgdfSTFGAM<- avgdfSTFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSTFAM = rbind(v.cloud202AMST,v.cloud208AMST, 
                     v.cloud215AMST, v.cloud222AMST)
#create fixed variogram fit based on averaged parameters 
vgmSTFAM = vgm(psill = 312578.2, nugget=1143882.7 , range= 118.575, model= "Sph")
#plot
plot(semidatSTFAM, vgmSTFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Stretch with Spherical Fit for Flight Group in AM Region")







###LOG WIDTH

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfWIFPL = rbind(fitWI202PL, fitWI208PL, fitWI215PL, 
                   fitWI222PL)

avgdfWIFGPL<- avgdfWIFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatWIFPL = rbind(v.cloud202PLWI,v.cloud208PLWI, 
                     v.cloud215PLWI, v.cloud222PLWI)
#create fixed variogram fit based on averaged parameters 
vgmWIFPL = vgm(psill = 0.03619624, nugget=0.02989208 , range= 133.5875, model= "Sph")
#plot
plot(semidatWIFPL, vgmWIFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Width with Spherical Fit for Flight Group in PL Region")


#FLIGHT, AM REGION
#get avgeraged parameters
avgdfWIFAM = rbind(fitWI202AM, fitWI208AM, fitWI215AM, 
                   fitWI222AM)

avgdfWIFGAM<- avgdfWIFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatWIFAM = rbind(v.cloud202AMWI,v.cloud208AMWI, 
                     v.cloud215AMWI, v.cloud222AMWI)
#create fixed variogram fit based on averaged parameters 
vgmWIFAM = vgm(psill = 0.02287550, nugget=0.03639531 , range= 151.425, model= "Sph")
#plot
plot(semidatWIFAM, vgmWIFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Width with Spherical Fit for Flight Group in AM Region")













###RANK FITTED ELLIPSOID VOLUME 

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfFEVFPL = rbind(fitFEV202PL, fitFEV208PL, fitFEV215PL, 
                   fitFEV222PL)

avgdfFEVFGPL<- avgdfFEVFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFEVFPL = rbind(v.cloud202PLFEV,v.cloud208PLFEV, 
                     v.cloud215PLFEV, v.cloud222PLFEV)
#create fixed variogram fit based on averaged parameters 
vgmFEVFPL = vgm(psill = 555518.2, nugget=745326.0 , range= 141.6573, model= "Sph")
#plot
plot(semidatFEVFPL, vgmFEVFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Fitted Elippsoid Volume with Spherical Fit for Flight Group in PL Region")




#FLIGHT, AM REGION
#get avgeraged parameters
avgdfFEVFAM = rbind(fitFEV202AM, fitFEV208AM, fitFEV215AM, 
                    fitFEV222AM)

avgdfFEVFGAM<- avgdfFEVFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFEVFAM = rbind(v.cloud202AMFEV,v.cloud208AMFEV, 
                      v.cloud215AMFEV, v.cloud222AMFEV)
#create fixed variogram fit based on averaged parameters 
vgmFEVFAM = vgm(psill = 598672.6, nugget=838468.4 , range= 144.9562, model= "Sph")
#plot
plot(semidatFEVFAM, vgmFEVFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Rank Fitted Elippsoid Volume with Spherical Fit for Flight Group in AM Region")







###LOG LENGTH

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfLEFPL = rbind(fitLE202PL, fitLE208PL, fitLE215PL, 
                    fitLE222PL)

avgdfLEFGPL<- avgdfLEFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatLEFPL = rbind(v.cloud202PLLE,v.cloud208PLLE, 
                      v.cloud215PLLE, v.cloud222PLLE)
#create fixed variogram fit based on averaged parameters 
vgmLEFPL = vgm(psill = 0.03000395, nugget=0.04246808 , range= 1110.51, model= "Sph")
#plot
plot(semidatLEFPL, vgmLEFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Length with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfLEFAM = rbind(fitLE202AM, fitLE208AM, fitLE215AM, 
                   fitLE222AM)

avgdfLEFGAM<- avgdfLEFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatLEFAM = rbind(v.cloud202AMLE,v.cloud208AMLE, 
                     v.cloud215AMLE, v.cloud222AMLE)
#create fixed variogram fit based on averaged parameters 
vgmLEFAM = vgm(psill = 0.008526518, nugget=0.049508395 , range= 191.0087, model= "Sph")
#plot
plot(semidatLEFAM, vgmLEFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Length with Spherical Fit for Flight Group in AM Region")








###LOG HEIGHT

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfHEFPL = rbind(fitHE202PL, fitHE208PL, fitHE215PL, 
                   fitHE222PL)

avgdfHEFGPL<- avgdfHEFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatHEFPL = rbind(v.cloud202PLHE,v.cloud208PLHE, 
                     v.cloud215PLHE, v.cloud222PLHE)
#create fixed variogram fit based on averaged parameters 
vgmHEFPL = vgm(psill = 0.008144364, nugget=0.0442588892 , range= 531.3212, model= "Sph")
#plot
plot(semidatHEFPL, vgmHEFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Height with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfHEFAM = rbind(fitHE202AM, fitHE208AM, fitHE215AM)#, 
                  # fitHE222AM)

avgdfHEFGAM<- avgdfHEFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatHEFAM = rbind(v.cloud202AMHE,v.cloud208AMHE, 
                     v.cloud215AMHE, v.cloud222AMHE)
#create fixed variogram fit based on averaged parameters 
vgmHEFAM = vgm(psill = 0.009410391, nugget=0.051489475 , range= 176.2612, model= "Sph")
#plot
plot(semidatHEFAM, vgmHEFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Height with Spherical Fit for Flight Group in AM Region")








###LOG FITTED ELLIPSOID SURFACE AREA 

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfFSAFPL = rbind(fitFSA202PL, fitFSA208PL, fitFSA215PL, 
                   fitFSA222PL)

avgdfFSAFGPL<- avgdfFSAFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFSAFPL = rbind(v.cloud202PLFSA,v.cloud208PLFSA, 
                     v.cloud215PLFSA, v.cloud222PLFSA)
#create fixed variogram fit based on averaged parameters 
vgmFSAFPL = vgm(psill = 0.05106002, nugget=0.056585509 , range= 790.8396, model= "Sph")
#plot
plot(semidatFSAFPL, vgmFSAFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Fitted Ellipsoid Surface Area with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfFSAFAM = rbind(fitFSA202AM, fitFSA208AM, fitFSA215AM, 
                    fitFSA222AM)

avgdfFSAFGAM<- avgdfFSAFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatFSAFAM = rbind(v.cloud202AMFSA,v.cloud208AMFSA, 
                      v.cloud215AMFSA, v.cloud222AMFSA)
#create fixed variogram fit based on averaged parameters 
vgmFSAFAM = vgm(psill =0.03029992, nugget=0.09030306 , range= 663.8104, model= "Sph")
#plot
plot(semidatFSAFAM, vgmFSAFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Log Fitted Ellipsoid Surface Area with Spherical Fit for Flight Group in AM Region")





###SPHERICITY

#FLIGHT, PL REGION
#get avgeraged parameters
avgdfSPFPL = rbind(fitSP202PL, fitSP208PL, fitSP215PL, 
                    fitSP222PL)

avgdfSPFGPL<- avgdfSPFPL %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSPFPL = rbind(v.cloud222PLFSP,v.cloud208PLFSP, 
                     v.cloud215PLFSP, v.cloud202PLFSP)
#create fixed variogram fit based on averaged parameters 
vgmSPFPL = vgm(psill = 0.000516366, nugget=0.001261995 , range= 370.4296, model= "Sph")
#plot
plot(semidatSPFPL, vgmSPFPL,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Sphericity with Spherical Fit for Flight Group in PL Region")



#FLIGHT, AM REGION
#get avgeraged parameters
avgdfSPFAM = rbind(fitSP202AM, fitSP208AM, fitSP215AM, 
                   fitSP222AM)

avgdfSPFGAM<- avgdfSPFAM %>% group_by(model) %>% 
  summarise_at(c("psill", "range","kappa",
                 "ang1",
                 "ang2",
                 "ang3", "anis1", "anis2"), mean, na.rm = TRUE)
#create semivariogram data 
semidatSPFAM = rbind(v.cloud222AMFSP,v.cloud208AMFSP, 
                     v.cloud215AMFSP, v.cloud202AMFSP)
#create fixed variogram fit based on averaged parameters 
vgmSPFAM = vgm(psill = 0.000437971, nugget=0.001602023 , range=201.7773, model= "Sph")
#plot
plot(semidatSPFAM, vgmSPFAM,xlab="Distance",
     ylab="Semivariance",
     main="Semivariogram for Sphericity with Spherical Fit for Flight Group in AM Region")



































################ SEMIVARIOGRAM BY LOADING, AND STRAIN ###############

#######LOG VOLUME
#######PL FLIGHT
#get PL region and flight group
lacundfPLF = lacundf%>%filter(group==1 &
                                location=="PL")
coordinates(lacundfPLF)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
v.cloudvolPLF<-variogram(log(Volume..µm..)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLF, main = "Variogram - Log Volume: PL Flight")

#######AM FLIGHT
#get AM region and flight group
lacundfAMF = lacundf%>%filter(group==1&
                                location=="AM")
coordinates(lacundfAMF)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
v.cloudvolAMF<-variogram(log(Volume..µm..)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMF, main = "Variogram - AM Flight")

#######PL GROUND
#get PL region and ground group
lacundfPLG = lacundf%>%filter(group==0 &
                                location=="PL")
coordinates(lacundfPLG)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
v.cloudvolPLG<-variogram(log(Volume..µm..)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLG, main = "Variogram - PL Ground")


#######AM GROUND
#get AM region and ground group
lacundfAMG = lacundf%>%filter(group==0 &
                                location=="AM")
coordinates(lacundfAMG)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
v.cloudvolAMG<-variogram(log(Volume..µm..)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMG, main = "Variogram - AM Ground")








#######RANK OBLATENESS

#######PL FLIGHT
v.cloudvolPLFOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFOB, main = "Variogram - Rank Oblateness: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFOB, main = "Variogram - Rank Oblateness: AM Flight")

#######PL GROUND
v.cloudvolPLGOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGOB, main = "Variogram - Rank Oblateness: PL Ground")

#######AM GROUND
v.cloudvolAMGOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGOB, main = "Variogram - Rank Oblateness: AM Ground")





#######LOG SURFACE AREA

#######PL FLIGHT
v.cloudvolPLFSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFSA, main = "Variogram - Log Surface Area: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFSA, main = "Variogram - Log Surface Area: AM Flight")

#######PL GROUND
v.cloudvolPLGSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGSA, main = "Variogram - Log Surface Area: PL Ground")

#######AM GROUND
v.cloudvolAMGSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGSA, main = "Variogram - Log Surface Area: AM Ground")






######LOG ASPECT RATIO

#######PL FLIGHT
v.cloudvolPLFAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFAR, main = "Variogram - Log Aspect Ratio: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFAR, main = "Variogram - Log Aspect Ratio: AM Flight")

#######PL GROUND
v.cloudvolPLGAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGAR, main = "Variogram - Log Aspect Ratio: PL Ground")

#######AM GROUND
v.cloudvolAMGAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGAR, main = "Variogram - Log Aspect Ratio: AM Ground")





#######RANK STRETCH

#######PL FLIGHT
v.cloudvolPLFST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFST, main = "Variogram - Rank Stretch: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFST, main = "Variogram - Rank Stretch: AM Flight")

#######PL GROUND
v.cloudvolPLGST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGST, main = "Variogram - Rank Stretch: PL Ground")

#######AM GROUND
v.cloudvolAMGST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGST, main = "Variogram - Rank Stretch: AM Ground")










#######LOG WIDTH

#######PL FLIGHT
v.cloudvolPLFWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFWI, main = "Variogram - Log Width: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFWI, main = "Variogram - Log Width: AM Flight")

#######PL GROUND
v.cloudvolPLGWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGWI, main = "Variogram - Log Width: PL Ground")

#######AM GROUND
v.cloudvolAMGWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGWI, main = "Variogram - Log Width: AM Ground")






##############RANK FITTED ELLIPSOID VOLUME
#######PL FLIGHT
v.cloudvolPLFFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Flight")

#######PL GROUND
v.cloudvolPLGFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: PL Ground")

#######AM GROUND
v.cloudvolAMGFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: AM Ground")






#######LOG LENGTH
#######PL FLIGHT
v.cloudvolPLFLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFLE, main = "Variogram - Log Length: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFLE, main = "Variogram - Log Length: AM Flight")

#######PL GROUND
v.cloudvolPLGLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGLE, main = "Variogram - Log Length: PL Ground")

#######AM GROUND
v.cloudvolAMGLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGLE, main = "Variogram - Log Length: AM Ground")







#######LOG FITTED ELLIPSOID SURFACE AREA
#######PL FLIGHT
v.cloudvolPLFFESA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFFESA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFFESA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFFESA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Flight")

#######PL GROUND
v.cloudvolPLGFESA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGFESA, main = "Variogram - Log Fitted Ellipsoid Surface Area: PL Ground")

#######AM GROUND
v.cloudvolAMGFESA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGFESA, main = "Variogram - Log Fitted Ellipsoid Surface Area: AM Ground")






#######LOG HEIGHT 

#######PL FLIGHT
v.cloudvolPLFHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFHE, main = "Variogram - Log Height: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFHE, main = "Variogram - Log Height: AM Flight")

#######PL GROUND
v.cloudvolPLGHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGHE, main = "Variogram - Log Height: PL Ground")

#######AM GROUND
v.cloudvolAMGHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGHE, main = "Variogram - Log Height: AM Ground")





#######SPHERICITY

#######PL FLIGHT
v.cloudvolPLFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfPLF, cloud = F)
plot(v.cloudvolPLFSP, main = "Variogram - Sphericity: PL Flight ")

#######AM FLIGHT
v.cloudvolAMFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfAMF, cloud = F)
plot(v.cloudvolAMFSP, main = "Variogram - Sphericity: AM Flight")

#######PL GROUND
v.cloudvolPLGSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfPLG, cloud = F)
plot(v.cloudvolPLGSP, main = "Variogram - Sphericity: PL Ground")

#######AM GROUND
v.cloudvolAMGSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfAMG, cloud = F)
plot(v.cloudvolAMGSP, main = "Variogram - Sphericity: AM Ground")









































################ SEMIVARIOGRAM BY LOADING ###############
lacundfF = lacundf%>%filter(group==1)
lacundfG = lacundf%>%filter(group==0)
coordinates(lacundfF)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")
coordinates(lacundfG)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")


#######LOG VOLUME
#######FLIGHT
v.cloudvolF<-variogram(log(Volume..µm..)~1, data = lacundfF, cloud = F)
plot(v.cloudvolF, main = "Variogram - Log Volume: Flight ")
#######GROUND
v.cloudvolG<-variogram(log(Volume..µm..)~1, data = lacundfG, cloud = F)
plot(v.cloudvolG, main = "Variogram - Log Volume: Ground ")


#######RANK OBLATENESS
#######FLIGHT
v.cloudFOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfF, cloud = F)
plot(v.cloudFOB, main = "Variogram - Rank Oblateness: Flight ")
#######GROUND
v.cloudGOB<-variogram(rank(Lacuna.Oblateness..Untested...NA.)~1, data = lacundfG, cloud = F)
plot(v.cloudGOB, main = "Variogram - Rank Oblateness: Ground ")




#######LOG SURFACE AREA
#######FLIGHT
v.cloudFSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfF, cloud = F)
plot(v.cloudFSA, main = "Variogram - Log Surface Area: Flight ")
#######GROUND
v.cloudGSA<-variogram(log(Surface.Area..Lindblad.2005...µm..)~1, data = lacundfG, cloud = F)
plot(v.cloudGSA, main = "Variogram - Log Surface Area: Ground ")





#######LOG ASPECT RATIO
#######FLIGHT
v.cloudFAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfF, cloud = F)
plot(v.cloudFAR, main = "Variogram - Log Aspect Ratio: Flight ")
#######GROUND
v.cloudGAR<-variogram(log(Aspect.Ratio..NA.)~1, data = lacundfG, cloud = F)
plot(v.cloudGAR, main = "Variogram - Log Aspect Ratio: Ground ")





#######RANK STRETCH
#######FLIGHT
v.cloudFST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfF, cloud = F)
plot(v.cloudFST, main = "Variogram - Rank Stretch: Flight ")
#######GROUND
v.cloudGST<-variogram(rank(Lacuna.Stretch..Untested...NA.)~1, data = lacundfG, cloud = F)
plot(v.cloudGST, main = "Variogram - Rank Stretch: Ground ")




#######LOG WIDTH
#######FLIGHT
v.cloudFWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfF, cloud = F)
plot(v.cloudFWI, main = "Variogram - Log Width: Flight ")
#######GROUND
v.cloudGWI<-variogram(log(Lacuna.Width..Untested...µm.)~1, data = lacundfG, cloud = F)
plot(v.cloudGWI, main = "Variogram - Log Width: Ground ")



#######RANK FITTED ELLIPSOID VOLUME 
#######FLIGHT
v.cloudFFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfF, cloud = F)
plot(v.cloudFFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: Flight ")
#######GROUND
v.cloudGFEV<-variogram(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~1, data = lacundfG, cloud = F)
plot(v.cloudGFEV, main = "Variogram - Rank Fitted Ellipsoid Volume: Ground ")





#######LOG LENGTH
#######FLIGHT
v.cloudFLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfF, cloud = F)
plot(v.cloudFLE, main = "Variogram - Log Length: Flight ")
#######GROUND
v.cloudGLE<-variogram(log(Lacuna.Length..Untested...µm.)~1, data = lacundfG, cloud = F)
plot(v.cloudGLE, main = "Variogram - Log Length: Ground ")





#######LOG HEIGHT
#######FLIGHT
v.cloudFHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfF, cloud = F)
plot(v.cloudFHE, main = "Variogram - Log Height: Flight ")
#######GROUND
v.cloudGHE<-variogram(log(Lacuna.Height..Untested...µm.)~1, data = lacundfG, cloud = F)
plot(v.cloudGHE, main = "Variogram - Log Height: Ground ")



#######LOG FITTED EPPIPSOID SURFACE AREA 
#######FLIGHT
v.cloudFESA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfF, cloud = F)
plot(v.cloudFESA, main = "Variogram - Log Fitted Ellipsoid Surface Area: Flight ")
#######GROUND
v.cloudGESA<-variogram(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~1, data = lacundfG, cloud = F)
plot(v.cloudGESA, main = "Variogram - Log Fitted Ellipsoid Surface Area: Ground ")



#######SPHERICITY 
#######FLIGHT
v.cloudFSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfF, cloud = F)
plot(v.cloudFSP, main = "Variogram - Sphericity: Flight ")
#######GROUND
v.cloudGSP<-variogram(Sphericity..Untested...NA.~1, data = lacundfG, cloud = F)
plot(v.cloudGSP, main = "Variogram - Sphericity: Ground ")


















































#############LINEAR MULTILEVEL MODELING FOR THE 4 GROUPS AS DISCUSSED WITH ERIC#################

#####LOG VOLUME
###one: varying intercept 
lmevolume1G = lmer(log(Volume..µm..)~Interaction+(1|id),data=lacundf)
summary(lmevolume1G)
display(lmevolume1G)
confint(lmevolume1G,level=0.95)
###two: varying intercept and slope , DIVERGENT... 
lmevolume2G = lmer(log(Volume..µm..)~Interaction+(1+Interaction|id),data=lacundf)
summary(lmevolume2G)
display(lmevolume2G)
confint(lmevolume2G,level=0.95)


#####RANK OBLATENESS
###one: varying intercept 
lmeRO = lmer(rank(Lacuna.Oblateness..Untested...NA.)~Interaction+(1|id),data=lacundf)
summary(lmeRO)
display(lmeRO)
confint(lmeRO,level=0.95)


#####LOG SURFACE AREA 
###one: varying intercept 
lmeSA = lmer(log(Surface.Area..Lindblad.2005...µm..)~Interaction+(1|id),data=lacundf)
summary(lmeSA)
display(lmeSA)
confint(lmeSA,level=0.95)


#####LOG ASPECT RATIO
###one: varying intercept 
lmeAR = lmer(log(Aspect.Ratio..NA.)~Interaction+(1|id),data=lacundf)
summary(lmeAR)
display(lmeAR)
confint(lmeAR,level=0.95)


#####RANK STRETCH
###one: varying intercept 
lmeST = lmer(rank(Lacuna.Stretch..Untested...NA.)~Interaction+(1|id),data=lacundf)
summary(lmeST)
display(lmeST)
confint(lmeST,level=0.95)


#####LOG WIDTH
###one: varying intercept 
lmeWI = lmer(log(Lacuna.Width..Untested...µm.)~Interaction+(1|id),data=lacundf)
summary(lmeWI)
display(lmeWI)
confint(lmeWI,level=0.95)


#####RANK FITTED ELLIPSOID VOLUME
###one: varying intercept 
lmeEV = lmer(rank(Fitted.Ellipsoid.Volume..Untested...µm..)~Interaction+(1|id),data=lacundf)
summary(lmeEV)
display(lmeEV)
confint(lmeEV,level=0.95)


#####LOG LENGTH
###one: varying intercept 
lmeLE = lmer(log(Lacuna.Length..Untested...µm.)~Interaction+(1|id),data=lacundf)
summary(lmeLE)
display(lmeLE)
confint(lmeLE,level=0.95)


#####LOG HEIGHT
###one: varying intercept 
lmeHE = lmer(log(Lacuna.Height..Untested...µm.)~Interaction+(1|id),data=lacundf)
summary(lmeHE)
display(lmeHE)
confint(lmeHE,level=0.95)


#####LOG FITTED ELLIPSOID SURFACE AREA 
###one: varying intercept 
lmeFSA = lmer(log(Fitted.Ellipsoid.Surface.Area..Untested...µm..)~Interaction+(1|id),data=lacundf)
summary(lmeFSA)
display(lmeFSA)
confint(lmeFSA,level=0.95)


#####SPHERICITY
###one: varying intercept 
lmeSP = lmer(Sphericity..Untested...NA.~Interaction+(1|id),data=lacundf)
summary(lmeSP)
display(lmeSP)
confint(lmeSP,level=0.95)





















##########messing around and truncating data for semi

######################MOUSE 222
lacundfGPL222= lacundf%>%filter(id=="222"&
                                  location =="PL")
range(lacundfGPL222$Center.Of.Mass.Z..µm.)
-673: -600
-3822:-3800
-769:-700
lacundfGPL222g = lacundfGPL222 %>% filter(between(Center.Of.Mass.X..µm., -673, -650))

coordinates(lacundfGPL222g)=c("Center.Of.Mass.X..µm.", "Center.Of.Mass.Y..µm.","Center.Of.Mass.Z..µm.")

####LOG VOLUME
#PL
v.cloudvol222PLmess<-variogram(log(Volume..µm..)~1, data = lacundfGPL222g, cloud = F)
plot(v.cloudvol222PLmess, main = "Variogram - Log Volume: PL Mouse 222 (FLIGHT)")
###fit empirical variogram to theoretical variogram
#now fitting 
#nugget= first x axis point
#sill = #first maximum 
#range = #distance at which sill occurs (x value)
modelvol222PLsph = vgm(psill = .14, nugget=15 , range= 200, model= "Sph")
#now fit
fitvol222PL = fit.variogram(object= v.cloudvol222PL, model =modelvol222PLsph )
plot(v.cloudvol222PL, model=fitvol222PL)


