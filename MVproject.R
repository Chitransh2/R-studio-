#*********** GROUP PROJECT************
#GROUP MEMBERS: 
#CHITRANSH JOSHI(02002659)
#GARIMA SHAH (01994997)
#MARISHA MAHADEVIA (01979594) 
#SHREYA PATEL (01988619)

#******** SEER Breast Cancer Dataset **************

#reading the dataset
cancerData <- read.csv("breast_cancer.csv")
head(cancerData)

#dimension of data
dim(cancerData)

#structure of data
str(cancerData)

#summary of all variables
summary(cancerData)

#creating categorical variable for Race
cancerData$raceType <- as.factor(ifelse(cancerData$Race =='Other (American Indian/AK Native, Asian/Pacific Islander)', 0,
                                 ifelse(cancerData$Race =='Black',1,2)))

#creating categorical variable for Tumor Stage
cancerData$TStageType <- as.factor(ifelse(cancerData$T.Stage == 'T1',1,
                                   ifelse(cancerData$T.Stage == 'T2',2,3)))

#finding missing values
is.na(cancerData)
#there are no missing values

#Data visualizations
library(ggplot2)

#Age and Tumor Size plot
ggplot(cancerData, aes(x=Age, y=Tumor.Size)) + geom_point() +ggtitle("Tumor sizes and Age plot")

#Race and Tumor Size plot
ggplot(cancerData, aes(x=raceType, y=Tumor.Size)) + geom_point() + ggtitle("Tumor sizes for other, Black and White people")

#Tumor Stage and Tumor Size plot
plot(cancerData$Tumor.Size ~ cancerData$TStageType, xlab = "Tumor Stage ", ylab = "Tumor Size", main="Tumor sizes based on Tumor Stage")

#simple linear regression model
#effect of age on Tumor size
library(lmtest)
library(sandwich)
linRegr1 <- lm(Tumor.Size~Age, data=cancerData)
coeftest(linRegr1, vcov. = vcovHC, type = "HC1")
#so if age increase by 1 unit then tumor size will decrease by 0.1819 millimeters
summary(linRegr1)

#change in effect on Tumor size after adding raceType factor
linRegr2 <- lm(Tumor.Size~Age+raceType, data=cancerData)
coeftest(linRegr2, vcov. = vcovHC, type = "HC1")
summary(linRegr2)
X#so if age increases by 1 unit then tumor size will decrease by 0.122 mm (not much difference from previous model)
#i.e.mean of tumor size for black people(raceType1) is 0.27 units less than the mean of tumor size of other people(raceType0)
#i.e.mean of tumor size for white people(raceType2) is 0.04 units less than the mean of tumor size of other people(raceType0)


#change in effect on Tumor size after adding Tumor Stage factor
linRegr3 <- lm(Tumor.Size~Age+raceType+TStageType, data=cancerData)
coeftest(linRegr3, vcov. = vcovHC, type = "HC1")
summary(linRegr3)
#if age increase by 1 unit then tumor size will decrease by 0.035 mm
#i.e mean of tumor size for stage 2(TStageType2) is 17.1 units more than the mean of tumor size of stage 1(TSTageType1)
#i.e. mean of tumor size for stage3(TStageType3) is 54.07 units more than the mean of tumor size of stage 1(TStageTyep1)

robustEE <- list(sqrt(diag(vcovHC(linRegr1, type = "HC1"))),
                 sqrt(diag(vcovHC(linRegr2, type = "HC1"))),
                 sqrt(diag(vcovHC(linRegr3, type = "HC1"))))

library(stargazer)
stargazer(linRegr1, linRegr2, linRegr3, 
          type="html",
          digits = 3,
          se = robustEE,
          dep.var.labels=c("Tumor.Size"),
          covariate.labels=c("Age","raceType1","raceType2","TStageType2","TStageType3"),
          ci = T,
          out="models1.htm")   

#concluding that:
#from all three models above, the model 3 has higher adjusted r square 
#So, we prefer regression model 3


#forming logarithmic models

#why
hist(cancerData$Tumor.Size)
hist(log(cancerData$Tumor.Size))
4#after the log transformation the histogram becomes more symmetric
#since we perform the analysis assuming normality, log transformation will help us in meeting this assumption
                    
#linear-log model
cancerData$lnAge = log(cancerData$Age)
linRegr4 <- lm(Tumor.Size~lnAge+raceType+TStageType, data=cancerData)
coeftest(linRegr4, vcov. = vcovHC, type = "HC1")
#if age increases by 1% then tumor size decreases by 0.019 millimeter
summary(linRegr4)

#log-linear model
cancerData$lnSize = log(cancerData$Tumor.Size)
linRegr5 <- lm(lnSize~Age+raceType+TStageType, data=cancerData)
coeftest(linRegr5, vcov. = vcovHC, type = "HC1")
#if the age increases by 1 unit then the tumor size will decrease by 0.13%
summary(linRegr5)

#log-log model
linRegr6 <- lm(lnSize~lnAge+raceType+TStageType, data=cancerData)
coeftest(linRegr6, vcov. = vcovHC, type = "HC1")
#if age changes by 1%, then tumor size is expected to change by 0.07%
summary(linRegr6)

#quadratic model
cancerData$AgeSq2 = cancerData$Age^2
linRegr7 <- lm(lnSize~Age+AgeSq2+raceType+TStageType, data=cancerData)
coeftest(linRegr7, vcov. = vcovHC, type = "HC1")
summary(linRegr7)

#building the stargazer table
robustEE <- list(sqrt(diag(vcovHC(linRegr4, type = "HC1"))),
                 sqrt(diag(vcovHC(linRegr5, type = "HC1"))),
                 sqrt(diag(vcovHC(linRegr6, type = "HC1"))),
                 sqrt(diag(vcovHC(linRegr7, type = "HC1"))))

library(stargazer)
stargazer(linRegr4, linRegr5, linRegr6, linRegr7,
          type="html",
          digits = 3,
          se = robustEE,
          dep.var.labels=c("Tumor.Size","lnSize"),
          covariate.labels=c("lnAge","Age","AgeSq2","raceType","TStageType"),
          ci = T,
          out="models2.htm")   


#so, we conclude that adjusted-R^2 for linear-log model is higher compared to others
#SO, we can prefer model 4 over others
#from the stargazer for model 4 we can see that Age is significant at 5* and 10% significance level
#TstageType2, TStageType3 are significant at both 5% and 10% significance level
#raceType is not significant at 5% and 10% significance level

#interaction term
library(fastDummies)
cancerData <- dummy_cols(cancerData, select_columns='T.Stage')
head(cancerData)
linRegr7 <- lm(Tumor.Size~Age+T.Stage_T2+T.Stage_T3+Age*T.Stage_T2+Age*T.Stage_T3, data=cancerData)
coeftest(linRegr7, vcov. = vcovHC, type = "HC1")
summary(linRegr7)
#so, -0.004 represents that there is 0.004 more effect of age on tumor size of Stage 2 than on Stage 1 
#so, -0/04 represents that there is 0.4 more effect of age on tumor size of Stage 3 then on Stage 1
#both the interaction terms are insignificant at 5% and 10% significance value

#so our final model suggests that age and whether the person has Tumor stage T2 and T3 might affect the increase or decrease in Tumor size
linRegr8 <- lm(Tumor.Size~lnAge+T.Stage_T2+T.Stage_T3, data=cancerData)
coeftest(linRegr8, vcov. = vcovHC, type = "HC1")

#hypothesis
#H0: There is no effect of age and Stage on tumor size
#H1: There is effect of age anf Stage on tumor size
#Conclusion;
#p-value for age is 0.07 which is significant at 10% 
#T.Stage_T2 and T.Stage_T3 are significant at 10% significance level
#Thus, we can reject the null hypothesis stating that atleast one of the coefficient might be zero
#thus there are relationships between them

#our final model 
#tumor size = 24.97 - 2.04*ln(Age) + 14.60*TStageType2 + 53.90*TStageType3

#interpretation for our final model
#Holding all other factors fixed, if age increase by 1% unit then the tumor size decreases by 2.04 mm
#Holding all other factors fixed, Stage 2 people have 14.60mm increase in tumor size on average 
#Holding all other factors fixed, Stage 3 people have 53.60 increase in tumor size on average
#the average tumor size will be 24.98 when age and stage 2 and 3 factor are zero, but this doesn't make sense as they will never be 0
