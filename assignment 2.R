library(dplyr)
library(tidyverse)
library(ggplot2)
assignment_data=read.csv("Assignment_dataset (2).csv")
attach(assignment_data)
##Pre and post no treatment arm division
  #IL1B
hist(IL1B_pre , prob=T)
curve(dnorm(x,mean=mean(IL1B_pre), sd=sd(IL1B_pre)), add=T, col="red")
qqnorm(IL1B_pre,main = 'IL1B_pre')
qqline(IL1B_pre, col='red')
shapiro.test(IL1B_pre)
?qqline
IL1B_post_Tras<-(IL1B_post[49:88])

plot(IL1B_post,IL1B_pre)

boxplot(IL1B_pre,IL1B_post,main='IL1B',names=c('IL1B Pre','IL1B Post'))

hist(IL1B_post, prob=T)
curve(dnorm(x,mean=mean(IL1B_post), sd=sd(IL1B_post)), add=T, col="red")
qqnorm(IL1B_post,main = 'IL1B_post')
qqline(IL1B_post, col='red')
shapiro.test(IL1B_post)
boxplot(IL1B_post)

wilcox.test(IL1B_pre,IL1B_post,paired=TRUE)

  #CX3Cl1
hist(CX3CL1_pre , prob=T)
curve(dnorm(x,mean=mean(CX3CL1_pre), sd=sd(CX3CL1_pre)), add=T, col="red")
qqnorm(CX3CL1_pre,main = 'CX3CL1_pre')
qqline(CX3CL1_pre, col='red')
shapiro.test(CX3CL1_pre)

hist(CX3CL1_post , prob=T)
curve(dnorm(x,mean=mean(CX3CL1_post), sd=sd(CX3CL1_post)), add=T, col="red")
qqnorm(CX3CL1_post,main = 'CX3CL1_post')
qqline(CX3CL1_post, col='red')
shapiro.test(CX3CL1_post)
boxplot(CX3CL1_post)
head(CX3CL1_post)
CX3CL1_post
?boxplot
count_higher<-assignment_data %>%
  count(Treatment_Arm,CX3CL1_post>3.5) 
count_respon

wilcox.test(CX3CL1_pre,CX3CL1_post,paired=TRUE)
boxplot(CX3CL1_pre,CX3CL1_post,main='CX3CL1',names=c('CX3CL1_pre','CX3CL1_post'))


  #TNFA
hist(TNFA_pre , prob=T)
curve(dnorm(x,mean=mean(TNFA_pre), sd=sd(TNFA_pre)), add=T, col="red")
qqnorm(TNFA_pre, main='TNFA_pre')
qqline(TNFA_pre, col='red')
shapiro.test(TNFA_pre)

hist(TNFA_post , prob=T)
curve(dnorm(x,mean=mean(TNFA_post), sd=sd(TNFA_post)), add=T, col="red")
qqnorm(TNFA_post, main='TNFA_post')
qqline(TNFA_post, col='red')
shapiro.test(TNFA_post)

wilcox.test(TNFA_pre,TNFA_post,paired=TRUE)
boxplot(TNFA_pre,TNFA_post,main='TNFA',names=c('TNFA_pre','TNFA_post'))

  #CCL20

hist(CCL20_pre , prob=T)
curve(dnorm(x,mean=mean(CCL20_pre), sd=sd(CCL20_pre)), add=T, col="red")
qqnorm(CCL20_pre, main='CCL20_pre')
qqline(CCL20_pre, col='red')
shapiro.test(CCL20_pre)

hist(CCL20_post , prob=T)
curve(dnorm(x,mean=mean(CCL20_post), sd=sd(CCL20_post)), add=T, col="red")
qqnorm(CCL20_post, main='CCL20_post')
qqline(CCL20_post, col='red')
shapiro.test(CCL20_post)

wilcox.test(CCL20_pre,CCL20_post,paired=TRUE)
boxplot(CCL20_pre,CCL20_post,main='CCL20',names=c('CCL20_pre','CCL20_post'))



##Q2
count_respon<-assignment_data %>%
  count(BMI_pre_treatment<20,BMI_post_treatment<20) 
count_respon

plot(count_respon)

shapiro.test(BMI_pre_treatment)

TreatmentTestData<- matrix(c(5,16,2,65), nr=2, dimnames= list("Underweight"= c("Before_yes", "Before_no") ,  "Underweight"=c("After_yes","After_no")))
TreatmentTestData

Underweight<-assignment_data %>%
  select(BMI_pre_treatment,BMI_post_treatment)%>%
filter(Underweight,BMI_post_treatment=<20)

?filter

mcnemar.test(TreatmentTestData)
table(BMI_post_treatment<20)
cor.test(BMI_pre_treatment,BMI_post_treatment, method = 'pearson')

library(RBioinf)
library(ape)
library(rcompanion)
library(PMCMRplus)
###Q3
stacked<-assignment_data%>%
  bind_cols(c(IL1B_pre,IL1B_during))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("reshape2")
a

data_long <- gather(assignment_data, variable, measurement, IL1B_pre:CCL20_post, factor_key=TRUE)
data_long

data_long2 <- gather(assignment_data, Pre_var, Pre_measurement, IL1B_pre:CCL20_pre, factor_key = TRUE)

data_long2 <- gather(data_long2, during_var, during_measurement, IL1B_during:CCL20_during, factor_key = TRUE)

data_long2<-  gather(data_long2,post_var,post_measure, IL1B_post:CCL20_post,factor_key = TRUE)

?gather

data_long3<-pivot_longer(assignment_data,IL1B_pre:CCL20_pre,IL1B_during:CCL20_during,IL1B_post:CCL20_post,names_to = Pre,pre_measure,During,during_measure,post,post_measure)

?vignette("pivot")

data_long3<-pivot_longer(assignment_data,7:10,11:14,15:18,names_to = 'Pre','During','post',values_to = 'Pre-measure','during_measure','post_measure',names_sep = ',')


friedman.test(data_long2$during_measurment,data_long2$during_var,data_long2$Patient_number)

friedman.test(during_measurement ~ during_var | Patient_number, data = data_long2)
 



Data_pre<-assignment_data%>%
  select(c(Patient_number,IL1B_pre:CCL20_pre))%>%
  pivot_longer(2:5,names_to = 'Pre_Var',values_to = 'Pre_measure')

friedman.test(Pre_measure ~ Pre_Var | Patient_number, data = Data_pre)

frdAllPairsConoverTest(Data_pre$Pre_measure,Data_pre$Pre_Var,Data_pre$Patient_number)

Data_during<-assignment_data%>%
  select(c(Patient_number,IL1B_during:CCL20_during))%>%
  pivot_longer(2:5,names_to = 'during_Var',values_to = 'during_measure')

friedman.test(during_measure ~ during_Var | Patient_number, data = Data_during)

frdAllPairsConoverTest(Data_during$during_measure,Data_during$during_Var,Data_during$Patient_number)




Data_post<-assignment_data%>%
  select(c(Patient_number,IL1B_post:CCL20_post))%>%
  pivot_longer(2:5,names_to = 'post_Var',values_to = 'post_measure')

friedman.test(post_measure ~ post_Var | Patient_number, data = Data_post)

frdAllPairsConoverTest(Data_post$post_measure,Data_post$post_Var,Data_post$Patient_number)




Data_CCL20<-assignment_data%>%
  select(c(Patient_number,CCL20_during,CCL20_post,CCL20_pre))%>%
  pivot_longer(2:4,names_to = 'CC_variable',values_to = 'CC_measure')


friedman.test(CC_measure ~ CC_variable | Patient_number, data = Data_CCL20)

frdAllPairsConoverTest(Data_CCL20$CC_measure,Data_CCL20$CC_variable,Data_CCL20$Patient_number)



?friedman.test              
              
maybe<-frdAllPairsConoverTest(data_long$measurement,data_long$variable,data_long$Patient_number)
