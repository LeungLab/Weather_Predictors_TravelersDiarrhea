#Load libraries-----------------------------------------------------------------
library(readr)
library(tidyverse)
library(ggplot2)
library(randomForest)
library(pROC)
library(cvAUC)
library(viridis)
# setwd


#Data - cleaned, unfiltered-----------------------------------------------------
dat1 <- read.csv("Data_demographics_symptoms_weather_micro.csv", header = TRUE)
#Please contact the authors for access to data if needed

sapply(dat1,function(x) sum(is.na(x)))

#Filtering data to compare bact only vs viral/non-detected ---------------------
datC = dat1 %>% 
  filter(Viral.infection.only == 1 | Bacterial.infection.only == 1 | No.bacteria.or.virus.identified == 1) %>%
  select(-c(Study.ID, Viral.infection.only, Other..Protozoa.parasite.infection, 
            Mixed.bacterial.viral.parasitic.infection, No.bacteria.or.virus.identified, Shigella.positive) )


#Identifying predictors and importance using randomForest ----------------------

RF=randomForest(x=(dplyr::select(datC,-Bacterial.infection.only,-contains("Age."),-contains("duration."),-contains("Diarrhea."))), y=(datC$Bacterial.infection.only), importance = TRUE,ntree=1000)

str(datC)
imps=importance(RF)
df_imps=data.frame(variables=row.names(imps),IncMSE=imps[,1]) %>%
  arrange(desc(IncMSE))

#Identify top 10 predictors in descending order
head(df_imps, 10)
 

#Cross Validation using Logistic regression------------------------------------- 
datC$index=1:dim(datC)[1]

niter=100  
var_opts=c(1:10, 15)
result=data.frame(No.Vars=NA,Iter=NA,Truth=NA,LR=NA,RF=NA,int_lr=NA,slope_lr=NA,int_rf=NA,slope_rf=NA)

#Start Loop
for(each in 1:niter){
  print(each)
  train=datC %>% sample_frac(.8,replace=F)
  test=datC[-train$index,]
  
  form = as.formula(paste("Bacterial.infection.only~",paste(colnames(train)[-(c(3:8,36:37))], collapse="+"), sep=""))
  RF=randomForest(Bacterial.infection.only~., importance = T,ntree=1000,data=train %>% select(-contains("Age."),-contains("duration."),-contains("Diarrhea."),-index))
  
  
  df_imps=data.frame(variables=row.names(importance(RF)),IncMSE=importance(RF)[,1]) %>%
    arrange(desc(IncMSE))
  
  for(each2 in var_opts){
    form = as.formula(paste("Bacterial.infection.only~",paste(df_imps$variables[1:each2], collapse="+"), sep=""))
    
    modtrain = glm(form, family = "binomial", data = train)
    modpred1 = predict(modtrain, test, type = "response")
    
    #Calibration code
    cali_num=predict(modtrain, test)
    int_lr=summary(glm(Bacterial.infection.only~1,offset=cali_num,data=test,family="binomial"))$coef[1,1]
    slope_lr=summary(glm(Bacterial.infection.only~cali_num,data=test,family="binomial"))$coef[2,1]
     
    
    formRF = as.formula(paste("(Bacterial.infection.only)~",paste(df_imps$variables[1:each2], collapse="+"), sep=""))
    modtrainRF = randomForest(formRF, train, ntree = 1000)
    modpredRF = predict(modtrainRF, test)
    modpredRF = ifelse(modpredRF==1,.999,modpredRF)
    logit_predRF=log(modpredRF/(1-modpredRF))
    
    int_rf=summary(glm(Bacterial.infection.only~1,offset=logit_predRF,data=test,family="binomial"))$coef[1,1]
    slope_rf=summary(glm(Bacterial.infection.only~logit_predRF,data=test,family="binomial"))$coef[2,1]
   
    result=rbind(result,data.frame(No.Vars=each2,Iter=each,Truth=test$Bacterial.infection.only,LR=modpred1,RF=modpredRF,int_lr=int_lr,slope_lr=slope_lr,
                                   int_rf=int_rf,slope_rf=slope_rf))
  }}
result=result[-1,]

with(result %>% filter(No.Vars==5,Iter==1),roc(Truth,RF,ci=T)$ci[2])

result

#calibration calc-------------------------------------------------------------
result %>% group_by(No.Vars,Iter) %>% summarize(int_lr=mean(int_lr),slope_lr=mean(slope_lr),int_rf=mean(int_rf),slope_rf=mean(slope_rf)) %>% 
  group_by(No.Vars) %>% summarize(int_lr=paste0(round(mean(int_lr),2), "(",round(sd(int_lr),2),")"),
                                  slope_lr=paste0(round(mean(slope_lr),2), "(",round(sd(slope_lr),2),")"),
                                  int_rf=paste0(round(mean(int_rf),2), "(",round(sd(int_rf),2),")"),
                                  slope_rf=paste0(round(mean(slope_rf),2), "(",round(sd(slope_rf),2),")"))


# Make AUC trend curves--------------------------------------------------------

dat1=result %>% 
  group_by(No.Vars,Iter) %>% dplyr::summarize(         LR=roc(.$Truth,.$LR,ci=T)$ci[2],
                                                       LR_lo=roc(.$Truth,.$LR,ci=T)$ci[1],
                                                       LR_hi=roc(.$Truth,.$LR,ci=T)$ci[3],
                                                       RF=roc(.$Truth,.$RF,ci=T)$ci[2],
                                                       RF_lo=roc(.$Truth,.$RF,ci=T)$ci[1],
                                                       RF_hi=roc(.$Truth,.$RF,ci=T)$ci[3]) %>%
  ungroup() %>% group_by(No.Vars) %>% 
  summarize(LR=mean(LR),LR_lo=mean(LR_lo),LR_hi=mean(LR_hi),RF=mean(RF),RF_lo=mean(RF_lo),RF_hi=mean(RF_hi))


AUCs=bind_rows(result) %>% split(.$No.Vars) %>% purrr::map(~ci.cvAUC(.$LR,.$Truth,folds=.$Iter))
AUCs2=bind_rows(result) %>% split(.$No.Vars) %>% purrr::map(~ci.cvAUC(.$RF,.$Truth,folds=.$Iter))

AUC_df=rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
             bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
AUC_df$nvar=var_opts


ggplot(AUC_df ,aes(x=nvar,y=AUC,group=Model,color=Model)) + geom_point() + geom_line() + 
  geom_ribbon(aes(x=nvar, y=AUC, ymax=upper, ymin=lower,group=Model,color=Model),alpha=0.2,linetype="dashed") +
  scale_x_continuous(breaks=seq(1,15,by=1),labels=seq(1,15,by=1)) +
  scale_y_continuous(breaks=seq(.5,.80,by=.1),labels=seq(.5,.80,by=.1),limits=c(.5,.80)) +
  scale_color_viridis(discrete=TRUE,end=.75,labels=c("LR","RF")) + theme_bw()+ xlab("No. of Variables") + ggtitle("Model 3: AUC for LR and RF models at different number of variables") + theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(colour = "transparent", fill = "transparent"))


ggplot(data=dat1 %>% gather(Type,AUC,-No.Vars),aes(x=No.Vars,y=AUC,color=Type,group=Type)) + geom_line(size=1.01) + geom_point(size=2) + 
  theme(legend.position="bottom",legend.title=element_blank()) + scale_y_continuous(breaks=seq(.56,.7,by=.02),labels=seq(.56,.7,by=.02),limits=c(.56,.7)) +
  xlab("Number of variables") +
  ggtitle("Model 3: AUC for LR and RF models at different number of variables") + theme_bw() #+ scale_color_viridis_d(begin=.15,end=.85)

# Makes ROC curves
tst=purrr::map(result %>% filter(No.Vars %in% c(5,10,15)) %>% split(.$No.Vars),~data.frame(sens=roc(.$Truth,.$LR)$sensitivities,spec=roc(.$Truth,.$LR)$specificities))
tst2=purrr::map(result %>% filter(No.Vars %in% c(5,10,15)) %>% split(.$No.Vars),~data.frame(sens=roc(.$Truth,.$RF)$sensitivities,spec=roc(.$Truth,.$RF)$specificities))

#tst=df2glm %>% filter(novars %in% c(5,10,20)) %>% split(.$novars)
approxs=purrr::map(tst,~approx(1-.$spec,.$sens,method="linear"))
plts=bind_rows(purrr::map(approxs,~bind_cols(.)),.id='No. of Variables')
#a=ggplot(data=plts,aes(x=x,y=y,group=NoVars,color=NoVars)) + geom_line() + xlim(c(0,1))+ylim(c(0,1)) +
#  ylab("Sensitivity") + xlab("1 - Specificity")
#tst2=df2rf %>% filter(novars %in% c(5,10,20)) %>% split(.$novars)
approxs=purrr::map(tst2,~approx(1-.$spec,.$sens,method="linear"))
plts2=bind_rows(purrr::map(approxs,~bind_cols(.)),.id='No. of Variables')
plt=rbind(plts,
          plts2)
plt$Model=c(rep("LR",150),rep("RF",150))
plt$`No. of Variables`=factor(plt$`No. of Variables`)
plt$`No. of Variables`=factor(plt$`No. of Variables`,levels(plt$`No. of Variables`)[c(3,1,2)])


ggplot(data=plt,aes(x=x,y=y,group=interaction(`No. of Variables`,Model),color=Model,linetype=`No. of Variables`)) + 
  geom_line() +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=.05),labels=seq(0,1,by=.05),expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=.05),labels=seq(0,1,by=.05),expand = c(0, 0)) +
  ylab("Sensitivity") + xlab("1 - Specificity") + theme_bw()  + 
  scale_linetype_manual(name="No. of Variables",values=c('solid','dashed','dotted')) + 
  theme(axis.text.x=element_text(angle=270),axis.title.y=element_blank(),axis.title.x=element_blank(),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),legend.position="bottom") + ggtitle("Model 3: AUC vs. number of variables")


mod=glm(as.formula(paste("Bacterial.infection.only ~",paste(df_imps$variables[1:11],collapse="+"),sep="")),data=datC);summary(mod)

mods=summary(mod)
df_table=data.frame(Variable=as.character(rownames(mods$coefficient)),
                    Estimate=mods$coefficients[,1],
                    SE=mods$coefficients[,2],
                    z=mods$coefficients[,3],
                    Pvalue=mods$coefficients[,4])

df_table2=df_table %>% #[1:8,] %>%
  mutate(OddsRatio=round(exp(Estimate),2),Low=round(exp(Estimate-qnorm(.975)*SE),2),High=round(exp(Estimate+qnorm(.975)*SE),2),
         Odds=paste(OddsRatio," (",Low," - ",High,")",sep=""),
         P.value=round(Pvalue,3),)  %>% select(Variable,Odds, P.value)


