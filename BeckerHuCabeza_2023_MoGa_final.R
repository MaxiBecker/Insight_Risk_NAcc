### R Script to analyse data for manuscript: MoGA project for Paper: 
### Becker, Yu & Cabeza, 2023

##  study 1: replication of Yu et al. - Solving problems with an Aha! biases the choice over monetary rewards
##  study 2: fMRI study (taken from Becker, Sommer & Cabeza, 2023)
## AHA response  (no AHA! moment) 1 <--- 4 ---> 7 (strong AHA! moment)

rm(list=ls())

### load packages #### 
library(Matrix)
library(lme4)
library(sjPlot)
library(ggstatsplot)
library(emmeans)
library(glmmTMB)
library(lmerTest)
#library(effectsize)
library(performance)
library(see)
library(ggplot2)
library(tidyverse)
library(ggeffects)
#library(effectsize)
library(dplyr)

### Functions ###############################################################
  plot_histogram <- function(df, feature) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)))) +
      geom_histogram(aes(y = ..density..), alpha=0.7, fill="#33AADE", color="black") +
      geom_density(alpha=0.3, fill="red") + theme_classic()+
      geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
      labs(x=feature, y = "Density")
    print(plt)
  }
  
  plot_multi_histogram <- function(df, feature, label_column) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
      geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
      geom_density(alpha=0.7) + theme_classic()+
      geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
      labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column))
  }


#### Study path ##########################
  
setwd('C:/Users/Maxi/Google Drive/studies/Maxi/MoGa/GitHub_4_publication')
  
  
### A) STUDY 1 ##############################################################

#### 1)  load data ####

#150 subjects for Exp 1
load("Study1_Behave_data.Rda")

# check demographics before excluding subjects
  s1_gender <- data %>% group_by(subject) %>% dplyr::count(sex) 
  library(stringr)
  s1_female_n<- sum(str_count(s1_gender$sex,"F"), na.rm = T)
  s1_male_n<- sum(str_count(s1_gender$sex,"M"),na.rm = T)
  
  max(data$age)
  min(data$age)
  
  mean(data[data$sex == 'F',]$age, na.rm = T)
  mean(data[data$sex == 'M',]$age, na.rm =T)

#### 2)  Excluding subjects ########################################################  

# based on zero variance in aha experience one subject had to be excluded from further analyses
data = data[data$subject != "A14NP6X071S7GK",]

#### filter out non-switchers between risk and no risk behavior 

    s1_dat1 <- data %>% group_by(subject) %>% dplyr::count(subject) 
    
    data_summary <- function(data, n, varname, groupnames){
      require(plyr)
      summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE),
          se = sd(x[[col]], na.rm=TRUE)/sqrt(n))}
      data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
      data_sum <- rename(data_sum, c("mean" = varname))
      return(data_sum)}
    
    n = length(unique(data$subject))
    s1_IDdata<- data_summary(data, n, varname="riskchoice", groupnames=c("subject"))
    
    s1_IDreinX= s1_IDdata[s1_IDdata$sd > 0.001,]
    s1_IDraus= s1_IDdata[s1_IDdata$sd < 0.001,]$subject # those 53 subjects additionally have to be excluded
    
    `%notin%` <- function(x,y) !(x %in% y) 
    newdata = data
    newdata = newdata[newdata$subject %notin% s1_IDraus ,]
    newdata$sex_fac = as.factor(newdata$sex)

    
###################################################################################-  

#### 3a) Descriptives & view variables #####
  summary(newdata)
  str(newdata)
  s1_behave_data0 = newdata[,c('subject', 'cor2_num', 'cor1', 'RT', 'RT_correct', 'riskchoice', 'gamble_RT',
                               'AHA_response', 'insight_mediansplit_num', 'insight_mediansplit_num_correct')]
  s1_behave_data1 = s1_behave_data0 %>% group_by(subject) %>% summarise_all(funs(mean), na.rm =T)   
  
  # solution related descriptive measures
  mean(s1_behave_data1$cor1)*100
  sd(s1_behave_data1$cor1)*100
  
  mean(s1_behave_data1$cor2_num)*100
  sd(s1_behave_data1$cor2_num)*100
  
  mean(s1_behave_data1$RT, na.rm =T)/1000 
  sd(s1_behave_data1$RT, na.rm = T)/1000 
  
  mean(s1_behave_data1$RT_correct, na.rm =T)/1000 #note for one subject RT data were not recorded for some reason, data seem legit though
  sd(s1_behave_data1$RT_correct, na.rm =T)/1000
  
  mean(s1_behave_data1$insight_mediansplit_num)*100
  sd(s1_behave_data1$insight_mediansplit_num)*100
  
  mean(s1_behave_data1$insight_mediansplit_num_correct)*100
  sd(s1_behave_data1$insight_mediansplit_num_correct)*100
  
  mean(s1_behave_data1$riskchoice)*100
  sd(s1_behave_data1$riskchoice)*100
  
  mean(s1_behave_data1$gamble_RT)
  sd(s1_behave_data1$gamble_RT)

#### 3b) Plotting data for study 1 #######################################

  newdata$insight_mediansplit_num_correct_fac = as.factor( newdata$insight_mediansplit_num_correct)
  s1_data1 = newdata[,c('subject', 'cor2', 'insight_mediansplit_num_correct_fac', 'riskchoice', 'gamble_RT',
                        'AHA_response', 'insight_mediansplit_num', 'insight_mediansplit_num_correct')]
  s1_data1$accuracy = s1_data1$cor2
  s1_data2 = s1_data1 %>% group_by(subject, insight_mediansplit_num_correct_fac) %>% summarise_all(funs(mean), na.rm =T)   
  s1_data3=s1_data2[!is.na(s1_data2$AHA_response) & !is.na(s1_data2$insight_mediansplit_num_correct_fac),]
  

# overall amount of solved trials divided by condition: HI-L LO-I and not solved
  ##rm(list=setdiff(ls(), "s2_AHA_amount_plot1"))
  s1_behave_data = newdata %>% group_by(subject,insight_mediansplit) %>% dplyr::count(insight_mediansplit)
  s1_behave_data$n = (s1_behave_data$n/80)*100
  
  s1_p <- ggplot(s1_behave_data, aes(insight_mediansplit,n,fill = insight_mediansplit))
  s1_AHA_amount_plot <- s1_p+ geom_bar( position = 'dodge', stat = 'summary', fun.y = 'mean', color = "black") +
    geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2, size = .95) +
    geom_jitter(position = position_jitter(width = 0.1, height = 0.1), color = "#808080")+
    labs(title = "Study 1",x = "",y = "overall amount in %")+theme_classic()+scale_y_continuous(limits = c(0, 85))+
  scale_fill_manual(breaks = c('no solve','LO-I','HI-I'),values=c('#babcba', "#37c0c6", "#f98578")) 
  s1_AHA_amount_plot
  

##################################################################################################- 

## control analysis whether accuracy differs between high and low insight trials
  fitmodel_cor0<- glmmTMB(cor2 ~  insight_mediansplit_fac + blocknum + +(1|subject) + (1|stimnumber), data = newdata[newdata$insight_mediansplit_fac != "no solve",] , family = binomial(link = "logit"), na.action  = na.omit)
  summary(fitmodel_cor0)
  tab_model(fitmodel_cor0)
  
#### 4a) Main analysis: Insight effect binary  ####
  # note had to use glmer instead of glmmTMB here because otherwise ggpredict is "out of bonds" (likely a bug in ggpredict) glmer & glmmTMB results are the same though
  ######- Does solving this mooney with insight lead to different choice behavior #####  AHA_response  + (1|moonie)
  #baseline model
  fitmodel_AHA0<- glmer(riskchoice_fac ~ cor2 + blocknum + (1|subject) , data = newdata[newdata$insight_mediansplit != "no solve",],family = binomial(link = "logit"), na.action  = na.omit)
  summary(fitmodel_AHA0)
  #qqnorm(residuals(fitmodel0)); qqline(residuals(fitmodel0)) #oops? bad model fit
  
  #main model1
  fitmodel_AHA1<- glmer(riskchoice_fac ~ cor2+ insight_mediansplit_fac+blocknum +(1|subject), data = newdata[newdata$insight_mediansplit != "no solve",], family = binomial(link = "logit"), na.action  = na.omit)
  summary(fitmodel_AHA1)
  #qqnorm(residuals(fitmodel1)); qqline(residuals(fitmodel1)) #oops? bad model fit
  
  # model comparison 
  anova(fitmodel_AHA0,fitmodel_AHA1) #fitmodel1 wins
  
  # Posthoc test
  emmeans(fitmodel_AHA1, list(pairwise ~ insight_mediansplit_fac  ), adjust = "tukey")
  # plot winning model
  #ggpredict(fitmodel_AHA1, c("insight_mediansplit_fac", "cor2" )) %>% plot()+ theme_classic() #note error bars are between not within subj -> that's why they are so large
  #  for display purproses (without RT)
  AHA_ggpredict = ggpredict(fitmodel_AHA1 , c('insight_mediansplit_fac')) 
  AHA_ggpredict_plot = ggplot(AHA_ggpredict, aes(x= x, y = predicted*100 , fill=  x)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = F, color = "black") +
    geom_errorbar(aes(ymin=predicted*100-std.error*10, ymax=predicted*100+std.error*10 ), width=.2,
                  position=position_dodge(.9)) + theme_classic() +labs(title = "", #Insight-related risk propensity
                  x = "AHA! experience",y = "Risk choice in %")+scale_y_continuous(limits = c(0, 60))+
                 scale_fill_manual(breaks = c('LO-I','HI-I') ,values=c( "#37c0c6", "#f98578"))
  # output correct model + effect size
  tab_model(fitmodel_AHA1)

##################################################################################################-  

#### 4b) Main analysis: Insight effect dimensional #####
####  Does solving this mooney with insight lead to different choice behavior #####  AHA_response #+ (1|moonie)

    #baseline model
    fitmodel_AHA0_dim<- glmmTMB(riskchoice_fac ~  cor2+ scale(blocknum) + (1|subject) , data = newdata[newdata$insight_mediansplit_fac != "no solve",],family = binomial, na.action  = na.omit)
    summary(fitmodel_AHA0_dim)
    #qqnorm(residuals(fitmodel0)); qqline(residuals(fitmodel0)) #oops? bad model fit
    
    #main model1
    fitmodel_AHA1_dim<- glmmTMB(riskchoice_fac ~ cor2 +  AHA_response + scale(blocknum) +(1|subject) , data = newdata[newdata$insight_mediansplit_fac != "no solve",], family = binomial, na.action  = na.omit)
    summary(fitmodel_AHA1_dim)
    #qqnorm(residuals(fitmodel1)); qqline(residuals(fitmodel1)) #oops? bad model fit
    
    # model comparison 
    anova(fitmodel_AHA0_dim,fitmodel_AHA1_dim) #fitmodel1 wins
    
    # plot winning model
    dat = newdata[,c('subject', 'AHA_response', 'riskchoice')]
    dat1 = dat %>% group_by(subject, AHA_response) %>% summarise_all(funs(mean), na.rm =T) 
    risk_choice_median = median(dat1$riskchoice, na.rm =T)
    hist(newdata$riskchoice)
    ggplot(data = dat1, aes(x=AHA_response, y=riskchoice*100))+
      geom_point(color='gray') +
      geom_smooth(method = "lm", se = T)+  theme_classic() +labs(title = "",
          x = "AHA! experience (dimensional)",y = "Risk choice in %")+scale_y_continuous(limits = c(0, 100))
    
    ggpredict(fitmodel_AHA1_dim, c("AHA_response" )) %>% plot()+ theme_classic() #note error bars are between not within subj -> that's why they are so large
    #  for display purproses (without RT)
    AHAdim_ggpredict =ggpredict(fitmodel_AHA1_dim , c(  'AHA_response'))
    AHAdim_ggpredict_plot = ggplot(AHAdim_ggpredict, aes(x= x, y = predicted*100 )) + 
      geom_bar(stat="identity",  position=position_dodge(),show.legend = F, color = "black", fill = "light gray") +
      geom_errorbar(aes(ymin=predicted*100-std.error*20, ymax=predicted*100+std.error*20 ), width=.2,
                    position=position_dodge(.9)) + theme_classic() +labs(title = "",
                    x = "AHA! experience",y = "Risk choice in %")+scale_y_continuous(limits = c(0, 60))
                     
    # estimate effect size of predictors
    tab_model(fitmodel_AHA1_dim)
    tab_model(fitmodel_AHA1, fitmodel_AHA1_dim)

##############################################################################-


    
    
### B) STUDY 2 ##############################################################

#### 1) load data ######################

load("Study2_Behave_data.Rda")
    
#### 2a) Descriptives & view variables #####
  summary(PPSS)
  str(PPSS)
  behave_PPSS = PPSS[,c('ID', 'cor2', 'cor1', 'RT', 'RT_correct', 'sex','age',
                      'AHA_response', 'insight_mediansplit_num', 'insight_mediansplit_correct_num')]
  behave_PPSS1 = behave_PPSS %>% group_by(ID) %>% summarise_all(funs(mean), na.rm =T)   

  #demographics
  nrow(behave_PPSS1)
  s2_amount_female = sum(behave_PPSS1$sex)
  s2_amount_male = nrow(behave_PPSS1)- s2_amount_female
    
  max(behave_PPSS1$age)
  min(behave_PPSS1$age)
  
  mean(behave_PPSS1[behave_PPSS1$sex == 1,]$age)
  mean(behave_PPSS1[behave_PPSS1$sex == 0,]$age)
  
  mean(behave_PPSS1$cor1)*100
  sd(behave_PPSS1$cor1)*100
  
  mean(behave_PPSS1$cor2)*100
  sd(behave_PPSS1$cor2)*100
  
  mean(behave_PPSS1$RT, na.rm =T)
  sd(behave_PPSS1$RT, na.rm = T) 
  
  mean(behave_PPSS1$RT_correct, na.rm =T)#note for one subject RT data were not recorded for some reason, data seem legit though
  sd(behave_PPSS1$RT_correct, na.rm =T)
  
  mean(behave_PPSS1$insight_mediansplit_num)*100
  sd(behave_PPSS1$insight_mediansplit_num)*100
  
  mean(behave_PPSS1$insight_mediansplit_correct_num)*100
  sd(behave_PPSS1$insight_mediansplit_correct_num)*100
  

#control check: difference in accuracy between insight conditions
  
  M0_cor2_study2 <- glmer(cor2 ~                        (1|ID) + (1|Item),data= PPSS[PPSS$insight_mediansplit != "no solve",],family = binomial(link ="logit"), na.action  = na.omit)
  M1_cor2_study2 <- glmer(cor2 ~ insight_mediansplit  + (1|ID) + (1|Item),data= PPSS[PPSS$insight_mediansplit != "no solve",],family = binomial(link ="logit"), na.action  = na.omit)
  anova(M0_cor2_study2, M1_cor2_study2)
  tab_model(M1_cor2_study2)
  ggpredict(M1_cor2_study2 , c(  'insight_mediansplit')) #%>% plot() + ggplot2::theme_classic()


#### 2b) plot distribution of overall amount of solved trials divided by condition: HI-L LO-I and not solved ####

  s2_behave_data = PPSS %>% group_by(ID,insight_mediansplit) %>% dplyr::count(insight_mediansplit)
  s2_behave_data = na.omit(s2_behave_data)
  s2_behave_data$n = (s2_behave_data$n/120)*100
  s2_p <- ggplot(s2_behave_data, aes(insight_mediansplit,n,fill = insight_mediansplit))
  s2_AHA_amount_plot <- s2_p+ geom_bar( position = 'dodge', stat = 'summary', fun.y = 'mean', color = "black") +
    geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2, size = .95) +
    geom_jitter(position = position_jitter(width = 0.1, height = 0.1), color = "#808080")+
    labs(title = "Study 2",x = "",y = "overall amount in %")+theme_classic()+scale_y_continuous(limits = c(0, 85))+
    scale_fill_manual(breaks = c('no solve','LO-I','HI-I'),values=c('#babcba', "#37c0c6", "#f98578")) 
  s2_AHA_amount_plot


#### 3) Main analysis: ROI NACC Analysis #############################################################
 r_NACC   <- read.table("Study2_s2_ALL_HIvsLO_Insight_AHAmsplit_SubjSp_r_NAcc_18-Feb-2023.dat", sep = ",", dec = ".", header=TRUE, na.strings=c("", " " , "NA", "NAN" ))    
 l_NACC   <- read.table("Study2_s2_ALL_HIvsLO_Insight_AHAmsplit_SubjSp_l_NAcc_18-Feb-2023.dat", sep = ",", dec = ".", header=TRUE, na.strings=c("", " " , "NA", "NAN" ))     

NACC= data.frame(matrix(nrow = 64, ncol = 1))
NACC$ROI = factor(rep( c('left','right'), each =32))
NACC$ID = l_NACC$ID
NACC$vis_HIInsight_1 <- c(l_NACC$vis_HIInsight_1 ,r_NACC$vis_HIInsight_1)
NACC$vis_HIInsight_2 <- c(l_NACC$vis_HIInsight_2 , r_NACC$vis_HIInsight_2)
NACC$vis_HIInsight_3 <- c(l_NACC$vis_HIInsight_3 , r_NACC$vis_HIInsight_3)
NACC$vis_HIInsight_4 <- c(l_NACC$vis_HIInsight_4 , r_NACC$vis_HIInsight_4)

NACC$vis_LOInsight_1 <- c(l_NACC$vis_LOInsight_1 , r_NACC$vis_LOInsight_1)
NACC$vis_LOInsight_2 <- c(l_NACC$vis_LOInsight_2 , r_NACC$vis_LOInsight_2)
NACC$vis_LOInsight_3 <- c(l_NACC$vis_LOInsight_3 , r_NACC$vis_LOInsight_3)
NACC$vis_LOInsight_4 <- c(l_NACC$vis_LOInsight_4 , r_NACC$vis_LOInsight_4)

data_HILOI <- rbind( NACC[,-1] ) 
#data_HILOI$roi=  rep(c( 'NACC'  ), each=)# 
library(tidyr)
data_HILOI_long <- data_HILOI %>% pivot_longer(cols = 'vis_HIInsight_1':'vis_LOInsight_4', names_to='condition', values_to="betas")
data_HILOI_long$session = rep(c(1,2,3,4), 2*64*1) #'3*conds*32subj*7ROIs'
#data_HILOI_long$condition <- gsub('vis_', '', as.character(data_HILOI_long$condition)) 
data_HILOI_long$condition <- gsub('_1', '', as.character(data_HILOI_long$condition)) 
data_HILOI_long$condition <- gsub('_2', '', as.character(data_HILOI_long$condition)) 
data_HILOI_long$condition <- gsub('_3', '', as.character(data_HILOI_long$condition)) 
data_HILOI_long$condition <- gsub('_4', '', as.character(data_HILOI_long$condition)) 
data_HILOI_long$condition = factor(data_HILOI_long$condition, 
                                   levels = c('vis_LOInsight','vis_HIInsight'), 
                                   labels=c('LO-I','HI-I'))
data_HILOI_long$sex = NA
data_HILOI_long$cor2 = NA

### including gender & accuracy variable (between subject)
PPSS1 = PPSS[c('ID', 'cor2', 'scan','insight_mediansplit')]
PPSS_agg <- PPSS1 %>% group_by(ID, insight_mediansplit,scan ) %>% summarise_all(funs(mean), na.rm =T)  
PPSS_agg1 = PPSS_agg[(PPSS_agg$insight_mediansplit != "no solve"),]
PPSS_agg1$ID <- gsub('PVI', 'sub-', as.character(PPSS_agg1$ID)) 

  #PPSS_agg1$beta = NA
  for(i in 1:nrow(PPSS_agg1) ){
    for(j in 1:nrow(data_HILOI_long)){
    if(is.na(PPSS_agg1$ID[i])){
    } else if((data_HILOI_long$ID[j] == PPSS_agg1$ID[i]) & (data_HILOI_long$session[j] == PPSS_agg1$scan[i]) &
         (data_HILOI_long$condition[j] == PPSS_agg1$insight_mediansplit[i]) ){
          #data_HILOI_long$sex[j] = PPSS_agg1$sex[i] 
          data_HILOI_long$cor2[j] = PPSS_agg1$cor2[i] 
      }}}

  
### mixed models for NAcc 
library(glmmTMB)
  HIvLOI_m0_NAcc <- glmmTMB(betas ~ session +  cor2 + ROI + (1|ID) ,data= data_HILOI_long, na.action  = na.omit) #
  HIvLOI_m1_NAcc <- glmmTMB(betas ~ condition   +session  +  cor2 +ROI + (1|ID) ,data= data_HILOI_long, na.action  = na.omit) #

  anova(HIvLOI_m0_NAcc, HIvLOI_m1_NAcc)
  tab_model(HIvLOI_m1_NAcc,show.std  = T)
  emmeans(HIvLOI_m1_NAcc, list(pairwise ~ condition), adjust = "none") 
  
  HIvLOI_NAcc_ggplot = ggpredict(HIvLOI_m1_NAcc , c('condition')) 
  
  HIvLOI_NAcc_finalplot= ggplot(HIvLOI_NAcc_ggplot, aes(x= x, y = predicted ,fill= x)) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T, color = "black") +
    geom_errorbar(aes(ymin=predicted-std.error , ymax=predicted+std.error ), width=.2,
                  position=position_dodge(.9)) + theme_classic() +labs(title = "",
                  x = "AHA! experience", y = "Beta estimate in NAcc")+theme(legend.position="none")+
  scale_fill_manual(breaks = c('LO-I','HI-I') ,values=c( "#37c0c6", "#f98578")) #, fill = ""
  


#### 4) do single trial analysis in subject space   #####

  load(file = "Study2_singletrial_fMRI_data.Rda") 

  AHAsum_m0_NAcc <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ cor2 +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  AHAsum_m1_NAcc <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ cor2+scale(insight_sum) +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #

  anova(AHAsum_m0_NAcc, AHAsum_m1_NAcc)
  tab_model(AHAsum_m1_NAcc,show.std  = T)

  HIvLOI_NAcc_ggplot_dim = ggpredict(AHAsum_m1_NAcc , c(  'insight_sum'))# %>%plot()+ ggplot2::theme_classic()
  
  HIvLOI_NAcc_finalplot_dim= ggplot(HIvLOI_NAcc_ggplot_dim, aes(x= x, y = predicted )) + 
    geom_bar(stat="identity",  position=position_dodge(),show.legend = T, color = "black",fill = "light gray") +
    geom_errorbar(aes(ymin=predicted-std.error , ymax=predicted+std.error ), width=.2,
                  position=position_dodge(.9)) + theme_classic() +labs(title = "",
                  x = "AHA! experience", y = "Single trial: Beta estimate in NAcc", fill = "")

  ### split by different AHA dimensions
  AHA_m0_NAccdim <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ scan+cor2 +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  AHA_m1_NAccdim <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ scan+cor2+scale(Certain) +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  AHA_m2_NAccdim <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ scan+cor2+scale(Certain) +scale(Sudden) +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  AHA_m3_NAccdim <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ scan+cor2+scale(Certain) +scale(Sudden) +scale(Aha) +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  AHA_m4_NAccdim <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ scan+cor2+scale(Certain)*scale(Aha) +scale(Sudden)  +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  AHA_m5_NAccdim <- glmmTMB(Beta_solve_NAcc ~ ROI+trial+ scan+cor2+scale(Certain)*scale(Aha) *scale(Sudden)  +(1|ID) + (1|Item) ,data= stPPSS1, na.action  = na.omit) #
  
  anova(AHA_m0_NAccdim, AHA_m1_NAccdim, AHA_m2_NAccdim, AHA_m3_NAccdim,AHA_m4_NAccdim,AHA_m5_NAccdim)
  tab_model( AHA_m3_NAccdim,show.std  = T)
  tab_model( AHA_m4_NAccdim,show.std  = T)
  emmeans(AHA_m4_NAccdim, list(pairwise ~ Aha | Certain), adjust = "tukey") 
  
  HIvLOI_NAcc_ggplot_2 = ggpredict(AHA_m4_NAccdim , c(  'Aha','Certain')) %>%plot()+ ggplot2::theme_classic()+labs(title = "",
                                                 x = "AHA: Positive Emotion", y = "Single trial: Beta estimate in NAcc", 
                                                 fill = "", colour = "AHA: Certainty")
  
  
############################################################################ - 

  
  
### C) Control analysis: combine study 1 und study 2 (alldata)for comparison ###############################
  
  #newdata$cor2 = newdata$cor2_num-1
  newdata$gender =NA
  newdata$gender[newdata$sex == "F"]=1
  newdata$gender[newdata$sex == "M"]=0
  newdata$gender[newdata$sex == "D"]=0
  newdata$AHA_response_scale = scale(newdata$AHA_response)
  newdata_neu = newdata[c("subject",  "moonie", "cor1", "cor2_num", "RT", "RT_correct", "AHA_response", "gender", "age", 
                          "insight_mediansplit_num", "insight_mediansplit_num_correct", "study", "AHA_response_scale")]
    colnames(newdata_neu)[1] <- 'ID'
    colnames(newdata_neu)[2] <- 'items'
    colnames(newdata_neu)[4] <- 'cor2'
    colnames(newdata_neu)[8] <- 'sex'
    colnames(newdata_neu)[11] <- 'insight_mediansplit_correct_num'
    newdata_neu$RT = newdata_neu$RT/1000
    newdata_neu$RT_correct = newdata_neu$RT_correct/1000
    
  PPSS$AHA_response_scale = scale(PPSS$AHA_response)
  PPSS_neu = PPSS[c('ID', "items", 'cor1', 'cor2', 'RT', 'RT_correct', 'AHA_response', "sex", "age",
                    "insight_mediansplit_num", "insight_mediansplit_correct_num","study", 'AHA_response_scale')]
      PPSS_neu$ID = as.factor(PPSS_neu$ID)
  
  alldata = rbind(newdata_neu, PPSS_neu)
  
  agg_data_all <- alldata %>% group_by(study) %>%  dplyr::summarise(across(c('cor2','cor1','RT', 'sex', 'RT_correct', 'insight_mediansplit_num',
                              'AHA_response_scale', 'age', 'AHA_response', 'insight_mediansplit_correct_num'), ~mean(., na.rm=T))) 
 
  
  ###### DOING THE STATISTICAL COMPARISON between study1 and 2
  library(sjstats)
  library(sjmisc)
  library(mediation)
#### COR1
  #baseline model
  basemodel_study12_cor1<- glmmTMB(cor1 ~  (1|ID) + (1|items), data = alldata,family = binomial, na.action  = na.omit)
  #main model1
  fitmodel_study12_cor1<- glmmTMB(cor1 ~  study  + (1|ID) + (1|items), data = alldata,family = binomial, na.action  = na.omit)
  # model comparison 
  anova(basemodel_study12_cor1,fitmodel_study12_cor1) #fitmodel1 wins

  
#### COR2
  #baseline model
  basemodel_study12_cor2<- glmmTMB(cor2 ~  (1|ID) + (1|items), data = alldata,family = binomial, na.action  = na.omit)
  #main model1
  fitmodel_study12_cor2<- glmmTMB(cor2 ~  study  + (1|ID) + (1|items), data = alldata,family = binomial, na.action  = na.omit)
  summary(fitmodel_study12_cor2)
  # model comparison 
  anova(basemodel_study12_cor2,fitmodel_study12_cor2) #fitmodel1 wins

#### Sex
  #baseline model
  basemodel_study12_sex<- glmmTMB(sex ~  (1|ID) + (1|items), data = alldata,family = binomial, na.action  = na.omit)
  #main model1
  fitmodel_study12_sex<- glmmTMB(sex ~  study  + (1|ID) + (1|items), data = alldata,family = binomial, na.action  = na.omit)
  summary(fitmodel_study12_sex)
  # model comparison 
  anova(basemodel_study12_sex,fitmodel_study12_sex) #fitmodel1 wins

#### RT
  #baseline model
  basemodel_study12_RT<- glmmTMB(log(RT+10) ~  (1|ID) + (1|items), data = alldata, na.action  = na.omit)
  summary(basemodel_study12_RT)
  qqnorm(residuals(basemodel_study12_RT)); qqline(residuals(basemodel_study12_RT)) #oops? bad model fit
  #main model1
  fitmodel_study12_RT<- glmmTMB(log(RT+10) ~  study  + (1|ID) + (1|items), data = alldata, na.action  = na.omit)
  summary(fitmodel_study12_RT)
  qqnorm(residuals(fitmodel_study12_RT)); qqline(residuals(fitmodel_study12_RT)) #oops? bad model fit
  # model comparison 
  anova(basemodel_study12_RT,fitmodel_study12_RT) #fitmodel1 wins  

#### RT correct
  #baseline model
  basemodel_study12_RT_cor<- glmmTMB(log(RT_correct+10) ~  (1|ID) + (1|items), data = alldata, na.action  = na.omit)
  summary(basemodel_study12_RT_cor)
  qqnorm(residuals(basemodel_study12_RT_cor)); qqline(residuals(basemodel_study12_RT_cor)) #oops? bad model fit
  #main model1
  fitmodel_study12_RT_cor<- glmmTMB(log(RT_correct+10) ~  study  + (1|ID) + (1|items), data = alldata, na.action  = na.omit)
  summary(fitmodel_study12_RT_cor)
  qqnorm(residuals(fitmodel_study12_RT_cor)); qqline(residuals(fitmodel_study12_RT_cor)) #oops? bad model fit
  # model comparison 
  anova(basemodel_study12_RT_cor,fitmodel_study12_RT_cor) #fitmodel1 wins    
  
#### AHA! mediansplit
  #baseline model
  basemodel_study12_AHA<- glmmTMB(insight_mediansplit_num ~  (1|ID) + (1|items), data = alldata, family = binomial,na.action  = na.omit)
  summary(basemodel_study12_AHA)
  #main model1
  fitmodel_study12_AHA<- glmmTMB(insight_mediansplit_num ~  study  + (1|ID) + (1|items), data = alldata, family = binomial,na.action  = na.omit)
  summary(fitmodel_study12_AHA)
  # model comparison 
  anova(basemodel_study12_AHA,fitmodel_study12_AHA) #fitmodel1 wins 
  
#### AHA! mediansplit (for correctly solved items)
  #baseline model
  basemodel_study12_AHA_cor<- glmmTMB(insight_mediansplit_correct_num ~  (1|ID) + (1|items), data = alldata, family = binomial,na.action  = na.omit)
  #main model1
  fitmodel_study12_AHA_cor<- glmmTMB(insight_mediansplit_correct_num ~  study  + (1|ID) + (1|items), data = alldata, family = binomial,na.action  = na.omit)
  # model comparison 
  anova(basemodel_study12_AHA_cor,fitmodel_study12_AHA_cor) #fitmodel1 wins

#### age
  #baseline model
  basemodel_study12_age<- glmmTMB(age ~  (1|ID) + (1|items), data = alldata, na.action  = na.omit)
  summary(basemodel_study12_age)
  qqnorm(residuals(basemodel_study12_age)); qqline(residuals(basemodel_study12_age)) #oops? bad model fit
  #main model1
  fitmodel_study12_age<- glmmTMB(age ~  study  + (1|ID) + (1|items), data = alldata, na.action  = na.omit)
  summary(fitmodel_study12_age)
  qqnorm(residuals(fitmodel_study12_age)); qqline(residuals(fitmodel_study12_age)) #oops? bad model fit
  # model comparison 
  anova(basemodel_study12_age,fitmodel_study12_age) #fitmodel1 wins    

  
  
### D) Plot figures for study 1 & 2 ##############  
  alldata1 <- alldata[,-2]
  agg_data_all1 <- alldata1 %>% group_by(study,ID) %>% dplyr::summarise_all(funs(mean), na.rm =T)   
  
  plot_multi_histogram_neu <- function(df, feature, label_column) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
      geom_density(alpha=0.6) + theme_classic()+
      geom_histogram(alpha=0.2, position="identity", aes(y = ..density..), color="gray") +
      labs(x=feature, y = "Density") +scale_fill_manual( values = c("#f7f34b","#b46f8f"))
    plt + guides(fill=guide_legend(title=label_column)) } 
  
 s12_RT<- plot_multi_histogram_neu(agg_data_all1[!is.na(agg_data_all1$RT),], 'RT', 'study')+labs(x = "Solution time in sec",y = "Density")+
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study1",]$RT, na.rm =T)), color="black", linetype="longdash", size=.8) +
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study2",]$RT, na.rm =T)), color="black", linetype="solid", size=.8)
           
 s12_RT_cor<- plot_multi_histogram_neu(agg_data_all1[!is.na(agg_data_all1$RT_correct),], 'RT_correct', 'study')+labs(x = "Solution time (correct) in sec",y = "Density")+
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study1",]$RT_correct, na.rm =T)), color="black", linetype="longdash", size=.8) +
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study2",]$RT_correct, na.rm =T)), color="black", linetype="solid", size=.8)
  
 s12_AHA <- plot_multi_histogram_neu(agg_data_all1[!is.na(agg_data_all1$AHA_response_scale),], 'AHA_response_scale', 'study')+labs(x = "AHA! experience (scaled)",y = "Density")+
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study1",]$AHA_response_scale, na.rm =T)), color="black", linetype="longdash", size=.8) +
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study2",]$AHA_response_scale, na.rm =T)), color="black", linetype="solid", size=.8)
  
 s12_cor2<- plot_multi_histogram_neu(agg_data_all1[!is.na(agg_data_all1$cor2),], 'cor2', 'study')+labs(x = "Accuracy",y = "Density")+
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study1",]$cor2, na.rm =T)), color="black", linetype="longdash", size=.8) +
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study2",]$cor2, na.rm =T)), color="black", linetype="solid", size=.8)
 
 s12_cor1<- plot_multi_histogram_neu(agg_data_all1[!is.na(agg_data_all1$cor1),], 'cor1', 'study')+labs(x = "Solution button",y = "Density")+
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study1",]$cor1, na.rm =T)), color="black", linetype="longdash", size=.8) +
    geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study2",]$cor1, na.rm =T)), color="black", linetype="solid", size=.8)
  
 s12_aha_ms<- plot_multi_histogram_neu(agg_data_all1[!is.na(agg_data_all1$insight_mediansplit_correct_num),], 'insight_mediansplit_correct_num', 'study')+labs(x = "AHA! experience (median split)",y = "Density")+
   geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study1",]$insight_mediansplit_correct_num, na.rm =T)), color="black", linetype="longdash", size=.8) +
   geom_vline(aes(xintercept=mean(agg_data_all1[agg_data_all1$study == "study2",]$insight_mediansplit_correct_num, na.rm =T)), color="black", linetype="solid", size=.8)
 
  library(ggpubr)

##### plot everything together   
 
 Figure2 <- ggarrange(AHA_ggpredict_plot, AHAdim_ggpredict_plot,HIvLOI_NAcc_finalplot, HIvLOI_NAcc_finalplot_dim, HIvLOI_NAcc_ggplot_2,# Fig1x1, 
                      ncol = 4, nrow = 2,  common.legend = T, legend = "bottom", labels = c('A','B','C','D', 'E'))
 
 Figure3 <- ggarrange( s12_cor1,  s12_RT, s12_AHA,s1_AHA_amount_plot,
                                  s12_cor2,s12_RT_cor, s12_aha_ms,s2_AHA_amount_plot,
                                 ncol = 4, nrow = 2, legend = "top",common.legend = TRUE, 
                                 labels = c('A','B','C','D', 'E','F','G', 'H'))

  