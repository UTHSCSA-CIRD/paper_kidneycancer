#' ---
#' title: "Kidney Cancer Disparities"
#' author: "Bokov, Michalek, Guerra, Shah, Kaushig, Rodriguez"
#' date: "02/19/2018"
#' ---
#'
#+ echo=FALSE, message=FALSE
knitr::opts_chunk$set(echo=F,warning = F,message=F);
#+ cache=FALSE
source('global.R');
#' Report date: `r date()`.
#'
# Revision: `r gitstamp()`.
#'
#' Data file: `r inputdata`.
#' 
#' Metadata file: `r inputmeta`
#' 
#' 
#+ source_run, cache=TRUE
source('run.R');
#' 
## Let's try printing out the tables of concordances and goodness-of-fit
# #' 
# #+ results='asis'
# if(!interactive()){
#   #cat('\nFrailty\n');
#   #stargazer(results_con_wald_frail,type = 'html');
#   cat('Cluster\n');
#   stargazer(results_con_wald_cluster,type = 'html');
#   cat('Interval\n');
#   stargazer(results_con_wald_t2,type = 'html');
#   cat('Numeric\n');
#   stargazer(results_con_wald_numeric,results_con_wald_cluster,type = 'html');
#   cat('Labs, Categoric and Numeric predictors side by side\n');
#   stargazer(data.frame(
#     Numeric=results_con_wald_numeric[rows_labs,]
#     ,Categoric=results_con_wald_t2[rows_labs,]),type = 'html',summary = F);
# };
# 
# Hispanic ethnicity, as  predictor in an interval-censored model
#' 
#+ startplots, results='asis'
#sapply(class_hisp_exact
#       ,function(xx) stargazer(cox_t2_demog[[xx]]
#                               ,type=if(interactive()) 'text' else 'html')) -> .junk;
#'
#' ---
#'
#autoplot(survfit(Surv(a_dxage,a_cens_1)~pred_hisp,d5),col=c('red','blue')) + 
#+ results='asis'
autoplot(update(sf0,.~pred_hisp)
         ,main='Ethntableone::CreateTableOne(data=setNames(mutate(d5[,class_mainvars_exact]
                                              ,a_age_at_stdx=a_age_at_stdx/365.25)
         ,mainvars_nicelabels)
         ,vars=setdiff(mainvars_nicelabels,hispanic_nicelabel)
         ,strata=hispanic_nicelabel) %>% 
         print(printToggle=F) %>% data.frame %>% 
         setNames(c('Non Hispanic','Hispanic','p-value','')) %>% 
         knitr::kable(format='markdown') %>% gsub('NA','-',.)icity as Risk Factor for Progression') +
  scale_color_discrete('Ethnicity\n(N=1162)'
                       ,labels=c('Non Hispanic','Hispanic'));
tidy(cox_t2_demog[[class_hisp_exact[1]]])[,2:5] %>% t %>% data.frame %>% 
  knitr::kable(format='markdown');

# We will go for the following:
feb2018_pres<- c("Platelet # Bld Auto (777-3)"
                 , "RDW RBC Auto-Rto (788-0)"
                 , "Body Mass Index");
feb2018_prescodes<-c('v038_Pltlt_At_GENERIC_KUH_COMPONENT_ID_5341_numnona'
                     ,'v050_RDW_RBC_At_Rt_GENERIC_KUH_COMPONENT_ID_5629_numnona'
                     ,'v002_Bd_Ms_Indx_numnona');

for(ii in seq_along(feb2018_pres)){
  iiname <- feb2018_pres[ii];
  iicode <- feb2018_prescodes[ii];
  print(plots_cph_numeric[[iiname]] + ggtitle(iiname) +
    scale_color_discrete('Lab Value\n(N=1162)',labels=c('Low','High')));
  tidy(cox_ph_models_numeric[[iicode]])[,2:5] %>% t %>% 
    knitr::kable(format='markdown') %>% print;
  cat('\n\n');
}
#multiplot(plotlist=plots_cph_numeric[feb2018_pres],cols=1);
#' 
#+ demog_table
class_mainvars_exact <- c('sex_cd','a_cens_1','a_age_at_stdx',class_hisp_exact
                          ,grep('_Pltlt_At_GENERIC_KUH_COMPONENT_ID_5341_num$|_RDW_RBC_At_Rt_GENERIC_KUH_COMPONENT_ID_5629_num$|_Bd_Ms_Indx_num$',names(d5),val=T));
                          
mainvars_nicelabels <- submulti(class_mainvars_exact
                                ,rbind(m0[,1:2]
                                       ,c('sex_cd','Sex')
                                       ,c('a_age_at_stdx','Diagnosis Age')
                                       ,c('a_cens_1','Metastasis')));
hispanic_nicelabel <- submulti(class_hisp_exact,m0[,1:2]);

#' ### Demographic Summary
#' 
#' These preliminary results suggest that Hispanic patients tend to 
#' be diagnosed at a _younger_ age and yet still progress to 
#' metastasis at a higher rate.
#' 
tableone::CreateTableOne(data=setNames(mutate(d5[,class_mainvars_exact]
                                              ,a_age_at_stdx=a_age_at_stdx/365.25)
                                       ,mainvars_nicelabels)
                         ,vars=setdiff(mainvars_nicelabels,hispanic_nicelabel)
                         ,strata=hispanic_nicelabel) %>% 
  print(printToggle=F) %>% data.frame %>% 
  setNames(c('Non Hispanic','Hispanic','p-value','')) %>% 
  knitr::kable(format='markdown') %>% gsub('NA','-',.);
#New table for coxph_mv1
#stargazer(coxph_mv1, covariate.labels = vars_mv1, 
#  dep.var.labels = 'Time in days until metastasis',star.cutoffs=c(0.05,0.01,1e-6),type='text');
#tidy(coxph_mv1)[,1:5] %>% knitr::kable(format = 'markdown');

#autoplot(update(sf0,.~pred_mv1>median(pred_mv1)));
