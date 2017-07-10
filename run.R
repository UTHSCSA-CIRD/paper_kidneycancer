#' ---
#' title: "Falls Among Older Women in Primary Care"
#' author: "Alex F. Bokov"
#' date: "05/09/2017"
#' output: html_document
#' ---
#' 
#' ## Load libraries
require(survival); library(magrittr); require("dplyr");
require('data.table'); require("ggplot2"); require('MASS'); 
require('Hmisc'); require('readr');
#' ...and local functions...
source('functions.R');
#' ## Set variables, 1st pass
session <- 'session.rdata';
inputdata <- 'falls.csv'; datadict <- 'df_dynsql.csv';
rebuild <- c();
if(file.exists(session)) load(session);
if(!exists('dd')) {
  rebuild <- c(rebuild,'dd');
  dd <- read_csv(datadict,na='');
}
cols2drop <- c('start_date','birth_date','sex_cd','v063_VTMN_D_info');
#' We are going to delete the inactive diagnoses for now, under the following 
#' reasoning: if the diagnosis is historical, it isn't a dynamic reflection of
#' the patient's state. If it is deleted or resolved, then it won't persist in
#' future visits anyway. Leaving the fall codes in for now, though, pending a 
#' closer look....
cols2drop <- c(cols2drop,
               subset(dd,rule='diag')$colname %>% 
                 grep('_inactive',.,val=T) %>% 
                 grep('_trpng_stmblng_|_ACCDNTL_FLLS_',.,inv=T,val=T));
demcols <- c('patient_num','race_cd','language_cd','age_at_visit_days')
#' ## Load data
if(session %in% list.files()) load(session);
if(!exists('d0')) {
  rebuild <- c(rebuild,'d0');
  d0 <- read_csv(inputdata,na='');
}
#' Is the data already arranged in order of increasing patient_num and age? The
#' when uncommented, following code should return TRUE if it is (and it does)
# all.equal(d0,arrange(d0,patient_num,age_at_visit_days))
#' ## Set variables, 2nd pass
demcols <- c(demcols,subset(dd,rule=='ethnicity')$colname);
syn_diag_active <- list(
  c('v065_trpng_stmblng','v041_ACCDNTL_FLLS'),
  c('v048_rsprtr_unspcfd','v014_rsprtr_unspcfd'),
  c('v049_Act_brnchts','v015_brnchts_brnchlts','v016_Act_brnchts'),
  c('v030_invlvng_rsprtr','v031_Cgh'),
  c('v050_Vsmtr_rhnts','v017_rhnts_unspcfd'),
  c('v039_lprtn_lpdms','v002_Dsrdrs_mtblsm'),
  c('v051_Gstr_esphgl','v018_Dss_esphgs','v019_dsrdrs_esphgs'),
  c('v054_dsrdrs_urnr','v022_infctn_spcfd'),
  c('v059_R_Plr','v033_Smptms_invlvng'),
  c('v052_Dvrtclr_intstn','v020_Dvrtcl_cln'),
  c('v057_invlvng_dgstv','v032_invlvng_dgstv'),
  c('v053_K_Cnstptn','v021_Cnstptn'),
  c('v040_elctrlt_acd_bs','v003_elctrlt_acd_bs'),
  c('v043_Unspcfd_dmnt','v005_Dmnt_unspcfd'),
  c('v046_Alzhmr_s_ds','v009_Alzhmr_s_ds','v008_crbrl_dgnrtns'),
  c('v062_Smptms_cncrng','v026_mtblsm_dvlpmnt','v027_Abnrml_undrwght'),
  c('v047_G_Slp_dsrdrs','v024_Slp_dstrbncs'),
  c('v044_dprsv_dsrdr','v007_Dprsv_clsfd'),
  c('v045_anxt_dsrdrs','v006_dsctv_smtfrm'),
  c('v038_Vtmn_dfcnc','v001_Vtmn_dfcnc'),
  c('v037_Dfcnc_vtmns','v000_Ddim(subset(d2,ncode>0|v055_Ofc_Vst))fcnc_cmpnts'),
  c('v061_Mls_and_ftg','v025_Mls_and_ftg')
);
#' ## Clean up data
if('d0' %in% rebuild) {
  # Remove impossible START_DATEs (2 of them)
  d0 <- subset(d0,age_at_visit_days > 0);
  # Drop non-informative columns
  dd$present <- T;
  d0[,cols2drop] <- NULL;
  dd[dd$colname%in%cols2drop,'present'] <- F;
  # English-speaker or not
  d0$language_cd <- cl_bintail(d0$language_cd,1);
  # Hispanic, Non-Hispanic, or other
  d0$v042_Ethnct <- cl_bintail(d0$v042_Ethnct,2);
  # Race: white, black, asian, other
  d0$race_cd <- cl_bintail(d0$race_cd,4,2);
  # replace all diagnoses with T/F values
  # DON'T USE: code = visit type
  # I think visit type has values of '' instead of NA, so always true
  # UNKNOWN_DATA_ELEMENT = medications
  cols2tf <- subset(
    dd,rule%in%c('diag','UNKNOWN_DATA_ELEMENT','code','codemod')&present)$colname;
  d0[,cols2tf]<-sapply(d0[,cols2tf],function(xx) !is.na(xx));
  # Flag aberrant values
  ifelse(grepl('\'TNP\'',d0$v064_VTMN_TTL_1990_1_info)
         ,NA,d0$v064_VTMN_TTL_1990_1_info) %>% 
    ifelse(grepl('\'L\'',.),'L',.) %>% 
    ifelse(grepl('\'H\'',.),'H',.) %>% 
    ifelse(is.na(.),ifelse(is.na(d0$v064_VTMN_TTL_1990_1_num),NA,'OK'),.) -> 
    d0$v064_VTMN_TTL_1990_1_info;
  for(ii in syn_diag_active) {
    d0[,ii[1]]<-apply(d0[,ii],1,any);
    d0[,ii[-1]]<-NULL;
    dd[dd$colname%in%ii[-1],'present'] <- F;
    };
  # create the event indicators
  group_by(d0,patient_num) %>% 
    mutate(tt=age_at_visit_days
           ,which_event=cumsum(v065_trpng_stmblng)
           ,cens=lead(which_event)
           ,first=c(1,rep_len(0,length(tt)-1))
           ,last=c(rep_len(0,length(tt)-1),1)
           # preserving for each individual their age at first visit
           ,agestart=min(age_at_visit_days)) -> d1; # 100072 x 56
  # d1 = dataset with cumulative counts and censoring indicators for all falls
  # d2 = dataset with observations _prior_ to the first event only (or where no
  # events have occurred)
  d2 <- subset(d1,which_event==0);  # 87110 x 56
  subset(summarise(d2
                   ,v065in=min(which(v065_trpng_stmblng_inactive))
                   ,v041in=min(which(v041_ACCDNTL_FLLS_inactive))
                   ,v065=min(which(v065_trpng_stmblng)))
         ,v065>v065in|v065>v041in|v065==1)$patient_num -> inactive_first;
  # The following command shows the first occurrences o TRUE in the two inactive 
  # falls columns and the 'real' falls column (v065). If the first inactive dates 
  # are earlier than the first real event, those are IDs to toss because they are
  # patients who may have come in after they alread fell for the first time. Also
  # removing patients whose very first visit is in connection with a fall
  d2 <- subset(d2,!patient_num%in%inactive_first);
  # to prevent individuals who only had one visit and it was not in connection
  # with a fracture from being treated as missing data
  # d2$cens <- ifelse(is.na(d2$cens),0,d2$cens)
  # instead, though, we'll remove them
  d2 <- subset(d2,!is.na(cens));
  rebuild<-c();
  save.image(session);
}
#' ## Exploration of possibly combinable variables
#' The active diagnoses
#' ## Univariate models
predictors <- grep('_inactive',names(d2)[sapply(d2,class)=='logical'],inv=T,val=T);
predictors <- grep('_Ofc_Vst$|_Apntmnt$|_Prcdr$|_trpng_stmblng$',predictors,inv=T,val=T);
d2$npred <- rowSums(d2[,predictors]);
d3 <- subset(d2,npred>0|v055_Ofc_Vst|v056_Prcdr|cens>0);
#' The `patient_num`s of the patients who had an event, for later use in validation
fallpats <- subset(d3,cens==1)$patient_num;
varclus(as.matrix(d3[,predictors])+0,similarity='bothpos') -> vc1;
plot(vc1$hclust);
cxm <- coxph(formula = Surv(tt, cens) ~ 1, data = d2);
wbm <- survreg(formula = Surv(tt, cens) ~ 1, data = d2);
paste0("update(cxm,.~",predictors,")") %>% 
  sapply(function(xx) parse(text=xx)) %>% 
  sapply(eval,simplify=F) -> unicox;
unicox3 <- sapply(unicox,update,data=d3,simplify=F);
paste0("update(wbm,.~",predictors,")") %>% 
  sapply(function(xx) parse(text=xx)) %>% 
  sapply(eval,simplify=F) -> uniwei;
uniwei3 <- sapply(uniwei,update,data=d3,simplify=F);
#' ## How good are they? 
c_unicox <- sapply(unicox,function(xx) summary(xx)$concordance[1]);
c_unicox3 <- sapply(unicox3,function(xx) summary(xx)$concordance[1]);
c_uniwei <- 1-sapply(uniwei,function(xx) survConcordance(Surv(tt,cens)~predict(xx),d2)$concord);
c_uniwei3 <- 1-sapply(uniwei3,function(xx) survConcordance(Surv(tt,cens)~predict(xx),d3)$concord);
names(unicox) <- names(uniwei) <- names(c_unicox) <- names(c_uniwei) <- names(c_unicox3) <- names(c_uniwei3) <- predictors;
#' There seems to be an inverse relationship between concordances for Weibull and Cox
#' I don't know why that is.
plot(c_unicox,c_uniwei); point(c_unicox3,c_uniwei3);
#' Overall, cox (red) seems to perform better.
plot(sort(c_unicox),type='s',col='red',ylim=range(c(c_unicox,c_uniwei,c_unicox3,c_uniwei3)),ylab='Concordance');
lines(sort(c_uniwei),type='s');
lines(sort(c_unicox3),type='s',lty=2,col='red');
lines(sort(c_uniwei3),type='s',lty=2);
#' ## Our starting model (for stepwise selection)
# sort(c_unicox) %>%  names %>% paste(collapse='+') %>% 
#   paste('coxph(Surv(tt,cens)~',.,'+agestart,data=d2)') %>% parse(text=.) %>% 
#   eval -> stepstart;
#' ## Interactions
#' There's no point in considering an interaction of binary variables if if its 
#' component terms hardly ever co-occur. So we need to narrow down the list of
#' candidate participants in such interactions. We step through each predictor
#' and count the total number of other predictors that co-occur with it.
ints<-list();
for(kk in predictors) 
  ints[[kk]]<-sum(d2[unlist(d2[,kk]),setdiff(predictors,kk)]+0);
ints<-sort(unlist(ints));
plot(ints,type='s');
intpredictors <- names(tail(ints,5));
# paste(intpredictors,collapse='+') %>% 
#   paste0('update.formula(stepstart$call$formula,.~.:(',.,'))') %>% 
#   formula -> fm_upper;
# paste(predictors,collapse='+') %>% 
#   paste0('stepAIC(update(cxm,.~agestart,data=d3),scope=list(lower=~1,upper=~.+',.,'),
#          direction="both")') %>% 
#   parse(text=.) %>% eval -> coxaic1;
# paste(intpredictors,collapse = '+') %>% 
#   paste0('stepAIC(update(coxaic1,.~.-agestart),
#          scope=list(lower=~1,upper=~(.+',.,')^2+agestart),direction="both")') %>% 
#   parse(text=.) %>% eval -> coxaic2;
#' Nevermind the above. We can do it in one step, if we run it on a server...
paste(predictors,collapse = '+') %>% 
  paste0('stepAIC(update(cxm,data=d3)
         ,scope=list(lower=.~1,upper=.~(agestart+',.,')^2),direction="both")') %>%
  parse(text=.) %>% eval -> coxaic3; save.image(session);
#' And here is the non-age version...
paste(predictors,collapse = '+') %>% 
  paste0('stepAIC(update(cxm,data=d3)
         ,scope=list(lower=.~1,upper=.~(',.,')^2),direction="both")') %>%
  parse(text=.) %>% eval -> coxaic4; save.image(session);


d3$coxaic2 <- predict(coxaic2,type = 'lp');
