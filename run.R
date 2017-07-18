#' ---
#' title: "Kidney Cancer"
#' author: "Guerra, Bokov, Wilson"
#' date: "07/10/2017"
#' ---
#' 
#' ## Load libraries
require(survival); require(dtplyr); require(magrittr); require("ggplot2"); require('MASS'); 
require('Hmisc'); require('readr');
#' ## Load local config file
source('./config.R');
source('./metadata.R');
#'
#' ...and local functions...
source('./functions.R');
#'
#' ## Set variables, 1st pass
session <- 'session.rdata';
#' The rebuild vector is going to c
rebuild <- c();
#cols2drop <- c('start_date','birth_date','sex_cd','v063_VTMN_D_info');
#' We are going to delete the inactive diagnoses for now, under the following 
#' reasoning: if the diagnosis is historical, it isn't a dynamic reflection of
#' the patient's state. If it is deleted or resolved, then it won't persist in
#' future visits anyway. Leaving the fall codes in for now, though, pending a 
#' closer look....
# cols2drop <- c(cols2drop,
#                subset(dd,rule='diag')$colname %>% 
#                  grep('_inactive',.,val=T) %>% 
#                  grep('_trpng_stmblng_|_ACCDNTL_FLLS_',.,inv=T,val=T));
#' ## Load data
if(session %in% list.files()) load(session);
#' Is the data already arranged in order of increasing patient_num and age? The
#' when uncommented, following code should return TRUE if it is (and it does)
# all.equal(d0,arrange(d0,patient_num,age_at_visit_days))
# Current return-value is TRUE
#' ## Read the input file
d0 <- read_csv(inputdata,na='');
#' ## Clean up data
if('d0' %in% rebuild) {
  # Remove impossible START_DATEs
  d0 <- subset(d0,age_at_visit_days > 0);
}
#' ## Convert columns
#' 
#' Create copy of original dataset
d1 <- d0;
#' Obtain the actual column names for the Yes/No columns in this dataset
class_yesno_tailgreps %>% paste0(collapse='|') %>% 
  grep(names(d0),val=T) -> class_yesno_exact;
#' Name of the variable marking the entry into our retrospective cohort
#' (i.e. in this case kidney cancer diagnosis)
item_starting_exact <- grep(item_starting_grep,names(d0),value = T);
#' Convert those columns to Yes/No values
d1[,class_yesno_exact] <- sapply(d1[,class_yesno_exact]
                                 ,function(xx){
                                   factor(is.na(xx)
                                          ,levels = c(F,T)
                                          ,labels = c('Yes','No'))}
                                 ,simplify = F);
#' Repeat for the T/F columns in this dataset
class_tf_tailgreps %>% paste0(collapse='|') %>% 
  grep(names(d0),val=T) -> class_tf_exact;
d1[,class_tf_exact] <- sapply(d1[,class_tf_exact],function(xx) !is.na(xx),simplify = F);
#' Create nominal values, binning the small groups into `other` using `cl_bintail()`
#' ...all as one command!
d1[,class_demog_exact] <- d1[,class_demog_exact] %>% 
  sapply(function(xx) cl_bintail(xx,topn=2),simplify=F);
#' Create a composite Hispanic column, `a_` prefix to signify 'analysis', the stage
#' during which this column gets created
# d1$a_hispanic <- d1$v020_Hspnc_or_Ltn | d1$language_cd=='spanish';
#' 
#' ## Scrap for later:
#' 
#' Finding valueflags for labs
# grep('\'vf\':\\[\'[LH]\'\\]',d0$v029_Prt_SrPl_mCnc_2885_2_info,val=T) %>% grepl('\'L\'',.) %>% ifelse(-1,1) %>% head
#' Looks like in this dataset the L and H vf's really are mutually exclusive on a per-visit basis
#'
#' ## TODO:
#' * Extract VF's for lab values and create flag columns
#' Create cancer and metastasis indicators (after re-running data-pull)
d1$a_metastasis <- d1[,class_diag_outcome_exact] %>% apply(1,any);
#' Create a copy of whatever the starting diagnosis column is called in
#' the current dataset, but this copy will always have the same name
d1$a_stdx <- d1[[item_starting_exact]];
#' * Look at fraction of patients covered by each variable
#' * Weed out patients who start out with an inactive cancer diag
unique(subset(d1,a_stdx=='Yes')$patient_num) -> pat_with_diag;
d2 <- subset(d1,patient_num%in%pat_with_diag);
#' * create the event indicators
d2 <- group_by(d2,patient_num) %>% 
  mutate(a_stdx1st = a_stdx=="Yes" & !duplicated(a_stdx=="Yes")
         ,a_metastasis1st = a_metastasis & !duplicated(a_metastasis)
         ,a_stdx_started = cumsum(a_stdx1st)
         ,a_cens_1 = lead(a_metastasis1st,1,default=0)
         ,a_age_at_stdx = age_at_visit_days[which(a_stdx1st)]
         ,a_dxage = age_at_visit_days - a_age_at_stdx
         ,a_metastasis_started = cumsum(a_metastasis1st));
# plot(survfit(Surv(a_dxage,a_cens_1)~I(a_age_at_stdx<21560),foo),xlim=c(0,4000),col=c('red','blue'))
# subset(summarise(d2
#                  ,v065in=min(which(v065_trpng_stmblng_inactive))
#                  ,v041in=min(which(v041_ACCDNTL_FLLS_inactive))
#                  ,v065=min(which(v065_trpng_stmblng)))
#        ,v065>v065in|v065>v041in|v065==1)$patient_num -> inactive_first;
# The following command shows the first occurrences o TRUE in the two inactive 
# falls columns and the 'real' falls column (v065). If the first inactive dates 
# are earlier than the first real event, those are IDs to toss because they are
# patients who may have come in after they alread fell for the first time. Also
# removing patients whose very first visit is in connection with a fall
# d2 <- subset(d2,!patient_num%in%inactive_first);
# to prevent individuals who only had one visit and it was not in connection
# with a fracture from being treated as missing data
# d2$cens <- ifelse(is.na(d2$cens),0,d2$cens)
# instead, though, we'll remove them
# d2 <- subset(d2,!is.na(cens));
#rebuild<-c();
#save.image(session);
#' ## Exploration of possibly combinable variables
#' The active diagnoses
#' ## Univariate models
#d3 <- subset(d2,npred>0|v055_Ofc_Vst|v056_Prcdr|cens>0);
#' The `patient_num`s of the patients who had an event, for later use in validation
