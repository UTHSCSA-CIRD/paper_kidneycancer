#' ---
#' title: "Kidney Cancer Main Analysis"
#' author: "Guerra, Bokov, Wilson"
#' date: "07/10/2017"
#' ---
#' 
#' ## Load libraries
require(survival); require(magrittr); require("dplyr");
require('data.table'); require("ggplot2"); require('MASS'); 
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
demcols <- c('patient_num','race_cd','language_cd','age_at_visit_days');
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
#' Convert those columns to Yes/No values
d1[,class_yesno_exact] <- sapply(d1[,class_yesno_exact]
                                 ,function(xx) 
                                   factor(is.na(xx)
                                          ,levels = c(F,T)
                                          ,labels = c('Yes','No'))
                                 ,simplify = F);
#' Repeat for the T/F columns in this dataset
class_tf_tailgreps %>% paste0(collapse='|') %>% 
  grep(names(d0),val=T) -> class_tf_exact;
d1[,class_tf_exact] <- sapply(d1[,class_tf_exact],function(xx) !is.na(xx),simplify = F);
#' Create a composite Hispanic column, `a_` prefix to signify 'analysis', the stage
#' during which this column gets created
d1$a_hispanic <- d1$v020_Hspnc_or_Ltn | d1$language_cd=='spanish';
#' 
#' ## Scrap for later:
#' 
#' Finding valueflags for labs
# grep('\'vf\':\\[\'[LH]\'\\]',d0$v029_Prt_SrPl_mCnc_2885_2_info,val=T) %>% grepl('\'L\'',.) %>% ifelse(-1,1) %>% head
#' Looks like in this dataset the L and H vf's really are mutually exclusive on a per-visit basis
#'
#' ## TODO:
#' * Extract VF's for lab values and create flag columns
#' * Create cancer and metastasis indicators (after re-running data-pull)
#' * Weed out patients who start out with an inactive cancer diag
#' * Create (in metadata.R) a list of predictor variables
#' * Look at fraction of patients covered by each variable
# create the event indicators
# group_by(d0,patient_num) %>% 
#   mutate(tt=age_at_visit_days
#          ,which_event=cumsum(v065_trpng_stmblng)
#          ,cens=lead(which_event)
#          ,first=c(1,rep_len(0,length(tt)-1))
#          ,last=c(rep_len(0,length(tt)-1),1)
#          # preserving for each individual their age at first visit
#          ,agestart=min(age_at_visit_days)) -> d1; # 100072 x 56
# d1 = dataset with cumulative counts and censoring indicators for all falls
# d2 = dataset with observations _prior_ to the first event only (or where no
# events have occurred)
# d2 <- subset(d1,which_event==0);  # 87110 x 56
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
