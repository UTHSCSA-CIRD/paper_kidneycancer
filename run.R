#' ---
#' title: "Kidney Cancer"
#' author: "Guerra, Bokov, Wilson"
#' date: "07/10/2017"
#' ---
#' 
#' ## Load libraries
require(survival); require(dtplyr); require(magrittr); require("ggplot2"); require('MASS'); 
require('Hmisc'); require('readr');require(dplyr)
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
#' ## Extract VF's for lab values and create flag columns
#' New code for sapply lab values
#' Preliminary step -- exact names of columns containing full lab info
class_labs_tailgreps %>% gsub('_num','_info',.) %>% paste0(collapse = "|") %>%
  grep(names(d1), value=T)-> class_lab_info_exact;
#' Exact names of columns containing the value flags only
class_lab_vf_exact <- gsub('_info$','_vf',class_lab_info_exact);  

d1 <-sapply(d1[ ,class_lab_info_exact], function(xx){
  ifelse(grepl(pattern ="\'[HL]\'",x = xx),xx,"NONE") %>% 
    factor(levels=c("NONE","'vf':['L']","'vf':['H']"),labels= c("None","Low","High")) %>% relevel(ref = 'None')
  }) %>% data.frame() %>% 
  setNames(class_lab_vf_exact) %>% cbind(d1,.);


#' Create cancer and metastasis indicators (after re-running data-pull)
# Changing them to column names
class_diag_outcome_exact <- class_diag_outcome_grep  %>% paste0(collapse = "|") %>%  
grep(names(d1), value=T)

d1$a_metastasis <- d1[,class_diag_outcome_exact] %>% apply(1,any);
#' # TODO: create a metastasis OR deceased indicator, the trick is we need to use
#' one of the grep patterns, not hardcode. Never hardcode.
#' 
#' # TODO: double-check that patients patients for whom an inactive entry 
#' (for KC or any meta) occurs before the first active KC diagnosis get eliminated
#' 
#' Create a copy of whatever the starting diagnosis column is called in
#' the current dataset, but this copy will always have the same name
d1$a_stdx <- d1[[item_starting_exact]];
#' * Look at fraction of patients covered by each variable
#' * Weed out patients who start out with an inactive cancer diag
unique(subset(d1,a_stdx=='Yes')$patient_num) -> pat_with_diag;
d2 <- subset(d1,patient_num%in%pat_with_diag);
#' Create the event indicators
d2 <- group_by(d2,patient_num) %>% 
 mutate(
   # first KC diagnosis
   a_stdx1st = a_stdx=="Yes" & !duplicated(a_stdx=="Yes")
   # first metastasis diagnosis
   ,a_metastasis1st = a_metastasis & !duplicated(a_metastasis)
   # if this variable is 1 then this patient is in the interval occuring after first diagnosis
   ,a_stdx_started = cumsum(a_stdx1st)
   # whether or not this is the visit BEFORE the first metastasis diagnosis
   ,a_cens_1 = lead(a_metastasis1st,1,default=0)
   # age at first KC diagnosis
   ,a_age_at_stdx = age_at_visit_days[which(a_stdx1st)]
   # The number of days since first KC diagnosis
   ,a_dxage = age_at_visit_days - a_age_at_stdx
   # if this variable is 1 then this patient is in the interval occuring after first diagnosis
   ,a_metastasis_started = cumsum(a_metastasis1st));
#' Dataset with only the intervals between primary diagnosis and metastasis if any
d3 <- subset(d2, a_dxage>=0&a_metastasis_started==0);
#' Patients who only have one visit on or after diagnosis
summarise(d3,a=n()) %>% subset(a<2,select=patient_num) %>% 
  unlist -> pat_single_event;
#' Remove patients with only one visit from d3
d3 <- subset(d3,!patient_num%in%pat_single_event);
#' Example survival plot
d3[,c('patient_num','a_dxage','a_cens_1','a_age_at_stdx')] %>%
  mutate(a_dxage=last(a_dxage)
         ,a_cens_1=last(a_cens_1)
         ,a_age_at_stdx=last(a_age_at_stdx)) %>% 
  unique %>% 
  survfit(Surv(a_dxage,a_cens_1)~I(a_age_at_stdx<21560),.) %>% 
  plot(xlim=c(0,2000),col=c('red','blue'),mark.time=T);
#' ## We create out univariate baseline model to update a whole bunch of times soon
cox_univar<-coxph(Surv(a_dxage,a_cens_1) ~ a_age_at_stdx + cluster(patient_num),d3);
#' Example code... AFTER you have gone over all the revisions, see if you
#' can turn this into a non-hardcoded `sapply()` call.
sprintf('update(cox_univar,.~.-a_age_at_stdx+%s)','v026_Hct_VFr_Bld_At_4544_3_vf') %>% 
  parse(text=.) %>% eval;
cox_ph_models<-sapply(class_lab_vf_exact, function(xx) sprintf('update(cox_univar,.~.-a_age_at_stdx+%s)',xx) %>% 
         parse(text=.) %>% eval);

#Applying summary function uto the coxph models and finding their concordance, wald test, and 
sapply(cox_ph_models,function(xx) c(summary(xx)[['concordance']], summary(xx)[['waldtest']])) %>% t -> Concordance_and_Wald_results;

#' Hint: you don't need to understand `sprintf()` in order to do this, just look and
#' think about what the one part of the above code should not be a static value.
#' However, even though you don't need to understand `sprintf()` you should do
#' `?sprintf` anyway, and do that with any function you don't understand as a
#' general rule.
#' 
#' ## Interlude: thought process for constructing sapply statements
#' 
#' Let's say you have a function `foo()` that takes `xx` as a variable
#' and you have an object `baz`, and you are wondering whether to run 
#' `sapply(baz,foo)`. If the following doesn't work...
# foo(baz[1])
#' ...then `sapply(baz,foo)` will definitely not work.
#' If you can find a `baz` for which `foo(baz[1])` does give a valid result, and
#' so do `foo(baz[2])` and `foo(baz[3])` then the next step in testing might be
# sapply(baz[1:3], foo)

#' 
#rebuild<-c();
#save.image(session);
