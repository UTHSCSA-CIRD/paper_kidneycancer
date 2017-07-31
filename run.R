#' ---
#' title: "Kidney Cancer Main Analysis"
#' author: "Guerra, Bokov, Wilson"
#' date: "07/10/2017"
#' ---
#' 
#' ## Load libraries
#+ warning=FALSE, message=FALSE
rq_libs <- c('survival','MASS','Hmisc','zoo'       # various analysis methods
             ,'readr','dplyr','stringr','magrittr' # data manipulation & piping
             ,'ggplot2','ggfortify','grid'         # plotting
             ,'stargazer');                        # table formatting
rq_installed <- sapply(rq_libs,require,character.only=T);
rq_need <- names(rq_installed[!rq_installed]);
if(length(rq_need)>0) install.packages(rq_need,repos='https://cran.rstudio.com/',dependencies = T);
sapply(rq_need,require,character.only=T);
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
if(file.exists(inputmeta)) m0 <- read_csv(inputmeta,na='');
#' ## Clean up data
if('d0' %in% rebuild) {
  # Remove impossible START_DATEs
  d0 <- subset(d0,age_at_visit_days > 0);
}
#' ## Create the groups of exact column names for this dataset
#' 
#' Obtain the actual column names for the Yes/No columns in this dataset
class_yesno_tailgreps %>% paste0(collapse='|') %>% 
  grep(names(d0),val=T) -> class_yesno_exact;
#' Repeat for the T/F columns in this dataset
class_tf_tailgreps %>% paste0(collapse='|') %>% 
  grep(names(d0),val=T) -> class_tf_exact;
#' Columns indicating Hispanic ethnicity
class_hisp_grep %>% paste0(collapse='|') %>% 
  grep(names(d0),val=T) -> class_hisp_exact;
#' Name of the variable marking the entry into our retrospective cohort
#' (i.e. in this case kidney cancer diagnosis)
item_starting_exact <- grep(item_starting_grep,names(d0),value = T);
#' names of lab value columns and again for vitals
class_labs_tailgreps %>% paste0(collapse = "|") %>%
  grep(names(d0), value=T)-> class_labs_exact;
class_vitals_tailgreps %>% paste0(collapse = "|") %>%
  grep(names(d0), value=T)-> class_vitals_exact;
#' Preliminary step -- exact names of columns containing full lab info
class_labs_tailgreps %>% gsub('_num','_info',.) %>% paste0(collapse = "|") %>%
  grep(names(d0), value=T)-> class_lab_info_exact;
#' Exact names of columns containing the value flags only
class_lab_vf_exact <- gsub('_info$','_vf',class_lab_info_exact);
#' Create cancer and metastasis indicators
# Changing them to column names
class_diag_outcome_exact <- class_diag_outcome_grep  %>% paste0(collapse = "|") %>%  
  grep(names(d0), value=T);
class_locf_exact <- c(class_labs_exact,class_vitals_exact);

#' ## Convert columns
#' 
#' Create copy of original dataset
d1 <- d0;
#' Convert those columns to Yes/No values
d1[,class_yesno_exact] <- sapply(d1[,class_yesno_exact]
                                 ,function(xx){
                                   factor(is.na(xx)
                                          ,levels = c(F,T)
                                          ,labels = c('Yes','No'))}
                                 ,simplify = F);
d1[,class_tf_exact] <- sapply(d1[,class_tf_exact],function(xx) !is.na(xx),simplify = F);
#' Create nominal values, binning the small groups into `other` using `cl_bintail()`
#' ...all as one command!
d1[,class_demog_exact] <- d1[,class_demog_exact] %>% 
  sapply(function(xx) cl_bintail(xx,topn=2),simplify=F);
#' ## Extract VF's for lab values and create flag columns
#' New code for sapply lab values

d1 <-sapply(d1[ ,class_lab_info_exact], function(xx){
  ifelse(grepl(pattern ="\'[HL]\'",x = xx),xx,"NONE") %>% 
    factor(levels=c("NONE","'vf':['L']","'vf':['H']"),labels= c("None","Low","High")) %>% relevel(ref = 'None')
  }) %>% data.frame() %>% 
  setNames(class_lab_vf_exact) %>% cbind(d1,.);

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
   ,a_dxage2 = lead(a_dxage,default=Inf)
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
#' ## Fill in missing values!!
d4 <- mutate_all(d3[,c('patient_num',class_locf_exact)],na.locf0) %>% 
  mutate_all(na.locf0,fromLast=T);
#' There are _still_ missing values, so fill them in with medians
d4[,-1] <- lapply(d4[,-1],na.aggregate,median);
#' Make the names not colide with existing ones
names(d4)[-1] <- paste0(names(d4)[-1],'nona');
#' Now mush these back into d3
d3[,names(d4)[-1]]<-d4[,-1];
#' ## We create out univariate baseline model to update a whole bunch of times soon
cox_univar<-coxph(Surv(a_dxage,a_cens_1) ~ a_age_at_stdx + cluster(patient_num),d3);

##Plotting the survival curve using the package "ggfortify"
survival_Curve<-d3[,c('patient_num','a_dxage','a_cens_1','a_age_at_stdx')] %>%mutate(a_dxage=last(a_dxage) ,a_cens_1=last(a_cens_1),a_age_at_stdx=last(a_age_at_stdx)) %>% unique %>% survfit(Surv(a_dxage,a_cens_1)~I(a_age_at_stdx<21560),.)
autoplot(survival_Curve);

#Cox_Univariate and Cox_Univariate_Frailty
cox_univar<-coxph(Surv(a_dxage,a_cens_1) ~ a_age_at_stdx + cluster(patient_num),d3);
cox_univar_numeric<-coxph(Surv(a_dxage,a_dxage2,a_cens_1) ~ a_age_at_stdx,d3);

#' Cox_ph_models with the vfs and yes and no
#+ warning=FALSE
cox_ph_models<-sapply(c(class_lab_vf_exact,class_yesno_exact)
                      ,function(xx) sprintf('update(cox_univar,.~.-a_age_at_stdx+%s)',xx) %>% 
                        parse(text=.) %>% eval);
#' maybe the interval formulation is a better
#' way to handle time-varying covariates than cluster. So let's try this one too.
#' 
#+ warning=FALSE
cox_t2_models<- lapply(cox_ph_models
                       ,update
                       ,Surv(a_dxage,a_dxage2,a_cens_1)~.-cluster(patient_num));

#+ warning=FALSE
cox_ph_models_fraility<-sapply(cox_ph_models,function(xx) update(xx,.~.-cluster(patient_num)+frailty(patient_num)));
#' Create univariate models with numeric predictors
cox_ph_models_numeric <- sapply(paste0(class_locf_exact,'nona')
                         ,function(xx) sprintf('update(cox_univar_numeric,.~%s)',xx) %>% 
                           parse(text=.) %>% eval,simplify = F);
#' demographics (non-time dependent)
cox_t2_demog <- sapply(c(class_demog_exact,class_hisp_exact)
                       ,function(xx) sprintf('update(cox_univar_numeric,.~%s)',xx) %>%
                         parse(text=.) %>% eval,simplify = F);
#' Some static models
#' The last visits only, for plotting
d5 <- summarise_all(d3,last);
d5 <- lapply(cox_ph_models_numeric,predict,d5,type='lp') %>% 
  lapply(function(xx) xx>median(xx)) %>% 
  setNames(.,gsub('nona','lp',names(.))) %>% data.frame %>% cbind(d5,.)

#' ## Results
#' 
#' frailty
results_con_wald_frail <- sapply(cox_ph_models_fraility
                                 ,function(xx) 
                                   with(summary(xx),c(concordance,logtest))) %>% t;
#' cluster
results_con_wald_cluster <- sapply(cox_ph_models
                                   ,function(xx) 
                                     with(summary(xx),c(concordance,logtest))) %>% t;
#' intervals
results_con_wald_t2 <- sapply(cox_t2_models
                                   ,function(xx) 
                                     with(summary(xx),c(concordance,logtest))) %>% t;
#' numeric
results_con_wald_numeric <- sapply(cox_ph_models_numeric
                              ,function(xx) 
                                with(summary(xx),c(concordance,logtest))) %>% t;

#' New Rownames for the fraility and cluster table
rownames(results_con_wald_t2) <- 
  rownames(results_con_wald_frail) <- 
  rownames(results_con_wald_cluster) <- 
  submulti(gsub('_vf$','_num'
                ,rownames(results_con_wald_cluster)),m0[,1:2],method = 'exact');
#' These specifically for the numeric models 
rownames(results_con_wald_numeric) <- 
  submulti(gsub('nona$',''
                ,rownames(results_con_wald_numeric)),m0[,1:2],method = 'exact'); 
#' Comparing the numeric to the categoric, for labs only
rows_labs <- intersect(rownames(results_con_wald_numeric)
                       ,rownames(results_con_wald_t2));
rows_vitals <- setdiff(rownames(results_con_wald_numeric),rows_labs);
rows_other <- setdiff(rownames(results_con_wald_cluster),rows_labs);
#' ## Do the Wilcoxon tests
#' 
#' Right-censored (and clustered) vs interval-censored, goodness-of-fit
wilcox.test(results_con_wald_cluster[,3]
            ,results_con_wald_t2[,3],paired = T,conf.int=T);
#' Right-censored (and clustered) vs interval-censored, concordance
wilcox.test(results_con_wald_cluster[,1]
            ,results_con_wald_t2[,1],paired = T,conf.int=T);
#' Right-censored (and clustered) vs interval-censored, goodness-of-fit
wilcox.test(results_con_wald_cluster[,3]
            ,results_con_wald_t2[,3],paired = T,conf.int=T);
#' discretized versus LOCF, goodness-of-fit
wilcox.test(results_con_wald_t2[rows_labs,3]
            ,results_con_wald_numeric[rows_labs,3],paired = T,conf.int=T);
#' discretized versus LOCF, concordance
wilcox.test(results_con_wald_t2[rows_labs,1]
            ,results_con_wald_numeric[rows_labs,1],paired = T,conf.int=T);
#' Let's try printing out the tables of concordances and goodness-of-fit
#+ results='asis'
if(!interactive()){
  cat('\nFrailty\n');
  stargazer(results_con_wald_frail,type = 'html');
  cat('Cluster\n');
  stargazer(results_con_wald_cluster,type = 'html');
  cat('Interval\n');
  stargazer(results_con_wald_t2,type = 'html');
  cat('Numeric\n');
  stargazer(results_con_wald_numeric,results_con_wald_cluster,type = 'html');
  cat('Labs, Categoric and Numeric predictors side by side\n');
  stargazer(data.frame(
    Numeric=results_con_wald_numeric[rows_labs,]
    ,Categoric=results_con_wald_t2[rows_labs,]),type = 'html',summary = F);
};
#' Also, here is hispanic ethnicity, as  predictor in an interval-censored model
#+ results='asis'
sapply(class_hisp_exact
       ,function(xx) stargazer(cox_t2_demog[[xx]]
                              ,type=if(interactive()) 'text' else 'html')) -> .junk;
#' Survival plot for Hispanic vs Non Hispanic
pred_hisp <- predict(cox_t2_demog[[class_hisp_exact[1]]],d5);
autoplot(survfit(Surv(a_dxage,a_cens_1)~pred_hisp,d5),col=c('red','blue')) + 
  scale_fill_discrete('Ethnicity',labels=c('Non Hispanic','Hispanic')) + 
  scale_color_discrete('Ethnicity',labels=c('Non Hispanic','Hispanic'));
#' ## Survival plots for numeric predictors
#+ warning=FALSE
plots_cph_numeric <- grep('lp$',names(d5),val=T) %>% sapply(function(xx) 
  sprintf("autoplot(survfit(Surv(a_dxage,a_cens_1)~%s,d5),col=c('red','blue'),mark.time = T,xlim=c(0,1500))",xx) %>% 
    parse(text=.) %>% eval,simplify = F) %>% 
  setNames(.,gsub('lp$','',names(.)) %>% submulti(m0[,1:2]))
plots_cph_numeric <- sapply(names(plots_cph_numeric)
                            ,function(xx) plots_cph_numeric[[xx]] + 
                              theme(legend.position = 'none')+ggtitle(xx),simplify=F);
multiplot(plotlist=plots_cph_numeric,cols=5);
#' # The Multivariable Model
#' 
#' Based on some ad-hoc exploration (picking the univariate predictors with the 
#' best concordances and a few facially reasonable demographic predictors,
#' fitting them all, and seeing which were significant) here is the _starting_ 
#' model. _Not_ the final one.
coxph_mv0 <- coxph(formula = Surv(a_dxage, a_dxage2, a_cens_1) ~ v029_Hspnc_or_Ltn + 
                     v037_CN_ANLGSCS + v050_RDW_RBC_At_Rt_GENERIC_KUH_COMPONENT_ID_5629_numnona, 
                   data = d3);
#' ## Set up the column names
#' 
#' The variables used on our initial multivariate model
class_mv0_tailgreps %>% paste0(collapse='|') %>% 
  grep(names(d3),val=T) -> class_mv0_exact;
#' Possible variables to consider for addition
class_mv1_candidates_exact <- c(demcols[1:3],'a_age_at_stdx'
                          ,paste0(class_locf_exact,'nona')
                          ,class_yesno_exact);
#' Remove the ones which are already in mv0
class_mv1_candidates_exact <- setdiff(class_mv1_candidates_exact,class_mv0_exact);
#' Fire up stepAIC and get a cup of coffee!
# We use sprintf to construct the expression out of a fixed string and 
# a pasted-together vector of variable names, because it's a very long
# expression and as elewhere in this script, the variable names will likely
# change from time to time
sprintf(
  # This is a complete stepAIC call, missing only the additiona candidate 
  # variables to try. Those will go where the %s currently is. The 'lower' part
  # of the 'scope' argument means "be willing to drop any variable completely if
  # it does not improve model fit". The 'upper' part means "consider adding these
  # other variables to the existing ones (where the other variables will replace
  # %s shortly) AND also consider every possible two-way interaction between
  # variables you keep. See what I mean when I say this will take a while?
  # After the 'list' argument there is a 'direction' argument, and 'both' means 
  # we will add and remove variables.
  'stepAIC(coxph_mv0,scope = list(lower=.~1,upper=.~(.+%s)^2),direction="both")'
  # This paste statement simply combines together the names of the additional 
  # candidate variables generated above with a '+' between them.
  ,paste0(class_mv1_candidates_exact,collapse='+')) %>% 
  # Now we pipe this string to parse which turns it into a call and eval which
  # attempts to execute that call. Good luck to us!
  parse(text=.) %>% eval -> coxph_mv1;

#' Note that you can also generate a big version of any of these manually
#' by doing e.g. `plots_cph_numeric[[10]]` or `plots_cph_numeric[["AST SerPl-cCnc (1920-8)"]]`
#' 
#' 
#' More hints:
#' 
#' * Most of the above are variations of stuff you have already done. The names 
#' of variables will be different, but that doesn't matter, the underlying logic 
#' of the code will be the same or almost the same.
#' * When you don't understand something, try breaking it up into smaller expressions
#' running them separately, and seeing what each one does.
#' * When you run into an error, try breaking the expression into smaller experssions
#' running them separately, and seeing if you can isolate the error to a specific one.
#' * The help operator, `?` is helpful (sometimes) for understanding what a function
#' is supposed to do. The help files are organized in sections, so try to develop
#' a feel for which section has answers to which kinds of questions you might have.
#' * The example section of a help file is a great way to figure out what a function 
#' does by starting with their working example and gradually substituting in your
#' own input until the function is doing what you need it to do (or failing and 
#' giving error messages you can copy-paste into google).
#' * You can capture the output of an unknown function into a variable and the 
#' following commands can help you investigate it...
#'  * class(foo)
#'  * names(foo)
#'  * sapply(foo,class)
#'  * methods(class=class(foo))
#'  * summary(foo)
#'  * ...and all the above for bar where bar is an element within foo
#' that you are interested in and can be indexed either as foo[['bar']]
#' or foo$bar (the two mean the same thing)
#' * Strings need to be quoted, the names of objects need to not be quoted.


#' 


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
