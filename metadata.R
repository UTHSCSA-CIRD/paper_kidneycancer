#' ---
#' title: "Kidney Cancer, Metadata"
#' author: "Guerra, Bokov, Wilson"
#' date: "07/11/2017"
#' ---
#' ## Non-analytic columns
#' ### Search term
class_nonan_tailgreps <- c('_date$','_info$','_unit$','_inactive$');
#' For the time being, do not treat these columns as analytic variables
class_complex_tailgreps <-c(
  '_Tbc_Usg$' # not a simple T/F, the binary indicator must be extracted from it
  ,'_Hstr_Dgns'
  ,'_Prvdr_Speclt$' # needs to be translated to human readable levels, also there
                    # are many of them, needs exploration
);
#' ## Classes of analytic columns
#' 
#' We will use the `item_` prefix to indicate matches for 
#' individual columns. The `_grep` suffix indicates that it
#' is a regexp rather than an exact match (because most 
#' columns have prefixes that change from one version of the
#' dataset to another).
item_spec_grep <- '_Prvdr_Spclt$';
#' We will use the `class_` prefix to indicate matches for
#' potentially multiple columns. The `_exact` suffix means 
#' that these are exact column names that do not change from
#' one version of the dataset to another (if you're using 
#' DataFinisher which this project is).
class_demog_exact <-c('sex_cd','language_cd','race_cd');
demcols <- c('patient_num','race_cd','language_cd','age_at_visit_days');

#' Living vs deceased
class_vital_grep <- c('_Dcsd_pr_SS$',"_Vtl_Sts$");
#' Tobacco usage
class_tobac_grep <- c('_Tbc_Usg$');
#class_hisp_grep <- c('_Hspnc_or_Ltn$|_Spnsh$');
class_hisp_grep <- '_Hspnc_or_Ltn$';
#' The `_tailgreps` suffix indicates that it is a vector of
#' regexps that can be combined into one regexp using 
#' `paste0(FOO,collapse='|')`
class_diag_tailgreps <- c('_elswhr_clsfd$','_Mls_and_ftg$'
                          ,"_Mlgnt_nplsm$");
#' Occurrence of malignant neoplasia using regex 
item_starting_grep <- "_Mlgnt_nplsm$";
class_diag_outcome_grep <-c("_Scndr_nrndcrn$","_mlgnt_unspcfd$"
, "_rsprtr_dgstv$","_unspcfd_mlgnt$");

class_vitals_tailgreps <- c(
  '_Bd_Ms_Indx_num$','_Tmprtr_F_num$'
  ,'_Sstlc_Prsr_num$','_Dstlc_Prsr_num$'
  ,'_Pls_num$','_Rsprtn_Rt_num$');
class_labs_tailgreps <- c(
  '_ALP_SrPl_cCnc_6768_6_num$','_ALT_SrPl_cCnc_1742_6_num$'
  ,'_AST_SrPl_cCnc_1920_8_num$','_Albmn_SrPl_mCnc_1751_7_num$'
  ,'_BN_SrPl_mCnc_3094_0_num$','_BN_Crt_SrPl_3097_3_num$'
  ,'_C_SrPl_sCnc_2028_9_num$','_Clcm_SrPl_mCnc_17861_6_num$'
  ,'_Chlrd_SrPl_sCnc_2075_0_num$','_Crt_SrPl_mCnc_2160_0_num$'
  ,'_Glcs_SrPl_mCnc_2345_7_num$','_Hct_VFr_Bld_At_4544_3_num$'
  ,'_Hgb_Bld_mCnc_GENERIC_KUH_COMPONENT_ID_3282_num$'
  ,'_MCH_RBC_Qn_At_GENERIC_KUH_COMPONENT_ID_4283_num$'
  ,'_MCHC_At_mCnc_GENERIC_KUH_COMPONENT_ID_4284_num$'
  ,'_MCV_RBC_At_GENERIC_KUH_COMPONENT_ID_4285_num$'
  ,'_Pltlt_At_GENERIC_KUH_COMPONENT_ID_5341_num$'
  ,'_Ptsm_SrPl_sCnc_2823_3_num$','_Prt_SrPl_mCnc_2885_2_num$'
  ,'_RBC__Bld_At_GENERIC_KUH_COMPONENT_ID_5638_num$'
  ,'_RDW_RBC_At_Rt_GENERIC_KUH_COMPONENT_ID_5629_num$'
  ,'_Sdm_SrPl_sCnc_2951_2_num$');
class_drugs_tailgreps <-c(
  '_EVRLMS_P_TBS$','_HDRCDN_ACTMNPHN$','_MTCLPRMD_TBS$'
  ,'_ONDNSTRN_TBS$','_PZPNB_HCL_TBS$','_SNTNB_MLT$'
  ,'_STNT_MG_P_CPS$','_AM_ANTMCRBLS$','_ANTNPLSTCS$'
  ,'_ANTNPLSTC_OTHR$','_NRVS_MDCTNS$','_CN_ANLGSCS$'
  ,'_OPD_ANLGSCS$','_GSTRNTSTNL$','_G_LXTVS$'
  ,'_LXTVS_OTHR$','_HRMNS_SNTHTCS_MDFRS$');
#' ## Factors and indicators
#' 
#' Columns that should have a YES/NO value (i.e. may directly be used in a model)
class_yesno_tailgreps <- c(class_drugs_tailgreps,class_diag_tailgreps);
#' Columns that should have a T/F values (i.e. used for logical operations)
class_tf_tailgreps <- c(class_vital_grep,class_hisp_grep,class_diag_outcome_grep);
#' Columns that should have multiple nomanal levels (and may need re-binning)
class_mult_exact <- c('language_cd','race_cd');
