.temp_coxph_mv1_call <-
quote(coxph(formula = Surv(a_dxage, a_dxage2, a_cens_1) ~ v029_Hspnc_or_Ltn + 
    race_cd + v014_Chlrd_SrPl_sCnc_2075_0_numnona + v016_Albmn_SrPl_mCnc_1751_7_numnona + 
    v053_Sdm_SrPl_sCnc_2951_2_numnona + v024_elswhr_clsfd + v028_Rsprtn_Rt_numnona + 
    v043_Ptsm_SrPl_sCnc_2823_3_numnona + v034_Tmprtr_F_numnona + 
    v019_ALP_SrPl_cCnc_6768_6_numnona + v010_ALT_SrPl_cCnc_1742_6_numnona + 
    v020_Clcm_SrPl_mCnc_17861_6_numnona + v033_MCH_RBC_Qn_At_GENERIC_KUH_COMPONENT_ID_4283_numnona, 
    data = d3))
