#' ---
#' title: "Kidney Cancer Main Analysis"
#' author: "Guerra, Bokov, Wilson"
#' date: "07/10/2017"
#' ---
#' 
#' ## Load libraries
require(magrittr); require("dtplyr");
require("ggplot2");require('readr');
#' ## Load local config file
source('./config.R');
#' ## Load data
#' 
#' This contains d0, which are counts by age, sex, and reference group
load(inci_inputdata);
d0$agebin <- cut(d0$AGE,c(0,10*(2:7),Inf),include.lowest = T);
#' ## First stage of tabulation
#' 
d0 %>% group_by(Reference,agebin) %>% 
  summarise(N=sum(N_REF),n=sum(N_GRP),pct=n/N) -> d1;
#' ## Age tables
#' 
cbind(
  subset(data.frame(d1),grepl('^NHW',Reference))
  ,subset(data.frame(d1),grepl('^HIS',Reference))[,-c(1:2)]
  );
#' ...and then I got lazy and copy-pasted to spreadsheet.
#' 
