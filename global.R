#' ---
#' title: "Kidney Cancer Project Settings"
#' author: "Bokov, Michalek, Shah, Kaushig, Rodriguez"
#' date: "02/19/2018"
#' ---
#' 
#' ## Load libraries
#+ warning=FALSE, message=FALSE
rq_libs <- c('compiler'                            # just-in-time compilation
             ,'survival','MASS','Hmisc','zoo'      # various analysis methods
             ,'readr','dplyr','stringr','magrittr' # data manipulation & piping
             ,'ggplot2','ggfortify','grid'         # plotting
             ,'stargazer','broom');                # table formatting
rq_installed <- sapply(rq_libs,require,character.only=T);
rq_need <- names(rq_installed[!rq_installed]);
if(length(rq_need)>0) install.packages(rq_need,repos='https://cran.rstudio.com/',dependencies = T);
sapply(rq_need,require,character.only=T);
#' Turn JIT to max: pre-compile all closures, `for`, `while`, and `repeat` loops
enableJIT(3);
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
#' possible non mutually exclusive values for `rebuild`: `d0`, `aic00`
rebuild <- c();
#' ## Load data
