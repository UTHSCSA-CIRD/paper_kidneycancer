gitstamp <- function(production=T) {
  if(production){
    if(length(gitdiff<-system("git update-index --refresh && git diff-index HEAD --",intern = T))!=0) stop(sprintf(
      "\ngit message: %s\n\nYou have uncommitted changes. Please do 'git commit' and then try again."
      ,gitdiff));
    system("git push && git log --pretty=format:'%h' -n 1",intern=T);
  } else system("git log --pretty=format:'%h' -n 1",intern=T);
}

#' Return a factor containing the top-N levels, and the rest binned
#' into the specified existing or new level. That level is listed last
#' in the resulting factor.
#' @param xx A \code{vector} (required).
#' @param topn \code{numeric} for the top how many levels to keep (optional, default =4)
#' @param binto \code{character} which new or existing level to dump the other values in
cl_bintail <- function(xx,topn=4,binto='other'){
  if(!is.factor(xx)) xx <- factor(xx);
  counts <- sort(table(xx),decreasing = T);
  if(is.numeric(binto)) binto <- names(counts)[binto];
  keep <- setdiff(names(counts)[1:min(length(counts),topn)],binto);
  droplevels(
    factor(
      ifelse(
        xx %in% keep, as.character(xx), binto
        ),levels=c(keep,binto)));
}

autoplot <- ggfortify:::autoplot.survfit;
formals(autoplot)$conf.int<-F;
formals(autoplot)$surv.size<-2;
formals(autoplot)$xlab<-'Days Since Diagnosis';
formals(autoplot)$ylab<-'% Metastasis Free';
formals(autoplot)$xlim<-c(0,1500);


#' Take a character vector and perform multiple search-replace 
#' operations on it.
#' @param xx A \code{vector} of type \code{character} (required)
#' @param searchrep A \code{matrix} with two columns of type \code{character} (required). The left column is the pattern and the right, the replacement.
#' @param method One of 'partial','full', or 'exact'. Controls whether to replace only the matching regexp, replace the entire value that contains a matching regexp, or replace the entire value if it's an exact match.
submulti <- function(xx,searchrep,method=c('partial','full','exact')){
  # if no method is specified by the user, this makes it take the first value
  # if a method is only partially written out, this completes it, and if the
  # method doesn't match any of the valid possibilities this gives an informativ
  # error message
  method<-match.arg(method);
  # rr is a sequence of integers spanning the rows of the searchrep matrix
  rr <- 1:nrow(searchrep);
  # oo will hold the result that this function will return
  oo <- xx;
  switch(method
         ,partial = {for(ii in rr)
           oo <- gsub(searchrep[ii,1],searchrep[ii,2],oo)}
         ,full =    {for(ii in rr)
           oo[grepl(searchrep[ii,1],oo)]<-searchrep[ii,2]}
         ,exact = {for(ii in rr)
           oo[grepl(searchrep[ii,1],oo,fixed=T)]<-searchrep[ii,2]}
           #oo <- gsub(searchrep[ii,1],searchrep[ii,2],oo,fixed = T)}
         );
  oo;
}

#' Take a data.frame and create a PL/SQL table definition for it, ready to paste into DBVis and such
#' @param xx A \code{data.frame} (required)
sql_create_table <- function(xx,tname,sqltemplate='CREATE TABLE %s (%s);'){
  mysearchrep <- cbind(c('integer','character','numeric','Date','logical','factor')
    ,c('NUMBER(32)','VARCHAR2(nn)','NUMBER(32)','DATE','NUMBER(1)','VARCHAR2(nn)'));
  if(missing(tname)) tname <- as.character(substitute(xx));
  thenames <- sapply(xx, class);
  thenames <- submulti(thenames, mysearchrep, 'exact');
  strlengthres<-apply(nsqip,2,function(xx){max(str_length(xx))});
  thenames <- apply(cbind(round(strlengthres*1.5), thenames),1,function(xx){gsub('nn',xx[1],xx[2])});
  thenames <- paste(' ',names(thenames), as.character(thenames), sep=' ', collapse= ',\n');
  thenames <- gsub('.', '_', thenames, fixed=TRUE);
  thenames <- gsub('_{1,}','_',thenames);
  thenames <- gsub('_ ',' ',thenames);
  thenames <- sprintf(sqltemplate, tname, thenames);
  cat(thenames);
  invisible(thenames);
  # TODO: collect the column names and classes as we did in console
  # TODO: use the above and submulti() to convert R classes to SQL data types
  # TODO: paste() and/or sprintf() to create the string that will be the output
}

#' Take a data.frame or character vector and a vector of grep targets and return
#' the values that match (for data.frame, column names that match). If no 
#' patterns given just returns the names
#' @param xx A \code{data.frame} or character vector (required)
#' @param patterns A character vector of regexp targets to be OR-ed
grepor <- function(xx,patterns='.') {
  if(is.list(xx)) xx <-names(xx);
  grep(paste0(patterns,collapse='|'),xx,val=T);
}

#' Delete all the junk in your environment, for testing
clearenv <- function(env=.GlobalEnv) rm(list=setdiff(ls(all=T,envir=env),'clearenv'),envir=env);

#' From ... http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  # Make the panel
  # ncol: Number of columns of plots
  # nrow: Number of rows needed, calculated from # of cols
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))}
  if (numPlots==1) print(plots[[1]]) else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} 

#' ## Functions unlikely to be useful in contexts other than this one
#'
#' Takes a vector of filenames, each of which is a saved R environment from a
#' session where stepAIC was run on resampled data and return a frequency table
#' as well as a table of parameter estimates
#' @param xx    character vector of filenames/paths of rdata files
#' @param type  the class that a non-error fit should have (default: 'coxph') 
summ_aicresamp <- function(xx,type='coxph',...){
  yy <- sapply(xx,function(xx) new.env());
  sapply(xx,function(xx) load(xx,envir=yy[[xx]]));
  # Which model-fits ran and returned results without error?
  sapply(yy,function(xx) sapply(xx$aic_resampled,class)=='coxph') -> goodones;
  # Combine the good fits into one list
  mapply(function(aa,bb) aa$aic_resampled[bb],yy,goodones) %>% do.call(c,.) -> fits;
  # Extract the term tables from this list
  sapply(fits, function(xx) tidy(xx),simplify=F) -> trms;
  # Tally up frequency of each term
  sapply(trms,`[`,,1) %>% unlist %>% table %>% sort(decr=T) -> trmcounts;
  # Create matrix of parameter estimates
  tmx <- matrix(0,nrow=length(fits),ncol=length(trmcounts));
  colnames(tmx) <- names(trmcounts);
  for(ii in seq(nrow(tmx))) tmx[ii,trms[[ii]]$term] <- trms[[ii]]$estimate;
  apply(tmx,2,function(xx) c(
    quantile(xx,c(.025,.25,.5,.75,.975),na.rm=T)
    ,mean=mean(xx,na.rm=T)
    ,se=sd(xx,na.rm=T)
    ,t=mean(xx,na.rm=T)/sd(xx,na.rm=T)
    ,npresent = sum(abs(xx)>0,na.rm=T)
    ,notna=sum(!is.na(xx)))) -> aggregates
  # Save the entries from rows_resampled for which models failed to fit, and the 
  # ones where they did not fail.
  do.call(c,mapply(function(aa,bb) aa$rows_resampled[!bb],yy,goodones)) -> badoffsets;
  do.call(c,mapply(function(aa,bb) aa$rows_resampled[bb],yy,goodones)) -> goodoffsets;
  out <- list(trmcounts=trmcounts,tmx=tmx,aggregates=aggregates
              ,badoffsets=badoffsets,goodoffsets=goodoffsets);
  invisible(out);
}
