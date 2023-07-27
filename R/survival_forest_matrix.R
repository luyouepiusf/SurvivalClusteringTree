##' Build a Survival Forest (Data Supplied as Matrices)
##' 
##' @title Build a Survival Forest (Data Supplied as Matrices)
##' @description The function 
##' \code{survival_forest_matrix} build a survival forest given the survival outcomes and predictors of numeric and factor variables.
##' @param time survival times, a numeric vector. 
##' \code{time[i]} is the survival time of the ith sample.
##' @param event survival events, a logical vector. 
##' \code{event[i]} is the survival event of the ith sample.
##' @param matrix_numeric numeric predictors, a numeric matrix. 
##' \code{matrix_numeric[i,j]} is the jth numeric predictor of the ith sample.
##' @param matrix_factor factor predictors, a character matrix. 
##' \code{matrix_factor[i,j]} is the jth predictor of the ith sample.
##' @param weights sample weights, a numeric vector. 
##' \code{weights[i]} is the weight of the ith sample.
##' @param significance significance threshold, a numeric value. 
##' Stop the splitting algorithm when no splits give a p-value smaller than \code{significance}.
##' @param min_weights minimum weight threshold, a numeric value. 
##' The weights in a node are greater than \code{min_weights}.
##' @param missing a character value that specifies the handling of missing data. 
##' If \code{missing=="omit"}, samples with missing values in the splitting variables will be discarded.
##' If \code{missing=="majority"}, samples with missing values in the splitting variables will be assigned to the majority node.
##' If \code{missing=="weighted"}, samples with missing values in the splitting variables will be weighted by the weights of branch nodes.
##' @param test_type a character value that specifies the type of statistical tests.
##' If \code{test_type=="univariate"}, then it performs a log-rank test without p-value adjustments.
##' If \code{test_type} is in \code{p.adjust.methods}, i.e., one of holm, hochberg, hommel, bonferroni, BH, BY, or fdr, 
##' then the p-values will be adjusted using the corresponding method.
##' @param nboot an integer value that specifies the number of bootstrap replications.
##' @param seed an integer value that specifies the seed.
survival_forest_matrix<-function(
  time,
  event,
  matrix_numeric,
  matrix_factor,
  weights=rep(1,length(time)),
  significance=0.05,
  min_weights=50,
  missing="omit",
  test_type="univariate",
  nboot=100,
  seed=0){
  
  # check [missing], [test_type]
  if(!missing%in%c("majority","omit","weighted"))stop("Invalid 'missing' argument.")
  if(!test_type%in%c("univariate",p.adjust.methods))stop("Invalid 'test_type' argument.")
  
  event<-as.logical(event)
  matrix_numeric<-as.matrix(matrix_numeric)
  matrix_factor<-as.matrix(matrix_factor)
  matrix_factor<-apply(matrix_factor,c(1,2),as.character)
  ndim_numeric<-ncol(matrix_numeric)
  ndim_factor<-ncol(matrix_factor)
  nind<-length(time)
  
  if(!is.numeric(time))stop("Invalid 'time' argument")
  if(!is.logical(event))stop("Invalid 'event' argument")
  if(!is.numeric(matrix_numeric))stop("Invalid 'matrix_numeric' formula")
  if(!is.character(matrix_factor))stop("Invalid 'matrix_factor' formula")
  
  if(ndim_numeric+ndim_factor<1)stop("There are no predictors in the model.")
  
  # check dimensions
  if(nrow(matrix_numeric)!=nind)stop("Dimension mismatch between 'matrix_numeric' and 'time'.")
  if(nrow(matrix_factor)!=nind)stop("Dimension mismatch between 'matrix_factor' and 'time'.")
  if(any(is.na(time))|any(is.na(event)))stop("Missing values in 'time' or 'event'.")
  
  # create names
  if(is.null(colnames(matrix_numeric))&ndim_numeric>=1)colnames(matrix_numeric)<-paste0("numeric",1:ncol(matrix_numeric),sep="")
  if(is.null(colnames(matrix_factor))&ndim_factor>=1)colnames(matrix_factor)<-paste0("factor",1:ncol(matrix_factor),sep="")
  variable_names<-c(colnames(matrix_numeric),colnames(matrix_factor))
  
  # convert matrix_factor to an integer matrix
  # factor_dictionary<-list()
  # matrix_factor_int<-matrix(NA,nind,ndim_factor)
  # if(ncol(matrix_factor)>0){
  #   colnames(matrix_factor_int)<-colnames(matrix_factor)
  #   for(idx in 1:ncol(matrix_factor)){
  #     aname<-colnames(matrix_factor)[idx]
  #     a_dictionary<-create_dictionary(matrix_factor[,idx])
  #     factor_dictionary[[aname]]<-a_dictionary
  #     matrix_factor_int[,idx]<-a_dictionary[matrix_factor[,aname]]
  #   }
  # }
  # matrix_factor<-matrix_factor_int
  
  a_survival_forest<-list()
  set.seed(seed)
  for(boot_idx in 1:nboot){
    cat(boot_idx)
    shuffle<-sample(1:nind,nind,replace=T)
    time_boot<-time[shuffle]
    event_boot<-event[shuffle]
    matrix_numeric_boot<-matrix_numeric[shuffle,,drop=FALSE]
    matrix_factor_boot<-matrix_factor[shuffle,,drop=FALSE]
    weights_boot<-weights[shuffle]
    
    a_survival_tree<-survival_tree_matrix(
      time=time_boot,
      event=event_boot,
      matrix_numeric=matrix_numeric_boot,
      matrix_factor=matrix_factor_boot,
      weights=weights_boot,
      significance=significance,
      min_weights=min_weights,
      missing=missing,
      test_type=test_type)
    cat(" - ")
    a_survival_forest[[boot_idx]]<-a_survival_tree
  }
  
  return(list(
    variable_names=variable_names,
    ndim_numeric=ndim_numeric,
    ndim_factor=ndim_factor,
    survival_forest=a_survival_forest))
}