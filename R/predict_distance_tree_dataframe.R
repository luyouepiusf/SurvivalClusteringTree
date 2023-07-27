##' Predict Distances Between Samples Based on a Survival Tree Fit (Data Supplied as a Dataframe)
##' 
##' @title Predict Distances Between Samples Based on a Survival Tree Fit (Data Supplied as a Dataframe)
##' @description The function 
##' \code{predict_distance_tree} predicts distances between samples based on a survival tree fit.
##' @param survival_tree a fitted survival tree
##' @param numeric_predictor a formula specifying the numeric predictors. 
##' As in \code{~x1+x2+x3}, the three numeric variables \code{x1}, \code{x2}, and \code{x3} are included as numeric predictors. 
##' \code{x1[i]}, \code{x2[i]}, and \code{x3[i]} are the predictors of the ith sample.
##' The best practice is to use the same variables names in the training and testing dataset.
##' @param factor_predictor a formula specifying the numeric predictors. 
##' As in \code{~z1+z2+z3}, the three character variables \code{z1}, \code{z2}, and \code{z3} are included as factor predictors. 
##' \code{z1[i]}, \code{z2[i]}, and \code{z3[i]} are the predictors of the ith sample.
##' The best practice is to use the same variables names in the training and testing dataset.
##' @param data the dataframe (test data) that stores the outcome and predictor variables.
##' Variables in the global environment will be used if \code{data} is missing.
##' @param missing a character value that specifies the handling of missing data. 
##' If \code{missing=="omit"}, samples with missing values in the splitting variables will be discarded.
##' If \code{missing=="majority"}, samples with missing values in the splitting variables will be assigned to the majority node.
##' If \code{missing=="weighted"}, samples with missing values in the splitting variables will be weighted by the weights of branch nodes.
##' The best practice is to use the same method as the trained random tree.
predict_distance_tree<-function(
  survival_tree,
  numeric_predictor,
  factor_predictor,
  data,
  missing="omit"){
  
  if(!formula.tools::is.one.sided(numeric_predictor))stop("Invalid 'numeric_predictor' formula.")
  if(!formula.tools::is.one.sided(factor_predictor))stop("Invalid 'factor_predictor' formula.")
  
  if(missing(data)){
    mf_numeric_predictor<-eval(substitute(model.frame(numeric_predictor,na.action="na.pass")))
    mf_factor_predictor<-eval(substitute(model.frame(factor_predictor,na.action="na.pass")))
  }else{
    mf_numeric_predictor<-eval(substitute(model.frame(numeric_predictor,data=data,na.action="na.pass")))
    mf_factor_predictor<-eval(substitute(model.frame(factor_predictor,data=data,na.action="na.pass")))
  }
  
  # if(!all(sapply(mf_numeric_predictor[[1]],class)%in%c("integer","numeric")))stop("Invalid 'numeric_predictor' formula")
  # if(!all(sapply(mf_factor_predictor[[1]],class)%in%c("factor","character")))stop("Invalid 'factor_predictor' formula")
  if(length(mf_numeric_predictor)>=1&!all(sapply(mf_numeric_predictor,class)%in%c("integer","numeric")))stop("Invalid 'numeric_predictor' formula")
  if(length(mf_factor_predictor)>=1&!all(sapply(mf_factor_predictor,class)%in%c("factor","character")))stop("Invalid 'factor_predictor' formula")
  
  matrix_numeric<-as.matrix(mf_numeric_predictor)
  matrix_factor<-as.matrix(mf_factor_predictor)
  matrix_factor<-apply(matrix_factor,c(1,2),as.character)
  ndim_numeric<-ncol(matrix_numeric)
  ndim_factor<-ncol(matrix_factor)
  nind_test<-nrow(matrix_numeric)
  
  # check dimensions
  if(nrow(matrix_numeric)!=nrow(matrix_factor))stop("'nrow(matrix_numeric)' and 'nrow(matrix_factor) are different.'")
  if(ndim_numeric!=survival_tree$ndim_numeric)stop("'ncol(matrix_numeric)' inconsistent with training data.'")
  if(ndim_factor!=survival_tree$ndim_factor)stop("'ncol(matrix_factor)' inconsistent with training data.'")
  
  # clean [matrix_numeric] and [matrix_factor]
  factor_dictionary<-survival_tree$factor_dictionary
  matrix_factor_int<-matrix(NA,nind_test,ndim_factor)
  if(ncol(matrix_factor)>0){
    colnames(matrix_factor_int)<-colnames(matrix_factor)
    for(idx in 1:ncol(matrix_factor)){
      aname<-colnames(matrix_factor)[idx]
      matrix_factor_int[,idx]<-(factor_dictionary[[aname]])[matrix_factor[,aname]]
    }
  }
  matrix_factor<-matrix_factor_int
  
  a_distance<-calculate_distance(
    survival_tree$survival_tree,
    matrix_numeric,
    matrix_factor,
    missing)
  
  return(a_distance)
}
