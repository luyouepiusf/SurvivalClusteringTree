##' Predict Distances Between Samples Based on a Survival Forest Fit (Data Supplied as a Dataframe)
##' 
##' @title Predict Distances Between Samples Based on a Survival Forest Fit (Data Supplied as a Dataframe)
##' @description The function 
##' \code{predict_distance_forest} predicts distances between samples based on a survival forest fit.
##' @param survival_forest a fitted survival forest
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
##' The best practice is to use the same method as the trained random forest.
predict_distance_forest<-function(
  survival_forest,
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
  if(ndim_numeric!=survival_forest$ndim_numeric)stop("'ncol(matrix_numeric)' inconsistent with training data.'")
  if(ndim_factor!=survival_forest$ndim_factor)stop("'ncol(matrix_factor)' inconsistent with training data.'")
  
  sum_distance<-matrix(0,nind_test,nind_test)
  sum_non_na<-matrix(0,nind_test,nind_test)
  
  for(boot_idx in 1:length(survival_forest$survival_forest)){
    a_distance<-predict_distance_tree_matrix(
      survival_forest$survival_forest[[boot_idx]],
      matrix_numeric,
      matrix_factor,
      missing=missing)$ind_distance
    sum_non_na<-sum_non_na+!is.na(a_distance)
    a_distance[is.na(a_distance)]<-0
    sum_distance<-sum_distance+a_distance
  }
  
  mean_distance<-sum_distance/sum_non_na
  return(list(
    mean_distance=mean_distance,
    sum_distance=sum_distance,
    sum_non_na=sum_non_na))
}
