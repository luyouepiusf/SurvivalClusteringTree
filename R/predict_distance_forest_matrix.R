##' Predict Distances Between Samples Based on a Survival Forest Fit (Data Supplied as Matrices)
##' (Works for raw matrices)
##' @title Predict Distances Between Samples Based on a Survival Forest Fit (Data Supplied as Matrices)
##' @description The function 
##' \code{predict_distance_forest_matrix} predicts distances between samples based on a survival forest fit.
##' @param survival_forest a fitted survival forest
##' @param matrix_numeric numeric predictors, a numeric matrix. 
##' \code{matrix_numeric[i,j]} is the jth numeric predictor of the ith sample.
##' The best practice is to have the same column names in the training and testing dataset.
##' @param matrix_factor factor predictors, a character matrix. 
##' \code{matrix_factor[i,j]} is the jth predictor of the ith sample.
##' The best practice is to have the same column names in the training and testing dataset.
##' @param missing a character value that specifies the handling of missing data. 
##' If \code{missing=="omit"}, samples with missing values in the splitting variables will be discarded.
##' If \code{missing=="majority"}, samples with missing values in the splitting variables will be assigned to the majority node.
##' If \code{missing=="weighted"}, samples with missing values in the splitting variables will be weighted by the weights of branch nodes.
##' The best practice is to use the same method as the trained random forest.
predict_distance_forest_matrix<-function(
  survival_forest,
  matrix_numeric,
  matrix_factor,
  missing="omit"){
  
  ndim_numeric<-ncol(matrix_numeric)
  ndim_factor<-ncol(matrix_factor)
  nind_test<-nrow(matrix_numeric)
  
  # check dimensions
  if(nrow(matrix_numeric)!=nrow(matrix_factor))stop("'nrow(matrix_numeric)' and 'nrow(matrix_factor) are different.'")
  if(ndim_numeric!=survival_forest$ndim_numeric)stop("'ncol(matrix_numeric)' inconsistent with training data.'")
  if(ndim_factor!=survival_forest$ndim_factor)stop("'ncol(matrix_factor)' inconsistent with training data.'")
  
  # clean [matrix_numeric] and [matrix_factor]
  # factor_dictionary<-survival_forest$factor_dictionary
  # matrix_factor_int<-matrix(NA,nind_test,ndim_factor)
  # if(ncol(matrix_factor)>0){
  #   colnames(matrix_factor_int)<-colnames(matrix_factor)
  #   for(idx in 1:ncol(matrix_factor)){
  #     aname<-colnames(matrix_factor)[idx]
  #     matrix_factor_int[,idx]<-(factor_dictionary[[aname]])[matrix_factor[,aname]]
  #   }
  # }
  # matrix_factor<-matrix_factor_int
  
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
