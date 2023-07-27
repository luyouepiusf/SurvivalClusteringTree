##' Predict Distances Between Samples Based on a Survival Tree Fit (Data Supplied as Matrices)
##' (Works for raw matrices)
##' @title Predict Distances Between Samples Based on a Survival Tree Fit (Data Supplied as Matrices)
##' @description The function 
##' \code{predict_distance_tree_matrix} predicts distances between samples based on a survival tree fit.
##' @param survival_tree a fitted survival tree
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
##' The best practice is to use the same method as the trained random tree.
predict_distance_tree_matrix<-function(
  survival_tree,
  matrix_numeric,
  matrix_factor,
  missing="omit"){
  
  # ndim_numeric<-survival_tree$ndim_numeric
  # ndim_factor<-survival_tree$ndim_factor
  # nind_test<-nrow(matrix_numeric)
  # 
  # # check dimensions
  # if(nrow(matrix_numeric)!=nrow(matrix_factor))stop("'nrow(matrix_numeric)' and 'nrow(matrix_factor) are different.'")
  # if(ncol(matrix_numeric)!=ndim_numeric)stop("'ncol(matrix_numeric)' inconsistent with training data.'")
  # if(ncol(matrix_factor)!=ndim_factor)stop("'ncol(matrix_factor)' inconsistent with training data.'")
  
  ndim_numeric<-ncol(matrix_numeric)
  ndim_factor<-ncol(matrix_factor)
  nind_test<-nrow(matrix_numeric)
  
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

