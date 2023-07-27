#' survival_tree
#' time: (mandatory) vector of survival times
#' event: (mandatory) vector of event indicators
#' xx: (mandatory) data frame of covariates
#' significance: (default=0.05) significance level
#' minsize: (default=20) minimum number of samples in the terminal
#' missing: (default="majority") one of "majority" or "drop"
#' test: (default="cox") one of "cox" or "logrank"
#' details: (default="complete") one of "complete" or "rule"
#' types: (optional) vector of characters "numeric", "ordinal", "factor", specifies the types of each column
#' 
#' RETURN: object survival_tree

#' predict.survival_tree
#' survival_tree: (mandatory) a survival_tree object
#' xx: (mandatory) data frame of covariates
#' type: (default="distance") one of "label", "distance" or "survival"
#' 
#' RETURN: predicted labels

#' plot.survival_tree
#' survival_tree: a survival_tree object
#' 
#' RETURN: a figure

#' survival_forest
#' (same as survival tree)
#' 
#' RETURN: object survival_forest, a bunch of survival trees 

#' predict.survival_forest
#' survival_forest
#' xx: (mandatory) data frame of covariates
#' type: (default="distance") probably just distance
#' 
#' 
#' 
#' 

#####################
# pacakge framework #
#####################

#' outcome: time,event
#' predictor: matrix_numeric,matrix_factor
#' 
#' grow_tree(outcome,predictor){
#'   1~split
#'   2~result<-list(
#'     left_node=grow_tree(...),
#'     right_node=grow_tree(...))
#' }
#' survival_tree_matrix(outcome,predictor){
#'   1~clean
#'   2~create dictionary, matrix_factor->(int)
#'   3~a_survival_tree<-grow_tree(outcome,predictor)
#'   4~return(list(
#'     ndim_numeric,
#'     ndim_factor,
#'     dictionary
#'     a_survival_tree))
#' }
#' survival_tree(outcome,predictor){
#'   1~clean
#'   2~create dictionary, matrix_factor->(int)
#'   3~a_survival_tree<-grow_tree(outcome,predictor)
#'   4~return(list(
#'     ndim_numeric,
#'     ndim_factor,
#'     dictionary,
#'     a_survival_tree))
#' }
#' survival_forest_matrix(outcome,predictor){
#'   1~clean
#'   2~a_survival_forest<-list()
#'   for(boot in 1:nboot){
#'     3~outcome_boot,predictor_boot
#'     4~a_survival_forest[[boot]]<-survival_tree(outcome_boot,predictor_boot)
#'   }
#'   5~return(list(
#'     ndim_numeric,
#'     ndim_factor,
#'     a_survival_forest))
#' }
#' survival_forest(outcome,predictor){
#'   1~clean
#'   2~a_survival_forest<-list()
#'   for(boot in 1:nboot){
#'     3~outcome_boot,predictor_boot
#'     4~a_survival_forest[[boot]]<-survival_tree(outcome_boot,predictor_boot)
#'   }
#'   5~return(list(
#'     ndim_numeric,
#'     ndim_factor,
#'     a_survival_forest))
#' }
#' 
#' tree_to_table<-function(node){
#'   1~a_table<-grow_table(node)
#'   2~return(table)
#' }
#' 
#' grow_table<-function(node){
#'   1~return(rbind(
#'     current_node_info
#'     grow_table(node$left_node),
#'     grow_table(node$right_node)))
#' }
#' 
#' predict_weights<-function(survival_tree,predictor<raw>){
#'   1~clean
#'   2~survival_tree$dictionary, matrix_factor->(int)
#'   3~a_table<-tree_to_table(survival_tree$survival_tree)
#'   4~weights<-calculate_weights_by_table(a_table,predictor<int>)
#'   5~return(weights)
#' }
#' 
#' ## calculate_weights_by_table?
#' ## calculate_weights_by_node?
#' calculate_weights_by_table<-function(a_table,predictor<int>){
#'   1~grow_weights<-function(idx,weights){...}
#'   2~weights<-grow_weights()
#'   3~return(weights)
#' }
#' 
#' ## calculate_distance_by_table?
#' ## calculate_distance_by_node?
#' calculate_distance<-function(node,predictor<int>){
#'   1~a_table<-tree_to_table(node)
#'   2~a_predict_weights<-calculate_weights_by_table(a_table,predictor<int>)
#'   3~root_to_node<-function(node_id){...}
#'   4~path_from_node_to_node<-function(node_id1,node_id2){...}
#'   5~node_distance_0<-node_distance
#'   6~ind_distance<-a_predict_weights%*%node_distance_0%*%t(a_predict_weights)
#'   7~return(list(node_distance,ind_distance,ind_weights=a_predict_weights))
#' }
#' 
#' predict_distance_tree_matrix<-function(survival_tree,predictor<raw>){
#'   1~clean
#'   2~survival_tree$dictionary, matrix_factor->(int)
#'   3~a_distance<-calculate_distance(survival_tree,predictor<int>)
#'   4~return(a_distance)
#' }
#' 
#' predict_distance_tree<-function(survival_tree,predictor<raw>){
#'   1~clean
#'   2~survival_tree$dictionary, matrix_factor->(int)
#'   3~a_distance<-calculate_distance(survival_tree,predictor<int>)
#'   4~return(a+distance)
#' }
#' 
#' predict_distance_forest_matrix<-function(survival_forest,predictor<raw>){
#'   1~clean
#'   2~sum_distance<-matrix(0,nind_test,nind_test);sum_non_na<-matrix(0,nind_test,nind_test)
#'   for(boot in 1:nboot){
#'     3~a_distance<-predict_distance_tree_matrix(
#'       survival_forest[[boot_idx]],
#'       predictor)
#'     4~sum_non_na<-sum_non_na+!is.na(a_distance);sum_distance<-sum_distance+a_distance
#'   }
#'   5~mean_distance<-sum_distance/sum_non_na
#'   6~return(list(
#'     mean_distance=mean_distance,
#'     sum_distance=sum_distance,
#'     sum_non_na=sum_non_na))
#' }
#' 
#' 

