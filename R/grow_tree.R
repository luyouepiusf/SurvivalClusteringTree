##' Grow a Tree From a Branch (Internal Function)
##' @noRd
grow_tree<-function(
  time,
  event,
  xx_numeric,
  xx_factor,
  weights,
  significance,
  min_weights,
  missing,
  test_type){
  
  xx_numeric<-data.matrix(xx_numeric)
  xx_factor<-data.matrix(xx_factor)
  
  result_find_best_split<-find_best_split_cox_bone(time,event,xx_numeric,xx_factor,weights,min_weights)
  result_find_best_split$sum_weights<-sum(weights)
  if(test_type%in%c("univariate","none")){
    result_find_best_split$best_pvalue_adjusted<-result_find_best_split$best_pvalue
  }else{
    all_pvalue<-result_find_best_split$all_pvalue
    all_pvalue[all_pvalue==1]<-NA
    all_pvalue_adjusted<-p.adjust(all_pvalue,method=test_type)
    result_find_best_split$best_pvalue_adjusted<-min(all_pvalue_adjusted)
  }
  
  if(result_find_best_split$best_pvalue_adjusted>significance|result_find_best_split$best_jj==0){
    result_find_best_split$best_zscore<-0
    result_find_best_split$best_chisq<-0
    result_find_best_split$best_pvalue<-1
    result<-list(
      left_node=NULL,
      right_node=NULL,
      terminal=TRUE,
      more_to_left=NA,
      time=time,
      event=event,
      weights=weights,
      ndim_numeric=ncol(xx_numeric),
      ndim_factor=ncol(xx_factor),
      split_info=result_find_best_split)
    return(result)
  }else{
    sum_weights_left<-result_find_best_split$sum_weights_left
    sum_weights_right<-result_find_best_split$sum_weights_right
    more_to_left<-sum_weights_left>=sum_weights_right
    
    if(missing=="majority"){
      left_weights<-dplyr::case_when(
        result_find_best_split$to_left~weights,
        result_find_best_split$to_right~0,
        result_find_best_split$to_unsure~ifelse(more_to_left,weights,0))
      right_weights<-dplyr::case_when(
        result_find_best_split$to_left~0,
        result_find_best_split$to_right~weights,
        result_find_best_split$to_unsure~ifelse(more_to_left,0,weights))
    }else if(missing=="omit"){
      left_weights<-dplyr::case_when(
        result_find_best_split$to_left~weights,
        result_find_best_split$to_right~0,
        result_find_best_split$to_unsure~0)
      right_weights<-dplyr::case_when(
        result_find_best_split$to_left~0,
        result_find_best_split$to_right~weights,
        result_find_best_split$to_unsure~0)
    }else if(missing=="weighted"){
      left_weights<-dplyr::case_when(
        result_find_best_split$to_left~weights,
        result_find_best_split$to_right~0,
        result_find_best_split$to_unsure~weights*sum_weights_left/(sum_weights_left+sum_weights_right))
      right_weights<-dplyr::case_when(
        result_find_best_split$to_left~0,
        result_find_best_split$to_right~weights,
        result_find_best_split$to_unsure~weights*sum_weights_right/(sum_weights_left+sum_weights_right))
    }else{
      stop("Wrong 'missing' argument")
    }
    
    left_idx<-left_weights>0
    right_idx<-right_weights>0
    
    time_left<-time[left_idx]
    event_left<-event[left_idx]
    weights_left<-left_weights[left_idx]
    xx_numeric_left<-xx_numeric[left_idx,]
    xx_factor_left<-xx_factor[left_idx,]
    
    time_right<-time[right_idx]
    event_right<-event[right_idx]
    weights_right<-right_weights[right_idx]
    xx_numeric_right<-xx_numeric[right_idx,]
    xx_factor_right<-xx_factor[right_idx,]
    
    result<-list(
      left_node=grow_tree(time_left,event_left,xx_numeric_left,xx_factor_left,weights_left,significance,min_weights,missing,test_type),
      right_node=grow_tree(time_right,event_right,xx_numeric_right,xx_factor_right,weights_right,significance,min_weights,missing,test_type),
      terminal=FALSE,
      more_to_left=more_to_left,
      time=time,
      event=event,
      weights=weights,
      ndim_numeric=ncol(xx_numeric),
      ndim_factor=ncol(xx_factor),
      split_info=result_find_best_split)
    return(result)
  }
}
