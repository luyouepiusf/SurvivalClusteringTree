##' Calculate Distances Between Samples (Internal Function)
##' (Works for coded factor matrices)
##' @noRd
calculate_distance<-function(
  node,
  matrix_numeric,
  matrix_factor,
  missing="omit"){
  
  a_table<-tree_to_table(node)
  a_predict_weights<-calculate_weights_by_table(a_table,matrix_numeric,matrix_factor,missing)
  
  root_to_node<-function(node_id){
    parent_node_id<-a_table$parent[node_id]
    if(is.na(parent_node_id)){
      return(node_id)
    }else{
      return(c(root_to_node(a_table$parent[node_id]),node_id))
    }
  }
  
  path_from_node_to_node<-function(node_id1,node_id2){
    path_root_to_node1<-root_to_node(node_id1)
    path_root_to_node2<-root_to_node(node_id2)
    
    minlength<-min(length(path_root_to_node1),length(path_root_to_node2))
    LCP_idx<-max(which(path_root_to_node1[1:minlength]==path_root_to_node2[1:minlength]))
    path_node1_to_node2<-
      c(rev(path_root_to_node1[LCP_idx:length(path_root_to_node1)][-1]),
        path_root_to_node1[LCP_idx],
        path_root_to_node2[LCP_idx:length(path_root_to_node2)][-1])
    
    path_zscore<-rep(NA,length(path_node1_to_node2))
    for(idx in 1:length(path_node1_to_node2)){
      temp_id<-path_node1_to_node2[idx]
      path_zscore[idx]<-a_table$z[temp_id]
    }
    return(list(
      path=path_node1_to_node2,
      path_zscore=path_zscore,
      path_total_zscore=sum(abs(path_zscore))))
  }
  
  node_distance<-matrix(NA,nrow(a_table),nrow(a_table))
  for(ii in 1:nrow(a_table)){
    if(!a_table$terminal[ii])next
    for(jj in 1:nrow(a_table)){
      if(!a_table$terminal[jj])next
      node_distance[ii,jj]<-path_from_node_to_node(ii,jj)$path_total_zscore
    }
  }
  
  node_distance_0<-node_distance
  node_distance_0[is.na(node_distance_0)]<-0
  ind_distance<-as.matrix(a_predict_weights%*%node_distance_0%*%t(a_predict_weights))
  
  if(missing=="omit"){
    non_terminal<-rowSums(a_predict_weights[,!a_table$terminal,drop=FALSE])>0
    ind_distance[non_terminal,]<-NA
    ind_distance[,non_terminal]<-NA
  }
  
  return(
    list(node_distance=node_distance,
         ind_distance=ind_distance,
         ind_weights=a_predict_weights))
}
