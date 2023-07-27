##' Visualize the Fitted Survival Tree
##' @title Visualize the Fitted Survival Tree
##' @param survival_tree a fitted survival tree object.
##' @param cex numeric character expansion factor.
plot_survival_tree<-function(survival_tree,cex=0.75){
  num_to_str_g<-function(x,digits=4)sprintf(paste0("%.",digits,"g"),x)
  num_to_str_f<-function(x,digits=2)sprintf(paste0("%.",digits,"f"),x)
  num_to_pvalue<-function(x,digits=3){
    dplyr::case_when(
      x<0.00001~"p<0.00001",
      x<0.0001~"p<0.0001",
      x<0.001~"p<0.001",
      TRUE~paste0("p=",sprintf(paste0("%.",digits,"f"),x)))
  }
  
  a_table<-tree_to_table(survival_tree$survival_tree)
  a_table<-a_table[order(a_table$plot_order),]
  
  maxtime<-max(unlist(lapply(a_table$survival,function(x)x$time)))
  variable_names<-survival_tree$variable_names
  
  n_node<-nrow(a_table)
  n_layer<-max(a_table$layer,na.rm=T)
  
  xx_node<-1:n_node
  yy_node<-n_layer-ifelse(a_table$terminal,n_layer-1,a_table$layer-1)
  
  grid.newpage()
  vp<-viewport(
    x=unit(0.5,"npc"),y=unit(0.5,"npc"),
    width=unit(1,"npc")-unit(2,"points"),
    height=unit(1,"npc")-unit(2,"points"),
    xscale=c(0,n_node+1),yscale=c(0,n_layer),gp=gpar(cex=cex))
  pushViewport(vp)
  
  ## lines
  for(ii in 1:nrow(a_table)){
    if(is.na(a_table$left_id[ii]))next
    if(is.na(a_table$right_id[ii]))next
    xx_node_parent<-xx_node[ii]
    yy_node_parent<-yy_node[ii]
    xx_node_left<-xx_node[which(a_table$id==a_table$left_id[ii])]
    yy_node_left<-yy_node[which(a_table$id==a_table$left_id[ii])]
    xx_node_right<-xx_node[which(a_table$id==a_table$right_id[ii])]
    yy_node_right<-yy_node[which(a_table$id==a_table$right_id[ii])]
    # yy_node_middle<-max((yy_node_parent+yy_node_left)/2,(yy_node_parent+yy_node_right)/2)
    yy_node_middle<-yy_node[ii]-0.5
    grid.lines(
      x=unit(c(xx_node_parent,xx_node_parent,xx_node_left,xx_node_left),"native"),
      y=unit(c(yy_node_parent,yy_node_middle,yy_node_middle,yy_node_left)-0.5,"native"),
      gp=gpar(col="gray"))
    grid.lines(
      x=unit(c(xx_node_parent,xx_node_parent,xx_node_right,xx_node_right),"native"),
      y=unit(c(yy_node_parent,yy_node_middle,yy_node_middle,yy_node_right)-0.5,"native"),
      gp=gpar(col="gray"))
    if(a_table$type[ii]=="numeric"){
      split_text<-num_to_str_g(a_table$split_numeric[[ii]],3)
      left_text<-paste0("<=",split_text)
      right_text<-paste0(">",split_text)
    }else{
      split_text<-paste0(num_to_str_g(unlist(a_table$split_factor[[ii]]),3),collapse=",")
      left_text<-paste0(" is ",split_text)
      right_text<-paste0(" not ",split_text)
    }
    grid.text(
      left_text,
      x=unit((xx_node_parent+xx_node_left)/2,"native"),y=unit(yy_node_middle-0.5,"native"))
    grid.text(
      right_text,
      x=unit((xx_node_parent+xx_node_right)/2,"native"),y=unit(yy_node_middle-0.5,"native"))
  }
  
  ## nodes
  for(ii in 1:nrow(a_table)){
    
    ot_sub<-a_table$survival[[ii]]$time
    delta_sub<-a_table$survival[[ii]]$event
    weights_sub<-a_table$survival[[ii]]$weights
    a_survfit<-survfit(Surv(ot_sub,delta_sub)~1,weights=weights_sub)
    a_dostep<-dostep(a_survfit$time,a_survfit$surv)
    xscale<-c(0,maxtime)
    yscale<-c(0,1)
    vp<-viewport(
      x=unit(xx_node[ii],"native"),y=unit(yy_node[ii]-0.5,"native"),
      width=unit(2,"native"),height=unit(1,"native")-unit(2,"line"),
      layout=grid.layout(nrow=3,heights=unit(c(1,1,2),c("line","null","line"))))
    pushViewport(vp)
    pushViewport(viewport(
      xscale=xscale+c(-0.01,0.01)*(xscale[2]-xscale[1]),
      yscale=yscale+c(-0.01,0.01)*(yscale[2]-yscale[1]),
      layout.pos.row=2,layout.pos.col=1,clip="off"))
    if(a_table$terminal[ii]){
      grid.xaxis(gp=gpar(cex=0.5))
    }
    grid.rect(gp=gpar(fill="white",col="gray"))
    popViewport()
    pushViewport(viewport(
      xscale=xscale+c(-0.01,0.01)*(xscale[2]-xscale[1]),
      yscale=yscale+c(-0.01,0.01)*(yscale[2]-yscale[1]),
      layout.pos.row=2,layout.pos.col=1,clip="on"))
    grid.lines(
      x=unit(c(0,a_dostep$x),"native"),
      y=unit(c(1,a_dostep$y),"native"))
    popViewport()
    pushViewport(viewport(layout.pos.row=1,layout.pos.col=1))
    grid.text(paste0("Node ",a_table$id[ii]," (n=",num_to_str_f(a_table$w[ii],1),")"))
    popViewport()
    pushViewport(viewport(layout.pos.row=3,layout.pos.col=1))
    if(!a_table$terminal[ii]){
      text_var<-variable_names[a_table$j[ii]]
      text_z<-paste0("z=",num_to_str_f(a_table$z[ii],1),", ",num_to_pvalue(a_table$p[ii],3))
      grid.text(paste0(text_var,"\n",text_z),gp=gpar(lineheight=0.8))
    }
    popViewport()
    popViewport()
  }
}

