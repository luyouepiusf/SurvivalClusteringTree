num_to_str_g<-function(x,digits=4)sprintf(paste0("%.",digits,"g"),x)
num_to_str_f<-function(x,digits=2)sprintf(paste0("%.",digits,"f"),x)
num_to_pvalue<-function(x,digits=3){
  dplyr::case_when(
    x<0.00001~"p<0.00001",
    x<0.0001~"p<0.0001",
    x<0.001~"p<0.001",
    TRUE~paste0("p=",sprintf(paste0("%.",digits,"f"),x)))
}
create_dictionary<-function(x,na.rm=T){
  x<-sort(na.omit(unique(x)))
  structure(1:length(x),names=x)
}
# dostep <- function(x, y) {
#   ### create a step function based on x, y coordinates
#   ### modified from `survival:print.survfit'
#   if (is.na(x[1] + y[1])) {
#     x <- x[-1]
#     y <- y[-1]
#   }
#   n <- length(x)
#   if (n > 2) {  
#     # replace verbose horizonal sequences like
#     # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
#     # with (1, .2), (3, .1).  They are slow, and can smear the looks
#     # of the line type.
#     dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
#     n2 <- sum(dupy)
#     #create a step function
#     xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
#     yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
#     RET <- list(x = xrep, y = yrep)
#   } else {
#     if (n == 1) {
#       RET <- list(x = x, y = y)
#     } else {
#       RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
#     }
#   }
#   return(RET)
# }
dostep<-function(x,y){
  n<-length(x)
  xrep<-rep(x,c(1,rep(2,n-1)))
  yrep<-rep(y,c(rep(2,n-1),1))
  result<-list(x=xrep,y=yrep)
  return(result)
}
