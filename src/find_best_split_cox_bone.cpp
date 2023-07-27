#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List find_best_split_cox_bone(
    NumericVector time, 
    LogicalVector event, 
    NumericMatrix xx_numeric,
    IntegerMatrix xx_factor,
    NumericVector weights,
    double min_weights=50.0) {
  
  const int nind=xx_numeric.nrow();
  const int ndim_numeric=xx_numeric.ncol();
  const int ndim_factor=xx_factor.ncol();
  NumericVector xx_numeric_column;
  IntegerVector xx_factor_column;
  
  NumericVector time_sub;
  LogicalVector event_sub;
  NumericVector weights_sub;
  NumericVector xx_numeric_sub;
  IntegerVector xx_factor_sub;
  
  LogicalVector bool_to_split(nind);
  LogicalVector bool_can_split(nind);
  LogicalVector to_unsure(nind);
  LogicalVector to_left(nind);
  LogicalVector to_right(nind);
  LogicalVector is_missing;
  
  NumericVector weights_left,weights_right;
  double sum_weights_left,sum_weights_right;
  
  NumericVector xx_numeric_unique;
  IntegerVector xx_factor_unique;
  
  int ii,ii2,jj,kk,nsplit,nsub,nfactor,splitidx;
  double zscore,a_split;
  double best_zscore=0.0,best_pvalue=1.0,best_chisq=0.0,best_split_numeric=0.0;
  int best_jj=-1;
  double expected_zscore,variance_zscore,denominator,numerator,fraction;
  LogicalMatrix factor_combinations;
  LogicalVector row_factor_combinations;
  IntegerVector factor_left;
  IntegerVector best_split_factor_left;
  
  NumericVector zscore_numeric(ndim_numeric);zscore_numeric.fill(0.0);
  NumericVector zscore_factor(ndim_factor);zscore_factor.fill(0.0);
  
  for(jj=0;jj<ndim_numeric;jj++){
    xx_numeric_column=xx_numeric(_,jj);
    is_missing=is_na(xx_numeric_column);
    
    bool_can_split=!is_missing;
    
    xx_numeric_sub=xx_numeric_column[bool_can_split];
    time_sub=time[bool_can_split];
    event_sub=event[bool_can_split];
    weights_sub=weights[bool_can_split];
    
    nsub=xx_numeric_sub.length();
    
    xx_numeric_unique=sort_unique(na_omit(xx_numeric_sub));
    nsplit=xx_numeric_unique.length();
    
    for(splitidx=0;splitidx<nsplit-1;splitidx++){
      a_split=xx_numeric_unique(splitidx);
      to_left=xx_numeric_sub<=a_split;
      to_right=xx_numeric_sub>a_split;
      weights_left=weights_sub[to_left];
      sum_weights_left=Rcpp::sum(weights_left);
      weights_right=weights_sub[to_right];
      sum_weights_right=Rcpp::sum(weights_right);
      if(sum_weights_left<min_weights)continue;
      if(sum_weights_right<min_weights)continue;
      expected_zscore=0.0;
      variance_zscore=0.0;
      for(ii=0;ii<nsub;ii++){
        if(!event_sub(ii))continue;
        if(to_right(ii))expected_zscore=expected_zscore+weights_sub(ii);
        denominator=0.0;
        numerator=0.0;
        for(ii2=0;ii2<nsub;ii2++){
          if(time_sub(ii2)<time_sub(ii))continue;
          denominator=denominator+weights_sub(ii2);
          if(to_right(ii2))numerator=numerator+weights_sub(ii2);
        }
        fraction=numerator/denominator;
        expected_zscore=expected_zscore-weights_sub(ii)*fraction;
        variance_zscore=variance_zscore+weights_sub(ii)*(fraction-fraction*fraction);
      }
      zscore=expected_zscore/std::sqrt(variance_zscore);
      if(std::abs(zscore)>std::abs(zscore_numeric(jj)))zscore_numeric(jj)=zscore;
      if(std::abs(zscore)>std::abs(best_zscore)){
        best_zscore=zscore;
        best_split_numeric=a_split;
        best_jj=jj;
      }
    }
  }
  
  for(jj=0;jj<ndim_factor;jj++){
    xx_factor_column=xx_factor(_,jj);
    is_missing=is_na(xx_factor_column);
    
    bool_can_split=!is_missing;
    
    xx_factor_sub=xx_factor_column[bool_can_split];
    time_sub=time[bool_can_split];
    event_sub=event[bool_can_split];
    weights_sub=weights[bool_can_split];
    
    nsub=xx_factor_sub.length();
    
    xx_factor_unique=Rcpp::sort_unique(Rcpp::na_omit(xx_factor_sub));
    nfactor=xx_factor_unique.length();
    if(nfactor<=1)continue;
    nsplit=std::pow(2,nfactor-1);
    factor_combinations=LogicalMatrix(nsplit,nfactor);
    for(ii=0;ii<nsplit;ii++){
      for(kk=0;kk<nfactor;kk++){
        factor_combinations(ii,kk)=(ii&(1<<kk))!=0;
      }
    }
    for(splitidx=1;splitidx<nsplit;splitidx++){
      
      row_factor_combinations=factor_combinations(splitidx,_);
      factor_left=xx_factor_unique[row_factor_combinations];
      to_left=LogicalVector(nsub);
      for(ii=0;ii<nsub;ii++){
        to_left(ii)=false;
        for(kk=0;kk<factor_left.length();kk++){
          if(factor_left(kk)==xx_factor_sub(ii)){
            to_left(ii)=true;
            break;
          }
        }
      }
      to_right=!to_left;
      
      weights_left=weights_sub[to_left];
      sum_weights_left=Rcpp::sum(weights_left);
      weights_right=weights_sub[to_right];
      sum_weights_right=Rcpp::sum(weights_right);
      if(sum_weights_left<min_weights)continue;
      if(sum_weights_right<min_weights)continue;
      expected_zscore=0.0;
      variance_zscore=0.0;
      for(ii=0;ii<nsub;ii++){
        if(!event_sub(ii))continue;
        if(to_right(ii))expected_zscore=expected_zscore+weights_sub(ii);
        denominator=0.0;
        numerator=0.0;
        for(ii2=0;ii2<nsub;ii2++){
          if(time_sub(ii2)<time_sub(ii))continue;
          denominator=denominator+weights_sub(ii2);
          if(to_right(ii2))numerator=numerator+weights_sub(ii2);
        }
        fraction=numerator/denominator;
        expected_zscore=expected_zscore-weights_sub(ii)*fraction;
        variance_zscore=variance_zscore+weights_sub(ii)*(fraction-fraction*fraction);
      }
      zscore=expected_zscore/std::sqrt(variance_zscore);
      if(std::abs(zscore)>std::abs(zscore_factor(jj)))zscore_factor(jj)=zscore;
      if(std::abs(zscore)>std::abs(best_zscore)){
        best_zscore=zscore;
        best_split_factor_left=Rcpp::clone(factor_left);
        best_jj=ndim_numeric+jj;
      }
    }
  }
  
  if(best_jj<0){
    
    best_split_numeric=NumericVector::get_na();
    best_split_factor_left=NumericVector::get_na();
    to_left=NumericVector::get_na();
    to_right=NumericVector::get_na();
    to_unsure=NumericVector::get_na();
    weights_left=NumericVector::get_na();
    weights_right=NumericVector::get_na();
    
  }else if(best_jj<ndim_numeric){
    xx_numeric_column=xx_numeric(_,best_jj);
    is_missing=is_na(xx_numeric_column);
    
    bool_to_split=LogicalVector(nind);bool_to_split.fill(false);
    bool_can_split=LogicalVector(nind);bool_can_split.fill(false);
    to_left=LogicalVector(nind);to_left.fill(false);
    to_right=LogicalVector(nind);to_right.fill(false);
    to_unsure=is_missing;
    
    for(ii=0;ii<nind;ii++){
      bool_to_split(ii)=true;
      if(to_unsure(ii)){
        continue;
      }else if(xx_numeric_column(ii)<=best_split_numeric){
        to_left(ii)=true;
        bool_can_split(ii)=true;
      }else{
        to_right(ii)=true;
        bool_can_split(ii)=true;
      }
    }
    
    weights_left=weights[to_left];
    sum_weights_left=Rcpp::sum(weights_left);
    weights_right=weights[to_right];
    sum_weights_right=Rcpp::sum(weights_right);
    
    best_chisq=best_zscore*best_zscore;
    best_pvalue=R::pchisq(best_chisq,1.0,false,false);
    
    best_split_factor_left=NumericVector::get_na();
  }else{
    xx_factor_column=xx_factor(_,best_jj-ndim_numeric);
    is_missing=is_na(xx_factor_column);
    
    bool_to_split=LogicalVector(nind);bool_to_split.fill(false);
    bool_can_split=LogicalVector(nind);bool_can_split.fill(false);
    to_left=LogicalVector(nind);to_left.fill(false);
    to_right=LogicalVector(nind);to_right.fill(false);
    to_unsure=is_missing;
    
    for(ii=0;ii<nind;ii++){
      bool_to_split(ii)=true;
      if(to_unsure(ii)){
        continue;
      }else{
        bool_can_split(ii)=true;
        to_left(ii)=false;
        for(kk=0;kk<best_split_factor_left.length();kk++){
          if(best_split_factor_left(kk)==xx_factor_column(ii)){
            to_left(ii)=true;
            break;
          }
        }
        to_right(ii)=!to_left(ii);
      }
    }
    
    weights_left=weights[to_left];
    sum_weights_left=Rcpp::sum(weights_left);
    weights_right=weights[to_right];
    sum_weights_right=Rcpp::sum(weights_right);
    
    best_chisq=best_zscore*best_zscore;
    best_pvalue=R::pchisq(best_chisq,1.0,false,false);
    
    best_split_numeric=NumericVector::get_na();
  }
  
  NumericVector all_zscore(ndim_numeric+ndim_factor);
  NumericVector all_chisq(ndim_numeric+ndim_factor);
  NumericVector all_pvalue(ndim_numeric+ndim_factor);
  for(jj=0;jj<ndim_numeric+ndim_factor;jj++){
    if(jj<ndim_numeric){
      all_zscore(jj)=zscore_numeric(jj);
    }else{
      all_zscore(jj)=zscore_factor(jj-ndim_numeric);
    }
    all_chisq(jj)=all_zscore(jj)*all_zscore(jj);
    all_pvalue(jj)=R::pchisq(all_chisq(jj),1.0,false,false);
  }
  
  List result=Rcpp::List::create(
    Rcpp::Named("best_zscore") = best_zscore,
    Rcpp::Named("best_chisq") = best_zscore*best_zscore,
    Rcpp::Named("best_pvalue") = best_pvalue,
    Rcpp::Named("best_split_numeric") = best_split_numeric,
    Rcpp::Named("best_split_factor_left") = best_split_factor_left,
    Rcpp::Named("best_jj") = best_jj+1,
    Rcpp::Named("to_left") = to_left,
    Rcpp::Named("to_right") = to_right,
    Rcpp::Named("to_unsure") = to_unsure,
    Rcpp::Named("sum_weights_left") = sum_weights_left,
    Rcpp::Named("sum_weights_right") = sum_weights_right,
    Rcpp::Named("all_zscore") = all_zscore,
    Rcpp::Named("all_pvalue") = all_pvalue);
  return(result);
}
