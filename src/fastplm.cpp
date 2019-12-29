#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]  

using namespace Rcpp;

// [[Rcpp::export()]]
List fastplm(arma::mat data,
             arma::mat FE,
             arma::colvec weight,
             int FEcoefs = 0L 
){
  
  // parse data
  int n = data.n_rows;
  int k = data.n_cols;
  int m = FE.n_cols;
  int p = k-1; // No. of covariates
  arma::mat data_bak = data;
  arma::mat data_wei = data;
  arma::mat data_old = arma::zeros(n, k);
  double diff = 100;
  
  /*if(weight.is_empty()){
   weight = arma::ones(n, 1);
  }*/
  
  // declare variables
  arma::colvec y;  // n*1
  arma::mat X; // n*p
  arma::colvec resid; // n*1
  arma::colvec e; // n*1 (with fixed effects)
  arma::colvec coeff; // coefficient (full)
  //arma::colvec se; // SE (full)
  arma::colvec coef; // coefficient
  //arma::colvec stderror; // SE
  //int df; // degrees of freedom
  //double sig2; // sigma2
  double mu = 0; // grand mean
  arma::colvec LHS; // group means 
  arma::mat W; // big weighting matrix to calculate fixed effects
  arma::colvec alphas; // fixed effect coefficients
  arma::rowvec FEobs;
  arma::colvec FEindex;
  std::map<std::string,int> FEindex_map;
  
  /* count total number of groups (loss of degrees of freedom) */
  arma::mat FEvalues;  
  for(int ii=0; ii<m; ii++){ // cluster
    arma::colvec fe = arma::unique(FE.col(ii));   
    int g=fe.n_rows; // No. of group values
    arma::colvec gp(g);
    gp.fill(ii);
    arma::mat FEvalue = join_rows(gp, fe);
    FEvalues = join_cols(FEvalues, FEvalue);
  }
  
  int gtot = FEvalues.n_rows;
  
  
  /* FWL-MAP iteration */
  int niter = 0;
  while ((diff > 1e-5) & (niter<50)) { 
    
    // demean Y and X
    for(int ii=0; ii<m; ii++){ // cluster
      // weighted data
      for(int i=0; i<k; i++){
        data_wei.col(i)=data_wei.col(i)%weight.col(0);
      }
      
      arma::colvec fe=arma::unique(FE.col(ii));   
      int g=fe.n_rows; // No. of group values
      arma::mat mean_value=arma::zeros(g,(p+1)); // store Y, X means
      arma::colvec fe_length=arma::zeros(g,1);    // store number
      // sum up
      for(int j=0; j<g; j++){  // each group value
        for(int i=0; i<n; i++){   //row
          if(FE(i,ii)==fe(j)){
            fe_length(j)=fe_length(j)+weight(i,0);  // get sum and number 
            for(int t=0; t<(p+1); t++){  //variables             
              mean_value(j,t)=mean_value(j,t) + data_wei(i,t);
            } 
          }  
        }    
      } 
      // take average
      for(int i=0; i<(p+1); i++){
        for(int j=0; j<g; j++){
          if (fe_length(j)!=0) {
            mean_value(j,i)=mean_value(j,i)/fe_length(j);
          } else{
            mean_value(j,i)=0;
          } 
        }
      }
      // demean
      for(int j=0; j<n; j++){ // each observation
        for(int t=0; t<g; t++){ // each group varlue
          if(FE(j,ii)==fe(t)){ 
            for(int i=0; i<(p+1); i++){ // variable
              data(j,i) = data(j,i)-mean_value(t,i); 
            }
          }  
        }
      }
      data_wei = data;
    }
    // check convergence
    diff = arma::accu(abs(data - data_old));
    data_old = data;
    niter++;
  }
  
  
  /* Estimation */
  y = data.col(0);  // n*1
  if (p>0) {
    X = data.cols(1, p); // n*p
    //store coefficents and check variation of X
    coeff =arma::zeros(p,1); //store coefficients
    //se =arma::zeros(p,1);
    
    // check X variation
    int cc = p;
    for(int i=0;i<cc;i++){
      arma::colvec var =arma::unique(X.col(i));
      if(var.n_rows==1){ // no variation
        coeff(i)=arma::datum::nan ;
        //se(i)=arma::datum::nan ;
        X.shed_col(i);
        i=i-1;
        cc=cc-1;
      }
    }
    // weighted
    y = y%sqrt(weight);
    for(int i=0;i<cc;i++){
      X.col(i)=X.col(i)%sqrt(weight);
    } 
    
    // OLS
    coef =  solve(X, y);    // fit model y ~ X
    //arma::colvec coef = (X.t() * X ).i() * X.t() * y ;
    resid = y - X*coef;    // residuals
  }
  else {
    resid = y;
  }
  
  
  // std.err.
  // df = n - gtot - p;
  //sig2 = arma::as_scalar(resid.t()*resid/df);
  // if (p>0) { 
  //   stderror = arma::sqrt(sig2 * arma::diagvec(arma::inv(arma::trans(X)*X))); 
  // } 
  
  // fill in non-NaN coefficients
  int count=0;
  //String label = xname(0); 
  for(int i=0; i<p; i++){
    if(coeff(i)==0){
      //  if(coeff(i)==0||se(i)==0){
      coeff(i)=coef(count);
      //  se(i)=stderror(count); 
      count=count+1;
    }
  }
  
  // Calculate fixed effects coefficients
  if (FEcoefs == 1) {
    data = data_bak;
    y = data.col(0);  // n*1
    
    // grand mean
    mu = arma::mean(y);
    if (p > 0) {
      X = data.cols(1, p); // n*p
      coef = coeff; 
      for (int i=0; i<p; i++) {
        if (coef(i) == arma::datum::nan) {
          coef(i) = 0;
        }
      }
      mu = mu - arma::as_scalar(arma::mean(X, 0) * coef);
    }
    
    // residuals (with fixed effects)
    e = y - mu;
    if (p > 0) {
      e = e - X * coef;
    }
    
    arma::colvec LHS(gtot, arma::fill::zeros);
    arma::mat W(gtot, gtot, arma::fill::zeros);
    
    for(int i=0; i<gtot; i++){
      std::string tempFE = std::to_string(int(FEvalues(i,0)))+"."+std::to_string(int(FEvalues(i,1)));
      
      FEindex_map[tempFE] = i;
    }
    
    for(int i=0; i<n; i++) {
      arma::colvec FEindex(m,arma::fill::zeros);
      arma::rowvec FEobs = FE.row(i);
      for(int j=0; j<m; j++) {
        std::string tempFE = std::to_string(j)+"."+std::to_string(int(FEobs(j)));
        //std::cout << tempFE << std::endl;
        FEindex(j) = FEindex_map[tempFE];
      }
      for(int j=0;j<m;j++) {
        int index1 = FEindex(j);
        LHS(index1) = LHS(index1) + e(i);
        for(int k=0;k<m;k++) {
          int index2 = FEindex(k);
          W(index1, index2) = W(index1, index2) + 1;
        }
      }
    }
    alphas = arma::solve(W,LHS);
    FEvalues = join_rows(FEvalues, alphas);
  }
  
  
  // storage
  List output;
  if (p > 0) {
    output["coefficients"] = coeff;
    //output["stderr"] = se;
  }
  output["residuals"] = resid;
  output["niter"] = niter;
  output["FEvalues"] = FEvalues;
  if (FEcoefs == 1) {
    output["mu"] = mu ;
    output["ngroups"] = gtot; 
  }
  return(output); 
  
}
