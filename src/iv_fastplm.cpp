#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]  

using namespace Rcpp;

// [[Rcpp::export()]]
List iv_fastplm(arma::mat Y, //outcome
                arma::mat X, //endogenous variable
                arma::mat Z, //included exogenous variables
                arma::mat IV, //excluded exogenous variables
                arma::mat FE, //fixed effects index
                arma::colvec weight,
                int FEcoefs = 0L 
){
  // parse data
  int n = Y.n_rows;
  int k_X = X.n_cols;
  int k_Z = Z.n_cols;
  int k_IV = IV.n_cols;
  arma::mat data =  join_rows(Y,X,Z,IV); //for demeaning
  int k = data.n_cols;
  int m = FE.n_cols;
  int p = k-1; 
  arma::mat data_bak = data;
  arma::mat data_wei = data;
  arma::mat data_old = arma::zeros(n, k);
  double diff = 100;

  // declare variables
  arma::colvec y;  // n*1
  arma::colvec resid; // n*1
  arma::colvec e; // n*1 (with fixed effects)
  arma::mat coeff; // coefficient (full)
  arma::mat coef; // coefficient
  double mu = 0; // grand mean
  arma::colvec LHS; // group means 
  arma::mat W; // big weighting matrix to calculate fixed effects
  arma::colvec alphas; // fixed effect coefficients
  arma::rowvec FEobs;
  arma::colvec FEindex;
  std::map<std::string,int> FEindex_map;
  

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

  //weights
  for(int i=0;i<k;i++){
      data.col(i)=data.col(i)%sqrt(weight);
  } 

  //recover Y X Z IV
  Y = data.col(0);  // n*1
  X = data.cols(1,k_X); //n*k_X
  if(k_Z>0){
    Z = data.cols(k_X+1,k_X+k_Z); //n*k_Z  
  }
  IV = data.cols(k_X+k_Z+1,k_X+k_Z+k_IV); //n*k_IV

  arma::mat Z_full = join_rows(Z,IV);
  arma::mat X_full = join_rows(X,Z);
  arma::mat inv_ZZ = inv_sympd(Z_full.t()*Z_full);
  arma::mat PzX = Z_full*inv_ZZ*(Z_full.t()*X_full);
  arma::mat PzY = Z_full*inv_ZZ*(Z_full.t()*Y);
  coeff = arma::solve(PzX,PzY);
  coef = coeff; 
  for (int i=0; i<k_X+k_Z; i++) {
    if (coef(i,0) == arma::datum::nan) {
        coef(i,0) = 0;
    }
  }
  arma::mat Y_hat = X_full*coef;
  arma::colvec residual = (Y - Y_hat).as_col();
  residual = residual/sqrt(weight);

  data = data_bak;
  y = data.col(0);  // n*1

  // grand mean
  mu = arma::sum(y%weight)/arma::sum(weight);
  X_full = data.cols(1, k_X+k_Z); //
  arma::mat X_wei = X_full;
  for(int i=0;i<k_X+k_Z;i++){
    X_wei.col(i)=X_wei.col(i)%weight;
  }
  X_wei = arma::sum(X_wei,0)/arma::sum(weight);
  mu = mu - arma::as_scalar(X_wei * coef); 

  // Calculate fixed effects coefficients
  if (FEcoefs == 1) {
    // residuals (with fixed effects)
    e = y - mu;
    if (p > 0) {
      e = e - X_full * coef;
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
      double sub_weight = arma::as_scalar(weight(i));
      for(int j=0; j<m; j++) {
        std::string tempFE = std::to_string(j)+"."+std::to_string(int(FEobs(j)));
        //std::cout << tempFE << std::endl;
        FEindex(j) = FEindex_map[tempFE];
      }
      for(int j=0;j<m;j++) {
        int index1 = FEindex(j);
        LHS(index1) = LHS(index1) + e(i)*sub_weight;
        for(int k=0;k<m;k++) {
          int index2 = FEindex(k);
          W(index1, index2) = W(index1, index2) + 1*sub_weight;
        }
      }
    }
    alphas = arma::solve(W,LHS);
    FEvalues = join_rows(FEvalues, alphas);
  }

  // storage
  List output;
  if (p > 0) {
    output["coefficients"] = coef;
    //output["stderr"] = se;
  }
  output["residuals"] = residual;
  output["niter"] = niter;
  output["mu"] = mu ;
  if (FEcoefs == 1) {
    output["ngroups"] = gtot; 
    output["FEvalues"] = FEvalues;
  }
  return(output); 

}
