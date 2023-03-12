// This script implement generalized factor model based on variational inference.
// Date: 2023-02-11

// Revised log:
// 2023-01-16: replace the for loop with matrix operation in updating variational parameters of Z of E-step!

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 
// diag(W0^t* Cki * W0)
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  




double calELBO(const field<mat>& Xf, const vec& typeID, const sp_mat& A,  const mat& Mu_y, const mat& S_y1,
               const vec& invLambda, const mat& B, const rowvec& mu, const mat& Mu_h,
               const mat& S_h, const mat& Sigma_h){
  
  int i, p, d = Xf.n_elem, n = Xf(0).n_rows, q= B.n_cols;
  double  logS, entropy=0.0, ELBO=0.0;
  double term1=0.0, term2=0.0;
  mat S_y(S_y1);
  vec p_vec(d+1, fill::zeros);
  for(i=0; i<d; ++i){
    p = Xf(i).n_cols;
    p_vec(i+1) = p_vec(i) + p;
    
    if(typeID(i)==2){
      mat tmp_mat1 = Xf(i) % Mu_y.cols(p_vec(i),p_vec(i+1)-1);
      tmp_mat1.replace(datum::nan, 0);
      mat tmp_mat2 = exp(Mu_y.cols(p_vec(i),p_vec(i+1)-1) + S_y.cols(p_vec(i),p_vec(i+1)-1)/2);
      tmp_mat2.replace(datum::nan, 0);
      term1 += accu(tmp_mat1); 
      term2 += -accu(tmp_mat2);
    }
    if(typeID(i)==3){
      rowvec nvec = max(Xf(i));
      mat tmp_mat1 = (Xf(i)-repmat(nvec, n, 1)) % Mu_y.cols(p_vec(i),p_vec(i+1)-1);
      tmp_mat1.replace(datum::nan, 0);
      mat tmp_mat2 = repmat(nvec, n, 1) %  exp(-Mu_y.cols(p_vec(i),p_vec(i+1)-1) + S_y.cols(p_vec(i),p_vec(i+1)-1)/2);
      tmp_mat2.replace(datum::nan, 0);
      term1 += accu(tmp_mat1);  // 
      term2 += -accu(tmp_mat2); // term2 is not written now!
    }
    
  }
  
  
  mat S_bar = n* S_h;
  logS = n* log_det_sympd(S_h);
  
  // if typeID = 1: Mu_y submatrix == x_1
  mat dX = (Mu_y - Mu_h * B.t() - repmat(mu, n, 1) - A) % repmat(sqrt(invLambda.t()), n, 1);
  dX.replace(datum::nan, 0);
  S_y.replace(datum::nan, 0);
  mat LamB = B % repmat(sqrt(invLambda), 1, q); 
  vec rowSum_BSB = decomp(S_bar, LamB);
  double dimreduction_term1 = -0.5*(accu(dX % dX)+
                                    accu(S_y % repmat(invLambda.t(), n, 1))+ accu(rowSum_BSB) - n* accu(log(invLambda)) );
  
  vec rosSum_MuSMu = decomp( inv_sympd(Sigma_h), Mu_h); // trace(Mu_h * inv_sympd(Sigma_h) * Mu_h.t())
  double dimreduction_term2 = -0.5*(n* log_det_sympd(Sigma_h) +  accu(rosSum_MuSMu)+
                                    trace(S_bar * inv_sympd(Sigma_h) ));
  
  
  entropy = 0.5*accu(log(S_y+ (S_y==0))) + 0.5 * logS;
  
  ELBO = term1 + term2 + dimreduction_term1 + dimreduction_term2 + entropy;
  return ELBO;
}


void VB_Estep(const field<mat>& Xf, const sp_mat& A, const vec& typeID, mat& Mu_y,  mat& S_y,
              const vec& invLambda, const mat& B, const rowvec& mu, mat& Mu_h,
              mat& S_h, const mat& Sigma_h){
  
  // Basic information:
  int i,  d = Xf.n_elem, n=Xf(0).n_rows, q= B.n_cols, p;
  vec p_vec(d+1, fill::zeros);
  // Mu_y is a n*p matrix
  mat tZ = Mu_h * B.t() + repmat(mu, n, 1) + A; // missing value, Mu_h still has normal value
  mat invLamMat = repmat(invLambda.t(), n, 1);
  
  //  ## VB E-step
  // update posterior variance of y: S_y
  // update posterior mean of y: Mu_y
  // double elbo0 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
  
  for(i=0; i<d; ++i){ // update Mu_y and S_y by the variable types, note the missing values.
      
      p = Xf(i).n_cols;
      p_vec(i+1) = p_vec(i) + p;
      //Rprintf("good 3, i= %d!\n", i);
      if(typeID(i)==2){ // Poisson
        mat Mu_y1 = Mu_y.cols(p_vec(i), p_vec(i+1)-1); // Mu_y has missing location filled with nan.
        Mu_y1 = (Xf(i) -  exp(Mu_y1) % (1-Mu_y1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1)% tZ.cols(p_vec(i), p_vec(i+1)-1)) / 
          (exp(Mu_y1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );
        Mu_y.cols(p_vec(i), p_vec(i+1)-1) = Mu_y1; // has nan values.
        S_y.cols(p_vec(i), p_vec(i+1)-1) = 1.0 / (exp(Mu_y1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );
       
      }
      
      //Rprintf("good 4, i= %d!\n", i);
      
      if(typeID(i)==3){ // Binomial
        rowvec nvec = max(Xf(i));
        mat Mu_y2 = Mu_y.cols(p_vec(i), p_vec(i+1)-1); // take out this submatrix
        Mu_y2 = (Xf(i)- (1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1)% tZ.cols(p_vec(i), p_vec(i+1)-1)) / 
          ((1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );
        Mu_y.cols(p_vec(i), p_vec(i+1)-1) = Mu_y2;
        S_y.cols(p_vec(i), p_vec(i+1)-1) = 1.0 / ( (1/(1+exp(-Mu_y2))) %(1- 1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );
        
        
      }
      // Rprintf("good 5, i= %d!\n", i);
      
  }
  
  mat tX = Mu_y - repmat(mu, n, 1) -  A; // has nan values.
  
  // double elbo1 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
  // Rprintf("VB dY= %4f\n", elbo1 - elbo0);
  tX.replace(datum::nan, 0); // replace nan with zero for tX.
  // ## update posterior variance of h_i, S_i
  // ## update posterior mean of h_i, Mu_h
  mat Si_inv;
  Si_inv = B.t() *  (B % repmat(invLambda, 1, q)) + inv_sympd(Sigma_h); 
  S_h = inv_sympd(Si_inv);
  Mu_h = tX * (B % repmat(invLambda, 1, q)) * S_h;
  //Rprintf("good 6!\n");
  // double elbo2 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
  // Rprintf("VB dH= %4f\n", elbo2 - elbo1);
}

rowvec mean_comp(const mat& Mu_y){
  
  int j, p = Mu_y.n_cols;
  rowvec mean_Mu(p, fill::zeros);
  vec tmpvec;
  for(j=0; j<p; ++j){
    tmpvec = Mu_y.col(j);
    uvec idx =  find_finite(tmpvec);
    mean_Mu(j) = accu(tmpvec(idx))/idx.n_elem;
    
  }
  return mean_Mu;
}

mat update_B(const mat& Mu_A, const mat& Mu_h, const mat& S_h){
  
  int j, p = Mu_A.n_cols, q= Mu_h.n_cols;
  mat B(p, q, fill::zeros);
  for(j=0; j<p; ++j){
    uvec idx = find_finite(Mu_A.col(j));
    B.row(j) = Mu_A.col(j).t()* Mu_h * inv_sympd(Mu_h.rows(idx).t()* Mu_h.rows(idx) + S_h);
  }
  
  return B;
}

rowvec get_number_obs(const mat& Mu_A){
  
  int j, p = Mu_A.n_cols;
  rowvec num_obs_vec(p, fill::zeros);
  for(j=0; j<p; ++j){
    uvec idx = find_finite(Mu_A.col(j));
    num_obs_vec(j) = idx.n_elem;
  }
  
  return num_obs_vec;
}



// [[Rcpp::export]]
Rcpp::List VB_OrMIGcpp(const Rcpp::List& XList, const arma::vec& typeID, const arma::sp_mat& A, 
                             const arma::mat& Mu_y_int, 
                             const arma::mat& S_y_int,
                             const arma::vec& invLambda_int, const arma::mat& B_int, 
                             const arma::rowvec& mu_int, const arma::mat& Mu_h_int,
                             const arma::mat& S_h_int, const arma::mat& Sigma_h_int, 
                             const double& epsELBO, const int& maxIter, const bool& verbose){
  // typeID: 1 means Gaussian; 2 means Poisson; 3 means Binomial;
  // whether add nvec for Binomial variables???? Now assume n_j = 1 for all.
  int i, d = XList.length();
  int d2 = typeID.n_elem;
  if(d != d2){
    stop("The length of XList must be equal to the length of typeID!");
  }
  
  field<mat> Xmisf(d);
  for(i=0; i<d; ++i){ // put the component of list into a field.
    mat Xtmp = XList[i];
    Xmisf(i) = Xtmp;
    
  }
  int n = Xmisf(0).n_rows; //  p = Mu_y_int.n_cols; //  find_finite(Xf(0).col(j));
  
  
  // Initialize: Mu_y, S_y include nan values.
  mat Mu_y(Mu_y_int), S_y(S_y_int),  Mu_h(Mu_h_int), S_h(S_h_int), B(B_int), Sigma_h(Sigma_h_int);
  vec invLambda(invLambda_int);
  rowvec mu(mu_int); // 
  
  vec ELBO_vec(maxIter), bsb, Lambda;
  ELBO_vec(0) = INT_MIN;
  mat S_bar, dX, Mu_A;
  int iter;
  
  
  for(iter = 1; iter < maxIter; ++iter){
    
    
    // Rprintf("E step starting!\n");
    // VB E-step
    VB_Estep(Xmisf, A, typeID, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    
    // Rprintf("Finish E step!\n");
    //VB M-step
    S_bar = n* S_h;
    
    // double elbo1 = calELBO(Xmisf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    //update B
    // Rprintf("update B\n");
    //Mu_y.replace(datum::nan, 0); // replace the nan with zeros.
    Mu_A = Mu_y-repmat(mu, n, 1)-A;
    Mu_A.replace(datum::nan, 0);
    // num_obs_vec = get_number_obs(Mu_y);
    // B = trans(Mu_A / repmat(num_obs_vec, n, 1))  * Mu_h * inv_sympd(Mu_h.t() * Mu_h/n + S_h);
    B = update_B( Mu_A, Mu_h, S_bar);
    // double elbo2 = calELBO(Xmisf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    // Rprintf("dB= %4f \n", elbo2 - elbo1);
    
    
    // update mu
    // Rprintf("update mu\n");
    mu = mean_comp(Mu_y - Mu_h * B.t()-A); // Mu_y for gaussian: x_1
    // double elbo3 =  calELBO(Xmisf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    // Rprintf("dmu= %4f \n", elbo3 - elbo2);
    
    
    // update Sigma_h
    // Rprintf("update Sigma_h\n");
    Sigma_h = (Mu_h.t() * Mu_h + S_bar)/n;
    // double elbo5 =  calELBO(Xmisf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    // Rprintf("dSigma= %4f\n", elbo5 - elbo3);
    
    // update Lambda
    // Rprintf("update Lambda\n");
    dX = Mu_y -  Mu_h * B.t() - repmat(mu, n, 1)-A;
    // vec bsb = diagvec(B * S_bar * B.t()); // this can be speeded up by using decomp() function.
    bsb = decomp(S_bar, B);
    Lambda =  trans(mean_comp(dX % dX + S_y)) + bsb/n; // typeID=1: S_y submatix must be zeros.
    invLambda = 1.0 / Lambda;
    // double elbo6 =   calELBO(Xmisf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    // Rprintf("dLambda= %4f\n", elbo6 - elbo5);
    
    
    ELBO_vec(iter) =  calELBO(Xmisf, typeID, A, Mu_y, S_y, invLambda, B, mu, Mu_h, S_h, Sigma_h);
    
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n", 
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }
  
  // output return value
  List resList = List::create(
    Rcpp::Named("H") = Mu_h,
    Rcpp::Named("B") = B,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("invLambda") = invLambda,
    Rcpp::Named("Sigma") = Sigma_h,
    Rcpp::Named("Mu_y") = Mu_y,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);
  
}