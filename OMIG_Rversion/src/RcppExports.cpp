// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// VB_OrMIGcpp
Rcpp::List VB_OrMIGcpp(const Rcpp::List& XList, const arma::vec& typeID, const arma::sp_mat& A, const arma::mat& Mu_y_int, const arma::mat& S_y_int, const arma::vec& invLambda_int, const arma::mat& B_int, const arma::rowvec& mu_int, const arma::mat& Mu_h_int, const arma::mat& S_h_int, const arma::mat& Sigma_h_int, const double& epsELBO, const int& maxIter, const bool& verbose);
RcppExport SEXP _OMIG_VB_OrMIGcpp(SEXP XListSEXP, SEXP typeIDSEXP, SEXP ASEXP, SEXP Mu_y_intSEXP, SEXP S_y_intSEXP, SEXP invLambda_intSEXP, SEXP B_intSEXP, SEXP mu_intSEXP, SEXP Mu_h_intSEXP, SEXP S_h_intSEXP, SEXP Sigma_h_intSEXP, SEXP epsELBOSEXP, SEXP maxIterSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type XList(XListSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type typeID(typeIDSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Mu_y_int(Mu_y_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S_y_int(S_y_intSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type invLambda_int(invLambda_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_int(B_intSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type mu_int(mu_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Mu_h_int(Mu_h_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S_h_int(S_h_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma_h_int(Sigma_h_intSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsELBO(epsELBOSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(VB_OrMIGcpp(XList, typeID, A, Mu_y_int, S_y_int, invLambda_int, B_int, mu_int, Mu_h_int, S_h_int, Sigma_h_int, epsELBO, maxIter, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OMIG_VB_OrMIGcpp", (DL_FUNC) &_OMIG_VB_OrMIGcpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_OMIG(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
