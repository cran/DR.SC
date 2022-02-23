// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getneighborhood_fast
arma::sp_umat getneighborhood_fast(const arma::mat x, double radius);
RcppExport SEXP _DR_SC_getneighborhood_fast(SEXP xSEXP, SEXP radiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    rcpp_result_gen = Rcpp::wrap(getneighborhood_fast(x, radius));
    return rcpp_result_gen;
END_RCPP
}
// sp_means_Rcpp
arma::vec sp_means_Rcpp(arma::sp_mat sp_data, bool rowMeans);
RcppExport SEXP _DR_SC_sp_means_Rcpp(SEXP sp_dataSEXP, SEXP rowMeansSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type sp_data(sp_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type rowMeans(rowMeansSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_means_Rcpp(sp_data, rowMeans));
    return rcpp_result_gen;
END_RCPP
}
// sp_sums_Rcpp
arma::vec sp_sums_Rcpp(arma::sp_mat sp_data, bool rowSums);
RcppExport SEXP _DR_SC_sp_sums_Rcpp(SEXP sp_dataSEXP, SEXP rowSumsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type sp_data(sp_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type rowSums(rowSumsSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_sums_Rcpp(sp_data, rowSums));
    return rcpp_result_gen;
END_RCPP
}
// getPairDist
arma::mat getPairDist(const arma::mat x);
RcppExport SEXP _DR_SC_getPairDist(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(getPairDist(x));
    return rcpp_result_gen;
END_RCPP
}
// calYenergy2D_sp
arma::mat calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta);
RcppExport SEXP _DR_SC_calYenergy2D_sp(SEXP ySEXP, SEXP AdjSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(calYenergy2D_sp(y, Adj, K, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// obj_beta
double obj_beta(const arma::ivec& y, const arma::mat& R, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta);
RcppExport SEXP _DR_SC_obj_beta(SEXP ySEXP, SEXP RSEXP, SEXP AdjSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_beta(y, R, Adj, K, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// icmem_heterCpp
Rcpp:: List icmem_heterCpp(const arma::mat& X, const arma::sp_mat& Adj, const arma::imat& y_int, Rcpp::List& Mu_intList, const arma::mat& W_int, Rcpp::List& Sigma_intList, arma::vec& Lam_vec_int, Rcpp::List& alphaList, const arma::vec& beta_int, const arma::vec& beta_grid, const int& maxIter_ICM, const int& maxIter, const double& epsLogLik, const int& verbose, const bool& homo, const bool& diagSigmak, const int maxK, const int minK, const int coreNum);
RcppExport SEXP _DR_SC_icmem_heterCpp(SEXP XSEXP, SEXP AdjSEXP, SEXP y_intSEXP, SEXP Mu_intListSEXP, SEXP W_intSEXP, SEXP Sigma_intListSEXP, SEXP Lam_vec_intSEXP, SEXP alphaListSEXP, SEXP beta_intSEXP, SEXP beta_gridSEXP, SEXP maxIter_ICMSEXP, SEXP maxIterSEXP, SEXP epsLogLikSEXP, SEXP verboseSEXP, SEXP homoSEXP, SEXP diagSigmakSEXP, SEXP maxKSEXP, SEXP minKSEXP, SEXP coreNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type y_int(y_intSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Mu_intList(Mu_intListSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W_int(W_intSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Sigma_intList(Sigma_intListSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Lam_vec_int(Lam_vec_intSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type alphaList(alphaListSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_int(beta_intSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_grid(beta_gridSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter_ICM(maxIter_ICMSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsLogLik(epsLogLikSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool& >::type homo(homoSEXP);
    Rcpp::traits::input_parameter< const bool& >::type diagSigmak(diagSigmakSEXP);
    Rcpp::traits::input_parameter< const int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< const int >::type minK(minKSEXP);
    Rcpp::traits::input_parameter< const int >::type coreNum(coreNumSEXP);
    rcpp_result_gen = Rcpp::wrap(icmem_heterCpp(X, Adj, y_int, Mu_intList, W_int, Sigma_intList, Lam_vec_int, alphaList, beta_int, beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose, homo, diagSigmak, maxK, minK, coreNum));
    return rcpp_result_gen;
END_RCPP
}
// EMmPCpp_heter
Rcpp:: List EMmPCpp_heter(const arma::mat& X, Rcpp::List& Pi_int, Rcpp::List& Mu_int, const arma::mat& W_int, Rcpp::List& Sigma_int, arma::vec& Lam_vec_int, const int& maxIter, const double& epsLogLik, const bool& verbose, const bool& homo, const bool& diagSigmak, const int& maxK, const int& minK, const int& coreNum);
RcppExport SEXP _DR_SC_EMmPCpp_heter(SEXP XSEXP, SEXP Pi_intSEXP, SEXP Mu_intSEXP, SEXP W_intSEXP, SEXP Sigma_intSEXP, SEXP Lam_vec_intSEXP, SEXP maxIterSEXP, SEXP epsLogLikSEXP, SEXP verboseSEXP, SEXP homoSEXP, SEXP diagSigmakSEXP, SEXP maxKSEXP, SEXP minKSEXP, SEXP coreNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Pi_int(Pi_intSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Mu_int(Mu_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W_int(W_intSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Sigma_int(Sigma_intSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Lam_vec_int(Lam_vec_intSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsLogLik(epsLogLikSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool& >::type homo(homoSEXP);
    Rcpp::traits::input_parameter< const bool& >::type diagSigmak(diagSigmakSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< const int& >::type minK(minKSEXP);
    Rcpp::traits::input_parameter< const int& >::type coreNum(coreNumSEXP);
    rcpp_result_gen = Rcpp::wrap(EMmPCpp_heter(X, Pi_int, Mu_int, W_int, Sigma_int, Lam_vec_int, maxIter, epsLogLik, verbose, homo, diagSigmak, maxK, minK, coreNum));
    return rcpp_result_gen;
END_RCPP
}
// calculateWeight
arma::vec calculateWeight(const arma::mat& X, const int& nPCs);
RcppExport SEXP _DR_SC_calculateWeight(SEXP XSEXP, SEXP nPCsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int& >::type nPCs(nPCsSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateWeight(X, nPCs));
    return rcpp_result_gen;
END_RCPP
}
// wpcaCpp
Rcpp::List wpcaCpp(const arma::mat& X, const int& nPCs, const bool& weighted);
RcppExport SEXP _DR_SC_wpcaCpp(SEXP XSEXP, SEXP nPCsSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int& >::type nPCs(nPCsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(wpcaCpp(X, nPCs, weighted));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DR_SC_getneighborhood_fast", (DL_FUNC) &_DR_SC_getneighborhood_fast, 2},
    {"_DR_SC_sp_means_Rcpp", (DL_FUNC) &_DR_SC_sp_means_Rcpp, 2},
    {"_DR_SC_sp_sums_Rcpp", (DL_FUNC) &_DR_SC_sp_sums_Rcpp, 2},
    {"_DR_SC_getPairDist", (DL_FUNC) &_DR_SC_getPairDist, 1},
    {"_DR_SC_calYenergy2D_sp", (DL_FUNC) &_DR_SC_calYenergy2D_sp, 5},
    {"_DR_SC_obj_beta", (DL_FUNC) &_DR_SC_obj_beta, 6},
    {"_DR_SC_icmem_heterCpp", (DL_FUNC) &_DR_SC_icmem_heterCpp, 19},
    {"_DR_SC_EMmPCpp_heter", (DL_FUNC) &_DR_SC_EMmPCpp_heter, 14},
    {"_DR_SC_calculateWeight", (DL_FUNC) &_DR_SC_calculateWeight, 2},
    {"_DR_SC_wpcaCpp", (DL_FUNC) &_DR_SC_wpcaCpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_DR_SC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
