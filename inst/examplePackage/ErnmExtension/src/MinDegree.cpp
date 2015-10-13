// [[Rcpp::depends(ernm)]]
#include "MinDegree.h"
#include "ernm.h"


/**
* This function registers the new MinDegree statistic so that
* it can be used in ernm formula.
*
* RcppExport means this function can be called from R using
* .C("registerMinDegree")
* see: .onLoad in zzz.R
*/
RcppExport void registerMinDegree(){
	using namespace ernm;
	Rcpp::XPtr< Stat<Undirected> > ps1(new MinDegree<Undirected>());
	REGISTER_UNDIRECTED_STATISTIC(ps1);
}
