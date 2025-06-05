// [[Rcpp::depends(ernm)]]
#include "ernm_hello_world.h"
#include "ernm.h"

//' An basic example of a function in C++ using ernm
//' @return a list of a character vector, a numeric vector, and an ernm DirectedNet
//' @param NULL
//' @examples
//' ernm_hello_world()
// [[Rcpp::export]]
Rcpp::List ernm_hello_world(){
    using namespace Rcpp;
    using namespace ernm;
    
    IntegerMatrix tmp(0,2);
    DirectedNet net(tmp,20); 
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y, wrap(net) ) ;
    
    return z ;
}
