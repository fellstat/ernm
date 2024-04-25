/*
 * tests.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: ianfellows
 */
#include "tests.h"
#include "BinaryNetTests.h"
#include "StatTests.h"
#include "ConstraintTests.h"
#include "ToggleTests.h"
#include "ReModelTests.h"
#include <Rcpp/iostream/Rstreambuf.h>

namespace ernm{
namespace tests{

std::string testContext;


RcppExport void runErnmTests(){
	Rcpp::Rcout << "\n\t";
	testBinaryNet();
	testStats();
	testConstraints();
	testToggles();
	testLatent();
	testReModel();
	Rcpp::Rcout << "All C++ Tests Complete\n";
}


}
}
