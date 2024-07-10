#include <tests.h>
#include <Rcpp/iostream/Rstreambuf.h>
#include <test_Toggle.h>
#include <test_TaperedModel.h>
#include <test_Stat.h>
#include <test_Constraint.h>
#include <test_BinaryNet.h>

namespace ernm{
namespace tests{

std::string testContext;
std::map< std::string, void(*)()> testFunctions;

void addTestFunction(std::string name, void(*test)() ){
    testFunctions.insert(std::make_pair(name, test));
}

void registerErnmTests(){
    addTestFunction("testBinaryNet", testBinaryNet);
    addTestFunction("testStats", testStats);
    addTestFunction("testConstraints", testConstraints);
    addTestFunction("testToggles", testToggles);
    addTestFunction("testTaperedModel", testTaperedModel);
}

void runErnmTests(){
#ifdef INSIDE
    Rcpp::Rcout << "\n\t";
#endif
    registerErnmTests();
    std::map< std::string, void(*)()>::iterator it;
    for ( it = testFunctions.begin(); it != testFunctions.end(); it++ ){
        testContext = it->first;
        it->second();
    }
#ifdef INSIDE
    Rcpp::Rcout << "All C++ Tests Complete\n";
#endif
}

}
}