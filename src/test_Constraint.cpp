/*
 * ContraintTests.cpp
 *
 *  Created on: Jan 9, 2014
 *      Author: ianfellows
 */


#include <Rcpp.h>
#include <BinaryNet.h>
#include <Stat.h>
#include <Stats.h>
#include <Offsets.h>
#include <Constraint.h>
#include <Constraints.h>
#include <Model.h>
#include <DyadToggles.h>
#include <MetropolisHastings.h>
#include <VarAttrib.h>
#include <VertexToggles.h>
#include <tests.h>
#include <test_Constraint.h>
namespace ernm {

namespace tests{


template<class Engine>
void testBoundedDegree(){
	GetRNGstate();
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,30);
	boost::shared_ptr< AbstractStat<Engine> > ed(new Stat<Engine, Edges<Engine> >());
    Rcpp::List ll;
    ll.push_back(2);
    ll.push_back(10);
    boost::shared_ptr< AbstractOffset<Engine> > off(
        		new Constraint<Engine,BoundedDegree<Engine> >(ll));
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addOffsetPtr(off);
    model.calculate();
    model.setThetas(std::vector<double>(1,0));

    EXPECT_TRUE(model.offset().at(0) < -100000)

    MetropolisHastings<Engine> mh(model);
    mh.initialize();
    mh.run(4000);
    EXPECT_NEAR(mh.getModel()->offset().at(0),0.0)
    for(int i=0;i<net.size();i++){
    	int deg = mh.getModel()->network()->degree(i);
    	EXPECT_TRUE(deg<=10 && deg>=2);
    }
    PutRNGstate();
}

template<class Engine>
void testREffect(){
	GetRNGstate();
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,30);
	std::vector<double>vals2;
	for(int i=0; i<30; i++)
		vals2.push_back(Rf_runif(1,29));
	ContinAttrib attr2;
	attr2.setName("var");
	attr2.setLowerBound(1);
	attr2.setUpperBound(29);
	net.addContinVariable(vals2,attr2);


	boost::shared_ptr< AbstractStat<Engine> > ed(new Stat<Engine, Edges<Engine> >());
    Rcpp::List ll;
    ll.push_back("var");
    boost::shared_ptr< AbstractOffset<Engine> > off(
        		new Offset<Engine,REffect<Engine> >(ll));
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addOffsetPtr(off);
    model.calculate();
    model.setThetas(std::vector<double>(1,0));
    std::vector<int> togVars(1,0);
    togVars.push_back(0);
    model.setRandomVariables(togVars,false);
    MetropolisHastings<Engine> mh(model);
    mh.setDyadProbability(.5);
    mh.initialize();
    mh.run(4000);
    double val = mh.getModel()->offset().at(0);
    model.calculate();
    EXPECT_NEAR(val,model.offset().at(0))

    PutRNGstate();
}

void testConstraints(){
	testContext = "Constraints";
	RUN_TEST(testBoundedDegree<Undirected>());
	RUN_TEST(testREffect<Undirected>());
}

}

} /* namespace ernm */
