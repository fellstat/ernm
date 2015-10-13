/*
 * ToggleTests.cpp
 *
 *  Created on: May 6, 2014
 *      Author: ianfellows
 */

#include "ToggleTests.h"

#ifdef INSIDE
#include "../BinaryNet.h"
#include "../Stat.h"
#include "../Stats.h"
#include "../Constraint.h"
#include "../Model.h"
#include "../DyadToggles.h"
#include "../MetropolisHastings.h"
#include "../VarAttrib.h"
#include "../VertexToggles.h"
#else
#include "BinaryNet.h"
#include "Stat.h"
#include "Stats.h"
#include "Constraint.h"
#include "Model.h"
#include "DyadToggles.h"
#include "MetropolisHastings.h"
#include "VarAttrib.h"
#include "VertexToggles.h"
#endif
#include "tests.h"

namespace ernm{

namespace tests{

template <class Engine>
void toggleTest(){
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,30);
	GetRNGstate();
	Language call1("set.seed",wrap(1.0));
	call1.eval();
	for(int i=0;i<30;i++){
		std::pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	std::vector<double>vals2;
	for(int i=0; i<30; i++)
		vals2.push_back(Rf_runif(-90,90));
	ContinAttrib attr2;
	attr2.setName("contin");
	attr2.setLowerBound(-90.0);
	attr2.setUpperBound(90.0);
	net.addContinVariable(vals2,attr2);

	boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed(new Stat<Engine,Edges<Engine> >());
    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    //model.addOffsetPtr(off);
    model.calculate();

    std::vector<int> togVars(1,0);
    model.setRandomVariables(togVars,false);

    // run mcmc
    DyadToggle<Engine,TieDyad<Engine> > tog(net);
    VertexToggle<Engine, DefaultVertex<Engine> > vtog(net);

    MetropolisHastings<Engine> mh(model,tog,vtog);
    mh.setDyadProbability(.5);
    mh.initialize();
    mh.run(500);
    for(int i=0;i<30;i++){
    	//std::cout << net.continVariableValue(0,i) << " " << vals2[i] << "\n";
    	EXPECT_TRUE(net.continVariableValue(0,i) <= 90.0);
    	EXPECT_TRUE(net.continVariableValue(0,i) >= -90.0);
    }

    //test vertex
    std::vector<double> th;
    int n = 1000;
    th.push_back(0);
    th.push_back(-.5);
    BinaryNet<Engine> net1(tmp,n);
    Model<Engine> model1(net1);
	Rcpp::List l;
	l.push_back("contin1");
	boost::shared_ptr< Stat<Engine,Gauss<Engine> > > gauss(new Stat<Engine,Gauss<Engine> >(l));
	ContinAttrib attr;
	attr.setName("contin1");
	vals2 = std::vector<double>(n,0);
	for(int i=0;i<n;i++)
		vals2.at(i) = 0;
	net1.addContinVariable(vals2,attr);
	model1.addStatPtr(gauss);
	model1.setRandomVariables(togVars,false);
	model1.calculate();
	model1.setThetas(th);
	//std::cout << net1.size();
	//for(int i=0;i<model1.statistics().size();i++)
	//	std::cout << " " << model1.statistics()[i];
	//std::cout << "\n";

    DyadToggle<Engine,TieDyad<Engine> > tog1(net1);
    VertexToggle<Engine, DefaultVertex<Engine> > vtog1(net1);//,vector<int>(),togVars);
    MetropolisHastings<Engine> mh1(model1,tog1,vtog1);
    mh1.setDyadProbability(.5);
    mh1.initialize();
    double init = model1.statistics()[1];
    mh1.run(100000);
    double s = 0;
    double ssq = 0;
    for(int i=0;i<n;i++){
    	s += net1.continVariableValue(0,i);
    	ssq += pow(net1.continVariableValue(0,i),2.0);
    }
    //std::cout << "mean: " << s/(double)n <<" var: " << ssq/(double)n - pow(s/(double)n,2.0) << "\n";
    EXPECT_TRUE(s/(double)n > -.1 && s/(double)n<.1);
    EXPECT_TRUE(ssq/(double)n - pow(s/(double)n,2.0) >.9 && ssq/(double)n - pow(s/(double)n,2.0)  < 1.1);
    //EXPECT_TRUE(model1.statistics()[1] != init)

    //missing
    for(int i=0;i<50;i++){
    	net1.setContinVariableObserved(0,i,FALSE);
    }
    for(int i=50;i<100;i++){
    	net1.setContinVariableValue(0,i,-1.0);
    }
    VertexToggle<Engine, VertexMissing<Engine> > vtog2(net1);
    MetropolisHastings<Engine> mh2(model1,tog1,vtog2);
    mh2.setDyadProbability(0);
    mh2.initialize();

    init = net1.continVariableValue(0,0);
    mh2.run(10000);
    for(int i=50;i<100;i++){
    	EXPECT_NEAR(net1.continVariableValue(0,i),-1.0);
    }
    s = 0.0;
    for(int i=0;i<50;i++){
    	s += net1.continVariableValue(0,i);
    }
    s /=50.0;
    EXPECT_TRUE(s < .5)
    EXPECT_TRUE(s > -.5)
    EXPECT_TRUE(net1.continVariableValue(0,0) != init)
    PutRNGstate();
}

void testToggles(){
	testContext = "Toggle";
	RUN_TEST(toggleTest<Directed>());
	RUN_TEST(toggleTest<Undirected>())

}

}

} /* namespace ernm */
