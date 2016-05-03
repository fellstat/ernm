/*
 * ReModelTests.cpp
 *
 *  Created on: Jan 15, 2016
 *      Author: ianfellows
 */




#include <Rcpp.h>
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



template<class Engine>
void re(){
	using namespace std;
	using namespace ernm;
	vector<int> vals(30,1);
	vals[3] = 3;
	vals[2] = 2;
	vals[4] = 2;
	vector<string> labels(3,"a");
	labels[1] = "b";
	labels[2] = "c";
	DiscreteAttrib attr;
	attr.setName("fact");
	attr.setLabels(labels);

	vector<int> vals1(30,1);
	vals1[3] = 2;
	vals1[2] = 2;
	vals1[4] = 2;
	vector<string> labels1(2,"a");
	labels1[1] = "b";
	DiscreteAttrib attr1;
	attr1.setName("out");
	attr1.setLabels(labels1);

	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,30);
	GetRNGstate();
	//Language call1("set.seed",wrap(1.0));
	//call1.eval();
	for(int i=0;i<30;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	net.addDiscreteVariable(vals,attr);
	net.addDiscreteVariable(vals1,attr1);

	vector<double>vals2;
	for(int i=0; i<30; i++)
		vals2.push_back(Rf_runif(-90,90));
	ContinAttrib attr2;
	attr2.setName("contin");
	attr2.setLowerBound(-90);
	attr2.setUpperBound(90);
	net.addContinVariable(vals2,attr2);

	vals2.clear();
	for(int i=0; i<30; i++)
		vals2.push_back(Rf_runif(-180,180));
	ContinAttrib attr3;
	attr3.setName("contin1");
	attr3.setLowerBound(-180);
	attr3.setUpperBound(180);
	net.addContinVariable(vals2,attr3);

	//calculate network statistics
    vector<int> q;
    q.push_back(3);
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(14);
    Rcpp::List deg;
    deg.push_back(degrees);
    boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed(new Stat<Engine,Edges<Engine> >());
    boost::shared_ptr< AbstractStat<Engine> > stat;
    Rcpp::List fact;
    fact.push_back("fact");
    //string fact = "fact";
    Rcpp::List log;
    log.push_back("out");
    log.push_back("fact");
    Rcpp::List hom;
    hom.push_back("fact");
    hom.push_back(0+UNDIRECTED);
    hom.push_back(true);
    hom.push_back(false);
    stat = boost::shared_ptr< Stat<Engine, Triangles<Engine> > >(
            		new Stat<Engine, Triangles<Engine> >());
    vector<int> tmp1;
    tmp1.push_back(4);
    tmp1.push_back(1);
    tmp1.push_back(0);

    //create model
    std::vector<double> centers(2,0.0);
    std::vector<double> betas(2,0.0);
    betas.at(0) = .1;
    ReModel<Engine> model(net);
    model.addStatPtr(ed);
    model.addStatPtr(stat);
    model.setCenters(centers);
    model.setBetas(betas);
    model.thetaDependent(false);
    model.calculate();
    stat->vTheta().at(0) = 0.0;


    vector<int> togVars(1,0);
    togVars.push_back(1);
    model.setRandomVariables(togVars,false);
    model.setRandomVariables(togVars,true);

    // run mcmc
    DyadToggle<Engine,TieDyad<Engine> > tog(net);
    VertexToggle<Engine, DefaultVertex<Engine> > vtog(net);//,vector<int>(),togVars);
    MetropolisHastings<Engine> mh(model,tog,vtog);
    mh.setDyadProbability(0.5);
    //cout<<net.discreteVariableValue(0,0)<<net.discreteVariableValue(0,1)
	//	<<net.discreteVariableValue(0,4)<<" "<<net.discreteVariableValue(0,2)<<"\n";
    mh.initialize();

    Language call3("print",wrap(model.betaParams()));
    call3.eval();
   // for(int i=0;i<30;i++){
    mh.run(10);
    	//cout <<"acceptance rate: "<< mh.run(200) << "\n";
    	//cout << "edges: " << net.nEdges()<<"\n";
   // }
    //cout<<net.discreteVariableValue(0,0)<<net.discreteVariableValue(0,1)
    //	<<net.discreteVariableValue(0,4)<<" "<<net.discreteVariableValue(0,2)<<"\n";
    vector<double> mcmcStats = mh.getModel()->statistics();
    model.calculateStatistics();
    vector<double> realStats = model.statistics();

    for(int i=0;i<realStats.size();i++){
    	//cout << i << " " << mcmcStats.at(i) << " " << realStats.at(i) << " ";
    	EXPECT_NEAR((mcmcStats.at(i) + .0001)/(realStats.at(i) + .0001), 1.0);
    }
    PutRNGstate();
}


void testReModel(){
	testContext = "ReModel";

	RUN_TEST(re<Directed>());
    RUN_TEST(re<Undirected>());
}


}
}
