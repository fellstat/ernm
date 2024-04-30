/*
 * main.cpp
 *
 *  Created on: Sep 5, 2011
 *      Author: ianfellows
 */
//#define INSIDE
#ifdef INSIDE

#include <vector>
#include <Rcpp.h>
#include <BinaryNet.h>
#include <RInside.h> // for the embedded R via RInside
#include <Stat.h>
#include <Stats.h>
#include <StatController.h>
#include <Offset.h>
#include <Offsets.h>
#include <Constraint.h>
#include <Constraints.h>
#include <Model.h>
#include <DyadToggles.h>
#include <MetropolisHastings.h>
#include <VarAttrib.h>
#include <VertexToggle.h>
#include <VertexToggles.h>
#include <ToggleController.h>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <assert.h>
#include <tests.h>
#include <Quasi.h>
#include <CdSampler.h>
#include <boost/container/flat_set.hpp>
#include <ctime>
#undef NDEBUG
using namespace Rcpp;
using namespace std;
using namespace ernm;


void quasiTest(RInside& R){
	typedef Directed Engine;

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
	for(int i=0;i<30;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	net.addDiscreteVariable(vals,attr);
	net.addDiscreteVariable(vals1,attr1);
	//calculate network statistics
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(14);
    Rcpp::List fact;
    fact.push_back("fact");
    boost::shared_ptr< AbstractStat<Engine> > ed(new Stat<Engine,Edges<Engine> >());
    boost::shared_ptr< AbstractStat<Engine> > act(new Stat<Engine,DiffActivity<Engine> >(fact));
    boost::shared_ptr< AbstractStat<Engine> > stat(
    		new Stat<Engine, Homophily<Engine> >(fact));
    vector<int> tmp1;
    tmp1.push_back(0);
    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addStatPtr(stat);
    model.setRandomVariables(tmp1,true);
    //model.addOffsetPtr(off);
    model.calculate();
    //stat->theta()[0] = -1000;
    Language call1("print",wrap(model.statistics()));
    call1.eval();

    QuasiRandomPatch<Engine> q = QuasiRandomPatch<Engine>(14);
    q.initialize(model);
    q.calculate(model);
    cout<<"\nlog lik:" << q.logLik(model.thetas())<<"\n";


    model.calculateStatistics();
    Language call3("print",wrap(model.statistics()));
    call3.eval();
    //Language call4("print",wrap(mh.generateSampleStatistics(100,100,100)));
    //  call4.eval();
    PutRNGstate();
}

void changeStatTest(RInside& R){
	typedef Directed Engine;

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
	for(int i=0;i<30;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	net.addDiscreteVariable(vals,attr);
	net.addDiscreteVariable(vals1,attr1);
	//calculate network statistics
    vector<int> q;
    q.push_back(3);
    //q.push_back(10);
//    q.push_back(15);
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(14);
    Rcpp::List fact;
    fact.push_back("fact");
    boost::shared_ptr< AbstractStat<Engine> > ed(new Stat<Engine, Edges<Engine> >());
    boost::shared_ptr< AbstractStat<Engine> > act(new Stat<Engine, DiffActivity<Engine> >(fact));
    //ed->theta()[0] = -1000;
//    boost::shared_ptr< NodeMatchFull<Engine> > stat(
//    		new NodeMatchFull<Engine>("fact"));
//    boost::shared_ptr< DiscreteNode<Engine> > stat(
//    		new DiscreteNode<Engine>("fact"));
//    boost::shared_ptr< DegreeDistribution<Engine> > stat(
//    		new DegreeDistribution<Engine>());
//    boost::shared_ptr< DegreeQuantiles<Engine> > stat(
//    		new DegreeQuantiles<Engine>(q));
//   boost::shared_ptr< Transitivity<Engine> > stat(
//    		new Transitivity<Engine>());
//    boost::shared_ptr< LogTriangle<Engine> > stat(
//    		new LogTriangle<Engine>());
//    boost::shared_ptr< Logistic<Engine> > stat(
//    		new Logistic<Engine>("out","fact"));
//    boost::shared_ptr< DiffActivity<Engine> > stat(new DiffActivity<Engine>("fact"));
    boost::shared_ptr< AbstractStat<Engine> > stat(
    		new Stat<Engine, Homophily<Engine> >(fact));
//    boost::shared_ptr< Degree<Engine> > stat(
//    		new Degree<Engine>(degrees));

    vector<int> tmp1;
    tmp1.push_back(4);
    tmp1.push_back(1);
    tmp1.push_back(0);
    //boost::shared_ptr< BiasedSeedOffset<Engine> > off(
    //		new BiasedSeedOffset<Engine>("fact",tmp1));
    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addStatPtr(stat);
    //model.addOffsetPtr(off);
    model.calculate();
    //stat->theta()[0] = -1000;
    Language call1("print",wrap(model.statistics()));
    call1.eval();
    vector<int> togVars(2,0);
    togVars[1] = 1;

    // run mcmc
    DyadToggle<Engine,TieDyad<Engine> > tog(net);
    VertexToggle<Engine, DefaultVertex<Engine> > vtog(net);
    vtog.vSetDiscreteVars(togVars);
    MetropolisHastings<Engine> mh(model,tog,vtog);
    mh.setDyadProbability(1.0);
    //cout<<net.discreteVariableValue(0,0)<<net.discreteVariableValue(0,1)
	//	<<net.discreteVariableValue(0,4)<<" "<<net.discreteVariableValue(0,2)<<"\n";
    mh.initialize();
   // for(int i=0;i<30;i++){
    	mh.run(400);
    	//cout <<"acceptance rate: "<< mh.run(200) << "\n";
    	//cout << "edges: " << net.nEdges()<<"\n";
   // }
    //cout<<net.discreteVariableValue(0,0)<<net.discreteVariableValue(0,1)
    //	<<net.discreteVariableValue(0,4)<<" "<<net.discreteVariableValue(0,2)<<"\n";
    //Language call0("print",wrap(model.offset()));
    //call0.eval();
    Language call2("print",wrap(mh.getModel()->statistics()));
    call2.eval();
    model.calculateStatistics();
    Language call3("print",wrap(model.statistics()));
    call3.eval();
    //Language call4("print",wrap(mh.generateSampleStatistics(100,100,100)));
    //  call4.eval();
    PutRNGstate();
}


void constraintTest(RInside& R){
	typedef Undirected Engine;

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
	for(int i=0;i<30;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	net.addDiscreteVariable(vals,attr);
	net.addDiscreteVariable(vals1,attr1);
	//calculate network statistics
    vector<int> q;
    q.push_back(3);
    //q.push_back(10);
//    q.push_back(15);
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(14);
    boost::shared_ptr< AbstractStat<Engine> > ed(new Stat<Engine, Edges<Engine> >());
 //   boost::shared_ptr< DiffActivity<Engine> > act(new DiffActivity<Engine>("fact"));
    //ed->theta()[0] = -1000;
//    boost::shared_ptr< NodeMatchFull<Engine> > stat(
//    		new NodeMatchFull<Engine>("fact"));
//    boost::shared_ptr< DiscreteNode<Engine> > stat(
//    		new DiscreteNode<Engine>("fact"));
//    boost::shared_ptr< DegreeDistribution<Engine> > stat(
//    		new DegreeDistribution<Engine>());
//    boost::shared_ptr< DegreeQuantiles<Engine> > stat(
//    		new DegreeQuantiles<Engine>(q));
//   boost::shared_ptr< Transitivity<Engine> > stat(
//    		new Transitivity<Engine>());
//    boost::shared_ptr< LogTriangle<Engine> > stat(
//    		new LogTriangle<Engine>());
//    boost::shared_ptr< Logistic<Engine> > stat(
//    		new Logistic<Engine>("out","fact"));
//    boost::shared_ptr< DiffActivity<Engine> > stat(new DiffActivity<Engine>("fact"));
//    boost::shared_ptr< NodalHomophily<Engine> > stat(
//    		new NodalHomophily<Engine>("fact"));
//    boost::shared_ptr< Degree<Engine> > stat(
//    		new Degree<Engine>(degrees));
    vector<int> tmp1;
    tmp1.push_back(4);
    tmp1.push_back(1);
    tmp1.push_back(0);
    //boost::shared_ptr< BiasedSeedOffset<Engine> > off(
    //		new BiasedSeedOffset<Engine>("fact",tmp1));
    //boost::shared_ptr< FixedNode<Engine> > off(
    //		new FixedNode<Engine>(tmp1));
    //boost::shared_ptr< BoundedDegree<Engine> > off(
    //    		new BoundedDegree<Engine>(5,10));
    Rcpp:List ll;
    ll.push_back(tmp1);
    ll.push_back(degrees);
    boost::shared_ptr< AbstractOffset<Engine> > off(
        		new Constraint<Engine,FixedDegree<Engine> >(ll));
    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    //model.addStatPtr(stat);
    model.addOffsetPtr(off);
    model.calculate();
    //stat->theta()[0] = -1000;
    Language call1("print",wrap(model.statistics()));
    call1.eval();
    vector<int> togVars(2,0);
    togVars[1] = 1;

    // run mcmc
    DyadToggle<Engine, TieDyad<Engine> > tog(net);
    VertexToggle<Engine, DefaultVertex<Engine> > vtog(net);
    vtog.vSetDiscreteVars(togVars);
    MetropolisHastings<Engine> mh(model,tog,vtog);
    //mh.setDyadProbability(0.0);


    //cout<<net.discreteVariableValue(0,0)<<net.discreteVariableValue(0,1)
	//	<<net.discreteVariableValue(0,4)<<" "<<net.discreteVariableValue(0,2)<<"\n";
    cout<<"degrees: ";
    for(int i=0;i<net.size();i++){
    	cout<<net.degree(i)<<" ";
    }
    cout<<"\n";



    mh.initialize();
    Language cal("print",wrap(mh.getModel()->offset()));
    cal.eval();
    cout << "offset: ";
    for(int i=0;i<30;i++){
    	mh.run(40);
    	//cout <<"acceptance rate: "<< mh.run(200) << "\n";
    	cout << mh.getModel()->offset().at(0) <<" ";
   }
    cout<<"\n";


    //cout<<"\n"<<net.discreteVariableValue(0,0)<<net.discreteVariableValue(0,1)
    //	<<net.discreteVariableValue(0,4)<<" "<<net.discreteVariableValue(0,2)<<"\n";
    cout<<"degrees: ";
    for(int i=0;i<net.size();i++){
    	cout<<mh.getModel()->network()->degree(i)<<" ";
    }
    cout<<"\n";




    Language call0("print",wrap(mh.getModel()->offset()));
    call0.eval();
    Language call2("print",wrap(mh.getModel()->statistics()));
    call2.eval();
    //model.calculateStatistics();
    //Language call3("print",wrap(model.statistics()));
    //call3.eval();
    Language call4("print",wrap(mh.generateSampleStatistics(100,100,100)));
      call4.eval();
    PutRNGstate();
}

void tetradTest(RInside& R){
	typedef Undirected Engine;
	typedef DyadToggle<Engine, Tetrad<Engine> > DyadToggler ;
	//typedef SimpleToggle<Engine> DyadToggler ;
	typedef VertexToggle<Engine, DefaultVertex<Engine> > VertexToggler ;

	//create network
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,40);
	GetRNGstate();
	for(int i=0;i<100;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	PutRNGstate();
	//Language call2("print",net.edgelistR());
	//call2.eval();
	//calculate network statistics
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(4);
    degrees.push_back(5);
    degrees.push_back(6);
    Rcpp::List deg;
    deg.push_back(degrees);
    boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed(new Stat<Engine,Edges<Engine> >());
    //ed->theta()[0] = 1;
    boost::shared_ptr< Stat<Engine,Degree<Engine> > > st(new Stat<Engine,Degree<Engine> >(deg));
    boost::shared_ptr< Stat<Engine ,Transitivity<Engine> > > clust(new Stat<Engine ,Transitivity<Engine> >());
    clust->theta()[0] = 1;
    //ed->calculate(net);
    cout << "Edges: " << ed->statistics().at(0) << "\n";

	boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed1(new Stat<Engine,Edges<Engine> >());
	ed1->theta()[0] = -1;
    Model<Engine> model1(net);
    model1.addStatPtr(ed1);
    model1.addStatPtr(clust);
	model1.calculateStatistics();

	Language call3("print",model1.statistics());
	call3.eval();

    // run mcmc
    DyadToggler tog1(net);
    MetropolisHastings<Engine> mh1(model1);
    mh1.setDyadProbability(1.0);
    mh1.initialize();
	GetRNGstate();
	double sum1 = 0.0;
	mh1.run(1000);
    for(int i=0;i<300;i++){
    	//mh1.run(10);
    	cout <<"acceptance rate: "<< mh1.run(10) << "\n";
    	//sum1 += net.nEdges();
    	cout << "edges: " << net.nEdges()<<"\n";
    	//assert(net.hasEdge(1,0));
    	//assert(!net.hasEdge(0,2));
    }
	PutRNGstate();
	//Language call4("print",net.edgelistR());
	//call4.eval();
	Language call4("print",mh1.getModel()->statistics());
	call4.eval();
	mh1.getModel()->calculateStatistics();
	Language call5("print",mh1.getModel()->statistics());
	call5.eval();
}

void rdsTogglerTest(RInside& R){
	typedef Undirected Engine;
	typedef DyadToggle<Engine, Rds<Engine> > DyadToggler ;
	//typedef SimpleToggle<Engine> DyadToggler ;
	typedef VertexToggle<Engine, DefaultVertex<Engine> > VertexToggler ;
	int n=100;
	//create network
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,n);
	GetRNGstate();
	for(int i=0;i<3*n;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	PutRNGstate();
	net.setAllDyadsMissing();
	net.setMissing(0,1,false);
	net.addEdge(0,1);
	net.setMissing(1,2,false);
	net.addEdge(1,2);
	net.setMissing(1,4,false);
	net.addEdge(1,4);

	//Language call2("print",net.edgelistR());
	//call2.eval();
	//calculate network statistics
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(4);
    degrees.push_back(5);
    degrees.push_back(6);
    boost::shared_ptr< Edges<Engine> > ed(new Edges<Engine>());
    //ed->theta()[0] = 1;
    boost::shared_ptr< Degree<Engine> > st(new Degree<Engine>(degrees));
    boost::shared_ptr< Transitivity<Engine> > clust(new Transitivity<Engine>());
    clust->theta()[0] = 1;
    //ed->calculate(net);
    cout << "Edges: " << ed->statistics().at(0) << "\n";

	boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed1(new Stat<Engine, Edges<Engine> >());
	ed1->theta()[0] = -3.5;
    Model<Engine> model1(net);
    model1.addStatPtr(ed1);
    //model1.addStatPtr(clust);
	model1.calculateStatistics();

	//Language call3("print",model1.statistics());
	//call3.eval();
	cout <<"degrees: ";
	for(int i=0;i<n;i++)
		cout<<net.degree(i)<<" ";
	cout<<"\n";

    // run mcmc
    DyadToggler tog1(net);
    MetropolisHastings<Engine> mh1(model1);
    mh1.setDyadProbability(1.0);
    mh1.initialize();
	GetRNGstate();
	double sum1 = 0.0;
	mh1.run(1000);
    for(int i=0;i<30;i++){
    	//mh1.run(10);
    	cout <<"acceptance rate: "<< mh1.run(10000) << "\n";
    	//sum1 += net.nEdges();
    	//cout << "edges: " << net.nEdges()<<"\n";
    	//assert(net.hasEdge(1,0));
    	//assert(!net.hasEdge(0,2));
    }
	PutRNGstate();
	cout <<"degrees: ";
	for(int i=0;i<n;i++)
		cout<<net.degree(i)<<" ";
	cout<<"\n";
	//Language call4("print",net.edgelistR());
	//call4.eval();
	/*Language call4("print",mh1.getModel()->statistics());
	call4.eval();
	mh1.getModel()->calculateStatistics();
	Language call5("print",mh1.getModel()->statistics());
	call5.eval();
	*/
}

void rdsOffsetTest(RInside& R){
	typedef Undirected Engine;
	typedef DyadToggle<Engine, Rds<Engine> > DyadToggler ;
	//typedef SimpleToggle<Engine> DyadToggler ;
	typedef VertexToggle<Engine, DefaultVertex<Engine> > VertexToggler ;
	int n=20;
	//create network
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,n);
	GetRNGstate();
	for(int i=0;i<3*n;i++){
		pair<int,int> dyad = net.randomDyad();
		net.addEdge(dyad.first,dyad.second);
	}
	PutRNGstate();
	net.setAllDyadsMissing();
	net.setMissing(0,1,false);
	net.addEdge(0,1);
	net.setMissing(1,2,false);
	net.addEdge(1,2);
	net.setMissing(1,4,false);
	net.addEdge(1,4);
	net.setMissing(1,10,false);
	net.addEdge(1,10);
	net.setMissing(4,5,false);
	net.addEdge(4,5);
	net.setMissing(5,6,false);
	net.addEdge(5,6);
	net.setMissing(5,7,false);
	net.addEdge(5,7);
	net.setMissing(5,8,false);
	net.addEdge(5,8);
	net.setMissing(8,9,false);
	net.addEdge(8,9);
	//Language call2("print",net.edgelistR());
	//call2.eval();
	//calculate network statistics
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(4);
    degrees.push_back(5);
    degrees.push_back(6);
    boost::shared_ptr< Edges<Engine> > ed(new Edges<Engine>());
    //ed->theta()[0] = 1;
    boost::shared_ptr< Degree<Engine> > st(new Degree<Engine>(degrees));
    boost::shared_ptr< Transitivity<Engine> > clust(new Transitivity<Engine>());
    clust->theta()[0] = 0;
    //ed->calculate(net);
    cout << "Edges: " << ed->statistics().at(0) << "\n";

	boost::shared_ptr< Stat<Engine, Edges<Engine> > > ed1(new Stat<Engine,Edges<Engine> >());
	ed1->theta()[0] = -1;

	vector<int> order = vector<int>(n,-1);
	order[0] = 1;
	order[1] = 2;
	order[4] = 3;
	order[2] = 4;
	order[5] = 5;
	order[6] = 6;
	order[7] = 7;
	order[8] = 8;
	order[9] = 9;
	order[10] = 10;
	Rcpp::List l;
	l.push_back(order);
	boost::shared_ptr< AbstractOffset<Engine> > rds(new Offset<Engine, RdsBias<Engine> >(l));

    Model<Engine> model1(net);
    model1.addStatPtr(ed1);
    model1.addOffsetPtr(rds);
    //model1.addStatPtr(clust);
	model1.calculateStatistics();
	model1.calculateOffsets();

	Language call3("print",model1.offset());
	call3.eval();
	cout <<"degrees: ";
	for(int i=0;i<n;i++)
		cout<<net.degree(i)<<" ";
	cout<<"\n";

    // run mcmc
    DyadToggler tog1(net);
    MetropolisHastings<Engine> mh1(model1);
    mh1.setDyadProbability(1.0);
    mh1.initialize();
	GetRNGstate();
	double sum1 = 0.0;
	mh1.run(1000);
    for(int i=0;i<3;i++){
    	//mh1.run(10);
    	cout <<"acceptance rate: "<< mh1.run(20000) << "\n";
    	//sum1 += net.nEdges();
    	//cout << "edges: " << net.nEdges()<<"\n";
    	//assert(net.hasEdge(1,0));
    	//assert(!net.hasEdge(0,2));
    }
	PutRNGstate();
	cout <<"degrees: ";
	for(int i=0;i<n;i++)
		cout<<net.degree(i)<<" ";
	cout<<"\n";
	//Language call4("print",net.edgelistR());
	//call4.eval();
	Language call4("print",mh1.getModel()->offset());
	call4.eval();
	mh1.getModel()->calculateStatistics();
	mh1.getModel()->calculateOffsets();
	Language call5("print",mh1.getModel()->offset());
	call5.eval();

}

void togglerTest(RInside& R){
	typedef Directed Engine;

	//create network
	IntegerMatrix tmp(0,2);
	BinaryNet<Engine> net(tmp,40);
	net.setDyad(0,1,true);
	net.setDyad(1,0,true);
	net.setDyad(0,3,true);
	net.setDyad(3,1,true);
	net.setMissing(0,3,true);
	net.setMissing(0,1,true);
	net.setMissing(1,3,true);
	net.setAllDyadsMissing();
	//check missing
	//assert(net.isMissing(0,1));
	//assert(net.isMissing(1,0));

	boost::shared_ptr< vector< pair<int,int> > > miss = net.missingDyads();
	//cout << "Missing dyads:\n";
	//for(int i=0;i<miss->size();i++)
	//	cout << miss->at(i).first << " " << miss->at(i).second<<"\n";

	//calculate network statistics
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed(new Stat<Engine, Edges<Engine> >());
    //ed->theta()[0] = 1;
    //boost::shared_ptr< Istar<Engine> > st(new Istar<Engine>(degrees));
    ed->calculate(net);
    cout << "Edges: " << ed->statistics().at(0) << "\n";

    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.calculateStatistics();





	boost::shared_ptr< Stat<Engine, Edges<Engine> > > ed1(new Stat<Engine,Edges<Engine> >());
	ed1->calculate(net);
	ed1->theta()[0] = -1;
    Model<Engine> model1(net);
    model1.addStatPtr(ed1);
	model1.calculateStatistics();

	//Language call3("print",net.edgelistR());
	//    call3.eval();

    // run mcmc
	Rcpp::List l;
	//CompoundToggle<NTDToggle<Engine>, NeighborhoodToggle<Engine>, Engine > tmp3;
	//cout << tmp3.vName();
	CdSampler<Engine> cd(model1);
    MetropolisHastings<Engine> mh1(model1);
    mh1.setDyadToggleType("TieDyad",l);
    mh1.setDyadProbability(.5);
    mh1.initialize();
	GetRNGstate();
	double sum1 = 0.0;
	mh1.run(20000);
	mh1.setDyadToggleType("Neighborhood",l);
	mh1.initialize();
    for(int i=0;i<1000;i++){
    	mh1.run(100);
    	//cout <<"acceptance rate: "<< mh1.run(10) << "\n";
    	sum1 += net.nEdges();
    	cout << "edges: " << net.nEdges()<<"\n";
    	//assert(net.hasEdge(1,0));
    	//assert(!net.hasEdge(0,2));
    }
	PutRNGstate();
	cout << "average nEdges Rand: " << (sum1/1000.0) << "\n";
//
//    // run mcmc
//    NTDNonObservedToggle<Engine> tog(net);
//    MetropolisHastings< NTDNonObservedToggle<Engine>, Model<Engine> > mh(model,tog);
//    mh.initialize();
//	GetRNGstate();
//	double sum = 0.0;
//    for(int i=0;i<30000;i++){
//    	mh.run(10);
//    	//cout <<"acceptance rate: "<< mh.run(199) << "\n";
//    	sum += net.nEdges();
//    	//cout << "edges: " << net.nEdges()<<"\n";
//    	//assert(net.hasEdge(1,0));
//    	//assert(!net.hasEdge(0,2));
//    }
//	PutRNGstate();
//	cout << "average nEdges NTD : " << (sum/30000.0) << "\n";
}

void flatSetTest(){
	GetRNGstate();
	std::clock_t start;
	double d1,d2;
	boost::container::flat_set<int> fs;
	std::set<int> ss;
	int n = 500000;
	int max = 100000000;
	for(int i = 1;i<14;i++){
		int s = round(pow(2.0,i));
		d1=0.0;
		d2=0.0;
		fs.clear();
		fs.reserve(s*2);
		ss.clear();
		for(int j=0;j<s;j++){
			int val = round(((j+1)/(double)s)*max);//floor(Rf_runif(0.0,max));
			fs.insert(val);
			ss.insert(val);
		}
		start = std::clock();
		for(int j=0;j<n;j++){
			int rand = floor(Rf_runif(0.0,max));
			bool ins = fs.insert(rand).second;
			if(ins)
				fs.erase(rand);
		}
		d1 += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		for(int j=0;j<n;j++){
			int rand = floor(Rf_runif(0.0,max));
			bool ins = ss.insert(rand).second;
			if(ins)
				ss.erase(rand);
		}
		d2 += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << s << ", " << d1 << ", " << d2 << "," << endl;
	}
	PutRNGstate();
}



/*!
 * Main entry point when called from outside R (with R embedded via RInside).
 */
int main(int argc, char *argv[]) {
    RInside R(argc, argv); // create an embedded R instance
    initToggles();
    initStats();
/**/
    RInside R(argc, argv);              // create an embedded R instance
    initToggles();
    initStats();
    /*
    std::string txt = "Hello, world!\n";// assign a standard C++ string to 'txt'
    R.assign( txt, "txt");              // assign string var to R variable 'txt'
    std::string evalstr = "cat(txt)";
    for (int i=0; i<1e1; i++) {
        R.parseEvalQ(evalstr);          // eval the init string, ignoring any returns
    }
    evalstr = "txt <- \"foo\\n\"";
    for (int i=0; i<1e6; i++) {
        R.parseEvalQ(evalstr);          // eval the init string, ignoring any returns
    }
    evalstr = "cat(txt)";
    for (int i=0; i<1e1; i++) {
        R.parseEvalQ(evalstr);          // eval the init string, ignoring any returns
    }*/
    //quasiTest(R);
    tests::runErnmTests();
    //flatSetTest();
    //togglerTest( R);
    //changeStatTest(R);
    //constraintTest(R);
    //togglerTest(R);
    //tetradTest(R);
    //rdsOffsetTest(R);
    //cout << expectedSqrtHypergeometric(14,0,14);
    exit(0);
    return 0;
}


#endif
