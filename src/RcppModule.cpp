
#ifndef INSIDE

#include <Rcpp.h>
#include "BinaryNet.h"
#include "DyadToggles.h"
#include "VertexToggles.h"
#include "MetropolisHastings.h"
#include "CdSampler.h"
#include "LatentOrderLikelihood.h"

/*
 * Handles all functions and methods exported to R.
 */

RCPP_MODULE(ernm){
	using namespace Rcpp ;
	using namespace ernm;

//	typedef MetropolisHastings< SimpleToggle<Directed>, Model<Directed> > SimpleMetropolisDirected;
//	typedef MetropolisHastings< TDBasicToggle<Directed>, Model<Directed> > TDBasicMetropolisDirected;
//	typedef MetropolisHastings< TDToggle<Directed>, Model<Directed> > TDMetropolisDirected;
/*
	typedef MetropolisHastingsVertex< TDToggle<Directed>, DefaultVertexToggle<Directed>,
			Model<Directed> > DirectedVertexMetropolis;
	typedef MetropolisHastingsVertex< TDToggle<Undirected>, DefaultVertexToggle<Undirected>,
				Model<Undirected> > UndirectedVertexMetropolis;

	typedef MetropolisHastingsVertex< NTDNonObservedToggle<Directed>,
				NonObservedVertexToggle<Directed>,
				Model<Directed> > DirectedNonObservedVertexMetropolis;
	typedef MetropolisHastingsVertex< NTDNonObservedToggle<Undirected>,
					NonObservedVertexToggle<Undirected>,
					Model<Undirected> > UndirectedNonObservedVertexMetropolis;
	typedef MetropolisHastingsVertex< RDSToggle<Undirected>,
						NonObservedVertexToggle<Undirected>,
						Model<Undirected> > RDSMetropolis;
	//NeighborhoodToggle<Undirected>
	typedef MetropolisHastingsVertex< UndirectedNTDNBRToggle,
			DefaultVertexToggle<Undirected>,
					Model<Undirected> > UndirectedNTDNBRMetropolis;
	typedef MetropolisHastingsVertex< DirectedNTDNBRToggle,
				DefaultVertexToggle<Directed>,
						Model<Directed> > DirectedNTDNBRMetropolis;

	typedef MetropolisHastingsVertex< UndirectedNTDNBRNonObservedToggle,
			NonObservedVertexToggle<Undirected>,
					Model<Undirected> > UndirectedNTDNBRNonObservedMetropolis;
	typedef MetropolisHastingsVertex< DirectedNTDNBRNonObservedToggle,
			NonObservedVertexToggle<Directed>,
						Model<Directed> > DirectedNTDNBRNonObservedMetropolis;

	typedef MetropolisHastingsVertex< RandomNonObservedToggle<Undirected>,
			NonObservedVertexToggle<Undirected>,
					Model<Undirected> > UndirectedRandomNonObservedMetropolis;
	typedef MetropolisHastingsVertex< RandomNonObservedToggle<Directed>,
			NonObservedVertexToggle<Directed>,
						Model<Directed> > DirectedRandomNonObservedMetropolis;
*/

	class_<DirectedNet >("DirectedNet")
			.constructor<Rcpp::IntegerMatrix,int>()
			.constructor<SEXP>()
			.method("clone",&DirectedNet::cloneR)
			.method("size",&DirectedNet::size)
			.method("isDirected",&DirectedNet::isDirected)
			.method("setDyads",&DirectedNet::setDyadsR)
			.method("getDyads",&DirectedNet::getDyadsR)
			.method("edges",&DirectedNet::edgelistR1)
			.method("edges",&DirectedNet::edgelistR2)
			.method("[",&DirectedNet::getDyadMatrixR)
			.method("[<-",&DirectedNet::setDyadMatrixR)
			.method("variableNames",&DirectedNet::getVariableNamesR1)
			.method("variableNames",&DirectedNet::getVariableNamesR2)
			.method("[[",&DirectedNet::getVariableR)
			.method("getVariable",&DirectedNet::getVariableR)
			.method("getVariable",&DirectedNet::getVariableR1)
			.method("[[<-",&DirectedNet::setVariableR)
			.method("nMissing",&DirectedNet::nMissingR)
			.method("nEdges",&DirectedNet::nEdgesR1)
			.method("nEdges",&DirectedNet::nEdgesR2)
			.method("inDegree",&DirectedNet::indegreeR)
			.method("outDegree",&DirectedNet::outdegreeR)
			.method("outNeighbors",&DirectedNet::outneighborsR)
			.method("inNeighbors",&DirectedNet::inneighborsR)
			.method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR1)
			.method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR2)
			.method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR3)
;

	class_<UndirectedNet >("UndirectedNet")
			.constructor<Rcpp::IntegerMatrix,int>()
			.constructor<SEXP>()
			.method("clone",&UndirectedNet::cloneR)
			.method("size",&UndirectedNet::size)
			.method("isDirected",&UndirectedNet::isDirected)
			.method("setDyads",&UndirectedNet::setDyadsR)
			.method("getDyads",&UndirectedNet::getDyadsR)
			.method("edges",&UndirectedNet::edgelistR1)
			.method("edges",&UndirectedNet::edgelistR2)
			.method("[",&UndirectedNet::getDyadMatrixR)
			.method("[<-",&UndirectedNet::setDyadMatrixR)
			.method("variableNames",&UndirectedNet::getVariableNamesR1)
			.method("variableNames",&UndirectedNet::getVariableNamesR2)
			.method("[[",&UndirectedNet::getVariableR)
			.method("getVariable",&UndirectedNet::getVariableR)
			.method("getVariable",&UndirectedNet::getVariableR1)
			.method("[[<-",&UndirectedNet::setVariableR)
			.method("nMissing",&UndirectedNet::nMissingR)
			.method("nEdges",&UndirectedNet::nEdgesR1)
			.method("nEdges",&UndirectedNet::nEdgesR2)
			.method("degree",&UndirectedNet::degreeR)
			.method("neighbors",&UndirectedNet::neighborsR)
			.method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR1)
			.method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR2)
			.method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR3)
;

	//function("rcppHelloWorld",&rcppHelloWorld);

	//function("rcppTester",&rcppTester);
/*
	class_<Edges<Directed> >("DirectedEdges")
		.constructor()
		.method("calculate",&Edges<Directed>::calculate)
		.method("setTheta",&Edges<Directed>::setTheta)
		.method("statistics",&Edges<Directed>::statistics)
		;
	class_<Istar<Directed> >("Istar")
			.constructor< vector<int> >()
			.method("calculate",&Istar<Directed>::calculate)
			.method("setTheta",&Istar<Directed>::setTheta)
			.method("statistics",&Istar<Directed>::statistics)
			;*/
//	class_<SimpleMetropolisDirected>("SimpleMetropolisDirected")
//		.constructor<Model<Directed> >()
//		.method("setModel",&SimpleMetropolisDirected::setModel)
//		.method("getModel",&SimpleMetropolisDirected::getModelR)
//		.method("generateSample", &SimpleMetropolisDirected::generateSample)
//		.method("generateSampleStatistics", &SimpleMetropolisDirected::generateSampleStatistics)
//		;
//	class_<TDBasicMetropolisDirected>("TDBasicMetropolisDirected")
//		.constructor<Model<Directed> >()
//		.method("setModel",&TDBasicMetropolisDirected::setModel)
//		.method("getModel",&TDBasicMetropolisDirected::getModelR)
//		.method("generateSample", &TDBasicMetropolisDirected::generateSample)
//		;
//	class_<TDMetropolisDirected>("TDMetropolisDirected")
//		.constructor<Model<Directed> >()
//		.method("setModel",&TDMetropolisDirected::setModel)
//		.method("getModel",&TDMetropolisDirected::getModelR)
//		.method("generateSample", &TDMetropolisDirected::generateSample)
//		.method("generateSampleStatistics", &TDMetropolisDirected::generateSampleStatistics)
//		;
/*	class_<DirectedVertexMetropolis>("DirectedVertexMetropolis")
		.constructor<Model<Directed> >()
		.constructor<Model<Directed>,double >()
		.method("setModel",&DirectedVertexMetropolis::setModel)
		.method("getModel",&DirectedVertexMetropolis::getModelR)
		//.method("toggleDiscreteVars",&DirectedVertexMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&DirectedVertexMetropolis::setDyadProbability)
		.method("generateSample", &DirectedVertexMetropolis::generateSample)
		.method("generateSampleStatistics", &DirectedVertexMetropolis::generateSampleStatistics)
				;
	class_<DirectedNonObservedVertexMetropolis>("DirectedNonObservedVertexMetropolis")
		.constructor<Model<Directed> >()
		.constructor<Model<Directed>,double >()
		.method("setModel",&DirectedNonObservedVertexMetropolis::setModel)
		.method("getModel",&DirectedNonObservedVertexMetropolis::getModelR)
		//.method("toggleDiscreteVars",&DirectedNonObservedVertexMetropolis::toggleDiscreteVars)
		.method("generateSample", &DirectedNonObservedVertexMetropolis::generateSample)
		.method("setDyadProbability",&DirectedNonObservedVertexMetropolis::setDyadProbability)
		.method("generateSampleStatistics", &DirectedNonObservedVertexMetropolis::generateSampleStatistics)
		.method("generateSampleStatisticsSupplimental", &DirectedNonObservedVertexMetropolis::generateSampleStatisticsSupplimental)
		;
	class_<UndirectedNonObservedVertexMetropolis>("UndirectedNonObservedVertexMetropolis")
		.constructor<Model<Undirected> >()
		.constructor<Model<Undirected>,double >()
		.method("setModel",&UndirectedNonObservedVertexMetropolis::setModel)
		.method("getModel",&UndirectedNonObservedVertexMetropolis::getModelR)
		//.method("toggleDiscreteVars",&UndirectedNonObservedVertexMetropolis::toggleDiscreteVars)
		.method("generateSample", &UndirectedNonObservedVertexMetropolis::generateSample)
		.method("setDyadProbability",&UndirectedNonObservedVertexMetropolis::setDyadProbability)
		.method("generateSampleStatistics", &UndirectedNonObservedVertexMetropolis::generateSampleStatistics)
		.method("generateSampleStatisticsSupplimental", &UndirectedNonObservedVertexMetropolis::generateSampleStatisticsSupplimental)
		;

	class_<UndirectedVertexMetropolis>("UndirectedVertexMetropolis")
		.constructor<Model<Undirected> >()
		.constructor<Model<Undirected>,double >()
		.method("setModel",&UndirectedVertexMetropolis::setModel)
		.method("getModel",&UndirectedVertexMetropolis::getModelR)
		//.method("toggleDiscreteVars",&UndirectedVertexMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&UndirectedVertexMetropolis::setDyadProbability)
		.method("generateSample", &UndirectedVertexMetropolis::generateSample)
		.method("generateSampleStatistics", &UndirectedVertexMetropolis::generateSampleStatistics)
	;

	class_<UndirectedNTDNBRMetropolis>("UndirectedNTDNBRMetropolis")
		.constructor<Model<Undirected> >()
		.constructor<Model<Undirected>,double >()
		.method("setModel",&UndirectedNTDNBRMetropolis::setModel)
		.method("getModel",&UndirectedNTDNBRMetropolis::getModelR)
		//.method("toggleDiscreteVars",&UndirectedNTDNBRMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&UndirectedNTDNBRMetropolis::setDyadProbability)
		.method("generateSample", &UndirectedNTDNBRMetropolis::generateSample)
		.method("generateSampleStatistics", &UndirectedNTDNBRMetropolis::generateSampleStatistics)
		;
	class_<DirectedNTDNBRMetropolis>("DirectedNTDNBRMetropolis")
		.constructor<Model<Directed> >()
		.constructor<Model<Directed>,double >()
		.method("setModel",&DirectedNTDNBRMetropolis::setModel)
		.method("getModel",&DirectedNTDNBRMetropolis::getModelR)
		//.method("toggleDiscreteVars",&DirectedNTDNBRMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&DirectedNTDNBRMetropolis::setDyadProbability)
		.method("generateSample", &DirectedNTDNBRMetropolis::generateSample)
		.method("generateSampleStatistics", &DirectedNTDNBRMetropolis::generateSampleStatistics)
		;

	class_<UndirectedNTDNBRNonObservedMetropolis>("UndirectedNTDNBRNonObservedMetropolis")
		.constructor<Model<Undirected> >()
		.constructor<Model<Undirected>,double >()
		.method("setModel",&UndirectedNTDNBRNonObservedMetropolis::setModel)
		.method("getModel",&UndirectedNTDNBRNonObservedMetropolis::getModelR)
		//.method("toggleDiscreteVars",&UndirectedNTDNBRNonObservedMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&UndirectedNTDNBRNonObservedMetropolis::setDyadProbability)
		.method("generateSample", &UndirectedNTDNBRNonObservedMetropolis::generateSample)
		.method("generateSampleStatistics", &UndirectedNTDNBRNonObservedMetropolis::generateSampleStatistics)
		;
	class_<RDSMetropolis>("UndirectedRDSMetropolis")
		.constructor<Model<Undirected> >()
		.constructor<Model<Undirected>,double >()
		.method("setModel",&RDSMetropolis::setModel)
		.method("getModel",&RDSMetropolis::getModelR)
		//.method("toggleDiscreteVars",&RDSMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&RDSMetropolis::setDyadProbability)
		.method("generateSample", &RDSMetropolis::generateSample)
		.method("generateSampleStatistics", &RDSMetropolis::generateSampleStatistics)
		;
	class_<DirectedNTDNBRNonObservedMetropolis>("DirectedNTDNBRNonObservedMetropolis")
		.constructor<Model<Directed> >()
		.constructor<Model<Directed>,double >()
		.method("setModel",&DirectedNTDNBRNonObservedMetropolis::setModel)
		.method("getModel",&DirectedNTDNBRNonObservedMetropolis::getModelR)
		//.method("toggleDiscreteVars",&DirectedNTDNBRNonObservedMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&DirectedNTDNBRNonObservedMetropolis::setDyadProbability)
		.method("generateSample", &DirectedNTDNBRNonObservedMetropolis::generateSample)
		.method("generateSampleStatistics", &DirectedNTDNBRNonObservedMetropolis::generateSampleStatistics)
		;


	class_<UndirectedRandomNonObservedMetropolis>("UndirectedRandomNonObservedMetropolis")
		.constructor<Model<Undirected> >()
		.constructor<Model<Undirected>,double >()
		.method("setModel",&UndirectedRandomNonObservedMetropolis::setModel)
		.method("getModel",&UndirectedRandomNonObservedMetropolis::getModelR)
		//.method("toggleDiscreteVars",&UndirectedRandomNonObservedMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&UndirectedRandomNonObservedMetropolis::setDyadProbability)
		.method("generateSample", &UndirectedRandomNonObservedMetropolis::generateSample)
		.method("generateSampleStatistics", &UndirectedRandomNonObservedMetropolis::generateSampleStatistics)
		;
	class_<DirectedRandomNonObservedMetropolis>("DirectedRandomNonObservedMetropolis")
		.constructor<Model<Directed> >()
		.constructor<Model<Directed>,double >()
		.method("setModel",&DirectedRandomNonObservedMetropolis::setModel)
		.method("getModel",&DirectedRandomNonObservedMetropolis::getModelR)
		//.method("toggleDiscreteVars",&DirectedRandomNonObservedMetropolis::toggleDiscreteVars)
		.method("setDyadProbability",&DirectedRandomNonObservedMetropolis::setDyadProbability)
		.method("generateSample", &DirectedRandomNonObservedMetropolis::generateSample)
		.method("generateSampleStatistics", &DirectedRandomNonObservedMetropolis::generateSampleStatistics)
		;
*/

	class_<MetropolisHastings<Directed> >("DirectedMetropolisHastings")
			//.constructor<ReModel<Directed> >()
			//.constructor<ReModel<Directed>,double >()
			.constructor()
			.constructor<Model<Directed> >()
			.constructor<Model<Directed>,double >()
			.method("setModel",&MetropolisHastings<Directed>::setModelR)
			.method("getModel",&MetropolisHastings<Directed>::getModelR)
			.method("setDyadToggleType",&MetropolisHastings<Directed>::setDyadToggleType)
			.method("setVertexToggleType",&MetropolisHastings<Directed>::setVertexToggleType)
			.method("setDyadProbability",&MetropolisHastings<Directed>::setDyadProbability)
			.method("generateSample", &MetropolisHastings<Directed>::generateSample)
			.method("generateSampleStatistics", &MetropolisHastings<Directed>::generateSampleStatistics)
			;
	class_<MetropolisHastings<Undirected> >("UndirectedMetropolisHastings")
			//.constructor<ReModel<Undirected> >()
			//.constructor<ReModel<Undirected>,double >()
		    .constructor()
			.constructor<Model<Undirected> >()
			.constructor<Model<Undirected>,double >()
			.method("setModel",&MetropolisHastings<Undirected>::setModelR)
			.method("getModel",&MetropolisHastings<Undirected>::getModelR)
			.method("setDyadToggleType",&MetropolisHastings<Undirected>::setDyadToggleType)
			.method("setVertexToggleType",&MetropolisHastings<Undirected>::setVertexToggleType)
			.method("setDyadProbability",&MetropolisHastings<Undirected>::setDyadProbability)
			.method("generateSample", &MetropolisHastings<Undirected>::generateSample)
			.method("generateSampleStatistics", &MetropolisHastings<Undirected>::generateSampleStatistics)
			;

	class_<CdSampler<Undirected> >("UndirectedCdSampler")
			.constructor<Model<Undirected> >()
			.constructor<Model<Undirected>,int >()
			.method("setModel", &CdSampler<Undirected>::setModel)
			.method("getModel",&CdSampler<Undirected>::getModelR)
			.method("setDyadToggleType",&CdSampler<Undirected>::setDyadToggleType)
			.method("setVertexToggleType",&CdSampler<Undirected>::setVertexToggleType)
			.method("setDyadProbability",&CdSampler<Undirected>::setDyadProbability)
			.method("generateSample", &CdSampler<Undirected>::generateSample)
			.method("generateSampleStatistics", &CdSampler<Undirected>::generateSampleStatistics)
			;

	class_<GibbsCdSampler<Undirected> >("UndirectedGibbsCdSampler")
			.constructor<Model<Undirected> >()
			.constructor<Model<Undirected>,double >()
			.method("setModel", &GibbsCdSampler<Undirected>::setModel)
			.method("getModel",&GibbsCdSampler<Undirected>::getModelR)
			.method("setDyadToggleType",&GibbsCdSampler<Undirected>::setDyadToggleType)
			.method("setVertexToggleType",&GibbsCdSampler<Undirected>::setVertexToggleType)
			.method("setDyadProbability",&GibbsCdSampler<Undirected>::setDyadProbability)
			.method("generateSample", &GibbsCdSampler<Undirected>::generateSample)
			.method("generateSampleStatistics", &GibbsCdSampler<Undirected>::generateSampleStatistics)
			;
	class_<GibbsCdSampler2<Undirected> >("UndirectedGibbsCdSampler2")
			.constructor<Model<Undirected> >()
			.constructor<Model<Undirected>,int >()
			.method("setModel", &GibbsCdSampler2<Undirected>::setModel)
			.method("getModel",&GibbsCdSampler2<Undirected>::getModelR)
			.method("setDyadToggleType",&GibbsCdSampler2<Undirected>::setDyadToggleType)
			.method("setVertexToggleType",&GibbsCdSampler2<Undirected>::setVertexToggleType)
			.method("setDyadProbability",&GibbsCdSampler2<Undirected>::setDyadProbability)
			.method("generateSample", &GibbsCdSampler2<Undirected>::generateSample)
			.method("generateSampleStatistics", &GibbsCdSampler2<Undirected>::generateSampleStatistics)
			;

	class_<Model<Undirected> >("UndirectedModel")
		.constructor()
		.constructor< Model<Undirected> >()
		.method("setNetwork",&Model<Undirected>::setNetworkR)
		.method("getNetwork",&Model<Undirected>::getNetworkR)
		.method("addStatistic",&Model<Undirected>::addStatistic)
		.method("addOffset",&Model<Undirected>::addOffset)
		.method("calculate",&Model<Undirected>::calculate)
		.method("statistics",&Model<Undirected>::statisticsR)
		.method("offset",&Model<Undirected>::offset)
		.method("thetas",&Model<Undirected>::thetasR)
		.method("setThetas",&Model<Undirected>::setThetas)
		.method("setRandomGraph",&Model<Undirected>::setRandomGraph)
		.method("hasRandomGraph",&Model<Undirected>::hasRandomGraph)
		.method("setRandomVariables",&Model<Undirected>::setRandomVariablesR)
		.method("getRandomVariables",&Model<Undirected>::getRandomVariablesR)
		;
	class_<Model<Directed> >("DirectedModel")
		.constructor()
		.constructor< Model<Directed> >()
		.method("setNetwork",&Model<Directed>::setNetworkR)
		.method("getNetwork",&Model<Directed>::getNetworkR)
		.method("addStatistic",&Model<Directed>::addStatistic)
		.method("addOffset",&Model<Directed>::addOffset)
		.method("calculate",&Model<Directed>::calculate)
		.method("statistics",&Model<Directed>::statisticsR)
		.method("offset",&Model<Directed>::offset)
		.method("thetas",&Model<Directed>::thetasR)
		.method("setThetas",&Model<Directed>::setThetas)
		.method("setRandomGraph",&Model<Directed>::setRandomGraph)
		.method("hasRandomGraph",&Model<Directed>::hasRandomGraph)
		.method("setRandomVariables",&Model<Directed>::setRandomVariablesR)
		.method("getRandomVariables",&Model<Directed>::getRandomVariablesR)
		;

	class_<ReModel<Undirected> >("UndirectedReModel")
		.derives< Model<Undirected> >("UndirectedModel")
		.constructor()
		.constructor< ReModel<Undirected> >()
		.method("setBetas",&ReModel<Undirected>::setBetas)
		.method("betas",&ReModel<Undirected>::betaParams)
		.method("setCenters",&ReModel<Undirected>::setCenters)
		.method("centers",&ReModel<Undirected>::centerParams)
		.method("isThetaDependent",&ReModel<Undirected>::isThetaDependent)
		.method("thetaDependent",&ReModel<Undirected>::thetaDependent)
		;

	class_<ReModel<Directed> >("DirectedReModel")
		.derives< Model<Directed> >("DirectedModel")
		.constructor()
		.constructor< ReModel<Directed> >()
		.method("setBetas",&ReModel<Directed>::setBetas)
		.method("betas",&ReModel<Directed>::betaParams)
		.method("setCenters",&ReModel<Directed>::setCenters)
		.method("centers",&ReModel<Directed>::centerParams)
		.method("isThetaDependent",&ReModel<Directed>::isThetaDependent)
		.method("thetaDependent",&ReModel<Directed>::thetaDependent)
		;

	class_<LatentOrderLikelihood<Undirected> >("UndirectedLatentOrderLikelihood")
		.constructor< Model<Undirected> >()
		.method("setModel",&LatentOrderLikelihood<Undirected>::setModel)
		.method("getModel",&LatentOrderLikelihood<Undirected>::getModelR)
		.method("setThetas",&LatentOrderLikelihood<Undirected>::setThetas)
		.method("fullLogLik",&LatentOrderLikelihood<Undirected>::fullLogLik)
		.method("fullLogLikWithFunc",&LatentOrderLikelihood<Undirected>::fullLogLikWithFunc)
		.method("generateNetwork",&LatentOrderLikelihood<Undirected>::generateNetwork)
		;

	class_<LatentOrderLikelihood<Directed> >("DirectedLatentOrderLikelihood")
		.constructor< Model<Directed> >()
		.method("setModel",&LatentOrderLikelihood<Directed>::setModel)
		.method("getModel",&LatentOrderLikelihood<Directed>::getModelR)
		.method("setThetas",&LatentOrderLikelihood<Directed>::setThetas)
		.method("fullLogLik",&LatentOrderLikelihood<Directed>::fullLogLik)
		.method("fullLogLikWithFunc",&LatentOrderLikelihood<Directed>::fullLogLikWithFunc)
		.method("generateNetwork",&LatentOrderLikelihood<Directed>::generateNetwork)
		;

	function("initErnmStatistics",&initStats);

	function("registerDirectedStatistic",&registerDirectedStatistic);
	function("registerUndirectedStatistic",&registerUndirectedStatistic);
	function("registerDirectedStatistic",&registerDirectedOffset);
	function("registerUndirectedStatistic",&registerUndirectedOffset);
	//TODO: have not figured out how to get MakeVars to compile/link the tests dir
	//function(".runCppTests",&tests::runTests);
}                     


#endif

