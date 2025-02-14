
#ifndef INSIDE

#include <Rcpp.h>
#include <BinaryNet.h>
#include <DyadToggles.h>
#include <VertexToggles.h>
#include <MetropolisHastings.h>
#include <CdSampler.h>
#include <tests.h>

/*
 * Handles all functions and methods exported to R.
 */

RCPP_MODULE(ernm){
	using namespace Rcpp;
	using namespace ernm;

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
			.method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR3);

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
			.method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR3);
	
	class_<MetropolisHastings<Directed> >("DirectedMetropolisHastings")
			.constructor()
			.constructor<Model<Directed> >()
			.constructor<Model<Directed>,double >()
			.method("setModel",&MetropolisHastings<Directed>::setModelR)
			.method("getModel",&MetropolisHastings<Directed>::getModelR)
			.method("setDyadToggleType",&MetropolisHastings<Directed>::setDyadToggleType)
			.method("setVertexToggleType",&MetropolisHastings<Directed>::setVertexToggleType)
			.method("setDyadProbability",&MetropolisHastings<Directed>::setDyadProbability)
			.method("generateSample", &MetropolisHastings<Directed>::generateSample)
			.method("generateSampleStatistics", &MetropolisHastings<Directed>::generateSampleStatistics);
	
	class_<MetropolisHastings<Undirected> >("UndirectedMetropolisHastings")
		    .constructor()
			.constructor<Model<Undirected> >()
			.constructor<Model<Undirected>,double >()
			.method("setModel",&MetropolisHastings<Undirected>::setModelR)
			.method("getModel",&MetropolisHastings<Undirected>::getModelR)
			.method("setDyadToggleType",&MetropolisHastings<Undirected>::setDyadToggleType)
			.method("setVertexToggleType",&MetropolisHastings<Undirected>::setVertexToggleType)
			.method("setDyadProbability",&MetropolisHastings<Undirected>::setDyadProbability)
			.method("generateSample", &MetropolisHastings<Undirected>::generateSample)
			.method("generateSampleStatistics", &MetropolisHastings<Undirected>::generateSampleStatistics);

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
        .method("dyadUpdate",&Model<Undirected>::dyadUpdateR)
        .method("discreteVertexUpdate",&Model<Undirected>::discreteVertexUpdateR)
        .method("continVertexUpdate",&Model<Undirected>::continVertexUpdateR)
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
        .method("dyadUpdate",&Model<Directed>::dyadUpdateR)
        .method("discreteVertexUpdate",&Model<Directed>::discreteVertexUpdateR)
        .method("continVertexUpdate",&Model<Directed>::continVertexUpdateR)
		;

	class_<TaperedModel<Undirected> >("UndirectedTaperedModel")
		.derives< Model<Undirected> >("UndirectedModel")
		.constructor()
		.constructor< TaperedModel<Undirected> >()
		.method("setTau",&TaperedModel<Undirected>::setTau)
		.method("tau",&TaperedModel<Undirected>::tauParams)
		.method("setCenters",&TaperedModel<Undirected>::setCenters)
		.method("centers",&TaperedModel<Undirected>::centerParams)
		;

	class_<TaperedModel<Directed> >("DirectedTaperedModel")
		.derives< Model<Directed> >("DirectedModel")
		.constructor()
		.constructor< TaperedModel<Directed> >()
		.method("setTau",&TaperedModel<Directed>::setTau)
		.method("tau",&TaperedModel<Directed>::tauParams)
		.method("setCenters",&TaperedModel<Directed>::setCenters)
		.method("centers",&TaperedModel<Directed>::centerParams)
		;
    
    // functions to register statistics with the controller
	function("registerDirectedStatistic",&registerDirectedStatistic);
	function("registerUndirectedStatistic",&registerUndirectedStatistic);
	function("registerDirectedStatistic",&registerDirectedOffset);
	function("registerUndirectedStatistic",&registerUndirectedOffset);
    // test functions to be exposed and then run by testthat
	function("runErnmCppTests",&ernm::tests::runErnmTests);
}                     


#endif

