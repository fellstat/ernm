/*
 * ToggleController.cpp
 *
 *  Created on: Jan 3, 2014
 *      Author: ianfellows
 */

#include <ToggleController.h>
#include <DyadToggle.h>
#include <VertexToggle.h>
#include <DyadToggles.h>
#include <VertexToggles.h>
#include <CdSampler.h>

namespace ernm {

typedef boost::shared_ptr< AbstractDyadToggle<Directed> > DirDyadTogglePtr;
typedef boost::shared_ptr< std::map< std::string, DirDyadTogglePtr > > DtMapPtr;
typedef boost::shared_ptr< AbstractVertexToggle<Directed> > DirVertexTogglePtr;
typedef boost::shared_ptr< std::map< std::string, DirVertexTogglePtr > > VtMapPtr;
template<> DtMapPtr ToggleController<Directed>::dyadMapPtr =
		DtMapPtr(new std::map< std::string, DirDyadTogglePtr >);
template<> VtMapPtr ToggleController<Directed>::vertexMapPtr =
		VtMapPtr(new std::map< std::string, DirVertexTogglePtr >);


typedef boost::shared_ptr< AbstractDyadToggle<Undirected> > UndirDyadTogglePtr;
typedef boost::shared_ptr< std::map< std::string, UndirDyadTogglePtr > > UnDtMapPtr;
typedef boost::shared_ptr< AbstractVertexToggle<Undirected> > UndirVertexTogglePtr;
typedef boost::shared_ptr< std::map< std::string, UndirVertexTogglePtr > > UnVtMapPtr;
template<> UnDtMapPtr ToggleController<Undirected>::dyadMapPtr =
		UnDtMapPtr(new std::map< std::string, UndirDyadTogglePtr >);
template<> UnVtMapPtr ToggleController<Undirected>::vertexMapPtr =
		UnVtMapPtr(new std::map< std::string, UndirVertexTogglePtr >);
}

//[[Rcpp::export(name=".initToggles")]]
void initToggles(){
	using namespace ernm;
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedRandomDyadToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedTieDyadToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedNeighborhoodToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedNodeTieDyadToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedTetradToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedRdsToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedDefaultCdToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedNtdNbrToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedRandomDyadMissingToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedNodeTieDyadMissingToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedNeighborhoodMissingToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Directed> >(new DirectedNtdNbrMissingToggle()));

	registerToggle(boost::shared_ptr< AbstractVertexToggle<Directed> >(new DirectedDefaultVertexToggle()));
	registerToggle(boost::shared_ptr< AbstractVertexToggle<Directed> >(new DirectedVertexMissingToggle()));
	//Make registration available outside ernm compilation unit
	R_RegisterCCallable("ernm",
			"registerDirectedDyadToggle",(DL_FUNC) &registerDirectedDyadToggle);
	R_RegisterCCallable("ernm",
			"registerDirectedVertexToggle",(DL_FUNC) &registerDirectedVertexToggle);

	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedRandomDyadToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedTieDyadToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedNeighborhoodToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedNodeTieDyadToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedTetradToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedRdsToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedDefaultCdToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedNtdNbrToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedRandomDyadMissingToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedNodeTieDyadMissingToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedNeighborhoodMissingToggle()));
	registerToggle(boost::shared_ptr< AbstractDyadToggle<Undirected> >(new UndirectedNtdNbrMissingToggle()));

	registerToggle(boost::shared_ptr< AbstractVertexToggle<Undirected> >(new UndirectedDefaultVertexToggle()));
	registerToggle(boost::shared_ptr< AbstractVertexToggle<Undirected> >(new UndirectedVertexMissingToggle()));
	//Make registration available outside ernm compilation unit
	R_RegisterCCallable("ernm",
			"registerUndirectedDyadToggle",(DL_FUNC) &registerUndirectedDyadToggle);
	R_RegisterCCallable("ernm",
			"registerUndirectedVertexToggle",(DL_FUNC) &registerUndirectedVertexToggle);


}

void registerDirectedDyadToggle(Rcpp::XPtr< ernm::AbstractDyadToggle<ernm::Directed> > ps){
	ernm::ToggleController<ernm::Directed>::addToggle(
			boost::shared_ptr< ernm::AbstractDyadToggle<ernm::Directed> >(ps->vCloneUnsafe()));
}
void registerDirectedVertexToggle(Rcpp::XPtr< ernm::AbstractVertexToggle<ernm::Directed> > ps){
	ernm::ToggleController<ernm::Directed>::addToggle(
			boost::shared_ptr< ernm::AbstractVertexToggle<ernm::Directed> >(ps->vCloneUnsafe()));
}

void registerUndirectedDyadToggle(Rcpp::XPtr< ernm::AbstractDyadToggle<ernm::Undirected> > ps){
	ernm::ToggleController<ernm::Undirected>::addToggle(
			boost::shared_ptr< ernm::AbstractDyadToggle<ernm::Undirected> >(ps->vCloneUnsafe()));
}
void registerUndirectedVertexToggle(Rcpp::XPtr< ernm::AbstractVertexToggle<ernm::Undirected> > ps){
	ernm::ToggleController<ernm::Undirected>::addToggle(
			boost::shared_ptr< ernm::AbstractVertexToggle<ernm::Undirected> >(ps->vCloneUnsafe()));
}




