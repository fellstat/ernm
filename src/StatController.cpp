#include "StatController.h"
#include "Stat.h"
#include "Stats.h"
#include "Offset.h"
#include "Offsets.h"
#include "Constraint.h"
#include "Constraints.h"
namespace ernm{



typedef boost::shared_ptr< AbstractStat<Directed> > DirStatPtr;
typedef boost::shared_ptr< std::map< std::string, DirStatPtr > > StatMapPtr;
typedef boost::shared_ptr< AbstractOffset<Directed> > DirOffsetPtr;
typedef boost::shared_ptr< std::map< std::string, DirOffsetPtr > > OffsetMapPtr;
template<> StatMapPtr StatController<Directed>::statMapPtr = StatMapPtr(new std::map< std::string, DirStatPtr >);
template<> OffsetMapPtr StatController<Directed>::offsetMapPtr =
		OffsetMapPtr(new std::map< std::string, DirOffsetPtr >);

typedef boost::shared_ptr< AbstractStat<Undirected> > UndirStatPtr;
typedef boost::shared_ptr< std::map< std::string, UndirStatPtr > > UndirStatMapPtr;
typedef boost::shared_ptr< AbstractOffset<Undirected> > UndirOffsetPtr;
typedef boost::shared_ptr< std::map< std::string, UndirOffsetPtr > > UndirOffsetMapPtr;
template<> UndirStatMapPtr StatController<Undirected>::statMapPtr =
		UndirStatMapPtr(new std::map< std::string, UndirStatPtr >);
template<> UndirOffsetMapPtr StatController<Undirected>::offsetMapPtr =
		UndirOffsetMapPtr(new std::map< std::string, UndirOffsetPtr >);


RcppExport void initStats(){
	/*
	 * Directed network statistics
	 */
	registerStatistic( DirStatPtr( new DirectedEdges() ) );
	registerStatistic( DirStatPtr( new DirectedTriangles() ) );
	registerStatistic( DirStatPtr( new DirectedReciprocity() ) );
	registerStatistic( DirStatPtr( new DirectedNodeMatch() ) );
	registerStatistic( DirStatPtr( new DirectedDegree() ) );
	registerStatistic( DirStatPtr( new DirectedStar() ) );
	registerStatistic( DirStatPtr( new DirectedNodeCov() ) );
	registerStatistic( DirStatPtr( new DirectedGwesp() ) );
	registerStatistic( DirStatPtr( new DirectedGeoDist() ) );
    registerStatistic( DirStatPtr( new DirectedGwDegree() ) );
    registerStatistic( DirStatPtr( new DirectedGwdsp() ) );
    registerStatistic( DirStatPtr( new DirectedEsp() ) );
	////////			Offsets				/////////

	//Make registration available outside ernm compilation unit
	R_RegisterCCallable("ernm",
			"registerDirectedStatistic",(DL_FUNC) &registerDirectedStatistic);
	R_RegisterCCallable("ernm",
			"registerDirectedOffset",(DL_FUNC) &registerDirectedOffset);

	/*
	 * Undirected network statistics
	 */
	registerStatistic( UndirStatPtr( new UndirectedEdges() ) );
	registerStatistic( UndirStatPtr( new UndirectedTriangles() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeMix() ) );
	registerStatistic( UndirStatPtr( new UndirectedDegree() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeMatch() ) );
	registerStatistic( UndirStatPtr( new UndirectedStar() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeCov() ) );
	registerStatistic( UndirStatPtr( new UndirectedGwesp() ) );
	registerStatistic( UndirStatPtr( new UndirectedGeoDist() ) );
    registerStatistic( UndirStatPtr( new UndirectedGwDegree() ) );
    registerStatistic( UndirStatPtr( new UndirectedGwdsp() ) );
    registerStatistic( UndirStatPtr( new UndirectedEsp() ) );
    registerStatistic( UndirStatPtr( new UndirectedDegreeCrossProd() ) );
    registerStatistic( UndirStatPtr( new UndirectedPreferentialAttachment() ) );
    registerStatistic( UndirStatPtr( new UndirectedSharedNbrs() ) );

	////////			Offsets				/////////
	registerOffset( UndirOffsetPtr( new UndirectedBoundedDegreeConstraint() ) );
	//Make registration available outside ernm compilation unit
	R_RegisterCCallable("ernm",
			"registerUndirectedStatistic",(DL_FUNC) &registerUndirectedStatistic);
	R_RegisterCCallable("ernm",
			"registerUndirectedOffset",(DL_FUNC) &registerUndirectedOffset);

}


}

void registerDirectedStatistic(Rcpp::XPtr< ernm::AbstractStat<ernm::Directed> > ps){
	ernm::StatController<ernm::Directed>::addStat(
			boost::shared_ptr< ernm::AbstractStat<ernm::Directed> >(ps->vCloneUnsafe()));
}

void registerUndirectedStatistic(Rcpp::XPtr< ernm::AbstractStat<ernm::Undirected> > ps){
	ernm::StatController<ernm::Undirected>::addStat(
			boost::shared_ptr< ernm::AbstractStat<ernm::Undirected> >(ps->vCloneUnsafe()));
}

void registerDirectedOffset(Rcpp::XPtr< ernm::AbstractOffset<ernm::Directed> > ps){
	ernm::StatController<ernm::Directed>::addOffset(
			boost::shared_ptr< ernm::AbstractOffset<ernm::Directed> >(ps->vCloneUnsafe()));
}

void registerUndirectedOffset(Rcpp::XPtr< ernm::AbstractOffset<ernm::Undirected> > ps){
	ernm::StatController<ernm::Undirected>::addOffset(
			boost::shared_ptr< ernm::AbstractOffset<ernm::Undirected> >(ps->vCloneUnsafe()));
}






