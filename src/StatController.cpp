#include <StatController.h>
#include <Stat.h>
#include <Stats.h>
#include <Offset.h>
#include <Offsets.h>
#include <Constraint.h>
#include <Constraints.h>
#include <Temporary.h>
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

}

//[[Rcpp::export(name=".initStats")]]
void initStats(){
    using namespace ernm;
	/*
	 * Directed network statistics
	 */
	registerStatistic( DirStatPtr( new DirectedEdges() ) );
	registerStatistic( DirStatPtr( new DirectedTriangles() ) );
	registerStatistic( DirStatPtr( new DirectedReciprocity() ) );
	registerStatistic( DirStatPtr( new DirectedNodeCount() ) );
	registerStatistic( DirStatPtr( new DirectedNodeMatch() ) );
//	registerStatistic( DirStatPtr( new DirectedDegree() ) );
	registerStatistic( DirStatPtr( new DirectedLogistic() ) );
	registerStatistic( DirStatPtr( new DirectedLogisticNeighbors() ) );
	registerStatistic( DirStatPtr( new DirectedDegreeDispersion() ) );
	registerStatistic( DirStatPtr( new DirectedDegreeSkew() ) );
	registerStatistic( DirStatPtr( new DirectedHomophily() ) );
	registerStatistic( DirStatPtr( new DirectedDegree() ) );
	registerStatistic( DirStatPtr( new DirectedStar() ) );
	registerStatistic( DirStatPtr( new DirectedNodeCov() ) );
	registerStatistic( DirStatPtr( new DirectedGwesp() ) );
	registerStatistic( DirStatPtr( new DirectedGeoDist() ) );
    registerStatistic( DirStatPtr( new DirectedGwDegree() ) );
    registerStatistic( DirStatPtr( new DirectedGwdsp() ) );
    registerStatistic( DirStatPtr( new DirectedEsp() ) );
    registerStatistic( DirStatPtr( new DirectedDiffActivity() ) );
    registerStatistic( DirStatPtr( new DirectedSumOfSquares() ) );
    registerStatistic( DirStatPtr( new DirectedGauss() ) );
    registerStatistic( DirStatPtr( new DirectedGamma() ) );
    registerStatistic( DirStatPtr( new DirectedLogDegreeMoment() ) );
    registerStatistic( DirStatPtr( new DirectedHamming() ) );
    registerStatistic( DirStatPtr( new DirectedGaussRegression() ) );
    registerStatistic( DirStatPtr( new DirectedAbsDiff() ) );
    registerStatistic( DirStatPtr( new DirectedRegressNeighbors() ) );
	////////			Offsets				/////////
    //registerOffset( DirOffsetPtr( new DirectedREffectOffset() ) );
	registerOffset( DirOffsetPtr( new DirectedBiasedSeedOffset() ) );

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
	registerStatistic( UndirStatPtr( new UndirectedDegreeChangeCounter() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeMix() ) );
	registerStatistic( UndirStatPtr( new UndirectedTransitivity() ) );
	registerStatistic( UndirStatPtr( new UndirectedHomophily() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeCount() ) );
	registerStatistic( UndirStatPtr( new UndirectedLogistic() ) );
	registerStatistic( UndirStatPtr( new UndirectedLogisticNeighbors() ) );
	registerStatistic( UndirStatPtr( new UndirectedDegree() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeMatch() ) );
	registerStatistic( UndirStatPtr( new UndirectedDegreeDispersion() ) );
	registerStatistic( UndirStatPtr( new UndirectedDegreeSkew() ) );
	registerStatistic( UndirStatPtr( new UndirectedDegreeSpread() ) );
	registerStatistic( UndirStatPtr( new UndirectedDiffActivity() ) );
	registerStatistic( UndirStatPtr( new UndirectedStar() ) );
	registerStatistic( UndirStatPtr( new UndirectedNodeCov() ) );
	registerStatistic( UndirStatPtr( new UndirectedGwesp() ) );
	registerStatistic( UndirStatPtr( new UndirectedGeoDist() ) );
  registerStatistic( UndirStatPtr( new UndirectedGwDegree() ) );
  registerStatistic( UndirStatPtr( new UndirectedGwdsp() ) );
  registerStatistic( UndirStatPtr( new UndirectedEsp() ) );
  registerStatistic( UndirStatPtr( new UndirectedSumOfSquares() ) );
  registerStatistic( UndirStatPtr( new UndirectedGauss() ) );
  registerStatistic( UndirStatPtr( new UndirectedGamma() ) );
  registerStatistic( UndirStatPtr( new UndirectedLogDegreeMoment() ) );
  registerStatistic( UndirStatPtr( new UndirectedDegreeCrossProd() ) );
  registerStatistic( UndirStatPtr( new UndirectedPreferentialAttachment() ) );
  registerStatistic( UndirStatPtr( new UndirectedHamming() ) );
  registerStatistic( UndirStatPtr( new UndirectedGaussRegression() ) );
  registerStatistic( UndirStatPtr( new UndirectedAbsDiff() ) );
  registerStatistic( UndirStatPtr( new UndirectedRegressNeighbors() ) );

	////////			Offsets				/////////
    registerOffset( UndirOffsetPtr( new UndirectedREffectOffset() ) );
	registerOffset( UndirOffsetPtr( new UndirectedBiasedSeedOffset() ) );
	registerOffset( UndirOffsetPtr( new UndirectedRdsBiasOffset() ) );
	registerOffset( UndirOffsetPtr( new UndirectedBoundedDegreeConstraint() ) );
	registerOffset( UndirOffsetPtr( new UndirectedFixedNodeConstraint() ) );
	registerOffset( UndirOffsetPtr( new UndirectedFixedDegreeConstraint() ) );
    registerOffset( UndirOffsetPtr( new UndirectedStarPenalty() ) );
    registerOffset( UndirOffsetPtr( new UndirectedHammingOffset() ) );
	//Make registration available outside ernm compilation unit
	R_RegisterCCallable("ernm",
			"registerUndirectedStatistic",(DL_FUNC) &registerUndirectedStatistic);
	R_RegisterCCallable("ernm",
			"registerUndirectedOffset",(DL_FUNC) &registerUndirectedOffset);

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






