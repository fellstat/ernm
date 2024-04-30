


#ifndef STATCONTROLLERH_
#define STATCONTROLLERH_

#include <string>
#include <map>
#include "Stat.h"
#include "Offset.h"

/*
 * Use these macros to register statistics compiled code
 * outside the ernm compilation unit. i.e. from a different package,
 * or from inline cpp code
 */
#define REGISTER_UNDIRECTED_STATISTIC(x) ((void(*)(Rcpp::XPtr< ernm::AbstractStat<ernm::Undirected> >))R_GetCCallable("ernm", "registerUndirectedStatistic"))(x)
#define REGISTER_DIRECTED_STATISTIC(x) ((void(*)(Rcpp::XPtr< ernm::AbstractStat<ernm::Directed> >))R_GetCCallable("ernm", "registerDirectedStatistic"))(x)
#define REGISTER_UNDIRECTED_OFFSET(x) ((void(*)(Rcpp::XPtr< ernm::AbstractOffset<ernm::Undirected> >))R_GetCCallable("ernm", "registerUndirectedOffset"))(x)
#define REGISTER_DIRECTED_OFFSET(x) ((void(*)(Rcpp::XPtr< ernm::AbstractOffset<ernm::Directed> >))R_GetCCallable("ernm", "registerDirectedOffset"))(x)

namespace ernm{

/*!
 * This class controls the stats available for all network types to models
 * on the R side.
 *
 * to add a new user statistic just call registerStatistic
 */
template<class Engine>
class StatController {
private:
	typedef boost::shared_ptr< AbstractStat<Engine> > StatPtr;
	typedef boost::shared_ptr< std::map< std::string, StatPtr > > StatMapPtr;
	typedef boost::shared_ptr< AbstractOffset<Engine> > OffsetPtr;
	typedef boost::shared_ptr< std::map< std::string, OffsetPtr > > OffsetMapPtr;
protected:
	static StatMapPtr statMapPtr;
	static OffsetMapPtr offsetMapPtr;
public:
	StatController(){}
	~StatController(){}
	static void init(){
		if(!statMapPtr)
			statMapPtr = StatMapPtr(new std::map< std::string, StatPtr >);
		if(!offsetMapPtr)
			offsetMapPtr = OffsetMapPtr(new std::map< std::string, OffsetPtr >);
	}
	static void addStat(StatPtr pS){
		init();
		statMapPtr->insert(std::make_pair(pS->vName(),pS));
	}
	static void addOffset(OffsetPtr pS){
		init();
		offsetMapPtr->insert(std::make_pair(pS->vName(),pS));
	}

	static AbstractStat<Engine>* getStat(std::string name, Rcpp::List params){
		StatPtr pS;
		try{
			pS = statMapPtr->at(name);

		}catch(...){
			Rf_error((std::string("Unknown statistic: ") + name).c_str());
			return NULL;
		}
		if(pS==NULL){
			Rf_error((std::string("Unknown statistic: ") + name).c_str());
			return NULL;
		}
		return pS->vCreateUnsafe(params);
	}

	static AbstractOffset<Engine>* getOffset(std::string name, Rcpp::List params){
		OffsetPtr pS;
		try{
			pS = offsetMapPtr->at(name);

		}catch(...){
			Rf_error((std::string("Unknown offset: ") + name).c_str());
			return NULL;
		}
		if(pS==NULL){
			Rf_error((std::string("Unknown offset: ") + name).c_str());
			return NULL;
		}
		return pS->vCreateUnsafe(params);
	}
};





template<class Engine>
void registerStatistic(boost::shared_ptr< AbstractStat<Engine> > ps){
	StatController<Engine>::addStat(ps);
}

template<class Engine>
void registerOffset(boost::shared_ptr< AbstractOffset<Engine> > ps){
	StatController<Engine>::addOffset(ps);
}

/*!
 * Called upon loading ERNM, registering default stats/offsets
 */
void initStats();

}

void registerDirectedStatistic(Rcpp::XPtr< ernm::AbstractStat<ernm::Directed> > ps);

void registerUndirectedStatistic(Rcpp::XPtr< ernm::AbstractStat<ernm::Undirected> > ps);

void registerDirectedOffset(Rcpp::XPtr< ernm::AbstractOffset<ernm::Directed> > ps);

void registerUndirectedOffset(Rcpp::XPtr< ernm::AbstractOffset<ernm::Undirected> > ps);


#endif /* STATCONTROLLERH_ */
