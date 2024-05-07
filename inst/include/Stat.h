/*
 * Stat.h
 *
 *  Created on: Jun 11, 2011
 *      Author: ianfellows
 */

#ifndef STATH_
#define STATH_

#include <vector>
#include <string>
#include "BinaryNet.h"
#include <RcppCommon.h>
#include <math.h>
#include "util.h"
#include "Offset.h"

namespace ernm{

using namespace Rcpp;

template<class Engine>
class AbstractStat{
public:
	virtual ~AbstractStat(){};

	/*!
	 * create a new statistic with an arbitrary set of parameters.
	 *
	 * this allows the Stat to act as a factory. see StatController.
	 *
	 * \param params a set of defining characteristics of the statistic.
	 */
	virtual AbstractStat* vCreateUnsafe(Rcpp::List params) const = 0;

	/*!
	 * \return the name of the statistic
	 */
	virtual std::string vName() = 0;

	/*!
	 * \return an identical un-aliased version of the Stat
	 */
	virtual boost::shared_ptr< AbstractStat<Engine> > vClone() = 0;

	virtual AbstractStat* vCloneUnsafe() = 0;


	/*!
	 * calculate the statistic based on the supplied network
	 */
	virtual void vCalculate(const BinaryNet<Engine>& net) = 0;

	/*!
	 * update statistics with a hypothetical edge toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses calculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param from toggled edge (from)
	 * \param to toggled edge (to)
	 */
	virtual void vDyadUpdate(const BinaryNet<Engine>& net, int from, int to) = 0;
	
	/*!
	 * calculate the change in the statistics from a hypothetical vertex toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses calculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param vert the index of the vertex change
	 * \param variable the id of the variable
	 * \param newValue the hypothetical new value
	 */
	virtual void vDiscreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue) = 0;

	/*!
	 * calculate the change in the statistics from a hypothetical vertex toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses calculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param vert the index of the vertex change
	 * \param variable the id of the variable
	 * \param newValue the hypothetical new value
	 */
	virtual void vContinVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, double newValue) = 0;

	/*!
	 * \return names for the statistics
	 */
	virtual std::vector<std::string> vStatNames() = 0;

	/*!
	 * number of statistics
	 */
	virtual int vSize() = 0;

	/*!
	 * returns the models statistics
	 */
	virtual std::vector<double>& vStatistics() = 0;

	/*!
	 * set the model parameter values
	 */
	virtual void vSetStatistics(const std::vector<double>&st) = 0;

	/*!
	 * set the model parameter values
	 */
	virtual void vSetTheta(const std::vector<double>&th) = 0;

	/*!
	 * the model parameter values
	 */
	virtual std::vector<double>& vTheta() = 0;

	/*!
	 * \return the terms
	 */
	virtual std::vector<double> vValues() = 0;

	/*!
	 * \return stats * thetas
	 */
	virtual double vLogLik() = 0;
	
	/*!
	 * tells you whether it is safe to used cached stats rather than rerunning dyadUpdate
	 */
	virtual bool vGetDyadUpdateSafe() =0;
	virtual bool vGetDiscreteVertexUpdateSafe() =0;
	virtual bool vGetContinVertexUpdateSafe() =0;

};

/*!
 * templated interface for model statistics.
 *
 * Allows for access via virtual functions or non-virtual functions,
 * which will be usefull if we ever want to do template meta programming.
 *
 * The StatEngine just needs to implement the non-virtual functions.
 *
 */
template<class NetworkEngine, class StatEngine>
class Stat : public AbstractStat<NetworkEngine>{
protected:
	StatEngine stat;

public:
	Stat() : stat(){}

	Stat(Rcpp::List params) : stat(params){}

	virtual ~Stat() {}

	/*!
	 * create a new statistic with an arbitrary set of parameters.
	 *
	 * this allows the Stat to act as a factory. see StatController.
	 *
	 * \param params a set of defining characteristics of the statistic.
	 */
	virtual AbstractStat<NetworkEngine>* vCreateUnsafe(Rcpp::List params) const{
		return createUnsafe(params);
	}

	inline AbstractStat<NetworkEngine>* createUnsafe(Rcpp::List params) const{
		return new Stat(params);
	}

	/*!
	 * \return the name of the statistic
	 */
	virtual std::string vName(){
		return name();
	}

	inline std::string name(){
		return stat.name();
	}

	/*!
	 * \return an identical un-aliased version of the Stat
	 */
	virtual boost::shared_ptr< AbstractStat<NetworkEngine> > vClone(){
		return clone();
	}

	inline boost::shared_ptr< AbstractStat<NetworkEngine> > clone(){
		return boost::shared_ptr< AbstractStat<NetworkEngine> >(cloneUnsafe());
	}

	virtual AbstractStat<NetworkEngine>* vCloneUnsafe(){
			return cloneUnsafe();
	}

	inline AbstractStat<NetworkEngine>* cloneUnsafe(){
		return new Stat<NetworkEngine,StatEngine>(*this);
	}

	/*!
	 * calculate the statistic based on the supplied network
	 */
	virtual void vCalculate(const BinaryNet<NetworkEngine>& net){
		calculate(net);
	}

	inline void calculate(const BinaryNet<NetworkEngine>& net){
		stat.calculate(net);
	}

	/*!
	 * update statistics with a hypothetical edge toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses calculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param from toggled edge (from)
	 * \param to toggled edge (to)
	 */
	virtual void vDyadUpdate(const BinaryNet<NetworkEngine>& net, int from, int to){
		dyadUpdate(net,from,to);
	}

	inline void dyadUpdate(const BinaryNet<NetworkEngine>& net, int from, int to){
		stat.dyadUpdate(net,from,to);
	}

	/*!
	 * calculate the change in the statistics from a hypothetical vertex toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses calculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param vert the index of the vertex change
	 * \param variable the id of the variable
	 * \param newValue the hypothetical new value
	 */
	virtual void vDiscreteVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, int newValue){
		discreteVertexUpdate(net,vert,variable,newValue);
	}

	inline void discreteVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, int newValue){
		stat.discreteVertexUpdate(net,vert,variable,newValue);
	}

	/*!
	 * calculate the change in the statistics from a hypothetical vertex toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses calculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param vert the index of the vertex change
	 * \param variable the id of the variable
	 * \param newValue the hypothetical new value
	 */
	virtual void vContinVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, double newValue){
		continVertexUpdate(net,vert,variable,newValue);
	}

	inline void continVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, double newValue){
		stat.continVertexUpdate(net,vert,variable,newValue);
	}

	/*!
	 * \return names for the statistics
	 */
	virtual std::vector<std::string> vStatNames(){
		return statNames();
	}

	inline std::vector<std::string> statNames(){
		std::vector<std::string> names = stat.statNames();
		if(names.size() == 0){
			names = std::vector<std::string>(this->vSize());
			std::vector<std::string> nm = this->vStatNames();
			/*
            if(this->vSize()==1){
				names[0] = nm[0];
				return names;
			}
			for(int i=1;i<=this->vSize();i++){
				names[i-1] = nm[0] + asString(i);
			}
            */
		}
		return names;
	}

	/*!
	 * number of statistics
	 */
	virtual int vSize(){
		return size();
	}

	inline int size(){
		return stat.size();
	}

	/*!
	 * returns the models statistics
	 */
	virtual std::vector<double>& vStatistics(){
		return statistics();
	}

	inline std::vector<double>& statistics(){
		return stat.statistics();
	}

	/*!
	 * set the model parameter values
	 */
	virtual void vSetStatistics(const std::vector<double>&st){
		setStatistics(st);
	}

	inline 	void setStatistics(const std::vector<double>&st){
		stat.setStatistics(st);
	}

	/*!
	 * set the model parameter values
	 */
	virtual void vSetTheta(const std::vector<double>&th){
		setTheta(th);
	}

	inline void setTheta(const std::vector<double>&th){
		stat.setTheta(th);
	}

	/*!
	 * the model parameter values
	 */
	virtual std::vector<double>& vTheta(){
		return theta();
	}

	inline std::vector<double>& theta(){
		return stat.theta();
	}


	/*!
	 * \return the terms theta * stats
	 */
	virtual std::vector<double> vValues(){
		return values();
	}

	std::vector<double> values(){
		return stat.values();
	}

	/*!
	 * \return stats * thetas
	 */
	virtual double vLogLik(){
		return logLik();
	}

	inline double logLik(){
		return stat.logLik();
	}
    
    /*!
     * update Cache safety methods
     */
    virtual bool vGetDyadUpdateSafe(){
        return getDyadUpdateSafe();
    }
    
    inline bool getDyadUpdateSafe(){
        return stat.getDyadUpdateSafe();
    }
    
    virtual bool vGetDiscreteVertexUpdateSafe(){
        return getDiscreteVertexUpdateSafe();
    }
    
    inline bool getDiscreteVertexUpdateSafe(){
        return stat.getDiscreteVertexUpdateSafe();
    }
    
    virtual bool vGetContinVertexUpdateSafe(){
        return getContinVertexUpdateSafe();
    }
    
    inline bool getContinVertexUpdateSafe(){
        return stat.getContinVertexUpdateSafe();
    }

};


/*!
 * a class representing model statistics.
 *
 * Extend this to create your own ernm statistic.
 *
 * Subclasses must either extend vCalculate to calculate the statistic in
 * this case each toggle will result in a call to vCalculate, or
 * implement calculate, dyadUpdate and discreteVertexUpdate which allows
 * for much faster toggling.
 */
template<class Engine>
class BaseStat : public BaseOffset<Engine>{
protected:
	std::vector<double> thetas;/*!< the parameter values */
    
public:

	BaseStat(){}

	virtual ~BaseStat(){}

	/*!
	 * set the model parameter values
	 */
	void setTheta(const std::vector<double>&th){
		thetas=th;
	}

	/*!
	 * the model parameter values
	 */
	std::vector<double>& theta(){
		return this->thetas;
	}

	/*!
	 * \return the terms
	 */
	std::vector<double> values(){
		std::vector<double> v(this->stats.size());
		for(int i=0;i<this->stats.size();i++)
			v[i] = this->stats[i]*this->thetas[i];
		return v;
	}

	/*!
	 * \return stats * thetas
	 */
	double logLik(){
		double ll = 0.0;
		for(int i=0;i<this->stats.size();i++)
			ll += this->stats[i]*this->thetas[i];
		return ll;
	}
};

}

#endif /* STATH_ */
