/*
 * Offset.h
 *
 *  Created on: Mar 10, 2012
 *      Author: ianfellows
 */

#ifndef OFFSETH_
#define OFFSETH_

#include <string>
#include "BinaryNet.h"
#include <memory>
#include <boost/shared_ptr.hpp>
#include <math.h>
#include <float.h>

namespace ernm {



template<class Engine>
class AbstractOffset{
public:
	virtual ~AbstractOffset(){};

	/*!
	 * create a new offset with an arbitrary set of parameters.
	 *
	 * this allows the Stat to act as a factory. see StatController.
	 *
	 * \param params a set of defining characteristics of the statistic.
	 */
	virtual AbstractOffset* vCreateUnsafe(Rcpp::List params) const = 0;

	/*!
	 * \return the name of the statistic
	 */
	virtual std::string vName() = 0;

	/*!
	 * \return an identical un-aliased version of the Stat
	 */
	virtual boost::shared_ptr< AbstractOffset<Engine> > vClone() = 0;
	virtual AbstractOffset* vCloneUnsafe() = 0;


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
	 * number of statistics
	 */
	virtual int vSize() = 0;

	/*!
	 * \return the terms
	 */
	virtual std::vector<double> vValues() = 0;

	/*!
	 * \return stats * thetas
	 */
	virtual double vLogLik() = 0;

};


/*!
 * templated interface for model statistics.
 *
 * Allows for access via virtual functions or non-virtual functions,
 * which will be useful if we ever want to do template meta programming.
 *
 * The StatEngine just needs to implement the non-virtual functions.
 *
 */
template<class NetworkEngine, class OffsetEngine>
class Offset : public AbstractOffset<NetworkEngine>{
protected:
	OffsetEngine off;

public:
	Offset() : off(){}

	Offset(Rcpp::List params) : off(params){}

	virtual ~Offset() {}

	/*!
	 * create a new statistic with an arbitrary set of parameters.
	 *
	 * this allows the Stat to act as a factory. see StatController.
	 *
	 * \param params a set of defining characteristics of the statistic.
	 */
	virtual AbstractOffset<NetworkEngine>* vCreateUnsafe(Rcpp::List params) const{
		return createUnsafe(params);
	}

	inline AbstractOffset<NetworkEngine>* createUnsafe(Rcpp::List params) const{
		return new Offset(params);
	}

	/*!
	 * \return the name of the statistic
	 */
	virtual std::string vName(){
		return name();
	}

	inline std::string name(){
		return off.name();
	}

	/*!
	 * \return an identical un-aliased version of the Stat
	 */
	virtual boost::shared_ptr< AbstractOffset<NetworkEngine> > vClone(){
		return clone();
	}

	inline boost::shared_ptr< AbstractOffset<NetworkEngine> > clone(){
		return boost::shared_ptr< AbstractOffset<NetworkEngine> >(cloneUnsafe());
	}

	virtual AbstractOffset<NetworkEngine>* vCloneUnsafe(){
			return cloneUnsafe();
	}

	inline AbstractOffset<NetworkEngine>* cloneUnsafe(){
		return new Offset<NetworkEngine,OffsetEngine>(*this);
	}

	/*!
	 * calculate the statistic based on the supplied network
	 */
	virtual void vCalculate(const BinaryNet<NetworkEngine>& net){
		calculate(net);
	}

	inline void calculate(const BinaryNet<NetworkEngine>& net){
		off.calculate(net);
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
		off.dyadUpdate(net,from,to);
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
		off.discreteVertexUpdate(net,vert,variable,newValue);
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
		off.continVertexUpdate(net,vert,variable,newValue);
	}

	/*!
	 * number of statistics
	 */
	virtual int vSize(){
		return size();
	}

	inline int size(){
		return off.size();
	}

	/*!
	 * \return the terms theta * stats
	 */
	virtual std::vector<double> vValues(){
		return values();
	}

	std::vector<double> values(){
		return off.values();
	}

	/*!
	 * \return stats * thetas
	 */
	virtual double vLogLik(){
		return logLik();
	}

	inline double logLik(){
		return off.logLik();
	}

};


/*!
 * Offsets are statistics with no assosiated parameter.
 *
 * extend this class to create a new offset.
 */
template<class Engine>
class BaseOffset {
protected:
	std::vector<double> stats; /*!< the statistics */

	/*!
	 * calculate the offset based on the supplied network
	 *
	 * This function is good for prototyping statistics. Simply have the subclass implement
	 * it, and the new stat can be used without the need to program dyad and vertex change update
	 * functions.
	 */
	virtual void vCalculate(const BinaryNet<Engine>& net){
		Rf_error("Either calculate, dyadUpdate, etc. _OR_ vCalculate must be implemented in the subclass.");
	}

public:

	BaseOffset() {};

	virtual ~BaseOffset(){};

	std::string name(){
		return "baseOffset";
	}

	/*!
	 * calculate the statistic based on the supplied network
	 *
	 * (Either this or vCalculate must be implemented by deriving classes)
	 */
	void calculate(const BinaryNet<Engine>& net){
		vCalculate(net);
	}

	/*!
	 * update offset with a hypothetical edge toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses vCalculate to compute the changes, but can be overridden to
	 * get speed gains.
	 *
	 * \param net the network
	 * \param from toggled edge (from)
	 * \param to toggled edge (to)
	 */
	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		BinaryNet<Engine>* pnet = const_cast< BinaryNet<Engine>* > (&net);
		std::vector<double> tmp = stats;
		pnet->toggle(from,to);
		this->calculate(net);
		pnet->toggle(from,to);

	}

	/*!
	 * calculate the change in the offset from a hypothetical vertex toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses vCalculate to compute the changes, but can be overridden to
	 * get speed gains.
	 *
	 * \param net the network
	 * \param vert the index of the vertex change
	 * \param variable the id of the variable
	 * \param newValue the hypothetical new value
	 */
	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
		BinaryNet<Engine>* pnet = const_cast< BinaryNet<Engine>* > (&net);
		int oldValue = net.discreteVariableValue(variable,vert);
		pnet->setDiscreteVariableValue(variable,vert,newValue);
		this->calculate(net);
		pnet->setDiscreteVariableValue(variable,vert,oldValue);
	}

	/*!
	 * calculate the change in the statistics from a hypothetical vertex toggle,
	 * assuming that the network has not changed since the statistic was last calculated.
	 *
	 * by default this uses vCalculate to compute the changes, but can be overridden to
	 * get speed gains
	 *
	 * \param net the network
	 * \param vert the index of the vertex change
	 * \param variable the id of the variable
	 * \param newValue the hypothetical new value
	 */
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, double newValue){
		BinaryNet<Engine>* pnet = const_cast< BinaryNet<Engine>* > (&net);
		double oldValue = net.continVariableValue(variable,vert);
		pnet->setContinVariableValue(variable,vert,newValue);
		this->calculate(net);
		pnet->setContinVariableValue(variable,vert,oldValue);
	}

	/*!
	 * number of statistics
	 */
	int size(){
		return stats.size();
	}

	/*!
	 * \return names for the statistics
	 */
	std::vector<std::string> statNames(){
		return std::vector<std::string>();
	}

	/*!
	 * returns the models statistics
	 */
	std::vector<double>& statistics(){
		return this->stats;
	}

	/*!
	 * set the model parameter values
	 */
	void setStatistics(const std::vector<double>&st){
		this->stats=st;
	}

	/*!
	 * \return the terms
	 */
	std::vector<double> values(){
		return stats;
	}

	/*!
	 * \return stats * thetas
	 */
	double logLik(){
		double ll = 0.0;
		for(int i=0;i<stats.size();i++)
			ll += stats[i];
		return ll;
	}

};


} /* namespace ernm */
#endif /* OFFSETH_ */
