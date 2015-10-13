/*
 * Constraint.h
 *
 *  Created on: Mar 21, 2012
 *      Author: ianfellows
 */

#ifndef CONSTRAINTH_
#define CONSTRAINTH_
#include <vector>
#include "Offset.h"
#include "limits.h"

namespace ernm {




/*!
 * Constraints are just an offset which is log(0) if satisfied, and
 * log(-Inf) if not.
 */
template<class NetworkEngine, class OffsetEngine>
class Constraint : public Offset<NetworkEngine,OffsetEngine>{
public:
	Constraint() : Offset<NetworkEngine,OffsetEngine>(){}

	Constraint(Rcpp::List params) : Offset<NetworkEngine,OffsetEngine>(params){}

	virtual ~Constraint() {}

	virtual AbstractOffset<NetworkEngine>* vCreateUnsafe(Rcpp::List params) const{
		return createUnsafe(params);
	}

	inline AbstractOffset<NetworkEngine>* createUnsafe(Rcpp::List params) const{
		return new Constraint(params);
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
		return new Constraint<NetworkEngine,OffsetEngine>(*this);
	}

	virtual void vCalculate(const BinaryNet<NetworkEngine>& net){
		calculate(net);
	}

	inline void calculate(const BinaryNet<NetworkEngine>& net){
		this->off.updateLogLik(this->off.initialize(net));
	}

	virtual void vDyadUpdate(const BinaryNet<NetworkEngine>& net, int from, int to){
		dyadUpdate(net,from,to);
	}

	inline void dyadUpdate(const BinaryNet<NetworkEngine>& net, int from, int to){
		this->off.updateLogLik(this->off.dyadUpdateDistance(net, from, to));
	}

	virtual void vDiscreteVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, int newValue){
		discreteVertexUpdate(net,vert,variable,newValue);
	}

	inline void discreteVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, int newValue){
		this->off.updateLogLik(this->off.discreteVertexUpdateDistance(net, vert, variable,newValue));
	}

	virtual void vContinVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, double newValue){
		continVertexUpdate(net,vert,variable,newValue);
	}

	inline void continVertexUpdate(const BinaryNet<NetworkEngine>& net, int vert,
			int variable, double newValue){
		this->off.updateLogLik(this->off.continVertexUpdateDistance(net, vert, variable, newValue));
	}

};


/*!
 * Constraints restrict the sample space of the model by setting it's log likelihood to
 * -infinity when outside the restricted space. If the network does not satisfy the
 * constraint, it will be snapped to it when MCMC is run.
 */
template<class Engine>
class BaseConstraint{
private:
	double logValue;
protected:

public:
	virtual ~BaseConstraint(){};

	/*!
	 * calculate how many steps away the constraint is from being satisfied
	 */
	double initialize(const BinaryNet<Engine>& net){
		Rf_error("initialize must be implemented");
	}

	/*!
	 * dyad update for how many steps away the constraint is from being satisfied
	 */
	double dyadUpdateDistance(const BinaryNet<Engine>& net, int& from, int& to){
		Rf_error("distanceFromSatisfaction (dyad) must be implemented");
	}

	/*!
	 * disc vertex update for how many steps away the constraint is from being satisfied
	 */
	double discreteVertexUpdateDistance(const BinaryNet<Engine>& net,
			int& vert, int& variable,int& newValue){
		Rf_error("distanceFromSatisfaction (vertex) must be implemented");
	}

	/*!
	 * cont vertex update for how many steps away the constraint is from being satisfied
	 */
	double continVertexUpdateDistance(const BinaryNet<Engine>& net, int vert,
			int variable, double newValue){
		Rf_error("distanceFromSatisfaction (vertex cont) must be implemented");
	}


	/*!
	 * takes the distance from constraint satisfaction, and sets the value of the
	 * constraint offset o(T).
	 *
	 */
	void updateLogLik(double dist){
		if(near(dist,0.0)){
			logValue = 0.0;
		}else{
			logValue = -100000000.0 - 100000.0*dist;
		}
	}

	/*!
	 * value of offset
	 */
	double logLik(){
		return logValue;
	};

	std::vector<double> values(){
		return std::vector<double>(1,logValue);
	}

	/*!
	 * number of terms
	 */
	int size(){
		return 1;
	}

	void calculate(const BinaryNet<Engine>& net){
		Rf_error("BaseConstraint calculate should not be called");
	}
	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		Rf_error("BaseConstraint dyadUpdate should not be called");
	}
	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, int newValue){
		Rf_error("BaseConstraint discreteVertexUpdate should not be called");
	}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){
		Rf_error("BaseConstraint continVertexUpdate should not be called");
	}
};


} /* namespace ernm */
#endif /* CONSTRAINTH_ */
