/*
 * Constraints.h
 *
 *  Created on: Jan 9, 2014
 *      Author: ianfellows
 */


#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_

#include "Constraint.h"

namespace ernm{

/*!
 * Restricts the sample space to the set of networks with all nodes having degree
 * >=lower and <=upper
 *
 */
template<class Engine>
class BoundedDegree : public BaseConstraint< Engine >{
protected:
	int upper;
	int lower;
	std::vector<int> scratch;
	double dist;
public:

	BoundedDegree() : upper(10000000), lower(0), dist(0.0){}

	BoundedDegree(int low, int up) : dist(0.0){
		lower=low;
		upper=up;
	}

	BoundedDegree(List params) : dist(0.0){
		if(params.size()<2){
			::Rf_error("BoundedDegree: two parameters required");
			return;
		}
		try{
			lower = as< int >(params[0]);
		}catch(...){
			::Rf_error("BoundedDegree: Invalid lower bound");
		}
		try{
			upper = as< int >(params[1]);
		}catch(...){
			::Rf_error("BoundedDegree: Invalid upper bound");
		}

	}

	std::string name(){
		return "boundedDegree";
	}

	/*!
	 * calculate how many steps away the constraint is from being satisfied
	 */
	double initialize(const BinaryNet<Engine>& net){
		dist = 0.0;
		for(int i=0;i<net.size();i++){
			int deg = net.degree(i);
			if(deg>upper)
				dist+= deg-upper;
			if(deg<lower)
				dist+= lower - deg;
		}
		return dist;
	}

	//dyad update
	double dyadUpdateDistance(const BinaryNet<Engine>& net, int& from, int&to){
		bool addingEdge = !net.hasEdge(from,to);
		int dfrom = net.degree(from);
		int dto = net.degree(to);
		//if(dfrom<lower || dto<lower || dfrom>upper || dto>upper)
		//	::Rf_error("Network degrees outside degree bounds");
		if(addingEdge){
			if(dfrom<lower)
				dist--;
			else if(dfrom>=upper)
				dist++;
			if(dto<lower)
				dist--;
			else if(dto>=upper)
				dist++;
		}else{
			if(dfrom<=lower)
				dist++;
			else if(dfrom>upper)
				dist--;
			if(dto<=lower)
				dist++;
			else if(dto>upper)
				dist--;
		}
		return dist;
	}

	//vertex update
	double discreteVertexUpdateDistance(const BinaryNet<Engine>& net,
			int& vert, int& variable,int& newValue){
		return dist;
	}

	double continVertexUpdateDistance(const BinaryNet<Engine>& net, int vert,
			int variable, double newValue){
		return dist;
	}

};

typedef Constraint<Directed, BoundedDegree<Directed> > DirectedBoundedDegreeConstraint;
typedef Constraint<Undirected, BoundedDegree<Undirected> > UndirectedBoundedDegreeConstraint;


}


#endif /* CONSTRAINTS_H_ */
