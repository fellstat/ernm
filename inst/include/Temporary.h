/*
 * Temporary.h
 *
 *  Created on: Dec 30, 2015
 *      Author: ianfellows
 */

#ifndef TEMPORARY_H_
#define TEMPORARY_H_

#include "Offset.h"
#include "Stats.h"

namespace ernm{

template<class Engine>
class StarPenalty : public BaseOffset< Engine > {
protected:
	double logValue;
	double penalty;
	double observedValue;
	Star<Engine> star;

public:

	StarPenalty(){}

	StarPenalty(List params){
		int starNum;
		BinaryNet<Engine> network;
		try{
			SEXP nwPtr = params(0);
			network = BinaryNet<Engine>(nwPtr);
		}catch(...){
			::Rf_error("network required");
		}
		try{
			starNum = as< int >(params(1));
		}catch(...){
			::Rf_error("star num required");
		}
		try{
			penalty = as< double >(params(2));
		}catch(...){
			::Rf_error("penalty required");
		}
		List l;
		l.push_back(starNum);
		l.push_back(1);
		star = Star<Engine>(l);
		star.calculate(network);
		observedValue = star.statistics().at(0);
	}

	std::string name(){
		return "starPenalty";
	}

	void calculate(const BinaryNet<Engine>& net){
		star.calculate(net);
		logLik();
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		star.dyadUpdate(net, from, to);
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}

	std::vector<double> values(){
		return std::vector<double>(1,logLik());
	}

	double logLik(){
		logValue = penalty * pow(observedValue - star.statistics().at(0), 2.0);
		return logValue;
	}

	int size(){
		return 1;
	}

};

typedef Offset<Directed, StarPenalty<Directed> > DirectedStarPenalty;
typedef Offset<Undirected, StarPenalty<Undirected> > UndirectedStarPenalty;

}


#endif /* TEMPORARY_H_ */
