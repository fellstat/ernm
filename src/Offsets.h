/*
 * Offsets.h
 *
 *  Created on: Jan 9, 2014
 *      Author: ianfellows
 */

#ifndef OFFSETS_H_
#define OFFSETS_H_

#include "Offset.h"

namespace ernm{



/*!
 * MNAR offset adjusting for biased seed link tracing
 */
template<class Engine>
class BiasedSeed : public BaseOffset< Engine > {
protected:
	double logValue;
	std::vector<int> nseeds;
	std::string variableName;
	std::vector<int> counts;
	BiasedSeed(std::string varName, double lvalue,  std::vector<int> numSeeds,
			std::vector<int> cnts){
		variableName = varName;
		logValue = lvalue;
		nseeds = numSeeds;
		counts = cnts;
	}
public:
	BiasedSeed(){
		logValue = 0.0;
	}

	BiasedSeed(std::string varName, std::vector<int> numSeeds){
		logValue = 0.0;
		variableName = varName;
		nseeds = numSeeds;
	}



	BiasedSeed(List params){
		if(params.size()<2){
			::Rf_error("BiasedSeedOffset: two parameters required");
			return;
		}
		try{
			variableName = as< std::string >(params[0]);
		}catch(...){
			::Rf_error("BiasedSeedOffset requires a nodal variable name");
		}
		try{
			nseeds = as< std::vector<int> >(params[1]);
		}catch(...){
			::Rf_error("BiasedSeedOffset requires the number of seeds");
		}
		logValue = 0.0;
	}

	std::string name(){
		return "biasedSeed";
	}

	void calculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.discreteVarNames();
		int varIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				varIndex = i;
			}
		}
		if(varIndex<0)
			::Rf_error("BiasedSeed::calculate nodal attribute not found in network");

		int nlevels = net.discreteVariableAttributes(varIndex).labels().size();
		if(nlevels!=nseeds.size())
			::Rf_error("length of seeds not equal to number of levels");

		counts = std::vector<int>(nlevels,0.0);
		for(int i=0;i<net.size();i++){
			int val = net.discreteVariableValue(varIndex,i);
			counts[val-1]++;
		}
		calcLogValue();
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
		std::vector<std::string> vars = net.discreteVarNames();
			int varIndex = -1;
			for(int i=0;i<vars.size();i++){
				if(vars[i] == variableName){
					varIndex = i;
				}
			}
		if(variable != varIndex)
			return;
		int val = net.discreteVariableValue(varIndex,vert);
		counts[val-1]--;
		counts[newValue-1]++;
		calcLogValue();
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}

	void calcLogValue(){
		logValue = 0.0;
		for(int i=0;i<counts.size();i++){
			if(counts[i]<nseeds[i]){
				logValue = -DBL_MAX;
				return;
			}
			for(int j=counts[i]-nseeds[i]+1;j<=counts[i];j++)
				logValue += -log((double)j);
			//cout<<counts[i] <<" "<<logValue<<"\n";
		}
	}

	std::vector<double> values(){
		return std::vector<double>(1,logValue);
	}

	double logLik(){
		return logValue;
	}

	int size(){
		return 1;
	}

};

typedef Offset<Directed, BiasedSeed<Directed> > DirectedBiasedSeedOffset;
typedef Offset<Undirected, BiasedSeed<Undirected> > UndirectedBiasedSeedOffset;

/*!
 * MNAR offset for RDS
 */
template<class Engine>
class RdsBias : public BaseOffset< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;

	double logValue;

	std::vector<int> parent;
	std::vector<int> nChildren;
	std::vector<int> effectiveDegrees;
	std::vector<int> order;
	std::vector< std::vector<double> > params;
	int ordLastParent;
/*
	RDSOffset(double lv, std::vector<int> par, std::vector<int> ed, std::vector<int> ord){
		logValue = lv;
		parent = par;
		effectiveDegrees = ed;
		order=ord;
	}
*/
public:

	RdsBias(){
		logValue = ordLastParent = 0;
	}

	RdsBias(std::vector<int> ord){
			order=ord;
			logValue=0.0;
			ordLastParent = 0;
	}

	RdsBias(List pars){
		int nargs = pars.size();
		if(pars.size()<1){
			::Rf_error("RDSOffset: a parameter required");
			return;
		}

		try{
			order = as< std::vector<int> >(pars(0));
		}catch(...){
			::Rf_error("RDSOffset requires a vector of node ids for the sample in order of recruitment");
		}

		if(nargs>1){
			try{
				std::vector<double> p = as< std::vector<double> >(pars(1));
				if(p.size()<2){
					Language call3("print",pars(1));
					call3.eval();
					Rf_error("Invalid sampling probs");
				}
				for(int i=0;i<p.size();i++){
					//cout<<"here ";
					params.push_back(std::vector<double>(i+1,0.0));
					for(int j=0;j<p.size();j++){
						params.at(i).at(std::min(j,i)) += p.at(j);
					}
				}
				//cout<<" there";
			}catch(...){
			}
		}
		ordLastParent = 0;
		logValue = 0.0;
	}



	std::string name(){
		return "rds";
	}


	void calculate(const BinaryNet<Engine>& net){
/*		params.push_back(std::vector<double>());
		params.push_back(std::vector<double>());
		params.push_back(std::vector<double>());
		params.push_back(std::vector<double>());
		params[0].push_back(1);
		params[1].push_back(0.125);
		params[1].push_back(0.875);
		params[2].push_back(0.125);
		params[2].push_back(0.5625);
		params[2].push_back(0.3125);
		params[3].push_back(0.125);
		params[3].push_back(0.5625);
		params[3].push_back(0.1875);
		params[3].push_back(0.125);
*/
		if(net.isDirected())
			Rf_error("RDSOffset can only be used on undirected networks");
		if(order.size()!=net.size())
			Rf_error("recruitment ordering must be the same size as network with -1 marking unobserved nodes");
		parent = std::vector< int >(net.size(),-1);
		nChildren = std::vector< int >(net.size(),0);
		effectiveDegrees = std::vector< int >(net.size(),-1);
		ordLastParent = -1;
		for(int i=0;i<order.size();i++){
			int ind = order[i];
			if(ind<0)
				continue;
			NeighborIterator it = net.begin(i);
			NeighborIterator end = net.end(i);
			//const Set* nbs = &net.neighbors(i);
			//Set::iterator it = nbs->begin();
			while(it!=end){
				if(order[*it]>ind && !net.isMissing(i,*it)){
					if(parent[*it]!=-1){
						Rf_error("Multiple parents: Improper RDS ordering");
					}
					parent[*it] = i;
					nChildren[i]++;
				}
				it++;
			}
			if(nChildren[i]>0 && ind > ordLastParent)
				ordLastParent = ind;
			//cout<<"node: "<<i<<"\n";
			//cout<<"children: "<<nChildren[i]<<"\n";
			//cout<<"efd: "<<effectiveDegrees[i]<<"\n";
		}
		//cout<<"\n";
		this->logValue = 0.0;
		for(int i=0;i<net.size();i++){
			int ord = order[i];
			if(ord<0)
				continue;
			int deg = net.degree(i);
			NeighborIterator it = net.begin(i);
			NeighborIterator end = net.end(i);
			while(it!=end){
				int tmpOrd;
				if(parent[*it]>=0)
					tmpOrd = order[parent[*it]];
				else
					tmpOrd = -1;
				if(order[*it]>0 && tmpOrd<ord)
					deg--;
				it++;
			}
			effectiveDegrees[i] = deg;
			//for(int j=0;j<nChildren[i];j++)
			//	this->logValue -= deg-j;
			this->logValue -= Rf_lchoose(deg,nChildren[i]);
					//nChildren[i] - deg;
			if(params.size()>0 && nChildren[i]>std::min((int)params.size() - 1,deg)){
				//cout<<"offending node: "<<i<<"\n";
				//cout<<"children: "<<nChildren[i]<<"\n";
				//cout<<"efd: "<<effectiveDegrees[i]<<"\n";
				//cout<<"max allowed children: "<<(((int)params.size()) - 1)<<"\n";
				Rf_error("RDSOffset: Too many children");
			}
			if(params.size()>0)
				if(ord<=ordLastParent)
					this->logValue += log(params.at(std::min((int)params.size()-1,deg)).at(nChildren[i]));
		}

	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		if(!net.isMissing(from,to)){
			//cout << from << " "<<to << "\n";
			Rf_error("RDSOffset: toggling observed variable");
		}

		if(order[from]>=0){
			int node = from;
			int nbr = to;
			int change = net.hasEdge(from,to) ? -1 : 1;
			int tmpOrd;
			if(parent.at(nbr)>=0)
				tmpOrd = order.at(parent.at(nbr));
			else
				tmpOrd = -1;
			if(order.at(nbr)>-1 && tmpOrd<order.at(node))
				change = 0.0;
			int deg = effectiveDegrees.at(node);
			int nch = nChildren.at(node);
			this->logValue -= Rf_lchoose(deg+change,nch) - Rf_lchoose(deg,nch);
			if(params.size()>0)
				if(order.at(node)<=ordLastParent)
					this->logValue += log(params.at(std::min((int)params.size()-1,deg+change)).at(nch)) -
						log(params.at(std::min((int)params.size()-1,deg)).at(nch));
			effectiveDegrees.at(node) += change;
		}

		if(order[to]>=0){
			int node = to;
			int nbr = from;
			int change = net.hasEdge(from,to) ? -1 : 1;
			int tmpOrd;
			if(parent.at(nbr)>=0)
				tmpOrd = order.at(parent.at(nbr));
			else
				tmpOrd = -1;
			if(order.at(nbr)>-1 && tmpOrd<order.at(node))
				change = 0.0;
			int deg = effectiveDegrees.at(node);
			int nch = nChildren.at(node);
			this->logValue -= Rf_lchoose(deg+change,nch) - Rf_lchoose(deg,nch);
			if(params.size()>0)
				if(order.at(node)<=ordLastParent)
					this->logValue += log(params.at(std::min((int)params.size()-1,deg+change)).at(nch)) -
						log(params.at(std::min((int)params.size()-1,deg)).at(nch));
			effectiveDegrees.at(node) += change;
		}
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}

	std::vector<double> values(){
		return std::vector<double>(1,logValue);
	}

	double logLik(){
		return logValue;
	}

	int size(){
		return 1;
	}

};

typedef Offset<Directed, RdsBias<Directed> > DirectedRdsBiasOffset;
typedef Offset<Undirected, RdsBias<Undirected> > UndirectedRdsBiasOffset;

///**
// * An offset term equivelent to NodeCov statistic.
// */
//template<class Engine>
//class REffect : public BaseOffset< Engine > {
//protected:
//	EdgeDirection direction;
//	std::string variableName;
//	int varIndex;
//	bool isDiscrete;
//	double expectedValue;
//public:
//
//	REffect(){
//		varIndex =  0;
//		direction = UNDIRECTED;
//		isDiscrete = false;
//		expectedValue = 0;
//	}
//
//	REffect(std::string name,EdgeDirection d){
//		varIndex = 0;
//		direction = d;
//		variableName = name;
//		isDiscrete = false;
//		expectedValue = 0;
//	}
//
//	REffect(std::string name){
//		varIndex = 0;
//		direction = UNDIRECTED;
//		variableName = name;
//		isDiscrete = false;
//		expectedValue = 0;
//	}
//
//
//	REffect(List params){
//		varIndex = 0;
//		isDiscrete=false;
//		expectedValue = 0;
//		try{
//			variableName = as< std::string >(params(0));
//		}catch(...){
//			::Rf_error("NodeCov requires a nodal variable name");
//		}
//
//		try{
//			int tmp = as< int >(params(1));
//			if(tmp==0)
//				direction = UNDIRECTED;
//			else if(tmp==1)
//				direction = IN;
//			else if(tmp==2)
//				direction = OUT;
//			else
//				::Rf_error("invalid direction");
//		}catch(...){
//			direction = UNDIRECTED;
//		}
//	}
//
//	std::string name(){
//		return "reffect";
//	}
//
//	double getValue(const BinaryNet<Engine>& net, int ind){
//		double val;
//		if(isDiscrete)
//			val = net.discreteVariableValue(varIndex,ind);
//		else
//			val = net.continVariableValue(varIndex,ind);
//		return val;
//	}
//
//	void calculate(const BinaryNet<Engine>& net){
//		isDiscrete = false;
//		std::vector<std::string> vars = net.continVarNames();
//		int variableIndex = -1;
//		for(int i=0;i<vars.size();i++){
//			if(vars[i] == variableName){
//				variableIndex = i;
//			}
//		}
//		if(variableIndex == -1){
//			isDiscrete = true;
//			vars = net.discreteVarNames();
//			for(int i=0;i<vars.size();i++){
//				if(vars[i] == variableName){
//					variableIndex = i;
//				}
//			}
//		}
//		if(variableIndex<0)
//			::Rf_error("nodal attribute not found in network");
//		varIndex = variableIndex;
//		this->stats = std::vector<double>(1,0.0);
//		this->stats[0] = 0;
//		expectedValue = 0.0;
//		for(int i=0;i<net.size();i++)
//			expectedValue +=  getValue(net,i);
//		expectedValue /= (double) net.size();
//		for(int i=0;i<net.size();i++){
//			double val = getValue(net,i);
//			if(net.isDirected()){
//				if(direction == IN || direction == UNDIRECTED)
//					this->stats[0] += (val - expectedValue) * net.indegree(i);
//				if(direction == OUT || direction == UNDIRECTED)
//					this->stats[0] += (val - expectedValue) * net.outdegree(i);
//			}else{
//				this->stats[0] += (val - expectedValue) * net.degree(i);
//			}
//		}
//	}
//
//
//	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
//		double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
//		if(net.isDirected()){
//			if(direction == IN || direction == UNDIRECTED)
//				this->stats[0] += change * (getValue(net,to) - expectedValue);
//			if(direction == OUT || direction == UNDIRECTED)
//				this->stats[0] += change * (getValue(net,from) - expectedValue);
//		}else{
//			this->stats[0] += change * (getValue(net,to)+getValue(net,from) - 2.0 * expectedValue);
//		}
//	}
//
//	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
//						int variable, int newValue){
//		if(isDiscrete && variable==varIndex){
//			double n = net.size();
//			double nv = newValue;
//			double ov = getValue(net,vert);
//			double deg = 0.0;
//			if(net.isDirected()){
//				if(direction == IN || direction == UNDIRECTED)
//					deg += net.indegree(vert);
//				if(direction == OUT || direction == UNDIRECTED)
//					deg += net.outdegree(vert);
//			}else
//				deg = net.degree(vert);
//			this->stats[0] += deg*(nv - ov);
//			double ev = expectedValue - ov / n + nv / n;
//			this->stats[0] -= (!net.isDirected() || direction == UNDIRECTED ? 2.0 : 1.0) * net.nEdges() * (ev - expectedValue);
//			expectedValue = ev;
//		}
//	}
//
//	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
//				int variable, double newValue){
//		if(!isDiscrete && variable==varIndex){
//			double n = net.size();
//			double nv = newValue;
//			double ov = getValue(net,vert);
//			double deg = 0.0;
//			if(net.isDirected()){
//				if(direction == IN || direction == UNDIRECTED)
//					deg += net.indegree(vert);
//				if(direction == OUT || direction == UNDIRECTED)
//					deg += net.outdegree(vert);
//			}else
//				deg = net.degree(vert);
//			this->stats[0] += deg*(nv - ov);
//			double ev = expectedValue - ov / n + nv / n;
//			this->stats[0] -= (!net.isDirected() || direction == UNDIRECTED ? 2.0 : 1.0) * net.nEdges() * (ev - expectedValue);
//			expectedValue = ev;
//		}
//	}
//
//};


/**
 * An offset term equivelent to NodeCov statistic.
 */
template<class Engine>
class REffect : public BaseOffset< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	EdgeDirection direction;
	std::string variableName;
	int varIndex;
	bool geometric;

	std::vector<double> lvals;
	double n;
	double lnm1;
public:

	REffect() : geometric(false){
		varIndex =  0;
		direction = UNDIRECTED;
	}

	REffect(std::string name,EdgeDirection d) : geometric(false){
		varIndex = 0;
		direction = d;
		variableName = name;
	}

	REffect(std::string name) : geometric(false){
		varIndex = 0;
		direction = UNDIRECTED;
		variableName = name;
	}


	REffect(List params) : geometric(false){
		varIndex = 0;
		try{
			variableName = as< std::string >(params(0));
		}catch(...){
			::Rf_error("reffect requires a nodal variable name");
		}

		try{
			int tmp = as< int >(params(1));
			if(tmp==0)
				direction = UNDIRECTED;
			else if(tmp==1)
				direction = IN;
			else if(tmp==2)
				direction = OUT;
			else
				::Rf_error("invalid direction");
		}catch(...){
			direction = UNDIRECTED;
		}
	}

	std::string name(){
		return "reffect";
	}



	virtual void vCalculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.continVarNames();
		varIndex = indexOf(variableName,vars);
		if(varIndex<0)
			::Rf_error("nodal attribute not found in network");
		this->stats = std::vector<double>(1,0.0);
		this->stats.at(0) = 0;
		n = net.size();
		lnm1 = log(n-1);
		if(!geometric)
			lvals = std::vector<double>(net.size(),0.0);
		for(int i=0;i<net.size();i++){
			double val = net.continVariableValue(varIndex,i);
			NeighborIterator it = net.begin(i);
			NeighborIterator end = net.end(i);
			while(it != end && *it < i){
				double nval = net.continVariableValue(varIndex,*it);
				if(nval<0 || nval >= (net.size() - 1))
					::Rf_error("reffect: value out of range");
				//double p = sqrt(val*nval);
				//this->stats.at(0) += log(p) - log(net.size()-1.0-p);
				double q;
				if(!geometric)
					q = std::max((n-1.0-val),(n-1.0-nval)) / (n-1.0); //q = (.5*(net.size()-1.0-val) + .5*(net.size()-1.0-nval)) / (n-1.0);
				else
					q = sqrt((net.size()-1.0-val) * (net.size()-1.0-nval)) / (n-1.0);
				double p = 1.0-q;
				this->stats.at(0) += log(p) - log(q);
				it++;
			}
			if(!geometric){
				lvals.at(i) = log(n - 1.0 - val);
				for(int j=i+1;j<net.size();j++){
					double nval = net.continVariableValue(varIndex,j);
					//double q = (.5*(net.size()-1.0-val) + .5*(net.size()-1.0-nval)) / (n-1.0);
					double q = std::max((net.size()-1.0-val),(net.size()-1.0-nval)) / (n-1.0);
					this->stats.at(0) += log(q);
				//	this->stats.at(0) += .5*log(net.size()-1.0-val);
				//	this->stats.at(0) += .5*log(net.size()-1.0-nval);
				//	this->stats.at(0) -= log(n-1.0);
					//double p = sqrt(val*nval);
					//this->stats.at(0) += log(net.size() - 1.0 - p);
				}
			}else
				this->stats.at(0) += .5*(n-1.0)*log(n-1.0-val);
		}
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
		double fval =  net.continVariableValue(varIndex,from);
		double tval =  net.continVariableValue(varIndex,to);
		double q;
		if(!geometric)
			q = std::max((net.size()-1.0-fval),(net.size()-1.0-tval)) / (net.size()-1.0);//q = (.5*(net.size()-1.0-fval) + .5*(net.size()-1.0-tval)) / (net.size()-1.0);
		else
			q = sqrt((net.size()-1.0-fval) * (net.size()-1.0-tval)) / (net.size()-1.0);
		double p = 1.0-q;
		this->stats.at(0) += change * (log(p) - log(q));
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
						int variable, int newValue){

	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){
		if(variable==varIndex){
			double n = net.size();
			double oval = net.continVariableValue(varIndex,vert);
			if(oval<0 || oval >= (net.size() - 1))
				::Rf_error("reffect update: old value out of range");
			NeighborIterator it = net.begin(vert);
			NeighborIterator end = net.end(vert);
			while(it!=end){
				double nval = net.continVariableValue(varIndex,*it);
				if(nval<0 || nval >= (net.size() - 1))
					::Rf_error("reffect update: old value out of range");
				double q;
				if(!geometric)
					q = std::max((net.size()-1.0-newValue),(net.size()-1.0-nval)) / (n-1.0);//q = (.5*(net.size()-1.0-newValue) + .5*(net.size()-1.0-nval)) / (net.size()-1.0);
				else
					q = sqrt((net.size()-1.0-newValue) * (net.size()-1.0-nval)) / (net.size()-1.0);
				double p = 1.0-q;
				double qold;
				if(!geometric)
					qold = std::max((net.size()-1.0-oval) ,(net.size()-1.0-nval)) / (net.size()-1.0);
				else
					qold = sqrt((net.size()-1.0-oval) * (net.size()-1.0-nval)) / (net.size()-1.0);
				double pold = 1.0-qold;
				this->stats[0] += log( (p/pold)*(qold/q) );//log(p) - log(q) - log(pold) + log(qold);
				it++;
			}
			if(!geometric){
				double lnv = log(n-1.0-newValue);
				double lov = lvals[vert];
				double l = std::min(lnv,lov);
				double u = std::max(lnv,lov);
				bool red = newValue <= oval;
				double c=0.0;
				for(int i=0;i<net.size();i++){
					if(i==vert)
						continue;
					double lnbr = lvals[i];
					if(lnbr > l && lnbr < u){
						this->stats[0] += red ? (lnv-lnbr) : (lnbr-lov);
					}else if(lnbr<u)
						c++;
				}
				this->stats[0] += c * (lnv - lov);
				lvals[vert] = lnv;

			}else
				this->stats[0] += .5*(n-1.0)*log(n-1.0-newValue) -
					.5*(n-1.0)*log(n-1.0-oval);
		}
	}

};
typedef Offset<Directed, REffect<Directed> > DirectedREffectOffset;
typedef Offset<Undirected, REffect<Undirected> > UndirectedREffectOffset;

// Hamming offset for tapering to MRF models
template<class Engine>
class HammingOffset : public BaseOffset< Engine > {
protected:
  //List edges;
  boost::shared_ptr< std::vector< std::pair<int,int> > > edges;
  double theta;
public:
  HammingOffset(){
    std::vector<double> v(1,0.0);
    this->stats=v;
    theta = 0;
  }
  
  HammingOffset(Rcpp::List params){
    
    std::vector<double> v(1,0.0);
    this->stats=v;
    theta = params[1];
    
    Rcpp::NumericMatrix edgeList = params[0];
    boost::shared_ptr< std::vector< std::pair<int,int> > > edges_tmp(new std::vector<std::pair<int,int> >());
    edges_tmp->reserve(edgeList.nrow());
    for(int i=0;i<edgeList.nrow();i++){
      // Do the minus one stuff since this is a R interface
      int from = edgeList(i,0)-1;
      int to = edgeList(i,1)-1;
      if(from < 0 || to<0){
        Rf_error("Edgelist indices out of range");
      }
      std::pair<int,int> p = std::make_pair(from,to);
      edges_tmp->push_back(p);
    }
    edges = edges_tmp;
  }
  
  std::string name(){
    return "hamming";
  }
  
  void calculate(const BinaryNet<Engine>& net){
    // Start by assuming they are completely the same
    // Step through the edge list to calculate the edges that are missing
    // Then use nEdges to get the additional edges
    std::vector<double> v(1,0.0);
    std::vector< std::pair<int,int> > ::iterator it = edges->begin();
    int shared = 0;
    while(it != edges->end()){
      if(!net.hasEdge(it->first,it->second)){
        v[0] += 1;
      }else{
        shared +=1;
      }
      it++;
    }
    
    // add the surplus edges to the differing count
    v[0] += (net.nEdges() - shared);
    v[0] = v[0]*theta;
    this->stats = v;
    return;
  }
  
  void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
    //Check if the toggle is in the edge list:
    if(this->stats[0]==0){
      this->stats[0] +=1;
      return;
    }
    
    bool net_ind = net.hasEdge(from,to);
    int is_in_net = net_ind?1:0;
    
    int is_in_list = 0;
    std::vector< std::pair<int,int> > ::iterator it = edges->begin();
    if(net.isDirected()){
      while(it != edges->end()){
        if(it->first == from and  it->second == to){
          is_in_list = 1;
          break;
        } 
        it++;
      }
    }else{
      while(it != edges->end()){
        if((it->first == from and it->second == to) or (it->first == to and it->second == from)){
          is_in_list = 1;
          break;
        }
        it++;
      }
    }
    
    this->stats[0] += (theta)*((is_in_list == is_in_net)?1:(-1));
    return;
  }
  
  // Don't do anything - purely edge statistic
  void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
                            int variable, int newValue){}
  
  void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
                          int variable, double newValue){}
};

typedef Offset<Directed, HammingOffset<Directed> > DirectedHammingOffset;
typedef Offset<Undirected, HammingOffset<Undirected> > UndirectedHammingOffset;

}


#endif /* OFFSETS_H_ */
