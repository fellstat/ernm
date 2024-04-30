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

/*!
 * Restricts the sample space so that the degrees of a set of nodes are kept constant
 *
 */
template<class Engine>
class FixedDegree : public BaseConstraint< Engine >{
protected:
	std::vector<int> nodes;
	std::vector<int> toDegree;
	std::vector<bool> isFixed;
	std::vector<int> fixedDegree;
	bool allNodes;
	double dist;
public:

	FixedDegree() : dist(0.0){
		allNodes = true;
	}

	FixedDegree(std::vector<int> nds) : dist(0.0){
		nodes = nds;
		allNodes=false;
	}

	FixedDegree(std::vector<int> nds,std::vector<int> degs) : dist(0.0){
		nodes = nds;
		toDegree = degs;
		allNodes=false;
	}

	FixedDegree(List params) : dist(0.0){
		if(params.size()<1){
			::Rf_error("FixedDegree: 1 parameters required");
			return;
		}
		try{
			nodes = as< std::vector<int> >(params[0]);
			allNodes = false;
		}catch(...){
			allNodes = true;
		}

		try{
			toDegree = as< std::vector<int> >(params[1]);
		}catch(...){
			toDegree = std::vector<int>();
		}


	}

	std::string name(){
		return "fixedDegree";
	}


	/*!
	 * calculate how many steps away the constraint is from being satisfied
	 */
	double initialize(const BinaryNet<Engine>& net){
		int deg;
		dist = 0.0;
		if(allNodes){
			nodes = std::vector<int>();
			for(int i=0;i<net.size();i++)
				nodes.push_back(i);
		}
		isFixed = std::vector<bool>(net.size(),false);
		fixedDegree = std::vector<int>(net.size(),0);
		for(int i=0;i<nodes.size();i++){
			if(nodes[i]>=net.size() || nodes[i]<0){
				Rf_error("FixedDegree: attempting to fix invalid node ids");
				return 0.0;
			}
			//cout <<nodes[i] << isFixed.size();;
			isFixed.at(nodes[i]) = true;
			deg = net.degree(nodes[i]);
			if(toDegree.size()>i){
				dist += abs(deg - toDegree[i]);
				deg = toDegree[i];
			}
			fixedDegree.at(nodes[i]) = deg;
		}
		return dist;
	}

	//dyad update
	double dyadUpdateDistance(const BinaryNet<Engine>& net, int& from, int&to){
		bool addingEdge = !net.hasEdge(from,to);
		//if(dfrom<lower || dto<lower || dfrom>upper || dto>upper)
		//	::Rf_error("Network degrees outside degree bounds");
		if(isFixed[from]){
			int dfrom = net.degree(from);
			dist += addingEdge ? (dfrom<fixedDegree[from] ? -1.0 : 1.0) :
					(dfrom>fixedDegree[from] ? -1.0 : 1.0);
		}
		if(isFixed[to]){
			int dto = net.degree(to);
			dist += addingEdge ? (dto<fixedDegree[to] ? -1.0 : 1.0) :
					(dto>fixedDegree[to] ? -1.0 : 1.0);
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

typedef Constraint<Directed, FixedDegree<Directed> > DirectedFixedDegreeConstraint;
typedef Constraint<Undirected, FixedDegree<Undirected> > UndirectedFixedDegreeConstraint;


/*!
 * Restricted the sample space such that the vertex variables of a set of nodes
 * are kept constant.
 *
 * TODO:continuous variable toggleing not implemented
 */
template<class Engine>
class FixedNode : public BaseConstraint< Engine >{
protected:
	std::set<int> nodes;
	std::vector<std::map<int,int> > variates;
public:

	FixedNode(){}
	FixedNode(std::vector<int> tmp){
		setNodes(tmp);
	}

	FixedNode(List params){
		if(params.size()<1){
			::Rf_error("FixedNode: two parameters required");
			return;
		}
		try{
			std::vector<int> tmp = as< std::vector<int> >(params[0]);
			setNodes(tmp);
		}catch(...){
			::Rf_error("FixedNode: Invalid node ids.");
		}
	}

	void setNodes(std::vector<int> tmp){
		nodes.insert(tmp.begin(),tmp.end());
	}


	std::string name(){
		return "fixedNode";
	}

	double initialize(const BinaryNet<Engine>& net){
		std::set<int>::iterator it = nodes.begin();
		int nvar = net.discreteVarNames().size();
		variates.clear();
		variates.resize(nvar,std::map<int,int>());
		for(;it!=nodes.end();it++){
			for(int i=0;i<nvar;i++){
				variates.at(i).insert(std::make_pair(*it,net.discreteVariableValue(i,*it)));
			}
		}
		return 0.0;
	}

	//dyad update
	double dyadUpdateDistance(const BinaryNet<Engine>& net, int& from, int&to){
		return 0.0;
	}
	//vertex update
	double discreteVertexUpdateDistance(const BinaryNet<Engine>& net,
			int& vert, int& variable,int& newValue){
		std::set<int>::iterator it = nodes.find(vert);
		bool isValid=true;
		if(it!=nodes.end()){
			isValid = variates.at(variable).at(vert) == newValue;
		}
		//cout<<"is valid: " << isValid<<" node: " << vert<<"\n";
		return !isValid;
	}


};

typedef Constraint<Directed, FixedNode<Directed> > DirectedFixedNodeConstraint;
typedef Constraint<Undirected, FixedNode<Undirected> > UndirectedFixedNodeConstraint;


}


#endif /* CONSTRAINTS_H_ */
