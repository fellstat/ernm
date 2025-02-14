/*
 * Stats.h
 *
 *  Created on: Jan 9, 2014
 *      Author: ianfellows
 */

#ifndef STATS_H_
#define STATS_H_

#include "Stat.h"
#include <map>
#include <vector>
#include <utility>
#include <boost/container/flat_map.hpp>
namespace ernm{


/*!
 * the number of edges in the network
 */
template<class Engine>
class Edges : public BaseStat<Engine>{
public:
	Edges(){
	}

	/*!
	 * constructor. params is unused
	 */
	Edges(List params){
	}


	std::string name(){
		return "edges";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"edges");
        return statnames;
	}
    
	void calculate(const BinaryNet<Engine>& net){
		std::vector<double> v(1,net.nEdges());
		this->stats=v;
		if(this->thetas.size()!=1){
			//this starts theta at a reasonable value assuming erdos-renyi
			double ne = net.nEdges();
			double nd = net.maxEdges();
			this->thetas = std::vector<double>(1,log(ne) - log(nd - ne));
		}
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		this->stats[0] += net.hasEdge(from,to) ? -1.0 : 1.0;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Edges<Directed> > DirectedEdges;
typedef Stat<Undirected, Edges<Undirected> > UndirectedEdges;


/*!
 * in/out/k-Stars
 */
template<class Engine>
class Star : public BaseStat< Engine > {
protected:
	std::vector<int> starDegrees; /*!< the star degrees */
	EdgeDirection direction;
public:

	Star(){
			std::vector<double> v(1,0.0);
			std::vector<double> t(1,0.0);
			this->stats=v;
			this->thetas = t;
			direction=IN;
	}

	/*!
	 * \param params 	a list
	 */
	Star(List params){
		try{
			starDegrees = as< std::vector<int> >(params(0));
		}catch(...){
			::Rf_error("Star requires a single integer vector as a parameter");
		}

		try{
			int tmp = as< int >(params(1));
			if(tmp==1)
				direction = IN;
			else if(tmp==2)
				direction = OUT;
			else
				::Rf_error("invalid direction");
		}catch(...){
			direction = IN;
		}
		std::vector<double> v(starDegrees.size(),0.0);
		std::vector<double> t(starDegrees.size(),0.0);
		this->stats=v;
		this->thetas = t;
	}


	std::string name(){
		return "star";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<starDegrees.size();i++){
            int d = starDegrees[i];
            std::string nm = "star."+asString(d);
            statnames.push_back(nm);
        }
        return statnames;
	}
    
	void calculate(const BinaryNet<Engine>& net){
		std::vector<double> v(starDegrees.size(),0.0);
		for(int i=0; i<net.size();i++){
			double nEd;
			if(!net.isDirected())
				nEd = net.degree(i);
			else{
				if(direction == IN)
					nEd = net.indegree(i);
				else
					nEd = net.outdegree(i);
			}
			for(int j=0;j<starDegrees.size();j++){
				v[j] += nchoosek(nEd,(double)starDegrees[j]);
				//cout << "n:"<<nEd<<" s:"<<starDegrees[j]<<" nck:"<<nchoosek(nEd,(double)starDegrees[j])<<"\n" ;
			}
		}
		this->stats=v;
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		int n;
		if(!net.isDirected())
			n = net.degree(to);
		else{
			if(direction==IN)
				n = net.indegree(to);
			else
				n = net.outdegree(from);
		}
		bool edge = net.hasEdge(from,to);
		for(int i=0;i<starDegrees.size();i++){
				if(edge){
					this->stats[i] +=  -nchoosek(n,starDegrees[i]) +
							nchoosek(n-1.0,starDegrees[i]);

				}else{
					this->stats[i] +=  nchoosek(n+1.0,starDegrees[i])-nchoosek(n,starDegrees[i]);
				}
		}
		if(!net.isDirected()){
			n = net.degree(from);
			edge = net.hasEdge(from,to);
			for(int i=0;i<starDegrees.size();i++){
				if(edge){
					this->stats[i] +=  -nchoosek(n,starDegrees[i]) +
							nchoosek(n-1.0,starDegrees[i]);

				}else{
					this->stats[i] +=  nchoosek(n+1.0,starDegrees[i])-nchoosek(n,starDegrees[i]);
				}
			}
		}
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Star<Directed> > DirectedStar;
typedef Stat<Undirected, Star<Undirected> > UndirectedStar;


/*!
 * The number of triangles in the network.
 */
template<class Engine>
class Triangles : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	double sumTri;
public:


	Triangles(){
		std::vector<double> v(1,0.0);
		std::vector<double> t(1,0.0);
		this->stats=v;
		this->thetas = t;
		sumTri = 0.0;
	}
	Triangles(List params){
		std::vector<double> v(1,0.0);
		std::vector<double> t(1,0.0);
		this->stats=v;
		this->thetas = t;
		sumTri = 0.0;
	}

	std::string name(){
		return "triangles";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"triangles");
        return statnames;
	}
    int sharedNbrs(const BinaryNet<Engine>& net, int from, int to){
    	if(net.isDirected()){
    		return directedSharedNbrs(net, from, to);
    	}
    	return undirectedSharedNbrs(net, from, to);
    }
	int undirectedSharedNbrs(const BinaryNet<Engine>& net, int from, int to){
		NeighborIterator fit = net.begin(from);
		NeighborIterator fend = net.end(from);
		NeighborIterator tit = net.begin(to);
		NeighborIterator tend = net.end(to);
		int shared = 0;
		while(tit!=tend && fit!=fend){
			if(*tit==*fit){
				shared++;
				tit++;
				fit++;
			}else if(*tit<*fit){
				tit++;
			}else
				fit++;
		}
		return shared;
	}

	int directedSharedNbrs(const BinaryNet<Engine>& net, int from, int to){
		NeighborIterator ifit = net.inBegin(from);
		NeighborIterator ifend = net.inEnd(from);
		NeighborIterator ofit = net.outBegin(from);
		NeighborIterator ofend = net.outEnd(from);
		int shared = 0;
		while(ifit != ifend){
			shared += net.hasEdge(*ifit, to);
			shared += net.hasEdge(to, *ifit);
			ifit++;
		}
		while(ofit != ofend){
			shared += net.hasEdge(*ofit, to);
			shared += net.hasEdge(to, *ofit);
			ofit++;
		}
		return shared;
	}


	void calculate(const BinaryNet<Engine>& net){

		std::vector<double> v(1,0.0);
		this->stats = v;
		if(this->thetas.size()!=1)
			this->thetas = v;
		sumTri = 0.0;

		boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net.edgelist();

		std::vector< std::pair<int,int> >::iterator it = edges->begin();
		while(it != edges->end()){
			int shared = sharedNbrs(net, (*it).first,(*it).second);
			sumTri += shared;
			it++;
		}
		sumTri = sumTri/3.0;
		this->stats[0] = sumTri;//sumSqrtTri - sumSqrtExpected;
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		int shared = sharedNbrs(net, from, to);
		bool hasEdge = net.hasEdge(from,to);
		if(hasEdge){
			sumTri -= shared;
		}else{
			sumTri += shared;
		}
		this->stats[0] = sumTri;//sumSqrtTri - sumSqrtExpected;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Triangles<Directed> > DirectedTriangles;
typedef Stat<Undirected, Triangles<Undirected> > UndirectedTriangles;






/*!
 * the sum over all nodes of the square root of the number of triangles incident
 * on the node minus what would be expected by chance given the degrees of the node's
 * neighbors.
 *
 * A robust transitivity statistic with almost no degeneracy. only currently implemented
 * for undirected nets
 */
template<class Engine>
class Transitivity : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	double sumSqrtTri;
	double sumSqrtExpected;
	std::vector<double> triadCounts;
	std::vector<double> sumNbrDegrees;;
public:


	Transitivity(){
		std::vector<double> v(1,0.0);
		std::vector<double> t(1,0.0);
		this->stats=v;
		this->thetas = t;
		sumSqrtTri = sumSqrtExpected = 0.0;
	}
	Transitivity(List params){
		std::vector<double> v(1,0.0);
		std::vector<double> t(1,0.0);
		this->stats=v;
		this->thetas = t;
		sumSqrtTri = sumSqrtExpected = 0.0;
	}

	std::string name(){
		return "transitivity";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"transitivity");
        return statnames;
	}

	void calcAtNode(const BinaryNet<Engine>& net, int& node, std::vector<double>& results){
		//const Set* nbs = &net.neighbors(node);
		//Set::iterator it = nbs->begin();
		NeighborIterator it = net.begin(node);
		NeighborIterator end = net.end(node);
		NeighborIterator nit;
		double tri = 0.0;
		double degSum = 0.0;
		while(it!=end){
			nit = it;
			nit++;
			for(;nit!=end;nit++){
				tri += net.hasEdge(*it,*nit);
			}
			degSum += net.degree(*it);
			it++;
		}
		//double deg = net.degree(node);
		//double density = degSum / ((double)nbs->size()) / (net.size()-1.0);
		//double expected = density*deg*(deg-1.0)/2.0;
		results.at(0) = tri;
		results.at(1) = degSum;
	}

	double trans(double val){return sqrt(val + 3.0/8.0);}

	void calculate(const BinaryNet<Engine>& net){
		triadCounts = std::vector<double>(net.size(),0.0);
		sumNbrDegrees = std::vector<double>(net.size(),0.0);
		std::vector<double> v(1,0.0);
		this->stats = v;
		if(this->thetas.size()!=1)
			this->thetas = v;
		sumSqrtTri = 0.0;
		sumSqrtExpected = 0.0;
		std::vector<double> tmp(2,0.0);
		for(int i=0; i<net.size();i++){
			double deg = net.degree(i);
			calcAtNode(net,i,tmp);
			triadCounts[i] = tmp[0];
			sumNbrDegrees[i] = tmp[1];
			sumSqrtTri += trans(tmp[0]);
			//double density = (tmp[1] / deg - 1.0) / (net.size()-2.0);
			double nEdgesBetweenNbrs = tmp[1] - tmp[0] - deg;
			//if(deg<.5)
			//	density = 0.0;
			int nPosTri = round(deg*(deg-1.0)/2.0);
			double nPosEdges = ( deg * (net.size() - 2.0) - nPosTri);
			double nExpected = nEdgesBetweenNbrs * nPosTri / nPosEdges;
			if(nPosEdges<.5)
				nExpected=0.0;
			//sumSqrtExpected += expectedAnscombe(nExpected,round(nPosEdges));
			sumSqrtExpected += expectedAnscombe2(round(nPosEdges), nPosTri, round(nEdgesBetweenNbrs));
			//		trans(density*deg*(deg-1.0)/2.0);
			//cout <<sqrt(tmp[0])<<","<<nExpected<<","<<nPosEdges<<","<<deg<<",\n";

		}

		this->stats[0] = sumSqrtTri - sumSqrtExpected;
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){

		BinaryNet<Engine>* pnet = const_cast< BinaryNet<Engine>* > (&net);
		bool edge = net.hasEdge(from,to);
		double change = edge ? -1.0 : 1.0;

		//const Set* nbs = &net.neighbors(from);
		//Set::iterator it = nbs->begin();
		double curTriadValue, curSumNbrDegrees, deg, newTriadValue, newSumNbrDegrees;
		NeighborIterator it = net.begin(from);
		NeighborIterator end = net.end(from);
		while(it!=end){
			if(*it == to){
				it++;
				continue;
			}
			curTriadValue = triadCounts[*it];
			curSumNbrDegrees = sumNbrDegrees[*it];
			if(net.hasEdge(to,*it)){
				sumSqrtTri += trans(curTriadValue + change) - trans(curTriadValue);
				triadCounts[*it] = curTriadValue + change;
				newSumNbrDegrees = sumNbrDegrees[*it] = curSumNbrDegrees + 2.0*change;
			}else
				newSumNbrDegrees = sumNbrDegrees[*it] = curSumNbrDegrees + change;
			newTriadValue = triadCounts[*it];

			deg = net.degree(*it);

			double curNEdgesBetweenNbrs = curSumNbrDegrees - curTriadValue - deg;
			double newNEdgesBetweenNbrs = newSumNbrDegrees - newTriadValue - deg;
			int nPosTri = round(deg*(deg-1.0)/2.0);
			double nPosEdges = ( deg * (net.size() - 2.0) - nPosTri);
			double curNExpected = curNEdgesBetweenNbrs * nPosTri / nPosEdges;
			double newNExpected = newNEdgesBetweenNbrs * nPosTri / nPosEdges;
			if(nPosEdges<.5)
				curNExpected = newNExpected = 0.0;
			//sumSqrtExpected += expectedAnscombe(newNExpected,round(nPosEdges)) -
			//		expectedAnscombe(curNExpected,round(nPosEdges));
			sumSqrtExpected += expectedAnscombe2(round(nPosEdges), nPosTri, round(newNEdgesBetweenNbrs)) -
				expectedAnscombe2(round(nPosEdges), nPosTri, round(curNEdgesBetweenNbrs));

			it++;
		}

		it = net.begin(to);
		end = net.end(to);
		while(it!=end){
			if(*it == from){
				it++;
				continue;
			}
			curTriadValue = triadCounts[*it];
			curSumNbrDegrees = sumNbrDegrees[*it];
			if(net.hasEdge(from,*it)){
				it++;
				continue;
				//sumSqrtTri += trans(curTriadValue + change) - trans(curTriadValue);
				//triadCounts[*it] = curTriadValue + change;
			}
			newSumNbrDegrees = sumNbrDegrees[*it] = curSumNbrDegrees + change;
			newTriadValue = triadCounts[*it];

			deg = net.degree(*it);

			double curNEdgesBetweenNbrs = curSumNbrDegrees - curTriadValue - deg;
			double newNEdgesBetweenNbrs = newSumNbrDegrees - newTriadValue - deg;
			int nPosTri = round(deg*(deg-1.0)/2.0);
			double nPosEdges = ( deg * (net.size() - 2.0) - nPosTri);
			double curNExpected = curNEdgesBetweenNbrs * nPosTri / nPosEdges;
			double newNExpected = newNEdgesBetweenNbrs * nPosTri / nPosEdges;
			if(nPosEdges<.5)
				curNExpected = newNExpected = 0.0;
			//sumSqrtExpected += expectedAnscombe(newNExpected,round(nPosEdges)) -
			//		expectedAnscombe(curNExpected,round(nPosEdges));
			sumSqrtExpected += expectedAnscombe2(round(nPosEdges), nPosTri, round(newNEdgesBetweenNbrs)) -
				expectedAnscombe2(round(nPosEdges), nPosTri, round(curNEdgesBetweenNbrs));
			it++;
		}

		std::vector<double> tmp(2,0.0);

		for(int i=0;i<2;i++){
			int node;
			if(i==0)
				node = from;
			else
				node = to;
			double curDeg = net.degree(node);
			pnet->toggle(from,to);
			double newDeg = net.degree(node);
			curTriadValue = triadCounts[node];
			curSumNbrDegrees = sumNbrDegrees[node];
			calcAtNode(*pnet,node,tmp);
			newTriadValue = triadCounts[node] = tmp[0];
			newSumNbrDegrees = sumNbrDegrees[node] = tmp[1];
			sumSqrtTri += trans(triadCounts[node]) - trans(curTriadValue);

			double curNEdgesBetweenNbrs = curSumNbrDegrees - curTriadValue - curDeg;
			double newNEdgesBetweenNbrs = newSumNbrDegrees - newTriadValue - newDeg;
			int newNPosTri = round(newDeg*(newDeg-1.0)/2.0);
			int curNPosTri = round(curDeg*(curDeg-1.0)/2.0);
			double newNPosEdges = ( newDeg * (net.size() - 2.0) - newNPosTri);
			double curNPosEdges = ( curDeg * (net.size() - 2.0) - curNPosTri);
			double curNExpected = curNEdgesBetweenNbrs * curNPosTri / curNPosEdges;
			double newNExpected = newNEdgesBetweenNbrs * newNPosTri / newNPosEdges;
			if(newNPosEdges<.5)
				newNExpected =  0.0;
			if(curNPosEdges<.5)
				curNExpected =  0.0;
			//sumSqrtExpected += expectedAnscombe(newNExpected,round(newNPosEdges)) -
			//		expectedAnscombe(curNExpected,round(curNPosEdges));
			sumSqrtExpected += expectedAnscombe2(round(newNPosEdges), newNPosTri, round(newNEdgesBetweenNbrs)) -
				expectedAnscombe2(round(curNPosEdges), curNPosTri, round(curNEdgesBetweenNbrs));
			pnet->toggle(from,to);
		}
		this->stats[0] = sumSqrtTri - sumSqrtExpected;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Transitivity<Directed> > DirectedTransitivity;
typedef Stat<Undirected, Transitivity<Undirected> > UndirectedTransitivity;




/*!
 * the number of reciprocal edges in the network
 */
template<class Engine>
class Reciprocity : public BaseStat< Engine > {
public:
	Reciprocity(){
		std::vector<double> v(1,0.0);
		std::vector<double> t(1,0.0);
		this->stats = v;
		this->thetas = t;
	}

	/*!
	 * constructor. params is unused
	 */
	Reciprocity(List params){
		std::vector<double> v(1,0.0);
		std::vector<double> t(1,0.0);
		this->stats = v;
		this->thetas = t;
	}

	std::string name(){
		return "reciprocity";
	}
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"reciprocity");
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		if(!net.isDirected())
			Rf_error("Reciprocity only make sense for directed networks");

		double rec = 0.0;
		int from, to;
		boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
		for(int i=0;i<edges->size();i++){
			from = (*edges)[i].first;
			to = (*edges)[i].second;
			if(from<to && net.hasEdge(to,from))
				rec++;
		}
		std::vector<double> v(1,rec);
		this->stats=v;
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		bool addingEdge = !net.hasEdge(from,to);
		bool hasReverse = net.hasEdge(to,from);
		double change;
		if(addingEdge && hasReverse)
			change = 1.0;
		else if(!addingEdge && hasReverse)
			change = -1.0;
		else
			change = 0.0;
		this->stats[0] += change;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Reciprocity<Directed> > DirectedReciprocity;
typedef Stat<Undirected, Reciprocity<Undirected> > UndirectedReciprocity;


/*!
 * Adds a statistic for each edge in which the 'to' and 'from' nodes
 * match on a variable
 */
template<class Engine>
class NodeMatch : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	std::string variableName; /*!< the name of the matching variable */
	int varIndex; /*!< the index of the variable in the network */
	int nstats; /*!< the number of stats generated (i.e. the number of levels squared) */
	int nlevels; /*!< the number of levels of the variable */
public:
	NodeMatch(){
		variableName="";
		nstats=nlevels=varIndex = -1;
	}

	NodeMatch(std::string name){
		variableName=name;
		nstats=nlevels=varIndex = -1;
	}

	NodeMatch(List params){
		nstats=nlevels=varIndex = -1;
		try{
			variableName = as< std::string >(params[0]);
		}catch(...){
			::Rf_error("NodeMatch requires a nodal variable name");
		}
	}


	std::string name(){
		return "nodeMatch";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"nodematch."+variableName);
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		int from,to;
		int value1, value2;
		std::vector<std::string> vars = net.discreteVarNames();
		int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
		}
		if(variableIndex<0)
			::Rf_error("NodeMatch::calculate nodal attribute not found in network");
		varIndex = variableIndex;
		//nlevels = net.discreteVariableAttributes(variableIndex).labels().size();
		//nstats = nlevels*nlevels;
		nstats = 1;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size() != nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
		for(int i=0;i<edges->size();i++){
			from = (*edges)[i].first;
			to = (*edges)[i].second;
			value1 = net.discreteVariableValue(varIndex,from) - 1;
			value2 = net.discreteVariableValue(varIndex,to) - 1;
			//this->stats[value1 + nlevels*value2]++;
			if(value1==value2)
				this->stats[0]++;
		}

	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		bool addingEdge = !net.hasEdge(from,to);
		int value1 = net.discreteVariableValue(varIndex,from) - 1;
		int value2 = net.discreteVariableValue(varIndex,to) - 1;
		if(value1==value2){
			if(addingEdge)
				this->stats[0]++;
			else
				this->stats[0]--;
		}
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
		if(variable != varIndex)
			return;
		int val = net.discreteVariableValue(varIndex,vert);
		if(net.isDirected()){
			NeighborIterator it = net.outBegin(vert);
			NeighborIterator end = net.outEnd(vert);
			while(it!=end){
				int val2 = net.discreteVariableValue(varIndex,*it);
				if(val2==val)
					this->stats[0]--;
				if(val2==newValue)
					this->stats[0]++;
				it++;
			}
			it = net.inBegin(vert);
			end = net.inEnd(vert);
			while(it!=end){
				int val2 = net.discreteVariableValue(varIndex,*it);
				if(val2==val)
					this->stats[0]--;
				if(val2==newValue)
					this->stats[0]++;
				it++;
			}
		}else{
			NeighborIterator it = net.begin(vert);
			NeighborIterator end = net.end(vert);
			while(it!=end){
				int val2 = net.discreteVariableValue(varIndex,*it);
				if(val2==val)
					this->stats[0]--;
				if(val2==newValue)
					this->stats[0]++;
				it++;
			}
		}
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}

};

typedef Stat<Directed, NodeMatch<Directed> > DirectedNodeMatch;
typedef Stat<Undirected, NodeMatch<Undirected> > UndirectedNodeMatch;



/*!
 * Adds a statistic for each edge in which the 'to' and 'from' nodes
 * match on a variable
 */
template<class Engine>
class NodeMix : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	std::string variableName; /*!< the name of the matching variable */
	int varIndex; /*!< the index of the variable in the network */
	int nstats; /*!< the number of stats generated (i.e. the number of levels squared) */
	int nlevels; /*!< the number of levels of the variable */
	std::vector<std::string> levels;
public:
	NodeMix(){
		variableName="";
		nstats=nlevels=varIndex = -1;
	}

	NodeMix(std::string name){
		variableName=name;
		nstats=nlevels=varIndex = -1;
	}

	NodeMix(List params){
		nstats=nlevels=varIndex = -1;
		try{
			variableName = as< std::string >(params[0]);
		}catch(...){
			::Rf_error("NodeMatch requires a nodal variable name");
		}
	}


	std::string name(){
		return "nodeMix";
	}

	int getIndex(int i,int j){
		int c;
		if(i>j){
			c=i;
			i=j;
			j=c;
		}
		c = 0;
		for(int k=1;k<=i;k++){
			c +=  nlevels - k;
		}
		return c + j;
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(nstats,"");
        for(int i=0;i<levels.size();i++){
        	for(int j=i;j<levels.size();j++){
        		std::string name = "nodemix." + levels.at(j) + "." + levels.at(i);
        		statnames.at(getIndex(i,j)) = name;
        	}
        }
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		int from,to;
		int value1, value2;
		std::vector<std::string> vars = net.discreteVarNames();
		int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
		}
		if(variableIndex<0)
			::Rf_error("NodeMatch::calculate nodal attribute not found in network");
		varIndex = variableIndex;
		levels = net.discreteVariableAttributes(varIndex).labels();
		nlevels = levels.size();
		nstats = nlevels * (nlevels + 1) / 2;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size() != nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
		for(int i=0;i<edges->size();i++){
			from = (*edges)[i].first;
			to = (*edges)[i].second;
			value1 = net.discreteVariableValue(varIndex,from) - 1;
			value2 = net.discreteVariableValue(varIndex,to) - 1;
			//this->stats[value1 + nlevels*value2]++;
			this->stats[getIndex(value1,value2)]++;

		}

	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		bool addingEdge = !net.hasEdge(from,to);
		double change = addingEdge ? 1.0 : -1.0;
		int value1 = net.discreteVariableValue(varIndex,from) - 1;
		int value2 = net.discreteVariableValue(varIndex,to) - 1;
		this->stats[getIndex(value1,value2)] += change;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
		if(variable != varIndex)
			return;
		Rf_error("NodeMix unimplemented");
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}

};

typedef Stat<Directed, NodeMix<Directed> > DirectedNodeMix;
typedef Stat<Undirected, NodeMix<Undirected> > UndirectedNodeMix;

/*!
 * The counts of the number of nodes in each category (except for the specified baseIndex)
 *  of variableName.
 */
template<class Engine>
class NodeCount : public BaseStat< Engine > {
protected:
	std::string variableName,baseValue; /*!< the name of the matching variable */
	int varIndex,baseIndex; /*!< the index of the variable in the network */
	int nstats; /*!< the number of stats generated (i.e. the number of levels) */
public:
	NodeCount(){
		variableName="";
		varIndex=nstats=0;
		baseIndex = 0;
		baseValue = "";
	}

	NodeCount(std::string name){
		variableName=name;
		varIndex=nstats=0;
		baseIndex = 0;
		baseValue = "";
	}

	NodeCount(List params){
		varIndex=nstats=0;
		try{
			variableName = as< std::string >(params(0));
		}catch(...){
			::Rf_error("NodeCount requires a nodal variable name");
		}
		if(params.size() > 1){
		    try{
		        baseValue = as<std::string>(params[1]);
		    }catch(...){
		        baseValue = "";
		    }
		}else{
		    baseValue = "";
		}
	}

	std::string name(){
		return "nodeCount";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<nstats;i++){
            std::string nm = "nodecount."+variableName+"."+asString(i+1);
            statnames.push_back(nm);
        }
        return statnames;
	}
    
        
	void calculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.discreteVarNames();
		int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
		}
		if(variableIndex<0)
			::Rf_error("nodal attribute not found in network");
		varIndex = variableIndex;
		// Find which level is the base level
		std::vector<std::string> levels = net.discreteVariableAttributes(varIndex).labels();
		baseIndex=-1;
		for(int i=0;i<levels.size();i++){
		    if(levels[i] == baseValue){
		        baseIndex = i;
		    }
		}
		if(baseIndex<0){
		    baseIndex = 0;
		}
		if(baseIndex<0)
		    Rf_error("invalid baseIndex");
		int nlevels = net.discreteVariableAttributes(variableIndex).labels().size();
		nstats = nlevels-1;
		this->stats = std::vector<double>(nstats,0.0);
		if(nlevels<2){
		    ::Rf_error("NodeCount::calculate: variable has only one level, you need to remove it from the network");
		}
		if(this->thetas.size() != nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		int val = 0;
		for(int i=0;i<net.size();i++){
			val = net.discreteVariableValue(varIndex,i)-1;
		    if(val > baseIndex)
		        this->stats.at((val-1))++;
		    if(val < baseIndex)
		        this->stats.at((val))++;
		}
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
		if(variable != varIndex)
			return;
		//Check if the new value is contained in the vector of current values
		std::vector<std::string> levels = net.discreteVariableAttributes(varIndex).labels();
		int nlevels = levels.size();
		if(newValue < 1 || newValue > nlevels){
		    ::Rf_error("NodeCount::discreteVertexUpdate: new value not in levels");
		}
		
		int val = net.discreteVariableValue(varIndex,vert)-1;
		newValue--;
		if(val > baseIndex)
		    this->stats.at((val-1))--;
		if(val < baseIndex)
		    this->stats.at((val))--;
		if(newValue > baseIndex)
		    this->stats.at((newValue-1))++;
		if(newValue < baseIndex)
		    this->stats.at((newValue))++;
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}

};

typedef Stat<Directed, NodeCount<Directed> > DirectedNodeCount;
typedef Stat<Undirected, NodeCount<Undirected> > UndirectedNodeCount;







//template<class Engine>
//class LogisticModel : public Stat< Engine > {
//protected:
//	int nstats; /*!< the number of stats generated */
//	int nlevels;
//	Formula formula;
//public:
//	LogisticModel(){}
//
//	LogisticModel(List params){
//
//
//		try{
//			formula = as<Formula>(params[0]);
//		}catch(...){
//			::Rf_error("LogisticModel requires a formula");
//		}
//
//	}
//
//	virtual Stat<Engine>* create(List params) const{
//		return new LogisticModel(params);
//	}
//
//	virtual std::string name(){
//		return "logisticModel";
//	}
//
//	virtual void calculate(const BinaryNet<Engine>& net){
//
//		Formula form = clone(formula);
//
//		int from,to;
//		int value1, value2;
//		std::vector<std::string> vars = net.discreteVarNames();
//
//		DataFrame df = DataFrame();
//
//		for(int i=0;i<vars.size();i++){
//			IntegerVector tmp = wrap(net.discreteVariableValues(i));
//			tmp.attr("levels") = wrap(net.discreteVariableAttributes(i).labels());
//			tmp.attr("class") = wrap("factor");
//			df.push_back(tmp);
//		}
//		df.attr("class") = wrap("data.frame");
//		df.attr("names") = wrap(vars);
//
//		Language cenv("environment",form);
//		SEXP env = cenv.eval();
//		Language call("model.frame",form,df);
//		SEXP mf = call.eval(env);
//
//		Language call1("model.matrix",mf,mf);
//		NumericMatrix mm = call1.eval(env);
//		Language call2("model.response",mf);
//		NumericVector response = call2.eval(env);
//
///*
//		Language call3("print",mm);
//		call3.eval();
//		Language call4("print",response);
//		call4.eval();
//*/
//		std::vector<double> sums = std::vector<double>(mm.ncol());
//
//		for(int i=0;i<mm.nrow();i++){
//			if(response(i)<1.5)
//				continue;
//			for(int j=0;j<sums.size();j++)
//				sums[j] += mm(i,j);
//		}
//
//		this->stats = sums;
//		if(this->thetas.size() != sums.size())
//			this->thetas = std::vector<double>(sums.size(),0.0);
//
//
//		this->changes.resize(this->stats.size());
//		fill(this->changes.begin(),this->changes.end(),0.0);
//	}
//
//	virtual void deltaEdgeToggle(const BinaryNet<Engine>& net, int from, int to,std::vector<double>* deltas){
//		fill(deltas->begin(),deltas->end(),0.0);
//	}
//
//};
//
//

/*!
 * A logistic regression statistic for variableName regressed upon regressorName.
 */
template<class Engine>
class Logistic : public BaseStat< Engine > {
protected:
	int nstats; /*!< the number of stats generated */
	int nlevels;
	int variableIndex, regIndex,baseIndex;
	std::string variableName, regressorName, baseValue;
public:
	Logistic(){
		nstats=nlevels=variableIndex=regIndex = 0;
	}

	Logistic(std::string out,std::string reg){
		nstats=nlevels=variableIndex=regIndex = 0;
		variableName = out;
		regressorName = reg;
		baseValue = "";
	}
    
    Logistic(std::string out,std::string reg,std::string base_val){
        nstats=nlevels=variableIndex=regIndex = 0;
        variableName = out;
        regressorName = reg;
        baseValue = base_val;
    }

	Logistic(List params){
		nstats=nlevels=variableIndex=regIndex = 0;
	    if(params.size() < 2){
	        ::Rf_error("LogisticModel requires at least two arguments passed");
	    }
		try{
			variableName = as<std::string>(params[0]);
		}catch(...){
			::Rf_error("LogisticModel requires a formula");
		}
		try{
			regressorName = as<std::string>(params[1]);
		}catch(...){
			::Rf_error("LogisticModel requires a formula");
		}
		if(params.size() > 2){
		    try{
		        baseValue = as<std::string>(params[2]);
		    }catch(...){
		        baseValue = "";
		    }
		}else{
		    baseValue = "";
		}
	}

	std::string name(){
		return "logistic";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"logistic");
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.discreteVarNames();
		variableIndex = -1;
		regIndex = -1;
		baseIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
			if(vars[i] == regressorName){
					regIndex = i;
			}
		}
		if(regIndex<0 || variableIndex<0)
			Rf_error("invalid variables");
		// Find which level is the base level
		std::vector<std::string> levels = net.discreteVariableAttributes(regIndex).labels();
		    for(int i=0;i<levels.size();i++){
		        if(levels[i] == baseValue){
		            baseIndex = i;
		        }
		    }
		if(baseIndex<0){
		    baseIndex = 0;
		}
        if(baseIndex<0)
            Rf_error("invalid baseIndex");
		int nlevels = net.discreteVariableAttributes(regIndex).labels().size();
		nstats = nlevels - 1;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size() != nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		int val,val1;
		for(int i=0;i<net.size();i++){
			val = net.discreteVariableValue(variableIndex,i)-1;
			val1 = net.discreteVariableValue(regIndex,i)-1;

        	// if the variable is !=0 and the val is not equal to the baselevel we add one on
        	// note that since the base level could be anywhere in the n-vector of stats
        	// need to be a little careful
        	if(val>0){
        	    if(val1 > baseIndex)
        	        this->stats.at((val1-1))++;
        	    if(val1 < baseIndex)
        	        this->stats.at((val1))++;
        	}
		}
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){
		if(variable != variableIndex && variable != regIndex)
			return;
		int varValue = net.discreteVariableValue(variableIndex,vert)-1;
		int regValue = net.discreteVariableValue(regIndex,vert)-1;
		newValue--;
		if(variable == regIndex){
			if(varValue>0){
			    // if the variable is !=0 and the old val is not equal to the baselevel we remove one
				if(regValue > baseIndex)
					this->stats.at((regValue-1))--;
				if(regValue < baseIndex)
				    this->stats.at(regValue)--;
				// if the variable is !=0 and the new val is not equal to the baselevel we add one on
				if(newValue > baseIndex)
				    this->stats.at((newValue-1))++;
				if(newValue < baseIndex)
				    this->stats.at(newValue)++;
			}
		}else{
		    // must mean variable == variableIndex
			if(varValue>0){
			    if(regValue > baseIndex)
			        this->stats.at((regValue-1))--;
			    if(regValue < baseIndex)
			        this->stats.at(regValue)--;
			}
			if(newValue>0){
			    if(regValue > baseIndex)
			        this->stats.at((regValue-1))++;
			    if(regValue < baseIndex)
			        this->stats.at(regValue)++;
			}
		}
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Logistic<Directed> > DirectedLogistic;
typedef Stat<Undirected, Logistic<Undirected> > UndirectedLogistic;

/*!
 * A logistic regression statistic for variableName regressed the sum of neighbours regressorName.
 */
template<class Engine>
class LogisticNeighbors : public BaseStat< Engine > {
protected:
    int nstats; /*!< the number of stats generated */
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    int nlevels;
    std::vector<std::string> levels;
    int variableIndex, regIndex, baseIndex;
    std::string variableName, regressorName, baseValue;
public:
    LogisticNeighbors(){
        nstats=nlevels=variableIndex=regIndex = 0;
        levels = std::vector<std::string>(0);
    }
    
    LogisticNeighbors(std::string out,std::string reg){
        nstats=nlevels=variableIndex=regIndex = 0;
        levels = std::vector<std::string>(0);
        variableName = out;
        regressorName = reg;
        baseValue = "";
    }
    
    LogisticNeighbors(std::string out,std::string reg,std::string base_val){
        nstats=nlevels=variableIndex=regIndex = 0;
        levels = std::vector<std::string>(0);
        variableName = out;
        regressorName = reg;
        baseValue = base_val;
    }
    
    LogisticNeighbors(List params){
        nstats=nlevels=variableIndex=regIndex = 0;
        levels = std::vector<std::string>(0);
        if(params.size() < 2){
            ::Rf_error("LogisticNeighbors requires at least two arguments passed");
        }
        try{
            variableName = as<std::string>(params[0]);
        }catch(...){
            ::Rf_error("LogisticNeighbors requires a formula");
        }
        try{
            regressorName = as<std::string>(params[1]);
        }catch(...){
            ::Rf_error("LogisticNeighbors requires a formula");
        }
        if(params.size() > 2){
            try{
                baseValue = as<std::string>(params[2]);
            }catch(...){
                baseValue = "";
            }
        }else{
            baseValue = "";
        }
    }
    
    std::string name(){
        return "logisticNeighbors";
    }
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<levels.size();i++){
            if(i != baseIndex){
                std::string nm = "logisticNeighbors.level."+levels.at(i);
                statnames.push_back(nm);
            }
        }
        return statnames;
    }
    
    void calculate(const BinaryNet<Engine>& net){
        std::vector<std::string> vars = net.discreteVarNames();
        variableIndex = -1;
        regIndex = -1;
        baseIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == variableName)
                variableIndex = i;
            if(vars[i] == regressorName)
                regIndex = i;
        }
        if(regIndex<0 || variableIndex<0)
            Rf_error("invalid variables");
        // Find which level is the base level
        levels = net.discreteVariableAttributes(regIndex).labels();
        for(int i=0;i<levels.size();i++){
            if(levels[i] == baseValue){
                baseIndex = i;
            }
        }
        if(baseIndex<0){
            baseIndex = 0;
        }
        if(baseIndex<0)
            Rf_error("invalid baseIndex");
        // Get the neighbours values for each node
        int nlevels = levels.size();
        nstats = nlevels - 1;
        this->stats = std::vector<double>(nstats,0.0);
        if(this->thetas.size() != nstats)
            this->thetas = std::vector<double>(nstats,0.0);
        int val,val1;
        
        for(int i=0;i<net.size();i++){
            NeighborIterator it = net.begin(i);
            NeighborIterator end = net.end(i);
            val = net.discreteVariableValue(variableIndex,i)-1;
            
            if(val>0){
                while(it != end){
                    val1 = net.discreteVariableValue(regIndex,*it)-1;
                    if(val1>baseIndex){
                        this->stats.at(val1-1)++;
                    }
                    if(val1<baseIndex){
                        this->stats.at(val1)++;
                    }
                    it++;
                }
            }
        }
    }
    
    void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
                              int variable, int newValue){
        
        if(variable != variableIndex && variable != regIndex)
            return;
        
        int varValue = net.discreteVariableValue(variableIndex,vert)-1;
        int regValue = net.discreteVariableValue(regIndex,vert)-1;
        newValue--;
        if(variable == regIndex){
            // need to step trough all neighbours and recalculate the effect
            NeighborIterator it = net.begin(vert);
            NeighborIterator end = net.end(vert);
            while(it != end){
                int neighbor_varVal = net.discreteVariableValue(variableIndex,*it)-1;
                if(neighbor_varVal>0){
                    if(regValue > baseIndex)
                        this->stats.at((regValue-1))--;
                    if(regValue < baseIndex)
                        this->stats.at(regValue)--;
                    // if the variable is !=0 and the new val is not equal to the baselevel we add one on
                    if(newValue > baseIndex)
                        this->stats.at((newValue-1))++;
                    if(newValue < baseIndex)
                        this->stats.at(newValue)++;
                }
                it++;
            }
        }
        
        if(variable == variableIndex){
            // need to step through all neighbours and recalculate the effect
            // remove old values, and add back in new values - thats why there are two if clauses
            NeighborIterator it = net.begin(vert);
            NeighborIterator end = net.end(vert);
            
            while(it != end){
                int neighbor_regVal = net.discreteVariableValue(regIndex,*it)-1;
                if(varValue > 0){
                    if(neighbor_regVal < baseIndex)
                        this->stats.at(neighbor_regVal)--;
                    if(neighbor_regVal > baseIndex)
                        this->stats.at((neighbor_regVal-1))--;
                }
                if(newValue > 0){
                    if(neighbor_regVal < baseIndex)
                        this->stats.at(neighbor_regVal)++;
                    if(neighbor_regVal > baseIndex)
                        this->stats.at((neighbor_regVal-1))++;
                }
                it++;
            }
        }
    }
    
    //Need  to ADD this in now !
    void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
        bool addingEdge = !net.hasEdge(from,to);
        //get the variables
        int fromVarValue = net.discreteVariableValue(variableIndex,from)-1;
        int fromRegValue = net.discreteVariableValue(regIndex,from)-1;
        int toVarValue = net.discreteVariableValue(variableIndex,to)-1;
        int toRegValue = net.discreteVariableValue(regIndex,to)-1;
        // If we are removing the edge may remove a value for both from  and to
        int add = net.hasEdge(from,to)?-1:1;
        if(fromVarValue>0){
           if(toRegValue>baseIndex)
               this->stats.at(toRegValue-1)+=add;
           if(toRegValue<baseIndex)
               this->stats.at(toRegValue)+=add;
        }
        if(toVarValue>0){
           if(fromRegValue>baseIndex)
               this->stats.at(fromRegValue-1)+=add;
           if(fromRegValue<baseIndex)
               this->stats.at(fromRegValue)+=add;
        }
    }
    
    void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
                            int variable, double newValue){}
};
typedef Stat<Directed, LogisticNeighbors<Directed> > DirectedLogisticNeighbors;
typedef Stat<Undirected, LogisticNeighbors<Undirected> > UndirectedLogisticNeighbors;

/**
 * Log variance of the degrees
 */
template<class Engine>
class DegreeDispersion : public BaseStat< Engine > {
protected:
	double slogfact;
	double ssq;
	double sum;
	double n;
public:

	DegreeDispersion(){
		slogfact = ssq = sum = n = 0.0;
	}

	/*!
	 * \param params
	 */
	DegreeDispersion(List params){
		slogfact = ssq = sum = n = 0.0;

	}

	std::string name(){
		return "degreeDispersion";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"degreeDispersion");
        return statnames;
	}

	double lfact(double d){return log(d + 1.0);}//R::lgammafn(d+1.0);}

	double dt(double d){return d;}

	void calculate(const BinaryNet<Engine>& net){
		int nstats = 1;

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		n = net.size();
		sum = 0.0;
		ssq = 0.0;
		slogfact = 0.0;
		double mean = 0.0;
		double var = 0.0;
		//double skew = 0.0;
		double deg = 0.0;

		for(int i=0;i<n;i++){
			//inVar += pow(net.indegree(i) - expectedDegree,2.0);
			//outVar += pow(net.outdegree(i) - expectedDegree,2.0);
			if(net.isDirected())
				deg = (net.outdegree(i) + net.indegree(i)) ;
			else
				deg = net.degree(i);
			sum += dt(deg);
			ssq += pow( dt(deg), 2.0);
			slogfact += lfact(deg);
		}

		mean = sum / n;
		var = ssq/n-pow(mean,2.0);

		this->stats[0] = log(var) - log(mean);//log(var);
	}


	void dyadUpdate(const BinaryNet<Engine>& net,
			int from, int to){
		double toDeg;
		double fromDeg;
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = 2.0*(addingEdge - 0.5);

		if(net.isDirected()){
			toDeg = (net.indegree(to) + net.outdegree(to));
			fromDeg = (net.indegree(from) + net.outdegree(from));
		}else{
			toDeg = net.degree(to);
			fromDeg = net.degree(from);

		}
		sum += dt(toDeg+edgeChange) +
				dt(fromDeg+edgeChange) -
				dt(toDeg) - dt(fromDeg);
		//sum += log(toDeg+edgeChange+1.0) +
		//		log(fromDeg+edgeChange+1.0) -
		//		log(toDeg+1.0) - log(fromDeg+1.0);
		ssq += pow(dt(toDeg+edgeChange),2.0) +
				pow(dt(fromDeg+edgeChange),2.0) -
				pow(dt(toDeg),2.0) - pow(dt(fromDeg),2.0);
		slogfact += lfact(toDeg+edgeChange) +
				lfact(fromDeg+edgeChange) -
				lfact(toDeg) - lfact(fromDeg);
		double mean = sum / n;
		double var = ssq/n-pow(mean,2.0);
		//double skew = (scube/n - 3.0*mean*var - mean*mean*mean) / pow(var,3.0/2.0);
		//double p = mean/(n-1.0);
		this->stats[0] = log(var) - log(mean);//log(var);
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, DegreeDispersion<Directed> > DirectedDegreeDispersion;
typedef Stat<Undirected, DegreeDispersion<Undirected> > UndirectedDegreeDispersion;


template<class Engine>
class DegreeSkew : public BaseStat< Engine > {
protected:
	double scube;
	double ssq;
	double sum;
	double n;
public:

	DegreeSkew(){
		scube = ssq = sum = n = 0.0;
	}

	/*!
	 * \param params
	 */
	DegreeSkew(List params){
		scube = ssq = sum = n = 0.0;

	}

	std::string name(){
		return "degreeSkew" ;
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"degreeSkew");
        return statnames;
	}


	void calculate(const BinaryNet<Engine>& net){
		int nstats = 1;

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		n = net.size();
		sum = 0.0;
		ssq = 0.0;
		scube = 0.0;
		double mean = 0.0;
		double var = 0.0;
		double skew = 0.0;
		double deg = 0.0;

		for(int i=0;i<n;i++){
			//inVar += pow(net.indegree(i) - expectedDegree,2.0);
			//outVar += pow(net.outdegree(i) - expectedDegree,2.0);
			if(net.isDirected())
				deg = (net.outdegree(i) + net.indegree(i)) ;
			else
				deg = net.degree(i);
			sum += deg;
			ssq += pow( deg, 0.5);
			scube += pow( deg, 3.0);
		}

		mean = sum / n;
		//var = ssq/n-pow(mean,2.0);
		//skew = (scube / n - 3.0 * mean * var - pow(mean,3.0));
		//this->stats[0] = log(skew) - log(mean);
		this->stats[0] = (ssq/n) - sqrt(mean);
	}


	void dyadUpdate(const BinaryNet<Engine>& net,
			int from, int to){
		double toDeg;
		double fromDeg;
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = 2.0*(addingEdge - 0.5);

		if(net.isDirected()){
			toDeg = (net.indegree(to) + net.outdegree(to));
			fromDeg = (net.indegree(from) + net.outdegree(from));
		}else{
			toDeg = net.degree(to);
			fromDeg = net.degree(from);

		}
		sum += toDeg+edgeChange +
				fromDeg+edgeChange -
				toDeg - fromDeg;
		//sum += log(toDeg+edgeChange+1.0) +
		//		log(fromDeg+edgeChange+1.0) -
		//		log(toDeg+1.0) - log(fromDeg+1.0);
		ssq += pow(toDeg+edgeChange,0.5) +
				pow(fromDeg+edgeChange,0.5) -
				pow(toDeg,0.5) - pow(fromDeg,0.5);
		scube += pow(toDeg+edgeChange,3.0) +
				pow(fromDeg+edgeChange,3.0) -
				pow(toDeg,3.0) - pow(fromDeg,3.0);
		double mean = sum / n;
		//double var = ssq/n-pow(mean,2.0);
		//double skew = (scube / n - 3.0 * mean * var - pow(mean,3.0));
		this->stats[0] = (ssq/n) - sqrt(mean);
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, DegreeSkew<Directed> > DirectedDegreeSkew;
typedef Stat<Undirected, DegreeSkew<Undirected> > UndirectedDegreeSkew;


/**
 * sqrt(E(degree)) - E(sqrt(degree))
 *
 * a robust stat for degree overdispersion.
 */
template<class Engine>
class DegreeSpread : public BaseStat< Engine > {
protected:
	double scube;
	double ssq;
	double sum;
	double n;
public:

	DegreeSpread(){
		scube = ssq = sum = n = 0.0;
	}

	/*!
	 * \param params
	 */
	DegreeSpread(List params){
		scube = ssq = sum = n = 0.0;
	}


	std::string name(){
		return "degreeSpread";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"degreeSpread");
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		int nstats = 1;

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		n = net.size();
		sum = 0.0;
		ssq = 0.0;
		scube = 0.0;
		double mean = 0.0;
		double var = 0.0;
		//double skew = 0.0;
		double deg = 0.0;

		for(int i=0;i<n;i++){
			//inVar += pow(net.indegree(i) - expectedDegree,2.0);
			//outVar += pow(net.outdegree(i) - expectedDegree,2.0);
			if(net.isDirected())
				deg = (net.outdegree(i) + net.indegree(i)) / 2.0;
			else
				deg = net.degree(i);
			deg = sqrt(deg);
			sum += deg;
			ssq += pow( deg, 2.0);
			scube += pow( deg, 3.0);
		}

		mean = sum / n;
		var = ssq/n-pow(mean,2.0);
		//skew = (scube/n - 3.0*mean*var - mean*mean*mean) / pow(var,3.0/2.0);
		//double p = mean/(n-1.0);
		this->stats[0] = log(sqrt(ssq/n)) - log(mean) ;
		//this->stats[0] = var;
	}


	void dyadUpdate(const BinaryNet<Engine>& net,
			int from, int to){
		double toDeg;
		double fromDeg;
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = 2.0*(addingEdge - 0.5);

		if(net.isDirected()){
			toDeg = (net.indegree(to) + net.outdegree(to))/2.0;
			fromDeg = (net.indegree(from) + net.outdegree(from))/2.0;
			edgeChange = edgeChange / 2.0;
		}else{
			toDeg = net.degree(to);
			fromDeg = net.degree(from);

		}
		sum += pow(toDeg+edgeChange,.5) +
				pow(fromDeg+edgeChange,.5) -
				pow(toDeg,.5) - pow(fromDeg,.5);
		//sum += log(toDeg+edgeChange+1.0) +
		//		log(fromDeg+edgeChange+1.0) -
		//		log(toDeg+1.0) - log(fromDeg+1.0);
		ssq += pow(toDeg+edgeChange,1.0) +
				pow(fromDeg+edgeChange,1.0) -
				pow(toDeg,1.0) - pow(fromDeg,1.0);
		scube += pow(toDeg+edgeChange,2.0) +
				pow(fromDeg+edgeChange,2.0) -
				pow(toDeg,2.0) - pow(fromDeg,2.0);
		double mean = sum / n;
		//double var = ssq/n-pow(mean,2.0);
		//double skew = (scube/n - 3.0*mean*var - mean*mean*mean) / pow(var,3.0/2.0);
		//double p = mean/(n-1.0);
		this->stats[0] = log(sqrt(ssq/n)) - log(mean) ;
		//this->stats[0] = var;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, DegreeSpread<Directed> > DirectedDegreeSpread;
typedef Stat<Undirected, DegreeSpread<Undirected> > UndirectedDegreeSpread;


template<class Engine>
class LogDegreeMoment : public BaseStat< Engine > {
protected:
	std::vector<int> moments;
	EdgeDirection direction;
public:

	LogDegreeMoment(){
		direction=UNDIRECTED;
	}
	LogDegreeMoment(std::vector<int> mom){
		direction = UNDIRECTED;
		moments = mom;
	}
	/*!
	 * \param params
	 */
	LogDegreeMoment(List params){
		try{
			moments = as< std::vector<int> >(params(0));
		}catch(...){
			::Rf_error("error");
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
		return "logDegreeMoment";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<moments.size();i++){
            int d = moments[i];
            std::string nm = "logDegreeMoment."+asString(d);
            statnames.push_back(nm);
        }
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		int nstats = moments.size();

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		double sum = 0.0;
		//double skew = 0.0;
		double deg = 0.0;
		int n = net.size();
		for(int i=0;i<n;i++){
			//inVar += pow(net.indegree(i) - expectedDegree,2.0);
			//outVar += pow(net.outdegree(i) - expectedDegree,2.0);
			if(net.isDirected())
				deg = (net.outdegree(i) + net.indegree(i)) / 2.0;
			else
				deg = net.degree(i);
			deg = log(deg + 1);
			sum += deg;
			for(int j=0;j<moments.size();j++){
				this->stats.at(j) += pow(deg,moments.at(j));
			}
		}
	}


	void dyadUpdate(const BinaryNet<Engine>& net,
			int from, int to){
		double toDeg;
		double fromDeg;
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = 2.0*(addingEdge - 0.5);

		if(net.isDirected()){
			toDeg = (net.indegree(to) + net.outdegree(to))/2.0;
			fromDeg = (net.indegree(from) + net.outdegree(from))/2.0;
			edgeChange = edgeChange / 2.0;
		}else{
			toDeg = net.degree(to);
			fromDeg = net.degree(from);
		}
		for(int j=0;j<moments.size();j++){
			this->stats.at(j) += pow(log(toDeg + edgeChange + 1), moments.at(j)) - pow(log(toDeg+1), moments.at(j));
			this->stats.at(j) += pow(log(fromDeg + edgeChange + 1), moments.at(j)) - pow(log(fromDeg+1), moments.at(j));
		}
		//this->stats[0] = var;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, LogDegreeMoment<Directed> > DirectedLogDegreeMoment;
typedef Stat<Undirected, LogDegreeMoment<Undirected> > UndirectedLogDegreeMoment;



/*!
 * Adds a stat for the counts of degrees
 */
template<class Engine>
class Degree : public BaseStat< Engine > {
protected:
	EdgeDirection direction;
	std::vector<int> degrees;
public:

	Degree(){
		direction = UNDIRECTED;
	}

	Degree(std::vector<int> deg){
		direction = UNDIRECTED;
		degrees = deg;
	}

	/*!
	 * \param params 	a list of length 1, the first element of which is an integer vector
	 * 					of degrees
	 */
	Degree(List params){
		try{
			degrees = as< std::vector<int> >(params(0));
		}catch(...){
			::Rf_error("error");
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
		return "degree";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<degrees.size();i++){
            int d = degrees[i];
            std::string nm = "degree."+asString(d);
            statnames.push_back(nm);
        }
        return statnames;
	}
    
    
	void calculate(const BinaryNet<Engine>& net){
		int nstats = degrees.size();
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		double n = net.size();
		for(int i=0;i<n;i++){
			for(int j=0;j<nstats;j++){
				if(net.isDirected()){
					if(direction==UNDIRECTED){
						this->stats[j] += (net.outdegree(i) + net.indegree(i)) == degrees[j];
					}else if(direction==OUT)
						this->stats[j] += net.outdegree(i) == degrees[j];
					else if(direction==IN)
						this->stats[j] += net.indegree(i) == degrees[j];
				}else{
					this->stats[j] += net.degree(i) == degrees[j];
				}
			}
		}
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		int change = !net.hasEdge(from,to) ? 1 : -1;
		int fromDegree = 0;
		int fromDegreeNew = 0;
		int toDegree = 0;
		int toDegreeNew = 0;
		if(net.isDirected()){
			if(direction==UNDIRECTED){
				fromDegree = net.outdegree(from) + net.indegree(from);
				toDegree = net.outdegree(to) + net.indegree(to);
				fromDegreeNew = change;
				toDegreeNew = change;
			}else if(direction==OUT){
				fromDegree = net.outdegree(from);
				toDegree = net.outdegree(to);
				fromDegreeNew = change;
			}else if(direction==IN){
				fromDegree += net.indegree(from);
				toDegree += net.indegree(to);
				toDegreeNew = change;
			}
		}else{
			fromDegree = net.degree(from);
			toDegree = net.degree(to);
			fromDegreeNew = change;
			toDegreeNew = change;
		}
		toDegreeNew += toDegree;
		fromDegreeNew += fromDegree;

		for(int j=0;j<degrees.size();j++){
			if(degrees[j] == fromDegree)
				this->stats[j]--;
			if(degrees[j] == toDegree)
				this->stats[j]--;
			if(degrees[j] == fromDegreeNew)
				this->stats[j]++;
			if(degrees[j] == toDegreeNew)
				this->stats[j]++;
		}

	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Degree<Directed> > DirectedDegree;
typedef Stat<Undirected, Degree<Undirected> > UndirectedDegree;


template<class Engine>
class DegreeCrossProd : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	double nEdges;
	double crossProd;
public:

	DegreeCrossProd(){
		crossProd = nEdges = 0.0;
	}

	/*!
	 * \param params
	 */
	DegreeCrossProd(List params){
		nEdges = crossProd = 0.0;
	}

	std::string name(){
		return "degreeCrossProd";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"degreeCrossProd");
        return statnames;
	}


	void calculate(const BinaryNet<Engine>& net){
		int nstats = 1;

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		nEdges = net.nEdges();
		crossProd = 0.0;
		boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net.edgelist();

		std::vector< std::pair<int,int> >::iterator it = edges->begin();
		while(it != edges->end()){
			crossProd += net.degree((*it).first) * net.degree((*it).second);
			it++;
		}
		if(nEdges==0)
			this->stats[0] = 0;
		else
			this->stats[0] = crossProd / nEdges;
	}


	void dyadUpdate(const BinaryNet<Engine>& net,
			int from, int to){
		double toDeg;
		double fromDeg;
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = 2.0*(addingEdge - 0.5);

		if(addingEdge)
			crossProd += (net.degree(from) + 1.0) * (net.degree(to) + 1.0);
		else
			crossProd -= net.degree(from) * net.degree(to);

		NeighborIterator it = net.begin(from);
		NeighborIterator end = net.end(from);
		double deg = net.degree(from);
		while(it!=end){
			double deg2 = net.degree(*it);
			if(addingEdge)
				crossProd += deg2;//(deg+1.0)*deg2 - deg*deg2
			else if(*it != to)
				crossProd -= deg2;// (deg-1.0)*deg - deg*deg2;

			it++;
		}

		it = net.begin(to);
		end = net.end(to);
		deg = net.degree(to);
		while(it!=end){
			double deg2 = net.degree(*it);
			if(addingEdge)
				crossProd += deg2;//(deg+1.0)*deg2 - deg*deg2
			else if(*it != from)
				crossProd -= deg2;// (deg-1.0)*deg - deg*deg2;

			it++;
		}
		nEdges += edgeChange;
		if(nEdges==0)
			this->stats[0] = 0;
		else
			this->stats[0] = crossProd / nEdges;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, DegreeCrossProd<Directed> > DirectedDegreeCrossProd;
typedef Stat<Undirected, DegreeCrossProd<Undirected> > UndirectedDegreeCrossProd;

template<class Engine>
class DegreeChangeCounter : public BaseStat< Engine > {
protected:
	double scube;
	double ssq;
	double sum;
	double n;
public:

	DegreeChangeCounter(){
		scube = ssq = sum = n = 0.0;
	}

	/*!
	 * \param params
	 */
	DegreeChangeCounter(List params){
		scube = ssq = sum = n = 0.0;
	}


	std::string name(){
		return "degreeChangeCounter";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"degreeChangeCounter");
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		int nstats = 1;

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		n = net.size();
		sum = 0.0;
		ssq = 0.0;
		scube = 0.0;
		double mean = 0.0;
		double var = 0.0;
		//double skew = 0.0;
		double deg = 0.0;

		for(int i=0;i<n;i++){
			//inVar += pow(net.indegree(i) - expectedDegree,2.0);
			//outVar += pow(net.outdegree(i) - expectedDegree,2.0);
			if(net.isDirected())
				deg = (net.outdegree(i) + net.indegree(i)) / 2.0;
			else
				deg = net.degree(i);
			deg = sqrt(deg);
			sum += deg;
			ssq += pow( deg, 2.0);
			scube += pow( deg, 3.0);
		}

		mean = sum / n;
		var = ssq/n-pow(mean,2.0);
		//skew = (scube/n - 3.0*mean*var - mean*mean*mean) / pow(var,3.0/2.0);
		//double p = mean/(n-1.0);
		this->stats[0] = ssq;//log(sqrt(ssq/n)) - log(mean) ;
		//this->stats[0] = var;
	}


	void dyadUpdate(const BinaryNet<Engine>& net,
			int from, int to){
		double toDeg;
		double fromDeg;
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = 2.0*(addingEdge - 0.5);

		if(net.isDirected()){
			toDeg = (net.indegree(to) + net.outdegree(to))/2.0;
			fromDeg = (net.indegree(from) + net.outdegree(from))/2.0;
			edgeChange = edgeChange / 2.0;
		}else{
			toDeg = net.degree(to);
			fromDeg = net.degree(from);

		}
		sum += pow(toDeg+edgeChange,.5) +
				pow(fromDeg+edgeChange,.5) -
				pow(toDeg,.5) - pow(fromDeg,.5);
		//sum += log(toDeg+edgeChange+1.0) +
		//		log(fromDeg+edgeChange+1.0) -
		//		log(toDeg+1.0) - log(fromDeg+1.0);
		ssq += pow(toDeg+edgeChange,1.0) +
				pow(fromDeg+edgeChange,1.0) -
				pow(toDeg,1.0) - pow(fromDeg,1.0);
		scube += pow(toDeg+edgeChange,2.0) +
				pow(fromDeg+edgeChange,2.0) -
				pow(toDeg,2.0) - pow(fromDeg,2.0);
		double mean = sum / n;
		//double var = ssq/n-pow(mean,2.0);
		//double skew = (scube/n - 3.0*mean*var - mean*mean*mean) / pow(var,3.0/2.0);
		//double p = mean/(n-1.0);
		this->stats[0] = ssq;//log(sqrt(ssq/n)) - log(mean) ;
		//this->stats[0] = var;
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
			int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, DegreeChangeCounter<Directed> > DirectedDegreeChangeCounter;
typedef Stat<Undirected, DegreeChangeCounter<Undirected> > UndirectedDegreeChangeCounter;

/*!
 * Hamming distance - the number of edges that the network differs by 
 */
template<class Engine>
class Hamming : public BaseStat< Engine > {
protected:
    //List edges;
    boost::shared_ptr< std::vector< std::pair<int,int> > > edges;
    boost::shared_ptr< BinaryNet<Engine>> compareNet;
public:
    Hamming(){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas =t;
    }
    
    Hamming(Rcpp::List params){
        if (params.size() < 2) {
            ::Rf_error("Insufficient parameters passed to HammingOffset constructor");
        }
        if (!Rcpp::is<Rcpp::NumericMatrix>(params(0)) and !Rcpp::is<Rcpp::IntegerMatrix>(params(0))) {
            ::Rf_error("Parameter should be an integer of numeric matrix, type passed was this: %s", Rcpp::type2name(params(0)));
            
        }
        
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas=t;
        
        // Keep the edge list and a network version of it
        Rcpp::NumericMatrix edgeList = params(0);
        this->compareNet.reset(new BinaryNet<Engine>(
                Rcpp::as<Rcpp::IntegerMatrix>(params(0)),
                Rcpp::as<int>(params(1))
        ));
        
        boost::shared_ptr< std::vector< std::pair<int,int> > > edges_tmp(new std::vector<std::pair<int,int> >());
        edges_tmp->reserve(edgeList.nrow());
        for(int i=0;i<edgeList.nrow();i++){
            // Do the minus one stuff since this is a R interface
            int from = edgeList(i,0)-1;
            int to = edgeList(i,1)-1;
            if(from < 0 || to<0)
                Rf_error("Edgelist indices out of range");
            std::pair<int,int> p = std::make_pair(from,to);
            edges_tmp->push_back(p);
        }
        this->edges = edges_tmp;
    }
    
    std::string name(){
        return "hamming";
    }
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"hamming");
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        // Start by assuming they are completely the same
        // Step through the edge list to calculate the edges that are missing
        // Then use nEdges to get the additional edges
        std::vector<double> v(1,0.0);
        std::vector< std::pair<int,int> > ::iterator it = this->edges->begin();
        int shared = 0;
        while(it != this->edges->end()){
            if(!net.hasEdge(it->first,it->second)){
                v[0] += 1;
            }else{
                shared +=1;
            }
            it++;
        }
        
        // add the surplus edges to the differing count
        v[0] += (net.nEdges() - shared);
        this->stats = v;
        return;
    }
    
    void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
    
        // Check if edge is in the compareNet
        int is_in_net = net.hasEdge(from,to)?1:0;
        int is_in_compare_net = this->compareNet->hasEdge(from,to)?1:0;
        this->stats[0] += (is_in_compare_net == is_in_net)?1:(-1);
        return;
    }
    
    // Don't do anything - purely edge statistic
    void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
                              int variable, int newValue){}
    
    void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
                            int variable, double newValue){}
};
typedef Stat<Directed, Hamming<Directed> > DirectedHamming;
typedef Stat<Undirected, Hamming<Undirected> > UndirectedHamming;


/*!
 * ERNM's native homophily statistic. see paper
 */
template<class Engine>
class Homophily : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	std::string variableName;
	EdgeDirection direction;
	bool includeMixing;
	bool collapseLevels;

	int varIndex;

	std::vector<double> sumMix;
	std::vector<double> sumDiff;
	std::vector<double> counts;
	std::vector< std::map<int,int> > degreeCounts;
	double n;

	int nlevels;
public:
	Homophily(){
		includeMixing = false;
		collapseLevels = true;
		nlevels = varIndex = 0;
		n = 0.0;
		direction = UNDIRECTED;
	}

	Homophily(std::string name){
		variableName = name;
		includeMixing = false;
		collapseLevels = true;
		nlevels = varIndex = 0;
		n = 0.0;
		direction = UNDIRECTED;
	}

	Homophily(std::string name, EdgeDirection dir, bool mix, bool collapse){
		variableName = name;
		direction = dir;
		includeMixing = mix;
		collapseLevels = collapse;
		nlevels = varIndex = 0;
		n = 0.0;
	}

	/*!
	 * \param params 	a list of length 1, the first element of which is an integer vector
	 * 					of degrees
	 */
	Homophily(List params){
		nlevels = varIndex = 0;
		n = 0.0;
		int nargs = params.size();
		if(nargs==0)
			Rf_error("Homophily requires a nodal variable name");
		try{
			variableName = as< std::string >(params[0]);
		}catch(...){
			::Rf_error("Homophily requires a nodal variable name");
		}
		try{
			int tmp;
			if(nargs<2)
				tmp=0;
			else
				tmp = as< int >(params[1]);
			if(tmp==0)
				direction = UNDIRECTED;
			else if(tmp==1)
				direction = IN;
			else if(tmp==2)
				direction = OUT;
			else
				::Rf_error("invalid direction");
		}catch(...){
			::Rf_error("Homophily requires a direction");
		}

		try{
			if(nargs<3)
				collapseLevels = true;
			else
				collapseLevels = as< bool >(params[2]);
		}catch(...){
			::Rf_error("Homophily: invalid collapse levels parameter");
		}

		try{
			if(nargs<4)
				includeMixing=false;
			else
				includeMixing = as< bool >(params[3]);
			//if(includeMixing && collapseLevels)
			//	::Rf_error("collapseLevels and includeMixing may not be used together");
		}catch(...){
			::Rf_error("Homophily: invalid mixing parameter");
		}
	}


	std::string name(){
		return "homophily";
	}
    
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        statnames.assign(1,"homophily."+variableName);
        return statnames;
	}

	double degree(const BinaryNet<Engine>& net,int vertex){
		if(net.isDirected()){
			if(direction == UNDIRECTED)
				return net.indegree(vertex) + net.outdegree(vertex);
			if(direction == IN)
				return net.indegree(vertex);
			if(direction == OUT)
				return net.outdegree(vertex);
			Rf_error("error");
			return -1.0;
		}else{
			return net.degree(vertex);
		}
	}

	void calculate(const BinaryNet<Engine>& net){
		if(!net.isDirected() && direction != UNDIRECTED)
			::Rf_error("Homophily: directed statistics can not be used with undirected networks");
		std::vector<std::string> vars = net.discreteVarNames();
		int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
		}
		if(variableIndex<0)
			::Rf_error("NodeMatch::calculate nodal attribute not found in network");
		varIndex = variableIndex;
		nlevels = net.discreteVariableAttributes(variableIndex).labels().size();
		int nstats;
		if(collapseLevels)
			nstats = 1;
		else
			nstats = nlevels;

		if(includeMixing)
			nstats = nstats==1 ? 2 : nstats*nstats;

		degreeCounts = std::vector< std::map<int,int> >();
		for(int i=0; i<nlevels*nlevels;i++)
			degreeCounts.push_back(std::map<int,int>());

		//if(collapseLevels && includeLogvar && includeMixing)
		//	nstats--;

		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		n = net.size();
		//double expectedDegree = net.nEdges() / n;
		sumMix = std::vector<double>(nlevels*nlevels,0.0);

		counts = std::vector<double>(nlevels,0.0);
		for(int i = 0;i<net.size();i++){
			int val = net.discreteVariableValue(varIndex,i) - 1;
			counts[val]++;
		}
		for(int i=0;i<n;i++){
			int nodeVal = net.discreteVariableValue(varIndex,i) - 1;
			double deg = degree(net,i);
			//Set inn = Set();
			//if(direction==IN || direction == UNDIRECTED)
			//	inn = net.inneighbors(i);
			//Set outn  = Set();
			//if(direction==OUT || direction == UNDIRECTED)
			//	outn = net.outneighbors(i);

			std::vector<double> mix = std::vector<double>(nlevels,0.0);
			if(net.isDirected()){
				NeighborIterator inIt = net.inBegin(i);
				NeighborIterator inEnd = net.inEnd(i);
				NeighborIterator outIt = net.outBegin(i);
				NeighborIterator outEnd = net.outEnd(i);
				//const Set* inn = &net.inneighbors(i);
				//const Set* outn = &net.outneighbors(i);
				if(direction==IN || direction == UNDIRECTED ){
					for(;inIt!=inEnd;inIt++){
						int val = net.discreteVariableValue(varIndex,*inIt) - 1;
						mix[val]++;
					}
				}

				if(direction==OUT || direction == UNDIRECTED ){
					for(;outIt!=outEnd;outIt++){
						int val = net.discreteVariableValue(varIndex,*outIt) - 1;
						mix[val]++;
					}
				}
			}else{
				NeighborIterator it = net.begin(i);
				NeighborIterator end = net.end(i);
				//const Set* nbs = &net.neighbors(i);
				for(;it!=end;it++){
					int val = net.discreteVariableValue(varIndex,*it) - 1;
					mix[val]++;
				}
			}


			for(int j=0;j<nlevels;j++){
				int ind = nlevels*nodeVal + j;
				std::map<int,int>::iterator it = degreeCounts[ind].find(deg);
				if(it == degreeCounts[ind].end())
					degreeCounts[ind].insert(std::make_pair(deg,1));
				else
					it->second++;
				double tmp = sqrt(mix[j]) ;
						//- expectedSqrtHypergeometric(counts[j], n-counts[j],
						//indegree+outdegree);
				//- sqrt((indegree+outdegree)*counts[j] / n);
				sumMix[ind] += tmp;
			}
		}

		sumDiff = subtractExpectedCounts(sumMix,counts,degreeCounts);

		this->stats = calculateStats(sumDiff);
	}

	double calculateExpectedCount(int which,std::vector<double>& cnts,
			std::vector< std::map<int,int> >& dc){
		double result = 0.0;
		std::map<int,int>::iterator it = dc[which].begin();
		std::map<int,int>::iterator end = dc[which].end();
		int j = which % nlevels;
		while(it!=end){
			if(it->second > 0.5)
				result += it->second *
					expectedSqrtHypergeometric(cnts[j], n-cnts[j],it->first);
			it++;
		}
		return result;
	}

	std::vector<double> subtractExpectedCounts(std::vector<double>& sm,std::vector<double>& cnts,
			std::vector< std::map<int,int> >& dc){
		std::vector<double> result = sm;
		for(int i=0;i<result.size();i++){
			result[i] = result[i] - calculateExpectedCount(i,cnts, dc);
		}
		return result;
	}


	std::vector<double> calculateStats(std::vector<double> sm){
			//std::vector<double>& sums, std::vector<double>& cnts,
			//std::vector< map<int,int> >& degCnts){
		//std::vector<double> sm = subtractExpectedCounts(sums,cnts,degCnts);
		std::vector<double> res = std::vector<double>(this->stats.size(),0.0);
		int ind = 0;
		if(collapseLevels){
			double match = 0.0;
			double mix = 0.0;
			//double match2 = 0.0;
			//double mix2 = 0.0;
			for(int i=0;i<nlevels*nlevels;i++){
				int nv =i / nlevels ;
				if( nlevels*nv + nv == i){
					match += sm[i];
				}else{
					mix += sm[i];
				}
			}
			res[ind] = match;
			ind++;
			if(includeMixing){
				this->stats[ind] =mix;
				ind++;
			}
		}else if(!includeMixing){
			double tmp = 0.0;
			//double tmp2 = 0.0;
			for(int i = 0;i<nlevels;i++){
				tmp = sm[i*nlevels + i] ;
				res[ind] = tmp;
				ind++;
			}
		}else{
			double tmp = 0.0;
			for(int i = 0;i<nlevels*nlevels;i++){
				tmp = sm[i] ;
				res[ind] = tmp;
				ind++;
			}
		}
		return res;
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		bool addingEdge = !net.hasEdge(from,to);
		double edgeChange = addingEdge*2-1;
		int fromVal =  net.discreteVariableValue(varIndex,from) - 1;
		int toVal =  net.discreteVariableValue(varIndex,to) - 1;

		//std::vector<double> counts = std::vector<double>(nlevels,0.0);
		//for(int i = 0;i<net.size();i++){
		//	int val = net.discreteVariableValue(varIndex,i) - 1;
		//	counts[val]++;
		//}

		if(direction==OUT || direction == UNDIRECTED){
			double deg = degree(net,from);
			//double mix = 0.0;
			//std::vector<double> mix = std::vector<double>(nlevels,0.0);

			std::vector<double> mix = std::vector<double>(nlevels,0.0);
			if(net.isDirected()){
				NeighborIterator inIt = net.inBegin(from);
				NeighborIterator inEnd = net.inEnd(from);
				NeighborIterator outIt = net.outBegin(from);
				NeighborIterator outEnd = net.outEnd(from);
			//const Set* inn = &net.inneighbors(from);
			//const Set* outn = &net.outneighbors(from);
				if(direction==IN || direction == UNDIRECTED){
					for(;inIt!=inEnd;inIt++){
						int val = net.discreteVariableValue(varIndex,*inIt) - 1;
						mix[val]++;
					}
				}
				if(direction==OUT || direction == UNDIRECTED){
					for(;outIt!=outEnd;outIt++){
						int val = net.discreteVariableValue(varIndex,*outIt) - 1;
						mix[val]++;
					}
				}
			}else{
				NeighborIterator it = net.begin(from);
				NeighborIterator end = net.end(from);
				//const Set* nbs = &net.neighbors(from);
				for(;it!=end;it++){
					int val = net.discreteVariableValue(varIndex,*it) - 1;
					mix[val]++;
				}
			}
			for(int j=0;j<nlevels;j++){
				int ind = nlevels*fromVal + j;
				double tmp = sqrt(mix[j]) ;
				double sub = expectedSqrtHypergeometric(counts[j], n-counts[j],
										deg);
				double tmpHypo= sqrt(mix[j] + edgeChange*(j==toVal)) ;
				double subHypo = expectedSqrtHypergeometric(counts[j], n-counts[j],
								deg + edgeChange);
						//sqrt(mix[j] + edgeChange*(j==toVal)) - sqrt((indegree+outdegree + edgeChange)*counts[j] / n);
				sumMix[ind] += tmpHypo - tmp;
				sumDiff[ind] += tmpHypo - subHypo - tmp + sub;
				std::map<int,int>::iterator it = degreeCounts[ind].find(deg+
						edgeChange);
				if(it == degreeCounts[ind].end())
					degreeCounts[ind].insert(std::make_pair(deg +
							edgeChange,1));
				else
					it->second++;
				it = degreeCounts[ind].find(deg);
				if(it == degreeCounts[ind].end())
					::Rf_error("Homophily deltaEdgeChange.");
				else
					it->second--;
				if(it->second < 0.5)
					degreeCounts[ind].erase(it);
			}
		}

		if(direction==IN || direction == UNDIRECTED){
			double deg = degree(net,to);
			std::vector<double> mix = std::vector<double>(nlevels,0.0);
			if(net.isDirected()){
				NeighborIterator inIt = net.inBegin(to);
				NeighborIterator inEnd = net.inEnd(to);
				NeighborIterator outIt = net.outBegin(to);
				NeighborIterator outEnd = net.outEnd(to);
				//const Set* inn = &net.inneighbors(to);
				//const Set* outn = &net.outneighbors(to);
				if(direction==IN || direction == UNDIRECTED){
					for(;inIt!=inEnd;inIt++){
						int val = net.discreteVariableValue(varIndex,*inIt) - 1;
						mix[val]++;
					}
				}
				if(direction==OUT || direction == UNDIRECTED){
					for(;outIt!=outEnd;outIt++){
						int val = net.discreteVariableValue(varIndex,*outIt) - 1;
						mix[val]++;
					}
				}
			}else{
				NeighborIterator it = net.begin(to);
				NeighborIterator end = net.end(to);
				//const Set* nbs = &net.neighbors(to);
				for(;it!=end;it++){
					int val = net.discreteVariableValue(varIndex,*it) - 1;
					mix[val]++;
				}
			}
			for(int j=0;j<nlevels;j++){
				int ind = nlevels*toVal + j;
				double tmp = sqrt(mix[j]) ;
				double sub = expectedSqrtHypergeometric(counts[j], n-counts[j],
						deg);
				double tmpHypo= sqrt(mix[j] + edgeChange*(j==fromVal));
				double subHypo = expectedSqrtHypergeometric(counts[j], n-counts[j],
								deg + edgeChange);
				//double tmp = sqrt(mix[j]) - sqrt((indegree+outdegree)*counts[j] / n);
				//double tmpHypo= sqrt(mix[j] + edgeChange*(j==fromVal)) - sqrt((indegree+outdegree+ edgeChange)*counts[j] / n);
				sumMix[ind] += tmpHypo - tmp;
				sumDiff[ind] += tmpHypo - subHypo - tmp + sub;
				std::map<int,int>::iterator it = degreeCounts[ind].find(deg+
						edgeChange);
				if(it == degreeCounts[ind].end())
					degreeCounts[ind].insert(std::make_pair(deg +
							edgeChange,1));
				else
					it->second++;
				it = degreeCounts[ind].find(deg);
				if(it == degreeCounts[ind].end())
					::Rf_error("Homophily deltaEdgeChange.");
				else
					it->second--;
				if(it->second < 0.5)
					degreeCounts[ind].erase(it);

			}
		}

		//sumDiffHypo = subtractExpectedCounts(sumMixHypo,countsHypo,
		//		degreeCountsHypo);

		this->stats = calculateStats(sumDiff);

	}


	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){
		if(variable != varIndex)
			return;
		int curValue = net.discreteVariableValue(varIndex,vert) - 1;
		newValue--;
		std::set<int> affectedNodes;
		if(net.isDirected()){
			if(direction==OUT || direction == UNDIRECTED)
				affectedNodes.insert(net.inBegin(vert),net.inEnd(vert));
				//= net.inneighbors(vert);
			if(direction==IN || direction == UNDIRECTED){
				//std::set<int> neig = net.outneighbors(vert);
				affectedNodes.insert(net.outBegin(vert),net.outEnd(vert));
			}
		}else{
			affectedNodes.insert(net.begin(vert),net.end(vert));
			//= net.neighbors(vert);
		}


		for(std::set<int>::iterator affIt = affectedNodes.begin(); affIt!=affectedNodes.end(); affIt++){
			int node = *affIt;
			int nodeValue = net.discreteVariableValue(varIndex,node) - 1;
			double change = 0.0;

			std::vector<double> mix = std::vector<double>(nlevels,0.0);
			if(net.isDirected()){
				NeighborIterator inIt = net.inBegin(node);
				NeighborIterator inEnd = net.inEnd(node);
				NeighborIterator outIt = net.outBegin(node);
				NeighborIterator outEnd = net.outEnd(node);
				//const Set* inn = &net.inneighbors(node);
				//const Set* outn = &net.outneighbors(node);
				if(direction==IN || direction == UNDIRECTED){
					for(;inIt!=inEnd;inIt++){
						int val = net.discreteVariableValue(varIndex,*inIt) - 1;
						mix[val]++;
						if(*inIt == vert)
							change++;
					}
				}
				if(direction==OUT || direction == UNDIRECTED){
					for(;outIt!=outEnd;outIt++){
						int val = net.discreteVariableValue(varIndex,*outIt) - 1;
						mix[val]++;
						if(*outIt == vert)
							change++;
					}
				}
			}else{
				NeighborIterator it = net.begin(node);
				NeighborIterator end = net.end(node);
				//const Set* nbs = &net.neighbors(node);
				for(;it!=end;it++){
					int val = net.discreteVariableValue(varIndex,*it) - 1;
					mix[val]++;
				}
				change = 1.0;
			}
			for(int j=0;j<nlevels;j++){
				int ind = nlevels*nodeValue + j;
				double tmp = sqrt(mix[j]) ;
				double tmpHypo= sqrt(mix[j] + change*(j==newValue) - change*(j==curValue));
				sumMix[ind] += tmpHypo - tmp;
			}
		}

		double deg = degree(net,vert);

		std::vector<double> mix = std::vector<double>(nlevels,0.0);
		if(net.isDirected()){
			NeighborIterator inIt = net.inBegin(vert);
			NeighborIterator inEnd = net.inEnd(vert);
			NeighborIterator outIt = net.outBegin(vert);
			NeighborIterator outEnd = net.outEnd(vert);
			//const Set* inn = &net.inneighbors(vert);
			//const Set* outn = &net.outneighbors(vert);
			if(direction==IN || direction == UNDIRECTED){
				for(;inIt!=inEnd;inIt++){
					int val = net.discreteVariableValue(varIndex,*inIt) - 1;
					mix[val]++;
				}
			}
			if(direction==OUT || direction == UNDIRECTED){
				for(;outIt!=outEnd;outIt++){
					int val = net.discreteVariableValue(varIndex,*outIt) - 1;
					mix[val]++;
				}
			}
		}else{
			NeighborIterator it = net.begin(vert);
			NeighborIterator end = net.end(vert);
			//const Set* nbs = &net.neighbors(vert);
			for(;it!=end;it++){
				int val = net.discreteVariableValue(varIndex,*it) - 1;
				mix[val]++;
			}
		}
		for(int j=0;j<nlevels;j++){
			int curInd = nlevels*curValue + j;
			int newInd = nlevels*newValue + j;
			double tmp = sqrt(mix[j]);
			sumMix[newInd] += tmp;
			sumMix[curInd] -= tmp;
		}
		for(int j=0;j<nlevels;j++){
			int curInd = nlevels*curValue + j;
			int newInd = nlevels*newValue + j;
			std::map<int,int>::iterator it = degreeCounts[newInd].find(deg);
			if(it == degreeCounts[newInd].end())
				degreeCounts[newInd].insert(std::make_pair(deg,1));
			else
				it->second++;
			it = degreeCounts[curInd].find(deg);
			if(it == degreeCounts[curInd].end())
				::Rf_error("Homophily deltaDiscreteVertex.");
			else
				it->second--;
			if(it->second < 0.5)
				degreeCounts[curInd].erase(it);
		}


		counts[newValue]++;
		counts[curValue]--;

		sumDiff = subtractExpectedCounts(sumMix,counts,
				degreeCounts);

		this->stats = calculateStats(sumDiff);
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Homophily<Directed> > DirectedHomophily;
typedef Stat<Undirected, Homophily<Undirected> > UndirectedHomophily;


/*!
 * Differential activity by group
 *
 * E(degree | group) - E(degree)*E(# in group)
 */
template<class Engine>
class DiffActivity : public BaseStat< Engine > {
protected:
	EdgeDirection direction;
	std::string variableName;
	int varIndex;
	int nstats;

	double aveDeg;
	std::vector<double> counts;
public:

	DiffActivity(){
		varIndex = nstats = 0;
		aveDeg = 0.0;
		direction = UNDIRECTED;
	}

	DiffActivity(std::string name,EdgeDirection d){
		varIndex = nstats = 0;
		aveDeg = 0.0;
		direction = d;
		variableName = name;
	}

	DiffActivity(std::string name){
		varIndex = nstats = 0;
		aveDeg = 0.0;
		direction = UNDIRECTED;
		variableName = name;
	}


	DiffActivity(List params){
		varIndex = nstats = 0;
		aveDeg = 0.0;
		try{
			variableName = as< std::string >(params(0));
		}catch(...){
			::Rf_error("NodeCount requires a nodal variable name");
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
		return "diffActivity";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<nstats;i++){
            std::string nm = "diffActivity."+variableName+"."+asString(i+1);
            statnames.push_back(nm);
        }
        return statnames;
	}
    
    
	inline int degree(const BinaryNet<Engine>& net,int i){
		int result = 0;
		if(net.isDirected()){
			if(direction==OUT || direction==UNDIRECTED)
				result += net.outdegree(i);
			if(direction==IN || direction==UNDIRECTED)
				result += net.indegree(i);
		}else
			result = net.degree(i);
		return result;
	}

	void calculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.discreteVarNames();
		int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
		}
		if(variableIndex<0)
			::Rf_error("nodal attribute not found in network");
		varIndex = variableIndex;
		int nlevels = net.discreteVariableAttributes(variableIndex).labels().size();
		nstats = nlevels-1;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		double n = net.size();
		double degSum = 0.0;
		double deg = 0.0;
		counts = std::vector<double>(nlevels,0.0);
		for(int i=0;i<n;i++){
			deg = degree(net,i);
			degSum += deg;
			int val = net.discreteVariableValue(varIndex,i) - 1;
			counts[val]++;
			if(val<nstats)
				this->stats[val] += deg;
		}
		aveDeg = degSum / n;
		for(int i=0;i<nstats;i++)
			this->stats[i] = this->stats[i] - counts[i] * aveDeg;
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		int fromVal = net.discreteVariableValue(varIndex,from)-1;
		//int fromDeg = degree(net,from);
		int toVal = net.discreteVariableValue(varIndex,to)-1;
		//int toDeg = degree(net,to);
		int change;
		if(direction==UNDIRECTED)
			change = !net.hasEdge(from,to) ? 2 : -2;
		else
			change = !net.hasEdge(from,to) ? 1 : -1;
		double n = net.size();

		for(int i=0;i<nstats;i++){
			this->stats[i] -= counts[i] * (change / n );
		}
		aveDeg = aveDeg + change / n;

		change = change>0 ? 1 : -1;

		if( (direction==UNDIRECTED || direction==OUT) && fromVal<nstats)
			this->stats[fromVal] += change;
		if( (direction==UNDIRECTED || direction==IN) && toVal<nstats)
			this->stats[toVal] += change;

	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){
		if(variable != varIndex)
			return;
		int val = net.discreteVariableValue(varIndex,vert)-1;
		newValue--;

		double degDiff = degree(net,vert) - aveDeg;

		if(val<nstats)
			this->stats[val] -= degDiff;
		counts[val]--;

		if(newValue<nstats)
			this->stats[newValue] += degDiff;
		counts[newValue]++;
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, DiffActivity<Directed> > DirectedDiffActivity;
typedef Stat<Undirected, DiffActivity<Undirected> > UndirectedDiffActivity;

/*!
 * Main effect of a continuous covariate.
 */
/*!
 * Main effect of a continuous covariate.
 */
template<class Engine>
class NodeCov : public BaseStat< Engine > {
protected:
	EdgeDirection direction;
	std::string variableName;
	int varIndex;
	bool isDiscrete;
public:

	NodeCov(){
		varIndex =  0;
		direction = UNDIRECTED;
		isDiscrete = false;
	}

	NodeCov(std::string name,EdgeDirection d){
		varIndex = 0;
		direction = d;
		variableName = name;
		isDiscrete = false;
	}

	NodeCov(std::string name){
		varIndex = 0;
		direction = UNDIRECTED;
		variableName = name;
		isDiscrete = false;
	}


	NodeCov(List params){
		varIndex = 0;
		isDiscrete=false;
		try{
			variableName = as< std::string >(params(0));
		}catch(...){
			::Rf_error("NodeCov requires a nodal variable name");
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
		return "nodeCov";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        statnames.assign(1,"nodecov."+variableName);
        return statnames;
	}

	double getValue(const BinaryNet<Engine>& net, int ind){
		double val;
		if(isDiscrete)
			val = net.discreteVariableValue(varIndex,ind);
		else
			val = net.continVariableValue(varIndex,ind);
		return val;
	}

	void calculate(const BinaryNet<Engine>& net){
		isDiscrete = false;
		std::vector<std::string> vars = net.continVarNames();
		int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == variableName){
				variableIndex = i;
			}
		}
		if(variableIndex == -1){
			isDiscrete = true;
			vars = net.discreteVarNames();
			for(int i=0;i<vars.size();i++){
				if(vars[i] == variableName){
					variableIndex = i;
				}
			}
		}
		if(variableIndex<0)
			::Rf_error("nodal attribute not found in network");
		varIndex = variableIndex;
		int nstats = 1;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		this->stats[0] = 0;
		for(int i=0;i<net.size();i++){
			double val = getValue(net,i);
			if(net.isDirected()){
				if(direction == IN || direction == UNDIRECTED)
					this->stats[0] += val * net.indegree(i);
				if(direction == OUT || direction == UNDIRECTED)
					this->stats[0] += val * net.outdegree(i);
			}else{
				this->stats[0] += val * net.degree(i);
			}
		}
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
		if(net.isDirected()){
			if(direction == IN || direction == UNDIRECTED)
				this->stats[0] += change * getValue(net,to);
			if(direction == OUT || direction == UNDIRECTED)
				this->stats[0] += change * getValue(net,from);
		}else{
			this->stats[0] += change * (getValue(net,to)+getValue(net,from));
		}
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){
		if(isDiscrete && variable==varIndex){
			double oldValue = getValue(net,vert);
			int deg = 0;
			if(net.isDirected()){
				if(direction == IN || direction == UNDIRECTED)
					deg += net.indegree(vert);
				if(direction == OUT || direction == UNDIRECTED)
					deg += net.outdegree(vert);
			}else
				deg = net.degree(vert);
			this->stats[0] += deg*(newValue - oldValue);
		}
	}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){
		if(!isDiscrete && variable==varIndex){
			double oldValue = getValue(net,vert);
			int deg = 0;
			if(net.isDirected()){
				if(direction == IN || direction == UNDIRECTED)
					deg += net.indegree(vert);
				if(direction == OUT || direction == UNDIRECTED)
					deg += net.outdegree(vert);
			}else
				deg = net.degree(vert);
			this->stats[0] += deg*(newValue - oldValue);
		}
	}

};



typedef Stat<Directed, NodeCov<Directed> > DirectedNodeCov;
typedef Stat<Undirected, NodeCov<Undirected> > UndirectedNodeCov;



template<class Engine>
class Gwesp : public BaseStat< Engine > {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	double alpha;
	double oneexpa;
	double expa;
	std::vector< boost::container::flat_map<int,int> > sharedValues;
public:

	Gwesp() : alpha(.5), oneexpa(1.0 - exp(-alpha)), expa(exp(alpha)),
	sharedValues(){
	}

	virtual ~Gwesp(){};

	Gwesp(List params) : sharedValues(){
		try{
			alpha = as< double >(params(0));
			oneexpa = 1.0 - exp(-alpha);
			expa = exp(alpha);
		}catch(...){
			::Rf_error("gwesp requires an alpha");
		}
	}

	std::string name(){
		return "gwesp";
	}
    
    std::vector<std::string> statNames(){
        std::string a = asString(alpha);
        std::string termname = "gwesp."+a;
        std::vector<std::string> statnames(1,termname);
        return statnames;
        
	}


	//counts the number of shared neighbors between f and t.
	//in directed networks this only counts | t --> f --> neighbor --> t | cycles.
	int sharedNbrs(const BinaryNet<Engine>& net, int f, int t){
		if(!net.isDirected()){
			int tmp = f;
			f = std::min(f,t);
			t = std::max(tmp,t);
		}
		boost::container::flat_map<int,int>::iterator it = sharedValues[f].find(t);
		if(it != sharedValues[f].end()){
			return it->second;
		}
		NeighborIterator fit1, fend1, tit1, tend1;
		if(!net.isDirected()){
			fit1 = net.begin(f);
			fend1 = net.end(f);
			tit1 = net.begin(t);
			tend1 = net.end(t);
		}else{
			fit1 = net.inBegin(f);
			fend1 = net.inEnd(f);
			tit1 = net.outBegin(t);
			tend1 = net.outEnd(t);
		}
		int sn = 0;
		while(fit1 != fend1 && tit1 != tend1){
			if(*tit1 == *fit1){
				sn++;
				tit1++;
				fit1++;
			}else if(*tit1 < *fit1)
				tit1 = std::lower_bound(tit1,tend1,*fit1);
			else
				fit1 = std::lower_bound(fit1,fend1,*tit1);
		}
		return sn;
	}

	void setSharedValue(const BinaryNet<Engine>& net, int f, int t, int nbrs){
		if(!net.isDirected()){
			int tmp = f;
			f = std::min(f,t);
			t = std::max(tmp,t);
		}
		sharedValues[f][t] = nbrs;
	}
	void eraseSharedValue(const BinaryNet<Engine>& net, int f, int t){
		if(!net.isDirected()){
			int tmp = f;
			f = std::min(f,t);
			t = std::max(tmp,t);
		}
		sharedValues[f].erase(t);
	}

	virtual void vCalculate(const BinaryNet<Engine>& net){
		this->stats = std::vector<double>(1,0.0);
		if(this->thetas.size()!=1)
			this->thetas = std::vector<double>(1,0.0);
		double result = 0.0;
		sharedValues = std::vector< boost::container::flat_map<int,int> >();
		for(int i = 0 ; i<net.size();i++)
			sharedValues.push_back(boost::container::flat_map<int,int>());
		boost::shared_ptr<std::vector< std::pair<int,int> > > el = net.edgelist();
		for(int i=0;i<el->size();i++){
			int from = el->at(i).first;
			int to = el->at(i).second;
			int sn = sharedNbrs(net, from, to);
			setSharedValue(net,from,to,sn);
			result += 1.0 - pow(oneexpa,sn);
		}
		this->stats[0] = expa * result;
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		NeighborIterator fit, fend, tit, tend;
		if(!net.isDirected()){
			fit = net.begin(from);
			fend = net.end(from);
			tit = net.begin(to);
			tend = net.end(to);
		}else{
			fit = net.inBegin(from);
			fend = net.inEnd(from);
			tit = net.outBegin(to);
			tend = net.outEnd(to);
		}
		bool add = !net.hasEdge(from,to);
		double change = 2.0 * (add - 0.5);
		double delta = 0.0;
		int sn = 0;
		double mult = (1.0 - (!add ? 1.0/oneexpa : oneexpa));
		while(fit != fend && tit != tend){
			if(*tit == *fit){
				sn++;
				//tie from to --> shared neighbor
				int tnsn = sharedNbrs(net,to,*tit);
				setSharedValue(net,to,*tit,tnsn + (add ? 1 : -1));
				delta += pow(oneexpa,tnsn) * mult;//pow(oneexpa,tnsn) - pow(oneexpa,tnsn + change);

				//tie from shared neighbor --> from
				int nfsn = sharedNbrs(net,*tit,from);
				setSharedValue(net,*tit,from,nfsn + (add ? 1 : -1));
				delta += pow(oneexpa,nfsn) * mult;
				tit++;
				fit++;
			}else if(*tit < *fit)
				tit = std::lower_bound(tit,tend,*fit);
			else
				fit = std::lower_bound(fit,fend,*tit);
		}
		if(add)
			setSharedValue(net,from,to,sn);
		else
			eraseSharedValue(net,from,to);
		this->stats[0] += expa * (delta + change * (1.0 - pow(oneexpa,sn)));
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Gwesp<Directed> > DirectedGwesp;
typedef Stat<Undirected, Gwesp<Undirected> > UndirectedGwesp;




    
    
/*!
 * geometrically weighted degree term - jc
 */

template<class Engine>
class GwDegree : public BaseStat< Engine > {
protected:
    double alpha;
    EdgeDirection direction; //same way of dealing with directionality as stars term. second param 1=IN , 2 =OUT
    double oneexpa;
    double expalpha;
public:
    
    GwDegree() : alpha(.5), direction(), oneexpa(0), expalpha(0){ //jc
    }
    
    virtual ~GwDegree(){};
    
    GwDegree(List params) : oneexpa(0), expalpha(0){
        try{
            alpha = as< double >(params(0));
        }catch(...){
            ::Rf_error("gwdegree requires an alpha");
        }
        
        try{
            int tmp = as< int >(params(1));
			if(tmp==1)
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
        return "gwdegree";
    }
    
    std::vector<std::string> statNames(){
        std::string a = asString(alpha);
        std::string termname = "gwdegree."+a;
        std::vector<std::string> statnames(1,termname);
        return statnames;
        
	}
    
    virtual void vCalculate(const BinaryNet<Engine>& net){
        oneexpa = 1.0 - exp(-alpha);
        expalpha = exp(alpha);
        this->stats = std::vector<double>(1,0.0);
        if(this->thetas.size()!=1)
            this->thetas = std::vector<double>(1,0.0);
        double result = 0.0;
        if (!net.isDirected()) {
            for(int i=0;i<net.size();i++){
                result += 1.0 - pow(oneexpa,net.degree(i));
            }
        }
        else if (net.isDirected() && direction==IN) {
            for(int i=0;i<net.size();i++){
                result += 1.0 - pow(oneexpa,net.indegree(i));
            }
        } else {
            for(int i=0;i<net.size();i++){
                result += 1.0 - pow(oneexpa,net.outdegree(i));
            }
        }
        
            
        this->stats[0] = expalpha * result;
    }
    
    
    void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
        //we'll toggle the dyad betwen from and to
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5); //change in edge value (1, -1)
        double delta1 = 0.0;
        double delta2 = 0.0;
        
        if(!net.isDirected()){ //does checking every time slow down the toggling?
            delta1 = pow(oneexpa,net.degree(from)) - pow(oneexpa,net.degree(from)+change);
            delta2 = pow(oneexpa,net.degree(to)) - pow(oneexpa,net.degree(to)+change);
        }
        
        else if(net.isDirected() && direction==IN){
            delta1 = pow(oneexpa,net.indegree(to)) - pow(oneexpa,net.indegree(to)+change);
            delta2 = 0;
        }else {
            delta1 = pow(oneexpa,net.outdegree(from)) - pow(oneexpa,net.outdegree(from)+change);
            delta2 = 0;
        }
 
        this->stats[0] += expalpha*(delta1 + delta2);
    }
    
    void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
                              int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, GwDegree<Directed> > DirectedGwDegree;
typedef Stat<Undirected, GwDegree<Undirected> > UndirectedGwDegree;


template<class Engine>
class Gwdsp : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double alpha;
public:
    
    Gwdsp() : alpha(.5){
    }
    
    virtual ~Gwdsp(){};
    
    Gwdsp(List params){
        try{
            alpha = as< double >(params(0));
        }catch(...){
            ::Rf_error("gwdsp requires an alpha");
        }
    }
    
    std::string name(){
        return "gwdsp";
    }
    
    std::vector<std::string> statNames(){
        std::string a = asString(alpha);
        std::string termname = "gwdsp."+a;
        std::vector<std::string> statnames(1,termname);
        return statnames;
        
	}
    
    //counts the number of shared neighbors between f and t.
    //in directed networks this only counts | t --> f --> neighbor --> t | cycles.
    
    int sharedNbrs(const BinaryNet<Engine>& net, int f, int t){
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(f);
            fend = net.end(f);
            tit = net.begin(t);
            tend = net.end(t);
        }else{
            fit = net.inBegin(f);
            fend = net.inEnd(f);
            tit = net.outBegin(t);
            tend = net.outEnd(t);
        }
        
        int sn = 0;
        while(fit != fend && tit != tend){
            if(*tit == *fit){
                sn++;
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit++;
            else
                fit++;
        }
        return sn;
    }
    
    virtual void vCalculate(const BinaryNet<Engine>& net){
        this->stats = std::vector<double>(1,0.0);
        if(this->thetas.size()!=1)
            this->thetas = std::vector<double>(1,0.0);
        
        //for each node, how many neighbors does its neighbor have? Where end index is greater than starting index, to avoid duplicates
        
        double result = 0.0;
        
        double oneexpa = 1 - exp(-alpha);
        int n = net.size();
        //std::vector<int> dp ;
        
        for (int f=0;f<n;f++){
            std::set<int> twoaways;
            NeighborIterator fit, fend;
            if(!net.isDirected()){
                fit = net.begin(f);
                fend = net.end(f);
            }else{
                fit = net.inBegin(f);
                fend = net.inEnd(f);
            }
            while(fit != fend){
                
                NeighborIterator tit, tend;
                
                if(!net.isDirected()){
                    tit = net.begin(*fit);
                    tend = net.end(*fit);
                }else{
                    tit = net.inBegin(*fit);
                    tend = net.inEnd(*fit);
                }
                while(tit != tend){
                    if(f < *tit){
                        //dp.push_back(sharedNbrs(net,f,*tit));
                        twoaways.insert(*tit);
                        tit++;
                    }else{
                        tit++;
                    }
                }
                fit++;
            }
            std::set<int>::iterator it;//set iterator
            for (it = twoaways.begin() ; it!=twoaways.end(); ++it)                 result += 1.0 - pow(oneexpa,sharedNbrs(net,f,*it));
        }
        
        this->stats[0] = exp(alpha) * result;
    }
    
    
    void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
        double oneexpa = 1.0 - exp(-alpha);
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(from);
            fend = net.end(from);
            tit = net.begin(to);
            tend = net.end(to);
        }else{
            fit = net.inBegin(from);
            fend = net.inEnd(from);
            tit = net.outBegin(to);
            tend = net.outEnd(to);
        }
        //double add = !net.hasEdge(from,to);
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        double delta = 0.0;
        //Would this be faster if I used the edgelist as in gwesp? -jc
        //boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
        
        while(fit != fend){
            if (*fit != to){ //iterate over the neighbors of from, except for to
                int tnsn = sharedNbrs(net,*fit,to);
                delta += pow(oneexpa,tnsn) - pow(oneexpa,tnsn + change);
            }
            fit++;
        }
        while(tit != tend){
            if(*tit != from){
                int tnsn = sharedNbrs(net,from,*tit);
                delta += (pow(oneexpa,tnsn) - pow(oneexpa,tnsn + change));
            }
            tit++;
        }
        this->stats[0] += (exp(alpha)*delta);
    }
    
    
    void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert, int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Gwdsp<Directed> > DirectedGwdsp;
typedef Stat<Undirected, Gwdsp<Undirected> > UndirectedGwdsp;

//Edgewise Shared Parnters. One stat for each user-generated value.
template<class Engine>
class Esp : public BaseStat< Engine > {
    //first part of code based on the code for Degree
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    EdgeDirection direction;
    std::vector<int> esps;
    std::string variableName = "";
    bool homogeneous = false;
public:
    
    Esp(){
        direction = UNDIRECTED;
    }
    
    virtual ~Esp(){}; //stick with vcalculate?
    
    Esp(std::vector<int> esps1){
        direction = UNDIRECTED;
        esps = esps1;
    }
    
    /*!
     * \param params 	a list of length 1, the first element of which is an integer vector of edgewise shared partners
     */
    
    Esp(List params){
        try{
            esps = as< std::vector<int> >(params(0));
        }catch(...){
            ::Rf_error("Submit a valid vector of counts for esp stats");
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
        try{
            homogeneous = as< bool >(params(2));
        }catch(...){
            homogeneous = false;
        }
        try{
            variableName = as< std::string >(params(3));
        }catch(...){
            variableName = "";
        }
    }
    
    std::string name(){
        return "esp";
    }
    
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<esps.size();i++){
            int e = esps[i];
            std::string nm = "esp."+asString(e);
            statnames.push_back(nm);
        }
        return statnames;
	}
    

    
    //counts the number of shared neighbors between f and t.
    //in directed networks this only counts | t --> f --> neighbor --> t | cycles.
    int sharedNbrs(const BinaryNet<Engine>& net, int f, int t, int varIndex =-1,int varValue = -1){
        if(varIndex >= 0){
            int value1 = net.discreteVariableValue(varIndex,f) - 1;
            int value2 = net.discreteVariableValue(varIndex,t) - 1;
            if(varValue <0){
                varValue = value1;
            }

            if(value2!=varValue){
                return 0;}
        }
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(f);
            fend = net.end(f);
            tit = net.begin(t);
            tend = net.end(t);
        }else{
            fit = net.inBegin(f);
            fend = net.inEnd(f);
            tit = net.outBegin(t);
            tend = net.outEnd(t);
        }
        
        int sn = 0;
        while(fit != fend && tit != tend){
            if(*tit == *fit){
                if(varIndex <0){
                    sn++;
                }else{
                    int value1 = net.discreteVariableValue(varIndex,f) - 1;
                    int value3 = net.discreteVariableValue(varIndex,*tit) - 1;
                    if(value3 == value1){
                        sn++;
                    }
                }
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit++;
            else
                fit++;
        }
        return sn;
    }
    
    virtual void vCalculate(const BinaryNet<Engine>& net){
        int varIndex = -1;
        if(homogeneous == true){
            std::vector<std::string> vars = net.discreteVarNames();
            int variableIndex = -1;
            for(int i=0;i<vars.size();i++){
                if(vars[i] == variableName){
                    variableIndex = i;
                }
            }
            if(variableIndex<0){
                Rcpp::Rcout<<variableName;
                ::Rf_error("NodeMatch::calculate nodal attribute not found in network");
            }
            varIndex = variableIndex;
        }
        
        int nstats = esps.size();
        this->stats = std::vector<double>(nstats,0.0);
        if(this->thetas.size()!=nstats)
            this->thetas = std::vector<double>(nstats,0.0);
        
        boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
        
        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            int espi = sharedNbrs(net, from, to,varIndex);
            for(int j=0;j<nstats;j++){
                this->stats[j] += espi==esps[j];
            }
        }
    }
    
    void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
        int varIndex = -1;
        if(homogeneous == true){
            std::vector<std::string> vars = net.discreteVarNames();
            int variableIndex = -1;
            for(int i=0;i<vars.size();i++){
                if(vars[i] == variableName){
                    variableIndex = i;
                }
            }
            if(variableIndex<0){
                Rcpp::Rcout<<variableName;
                ::Rf_error("NodeMatch::calculate nodal attribute not found in network");
            }
            varIndex = variableIndex;
        }
        
        int nstats = esps.size();
        int espi = sharedNbrs(net, from, to);
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        for(int j=0;j<nstats;j++){ //main edge change from-to
            this->stats[j] += change*(espi==esps[j]);
        }
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(from);
            fend = net.end(from);
            tit = net.begin(to);
            tend = net.end(to);
        }else{
            fit = net.inBegin(from);
            fend = net.inEnd(from);
            tit = net.outBegin(to);
            tend = net.outEnd(to);
        }
        while(fit != fend && tit != tend){
            if(*tit == *fit){ //it's a shared neighbor
                int fnsn = sharedNbrs(net,from,*fit,varIndex);
                for(int j=0;j<nstats;j++){ // side edge change +/-1
                    this->stats[j] += (fnsn+change)==esps[j];
                    this->stats[j] -= fnsn==esps[j];
                }
                int tnsn = sharedNbrs(net,*fit,to,varIndex);
                for(int j=0;j<nstats;j++){ // side edge change +/-1
                    this->stats[j] += (tnsn+change)==esps[j];
                    this->stats[j] -= tnsn==esps[j];
                }
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit++;
            else
                fit++;
        }
        
    }
    
    void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
                              int variable, int newValue){
        int nstats = esps.size();
        int varIndex = -1;
        if(homogeneous == true){
            std::vector<std::string> vars = net.discreteVarNames();
            int variableIndex = -1;
            for(int i=0;i<vars.size();i++){
                if(vars[i] == variableName){
                    variableIndex = i;
                }
            }
            if(variableIndex<0){
                Rcpp::Rcout<<variableName;
                ::Rf_error("NodeMatch::calculate nodal attribute not found in network");
            }
            varIndex = variableIndex;
        }else{
            return;
        }
        
        int oldValue = net.discreteVariableValue(varIndex,vert)-1;
        newValue = newValue - 1; // to use in shared neighbor function
        boost::shared_ptr<std::vector< std::pair<int,int> > > el = net.edgelist();
        
        // Find the homogenous esp terms that it turns off:
        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            if(to==vert || from==vert){
                int old_sn = sharedNbrs(net,from,to,varIndex);
                for(int j=0;j<nstats;j++){
                    this->stats[j] -= old_sn==esps[j];
                }
                int new_sn = sharedNbrs(net,from,to,varIndex,newValue);
                for(int j=0;j<nstats;j++){
                    this->stats[j] += new_sn==esps[j];
                }
            }
            //check if vert is a shared neighbor of both from and to:
            if(net.hasEdge(from,vert) && net.hasEdge(to,vert)){
                int value_to = net.discreteVariableValue(varIndex,to)-1;
                int value_from = net.discreteVariableValue(varIndex,from)-1;
                
                if(value_to == value_from && value_to == oldValue){
                    //removing old
                    //Rcpp::Rcout<<"\n removing a ESP that vert was in homogenous shared neighbor";
                    int old_sn = sharedNbrs(net, from, to,varIndex,oldValue);
                    for(int j=0;j<nstats;j++){
                        this->stats[j] -= old_sn==esps[j];
                    }
                }
                if(value_to == value_from && value_to == newValue){
                    //adding new
                    //Rcpp::Rcout<<"\n adding a ESP that vert was in homogenous shared neighbor";
                    int new_sn = sharedNbrs(net, from, to,varIndex,newValue);
                    for(int j=0;j<nstats;j++){
                        this->stats[j] += new_sn==esps[j];
                    }
                }
            }
        }
        
        
        
        
        
    }
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, Esp<Directed> > DirectedEsp;
typedef Stat<Undirected, Esp<Undirected> > UndirectedEsp;

/*!
 * Great circle distance between two long-lat points
 */
template<class Engine>
class GeoDist : public BaseStat< Engine > {
protected:
	EdgeDirection direction;
	std::string latVarName;
	int latIndex;
	std::string longVarName;
	int longIndex;
public:

	GeoDist() : direction(UNDIRECTED), latIndex(-1), longIndex(-1){}

	virtual ~GeoDist(){};

	GeoDist(List params) : latIndex(-1), longIndex(-1){
		try{
			longVarName = as< std::string >(params(0));
		}catch(...){
			::Rf_error("The first parameter of geoDist should be the longitude variable");
		}
		try{
			latVarName = as< std::string >(params(1));
		}catch(...){
			::Rf_error("The first parameter of geoDist should be the latitude variable");
		}
		try{
			int tmp = as< int >(params(2));
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
		return "geoDist";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"geoDist");
        return statnames;
	}
    
	double dist(double th1, double ph1, double th2, double ph2)
	{
		double dx, dy, dz;
		ph1 -= ph2;
		ph1 *= (3.1415926536 / 180.0), th1 *= (3.1415926536 / 180.0), th2 *= (3.1415926536 / 180.0);

		dz = sin(th1) - sin(th2);
		dx = cos(ph1) * cos(th1) - cos(th2);
		dy = sin(ph1) * cos(th1);
		return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * 6371.0;
	}

	virtual void vCalculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.continVarNames();
		//int variableIndex = -1;
		for(int i=0;i<vars.size();i++){
			if(vars[i] == longVarName){
				longIndex = i;
			}
			if(vars[i] == latVarName){
				latIndex = i;
			}
		}
		if(latIndex<0)
			::Rf_error("latitude attribute not found in network");
		for(int i=0;i<net.size();i++){
			double deg = net.continVariableValue(latIndex,i);
			if(deg<-90 || deg>90)
				Rf_error("Latitude values out of range.");
		}

		if(longIndex<0)
			::Rf_error("longitude attribute not found in network");
		for(int i=0;i<net.size();i++){
			double deg = net.continVariableValue(longIndex,i);
			if(deg<-180 || deg>180)
				Rf_error("Longitude values out of range.");
		}

		int nstats = 1;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);
		this->stats[0] = 0;

		boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
		double result = 0.0;
		for(int i=0;i<el->size();i++){
			int from = el->at(i).first;
			int to = el->at(i).second;
			result += dist(
					net.continVariableValue(latIndex,from),
					net.continVariableValue(longIndex,from),
					net.continVariableValue(latIndex,to),
					net.continVariableValue(longIndex,to)
					);
		}
		this->stats[0] = result;
		//this->stats[0] = result / (double) net.nEdges();
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
		//double ne =net.nEdges();
/*		this->stats[0] = this->stats[0] * (ne / (ne + change)) + change * dist(
				net.continVariableValue(latIndex,from),
				net.continVariableValue(longIndex,from),
				net.continVariableValue(latIndex,to),
				net.continVariableValue(longIndex,to)
				) / ne;
*/
		this->stats[0] = this->stats[0] + change * dist(
						net.continVariableValue(latIndex,from),
						net.continVariableValue(longIndex,from),
						net.continVariableValue(latIndex,to),
						net.continVariableValue(longIndex,to)
						);
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}

};

typedef Stat<Directed, GeoDist<Directed> > DirectedGeoDist;
typedef Stat<Undirected, GeoDist<Undirected> > UndirectedGeoDist;


/*!
 * independent gaussian sufficient statistics
 */
template<class Engine>
class Gauss : public BaseStat< Engine > {
protected:
	std::vector< std::string > varNames;
	std::vector<int> indices;

public:

	Gauss() {}

	virtual ~Gauss(){};

	Gauss(List params) {
		try{
			varNames = as< std::vector<std::string> >(params(0));
		}catch(...){
			::Rf_error("The first parameter of guass should be a character vector of variable names");
		}
	}

	std::string name(){
		return "gauss";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"gauss");
        return statnames;
	}

	virtual void vCalculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.continVarNames();

		indices = std::vector<int>(varNames.size(),-1);
		for(int i=0;i<vars.size();i++){
			for(int j=0;j<varNames.size();j++){
				if(vars[i] == varNames[j]){
					indices[j] = i;
				}
			}
		}
		for(int i=0;i<varNames.size();i++)
			if(indices[i] < 0)
				::Rf_error("gauss: variable not found in network");
		int nstats = indices.size()*2;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats){
			this->thetas = std::vector<double>(nstats,-.5);
			for(int i=0;i<indices.size();i++)
				this->thetas[i] = 0.0;
		}

		for(int i=0;i<indices.size();i++){
			double s=0.0;
			double ssq=0.0;
			for(int j=0;j<net.size();j++){
				s += net.continVariableValue(indices[i], j);
				ssq += pow(net.continVariableValue(indices[i], j), 2.0);
			}
			this->stats[i] = s;// / (double) net.size();
			this->stats[indices.size() + i] = ssq;// / (double) net.size();// - pow(this->stats[i], 2.0);
		}
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){
		for(int i=0;i<indices.size();i++){
			if(indices[i] == variable){
				this->stats[i] += newValue - net.continVariableValue(variable, vert);
				this->stats[indices.size() + i] += pow(newValue, 2.0) -
					pow(net.continVariableValue(variable, vert), 2.0);
			}
		}
	}
};

typedef Stat<Directed, Gauss<Directed> > DirectedGauss;
typedef Stat<Undirected, Gauss<Undirected> > UndirectedGauss;


/*!
 * Gamma sufficient statistics
 */
template<class Engine>
class Gamma : public BaseStat< Engine > {
protected:
	std::string varName;
	int index;
	double epsilon;
public:

	Gamma() :index(-1), epsilon(1e-300){}

	virtual ~Gamma(){};

	Gamma(List params) : index(-1), epsilon(1e-300){
		try{
			varName = as< std::string >(params(0));
		}catch(...){
			::Rf_error("The first parameter of Gamma should be a variable name");
		}
	}

	std::string name(){
		return "gamma";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"gamma");
        return statnames;
	}

	virtual void vCalculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.continVarNames();
		index = indexOf(varName,vars);
		if(index < 0)
			::Rf_error("gamma: variable not found in network");

		int nstats = 2;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats){
			this->thetas = std::vector<double>(nstats,0);
		}

		double s=0.0;
		double slog=0.0;
		for(int j=0;j<net.size();j++){
			double val = net.continVariableValue(index, j);
			if(val<0)
				::Rf_error("gamma: Only defined for positive variables");
			s += val;
			slog += log(val+epsilon);
		}
		this->stats[0] = s;
		this->stats[1] = slog;

	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){
		if(index == variable){
			if(newValue < 0)
				::Rf_error("gamma update: Only defined for positive variables");
			double oldValue = net.continVariableValue(variable, vert);
			this->stats[0] += newValue - oldValue;
			this->stats[1] += log(newValue+epsilon) - log(oldValue+epsilon);
		}
	}
};

typedef Stat<Directed, Gamma<Directed> > DirectedGamma;
typedef Stat<Undirected, Gamma<Undirected> > UndirectedGamma;


template<class Engine>
class SumOfSquares : public BaseStat< Engine > {
protected:
	std::vector< std::string > varNames;
	std::vector<int> indices;

public:

	SumOfSquares() {}

	virtual ~SumOfSquares(){};

	SumOfSquares(List params) {
		try{
			std::string str = as< std::string >(params(0));
			varNames.push_back(str);
		}catch(...){
			::Rf_error("The first parameter of sumOfSquares should be a variable name");
		}
	}

	std::string name(){
		return "sumOfSquares";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"sumOfSquares");
        return statnames;
	}

	virtual void vCalculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.continVarNames();

		indices = std::vector<int>(varNames.size(),-1);
		for(int i=0;i<vars.size();i++){
			for(int j=0;j<varNames.size();j++){
				if(vars[i] == varNames[j]){
					indices[j] = i;
				}
			}
		}
		for(int i=0;i<varNames.size();i++)
			if(indices[i] < 0)
				::Rf_error("sumOfSquares: variable not found in network");
		int nstats = indices.size();
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,-.5);

		for(int i=0;i<indices.size();i++){
			double s=0.0;
			double ssq=0.0;
			for(int j=0;j<net.size();j++){
				s += net.continVariableValue(indices[i], j);
				ssq += pow(net.continVariableValue(indices[i], j), 2.0);
			}
			this->stats[i] = ssq;// / (double) net.size();// - pow(this->stats[i], 2.0);
		}
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}

	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){
		for(int i=0;i<indices.size();i++){
			if(indices[i] == variable){
				this->stats[i] += pow(newValue, 2.0) -
					pow(net.continVariableValue(variable, vert), 2.0);
			}
		}
	}
};

typedef Stat<Directed, SumOfSquares<Directed> > DirectedSumOfSquares;
typedef Stat<Undirected, SumOfSquares<Undirected> > UndirectedSumOfSquares;

/*!
 * distance measure
 */
template<class Engine>
class Dist : public BaseStat< Engine > {
protected:
	EdgeDirection direction;
	std::vector< std::string > varNames;
	std::vector<int> indices;
public:

	Dist() : direction(UNDIRECTED){}

	virtual ~Dist(){};

	Dist(List params){
		try{
			varNames = as< std::vector<std::string> >(params(0));
		}catch(...){
			::Rf_error("The first parameter of geoDist should be the longitude variable");
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
		return "dist";
	}
    
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"dist");
        return statnames;
	}

	double dist(const BinaryNet<Engine>& net, int from, int to){
		double ssq = 0.0;
		for(int j=0;j<indices.size();j++){
			ssq += pow(net.continVariableValue(indices[j],from) -
					net.continVariableValue(indices[j],to), 2.0);
		}
		return sqrt(ssq);
	}

	virtual void vCalculate(const BinaryNet<Engine>& net){
		std::vector<std::string> vars = net.continVarNames();
		//int variableIndex = -1;
		indices = std::vector<int>(varNames.size(),-1);
		for(int i=0;i<vars.size();i++){
			for(int j=0;j<varNames.size();j++){
				if(vars[i] == varNames[j]){
					indices[j] = i;
				}
			}
		}
		for(int i=0;i<varNames.size();i++)
			if(indices[i] < 0)
				::Rf_error("dist: variable not found in network");

		int nstats = 1;
		this->stats = std::vector<double>(nstats,0.0);
		if(this->thetas.size()!=nstats)
			this->thetas = std::vector<double>(nstats,0.0);


		boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
		double result = 0.0;
		for(int i=0;i<el->size();i++){
			int from = el->at(i).first;
			int to = el->at(i).second;
			result += dist(net, from,to);
		}
		this->stats[0] = result;
		//this->stats[0] = result / (double) net.nEdges();
	}


	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
		this->stats[0] = this->stats[0] + change * dist(net, from,to);
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
					int variable, int newValue){}

};

typedef Stat<Directed, Dist<Directed> > DirectedDist;
typedef Stat<Undirected, Dist<Undirected> > UndirectedDist;


/*!
 * barabasi-albert type change statistic. Order dependent
 */
template<class Engine>
class PreferentialAttachment : public BaseStat<Engine>{
public:
	PreferentialAttachment(){
	}

	/*!
	 * constructor. params is unused
	 */
	PreferentialAttachment(List params){
	}


	std::string name(){
		return "preferentialAttachment";
	}

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"preferentialAttachment");
        return statnames;
	}

	void calculate(const BinaryNet<Engine>& net){
		std::vector<double> v(1, 0.0); //should be NA
		this->stats=v;
		if(this->thetas.size()!=1){
			this->thetas = std::vector<double>(1,0.0);
		}
	}

	void dyadUpdate(const BinaryNet<Engine>& net, int from, int to){
		double direction = net.hasEdge(from,to) ? -1.0 : 1.0;
		double totDegree = net.nEdges() * 2.0;
		double deg = net.degree(to);
		if(deg < .5)
			deg = 1.0;
		//double toValue = net.degree(to) == 0 ? 0.0 : net.degree(to) / totDegree;
		if(totDegree <.5)
			totDegree = 1.0;
		this->stats[0] += direction * log( deg / totDegree);

		//this->stats[0] += direction * log((1.0 + net.degree(to)) / (1.0*net.size() + totDegree));
	}

	void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, int newValue){}
	void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
				int variable, double newValue){}
};

typedef Stat<Directed, PreferentialAttachment<Directed> > DirectedPreferentialAttachment;
typedef Stat<Undirected, PreferentialAttachment<Undirected> > UndirectedPreferentialAttachment;
    
#include <Rcpp.h>


}


#endif /* STATS_H_ */
