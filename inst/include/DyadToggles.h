/*
 * DyadToggles.h
 *
 *  Created on: Jan 9, 2014
 *      Author: ianfellows
 */

#ifndef DYADTOGGLES_H_
#define DYADTOGGLES_H_

#include "DyadToggle.h"

namespace ernm{


/*!
 * simple dyad toggle
 */
template<class  Engine >
class RandomDyad {
protected:
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
public:

	RandomDyad(){}

	RandomDyad(Rcpp::List l){}

	RandomDyad( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
	}

	virtual ~RandomDyad(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
	}
	void togglesAccepted(bool apply){}

	inline void generate(){
		//Rcpp::RNGScope scope;
		net->randomDyad(toggle[0]);
		//= std::pair<int,int>(i1,i2);
	}


	inline double logRatio(){return 0.0;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "RandomDyad";
	}
};

typedef DyadToggle<Directed, RandomDyad<Directed> > DirectedRandomDyadToggle;
typedef DyadToggle<Undirected, RandomDyad<Undirected> > UndirectedRandomDyadToggle;


/*!
 * tie - dyad toggle ergm Algorithm (Linear Time)
 */
template<class  Engine >
class TieDyadBasic{
protected:
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	double logProbRatio;
public:

	TieDyadBasic(Rcpp::List l){
		logProbRatio = 0.0;
	}

	TieDyadBasic( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		logProbRatio=0.0;
	}

	virtual ~TieDyadBasic(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
	}
	void togglesAccepted(bool apply){}

	inline void generate(){
		bool pickEdge = Rf_runif(0.0,1.0)>.5;
		double nEdges = net->nEdges();
		double ndyads = net->maxEdges();
		if(nEdges==0 && pickEdge){
			pickEdge=false;
		}

		if(pickEdge){
			toggle[0] = net->randomEdge();
			if(nEdges==1){
				this->logProbRatio = log(1.0/(ndyads + .5));
			}else{
				this->logProbRatio = log(nEdges / (ndyads + nEdges));
			}
		}else{
			net->randomDyad(toggle[0]);
			if(net->hasEdge(toggle[0].first,toggle[0].second)){
				this->logProbRatio = log((nEdges==1 ? 1.0/(ndyads + .5) :
									nEdges / (ndyads + nEdges)));
			}else{
				this->logProbRatio = log((nEdges==0 ? .5*(ndyads + 1.0) :
				        1.0 + (ndyads)/(nEdges + 1)));
			}
		}
	}


	inline double logRatio(){return this->logProbRatio;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "TieDyadBasic";
	}
};

typedef DyadToggle<Directed, TieDyadBasic<Directed> > DirectedTieDyadBasicToggle;
typedef DyadToggle<Undirected, TieDyadBasic<Undirected> > UndirectedTieDyadBasicToggle;

/*!
 * Tie - dyad toggle. constant time method.
 */
template<class  Engine >
class TieDyad {
protected:
	typedef  boost::shared_ptr< std::vector<int> > VecPtr;
	std::vector< std::pair<int,int> > toggle;		/*!< The toggle */
	boost::shared_ptr< BinaryNet<Engine> > net; 	/*!< The network */
	VecPtr edgesFrom; 								/*!< the from part of the edgelist */
	VecPtr edgesTo;								/*!< the to part of the edgelist */
	double logProbRatio; 							/*!< MCMC log likelihood ratio */
	int lastIndex; 								/*!< The index (in edgesForm/edgesTo) of the last toggle. -1 if non-edge toggle */
public:

	TieDyad(){
		logProbRatio=0.0;
		lastIndex = -1;
	};

	TieDyad(Rcpp::List l){
		logProbRatio=0.0;
		lastIndex = -1;
	}

	TieDyad( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		edgesFrom = VecPtr(new std::vector<int>());
		edgesTo = VecPtr(new std::vector<int>());
		logProbRatio=0.0;
		lastIndex = -1;
	}

	virtual ~TieDyad(){}

	/*!
	 * set the network
	 */
	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	/*!
	 * start-up the toggler
	 */
	void initialize(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net->edgelist();
		edgesFrom = VecPtr(new std::vector<int>());
		edgesTo = VecPtr(new std::vector<int>());
		edgesFrom->reserve(2*(*edges).size());
		edgesTo->reserve(2*(*edges).size());
		for(int i=0;i<edges->size();i++){
			edgesFrom->push_back((*edges)[i].first);
			edgesTo->push_back((*edges)[i].second);
		}
	}

	/*!
	 * generate a set of toggles
	 */
	inline void generate(){
		bool pickEdge = Rf_runif(0.0,1.0)>0.5;
		double nEdges = net->nEdges();
		double ndyads = net->maxEdges();
		if(nEdges==0 && pickEdge){
			pickEdge=false;
		}
		//cout << nEdges << " " << edgesFrom->size() << " " << edgesTo->size() << "\n";
		assert( (nEdges == edgesFrom->size()) && (nEdges == edgesTo->size()));
		if(pickEdge){
			int index = floor(Rf_runif(0.0,(double)edgesFrom->size()));
			toggle[0].first = (*edgesFrom)[index];
			toggle[0].second = (*edgesTo)[index];
			if(nEdges==1){
				this->logProbRatio = log(1.0/(ndyads + .5));
			}else{
				this->logProbRatio = log(nEdges / (ndyads + nEdges));
				// << nEdges<<" "<<ndyads<<" "<< this->logProbRatio<<" "<<log(nEdges / (ndyads + nEdges));
			}
			this->lastIndex=index;
		}else{
			//toggle[0] = net->randomNonEdge();
			net->randomDyad(toggle[0]);

			bool hasEdge = net->hasEdge(toggle[0].first,toggle[0].second);
			if(hasEdge){
				//just re-picking randomly is faster than find
				int index = floor(Rf_runif(0.0,(double)edgesFrom->size()));
				toggle[0].first = (*edgesFrom)[index];
				toggle[0].second = (*edgesTo)[index];
				this->lastIndex=index;
			}else{
				this->lastIndex=-1;
			}

			if(hasEdge){
				this->logProbRatio = log((nEdges==1 ? 2.0/(ndyads + 1.0) :
									nEdges / (ndyads + nEdges)));

			}else{
				this->logProbRatio = log((nEdges==0 ? .5*(ndyads + 1.0) :
				        1.0 + (ndyads)/(nEdges + 1.0)));
			}
		}
		//= std::pair<int,int>(i1,i2);
	}


	/*!
	 * log(p1/p2)
	 *
	 * \return the log likelihood ratio for the current toggles
	 */
	inline double logRatio(){return this->logProbRatio;}

	/*!
	 * get toggles
	 *
	 * \return the current set of toggles
	 */
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}
	/*!
	 * called if the dyad toggles are accepted
	 */
	void togglesAccepted(bool apply){
		if(!apply)
			return;
		if(lastIndex<0){
			// edge has been added
			edgesFrom->push_back(toggle[0].first);
			edgesTo->push_back(toggle[0].second);
		}else{
			//edge has been removed
			int size = edgesFrom->size();
			(*edgesFrom)[lastIndex] = (*edgesFrom)[size-1];
			(*edgesTo)[lastIndex] = (*edgesTo)[size-1];
			edgesFrom->pop_back();
			edgesTo->pop_back();
		}
	}

	inline std::string name(){
		return  "TieDyad";
	}
};

typedef DyadToggle<Directed, TieDyad<Directed> > DirectedTieDyadToggle;
typedef DyadToggle<Undirected, TieDyad<Undirected> > UndirectedTieDyadToggle;


/*!
 * toggle dyads connecting two nodes near each other in the graph.
 * Provides good mixing for very highly transitive graphs.
 *
 *   o
 *	|  \
 *	|    b
 *	|      \
 *	o---a---o
 *	 \      /
 *	  \    /
 *	   \  /
 *	    o   <--- a randomly chosen node
 *
 * 1. choose a random node
 * 2. choose two neighbors at random. the dyad between these two is denoted as "a"
 * 3. choose a random neighbor of the first node chosen in step 2. The connection
 *    between this node and the second node is denoted as "b"
 * 4. Toggle "a" with probability .5 and "b" with probability .5
 */
template<class  Engine >
class Neighborhood {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;

	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	bool twoSteps;
public:

	Neighborhood() : twoSteps(0){}

	Neighborhood(Rcpp::List l) : twoSteps(0){}

	Neighborhood( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		twoSteps = false;
	}

	virtual ~Neighborhood(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle.resize(1,std::make_pair(-1,-1));
	}
	void togglesAccepted(bool apply){

	}

	inline void generate(){
		int node = floor(Rf_runif(0.0,net->size()));
		generate(node);
	}

	inline void generate(int node){

		int from, to;
		//int node = floor(Rf_runif(0.0,net->size()));
		NeighborIterator it;
		NeighborIterator end;
		int degree;
		//const Set* nbs;
		if(net->isDirected()){
			//nbs = &net->outneighbors(node);
			it = net->outBegin(node);
			end = net->outEnd(node);
			degree = net->outdegree(node);
		}else{
			it = net->begin(node);
			end = net->end(node);
			degree = net->degree(node);
			//nbs = &net->neighbors(node);
		}
		if(degree<2){
			from = floor(Rf_runif(0.0,net->size()-1));
			to = floor(Rf_runif(0.0,net->size()-2));
			if(from>=node)
				node++;
			if(to>=std::min(node,from))
				to++;
			if(to>=std::max(node,from))
				to++;
		}else{
			int fromIndex = floor(Rf_runif(0.0,degree));
			std::advance(it,fromIndex);
			from = *it;
			std::advance(it,-fromIndex);
			//from = ithElement(*nbs,fromIndex);
			int toIndex = floor(Rf_runif(0.0,degree-1));
			if(toIndex>=fromIndex)
				toIndex++;
			std::advance(it,toIndex);
			to = *it;
			std::advance(it,-toIndex);
			//to = ithElement(*nbs,toIndex);
		}
		int first = from;
		int second = to;
		if(twoSteps){
			std::vector<int> v(3,0);
			v[0] = node;
			v[1] = to;
			v[2] = from;
			std::sort(v.begin(),v.end());
			//const Set* nbs1;
			NeighborIterator it1;
			NeighborIterator end1;
			int degree1;
			if(net->isDirected()){
				it1 = net->outBegin(from);
				end1 = net->outEnd(from);
				degree1 = net->outdegree(from);
				//nbs1 = &net->outneighbors(from);
			}else{
				it1 = net->begin(from);
				end1 = net->end(from);
				degree1 = net->degree(from);
				//nbs1 = &net->neighbors(from);
			}
			int reqSize = 1;
			if(net->hasEdge(from,node))
				reqSize++;
			if(net->hasEdge(from,to))
				reqSize++;
			if(degree1<reqSize){
				first = floor(Rf_runif(0.0,net->size()-3));
				if(first>=v[0])
					first++;
				if(first>=v[1])
					first++;
				if(first>=v[2])
					first++;
			}else{
				int ind = floor(Rf_runif(0.0,degree1+1-reqSize));
				//Set::iterator it = nbs1->begin();
				std::advance(it1,ind);
				if(net->hasEdge(from,std::min(to,node)) && *it1>=std::min(to,node))
					it1++;
				if(net->hasEdge(from,std::max(to,node)) && *it1>=std::max(to,node))
					it1++;
				first = *it1;
			}
		}

		toggle[0].first = first;
		toggle[0].second = second;
		twoSteps = !twoSteps;
	}


	inline double logRatio(){return 0.0;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "Neighborhood";
	}

};

typedef DyadToggle<Directed, Neighborhood<Directed> > DirectedNeighborhoodToggle;
typedef DyadToggle<Undirected, Neighborhood<Undirected> > UndirectedNeighborhoodToggle;


/*!
 * Nodal Tie-dyad toggling.
 *
 * pick a random node, and do tie-dyad toggling among (out)neighbors
 */
template<class  Engine >
class NodeTieDyad {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	boost::shared_ptr< BinaryNet<Engine> > net;
	std::vector< std::pair<int,int> > toggle;
	double logProbRatio;
public:

	NodeTieDyad() : logProbRatio(0.0){}

	NodeTieDyad(Rcpp::List l) : logProbRatio(0.0){}

	NodeTieDyad( BinaryNet<Engine> & network) : logProbRatio(0.0){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
	}

	virtual ~NodeTieDyad(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net = n;
	}

	void initialize(){
		toggle.resize(1,std::make_pair(-1,-1));
	}
	void togglesAccepted(bool apply){}

	inline void generate(){
		int node = floor(Rf_runif(0.0,net->size()));
		generate(node);
	}

	inline void generate(int node){

		int neighbor;
		NeighborIterator it;
		NeighborIterator end;
		int degree;
		if(!net->isDirected()){
			it = net->begin(node);
			end = net->end(node);
			degree = net->degree(node);
		}else{
			it = net->outBegin(node);
			end = net->outEnd(node);
			degree = net->outdegree(node);
		}
		double nedges = degree;
		double ndyads = net->size() - 1.0;
		bool pickEdge = Rf_runif(0.0,1.0)>0.5;
		if(nedges==0)
			pickEdge=false;
		bool hasEdge;
		if(pickEdge){
			hasEdge=true;
			int nbrIndex = floor(Rf_runif(0.0,(double)nedges));
			std::advance(it,nbrIndex);
			neighbor = *it;
			toggle[0].first = node;
			toggle[0].second = neighbor;
			assert(node != neighbor);

		}else{
			neighbor = floor(Rf_runif(0.0,net->size() - 1.0));
			if(neighbor>=node)
				neighbor++;
			toggle[0].first = node;
			toggle[0].second = neighbor;
			hasEdge = net->hasEdge(node,neighbor);
			assert(node != neighbor);
		}

		double tForward, tReverse;
		if(hasEdge){
			if(nedges<1.5)
				tReverse = 1.0/ndyads;
			else
				tReverse = 0.5/ndyads;
			tForward = 0.5/nedges + 0.5/ndyads;

		}else{
			if(nedges<.5)
				tForward = 1.0/ndyads;
			else
				tForward = 0.5/ndyads;
			tReverse = 0.5/(nedges + 1.0) + 0.5/ndyads;

		}
		this->logProbRatio = log( tReverse/tForward );

		//calculate log ratio for neighbor
		//do we need this?
		/*if(!net->isDirected()){
			nedges = net->degree(neighbor);
			ndyads = net->size() - 1.0;
			double tForwardNbr, tReverseNbr;
			if(hasEdge){
				if(nedges<1.5)
					tReverseNbr = 1.0/ndyads;
				else
					tReverseNbr = 0.5/ndyads;
				tForwardNbr = 0.5/nedges + 0.5/ndyads;
			}else{
				if(nedges<.5)
					tForwardNbr = 1.0/ndyads;
				else
					tForwardNbr = 0.5/ndyads;
				tReverseNbr = 0.5/(nedges + 1.0) + 0.5/ndyads;
			}
			this->logProbRatio = log( (tReverse + tReverseNbr) / (tForward + tForwardNbr) );
		}*/

	}


	inline double logRatio(){
		return logProbRatio;
	}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "NodeTieDyad";
	}

};

typedef DyadToggle<Directed, NodeTieDyad<Directed> > DirectedNodeTieDyadToggle;
typedef DyadToggle<Undirected, NodeTieDyad<Undirected> > UndirectedNodeTieDyadToggle;


/*!
 * Tetrad toggling
 *
 * Select two edges with no nodes in common, A1-A2 and B1-B2, s.t. A1-B2 and B1-A2 are not edges,
 * and propose to replace the former two by the latter two.
 */
template<class  Engine >
class Tetrad {
protected:
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	boost::shared_ptr<std::vector< std::pair<int,int> > > edges;
	int e1Index, e2Index;
public:

	Tetrad() : e1Index(0), e2Index(0){}

	Tetrad(Rcpp::List l) : e1Index(0), e2Index(0){}

	Tetrad( BinaryNet<Engine> & network) : e1Index(0), e2Index(0){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
	}

	virtual ~Tetrad(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		edges = net->edgelist();
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
	}
	void togglesAccepted(bool apply){
		if(!apply)
			return;
		edges->at(e1Index) = toggle[0];
		edges->at(e2Index) = toggle[1];
	}

	inline void generate(){
		std::pair<int,int>  e1, e2;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		int i = 0;
		double n = edges->size();
		while(i<100000){
			i++;
			e1Index = floor(Rf_runif(0,n));
			e2Index = floor(Rf_runif(0,n-1.0));
			if(e2Index>=e1Index)
				e2Index++;
			e1 = edges->at(e1Index);
			e2 = edges->at(e2Index);
			if(e1.first!=e2.first && e1.first != e2.second &&
					e1.second!=e2.first && e1.second!=e2.second &&
					!net->hasEdge(e1.first,e2.second) &&
					!net->hasEdge(e2.first,e1.second))
				break;
		}
		if(i>=100000)
			Rf_error("TetradToggle: could not find tetrad");

		toggle[0].first = e1.first;
		toggle[0].second = e2.second;
		toggle[1].first = e2.first;
		toggle[1].second = e1.second;
		toggle[2] = e1;
		toggle[3] = e2;
	}


	inline double logRatio(){return 0.0;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "Tetrad";
	}
};

typedef DyadToggle<Directed, Tetrad<Directed> > DirectedTetradToggle;
typedef DyadToggle<Undirected, Tetrad<Undirected> > UndirectedTetradToggle;

/*!
 * Respondent Driven Sampling toggling
 *
 * Does tetrad toggling among unobserved edges with prob .5 and simple toggling among
 * dyads connecting unobserved nodes.
 */
/*template<class  Engine >
class RDSBasicToggle : public AbstractDyadToggle<Engine> {
protected:
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	boost::shared_ptr<std::vector< std::pair<int,int> > > edges;
	std::vector<bool> observed;
	std::vector<int> unobservedNodes;

	bool wasTetrad;
	int e1Index, e2Index;
	int lastIndex;
public:
	DYAD_TOGGLE_BOILERPLATE(RDSBasicToggle)

	RDSBasicToggle(Rcpp::List l){}

	RDSBasicToggle( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
	}

	virtual ~RDSBasicToggle(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		if(net->isDirected())
			Rf_error("RDS is only applicable for undirected networks");
		edges = net->edgelist();

		observed = std::vector<bool>(net->size(),false);
		for(int i=0;i<edges->size();i++){
			if(!net->isMissing(edges->at(i).first,edges->at(i).second)){
				observed.at(edges->at(i).first) = true;
				observed.at(edges->at(i).second) = true;
				edges->at(i) = edges->at(edges->size()-1);
				edges->pop_back();
				i--;
			}
		}
		for(int i=0;i<net->size();i++){
			for(int j=0;j<net->discreteVarNames().size();j++){
				if(net->discreteVariableObserved(j,i))
					observed.at(i) = true;
			}
		}

		unobservedNodes.clear();
		for(int i=0;i<observed.size();i++){
			if(!observed[i])
				unobservedNodes.push_back(i);
		}
		//cout<<unobservedNodes.size()<<"\n";
	}
	void togglesAccepted(bool apply){
		if(!apply)
			return;
		if(wasTetrad){
			edges->at(e1Index) = toggle.at(0);
			edges->at(e2Index) = toggle.at(1);
		}else{
			if(lastIndex<0){
				// edge has been added
				edges->push_back(toggle[0]);

			}else{
				//edge has been removed
				int size = edges->size();
				(*edges)[lastIndex] = (*edges)[size-1];
				edges->pop_back();

			}
		}
	}

	inline void generate(){
		if(Rf_runif(0.0,1.0)<.5){
			generateUnobservedDyad();
			wasTetrad=false;
		}else{
			generateTetrad();
			wasTetrad=true;
		}
	}



	inline void generateUnobservedDyad(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		int index1 = floor(Rf_runif(0.0,(double)unobservedNodes.size()));
		int index2 = floor(Rf_runif(0.0,(double)unobservedNodes.size()-1.0));
		if(index2>=index1)
			index2++;
		int i1 = unobservedNodes[index1];
		int i2 = unobservedNodes[index2];
		toggle[0].first = i1;
		toggle[0].second = i2;
		bool hasEdge = net->hasEdge(i1,i2);
		if(hasEdge){
			double n = edges->size();
			bool found = false;
			if(!found){
				lastIndex=-1;
				for(int i=0;i<edges->size();i++){
					if((*edges)[i].first==i1 && (*edges)[i].second==i2){
						lastIndex=i;
						break;
					}
					if((*edges)[i].first==i2 && (*edges)[i].second==i1){
						lastIndex=i;
						break;
					}
					//cout<<(*edges)[i].first<<" "<<(*edges)[i].second<<"\n";
				}
				if(lastIndex==-1){
					//cout<<"\n"<<i1<<" "<<i2;
					Rf_error("RDSToggle: edge not found in edgelist");
				}
			}

			toggle[0]= edges->at(lastIndex);

		}else{
			lastIndex=-1;
		}

	}

	inline void generateTetrad(){
		std::pair<int,int>  e1, e2;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		int i = 0;
		double n = edges->size();
		while(i<100000000){
			i++;
			e1Index = floor(Rf_runif(0,n));
			e2Index = floor(Rf_runif(0,n-1.0));
			if(e2Index>e1Index)
				e2Index++;
			e1 = edges->at(e1Index);
			e2 = edges->at(e2Index);
			if(e1.first!=e2.first && e1.first != e2.second &&
					e1.second!=e2.first && e1.second!=e2.second &&
					!net->hasEdge(e1.first,e2.second) &&
					!net->hasEdge(e2.first,e1.second))
				break;
		}
		if(i>=100000)
			Rf_error("TetradToggle: could not find tetrad");

		toggle[0].first = e1.first;
		toggle[0].second = e2.second;
		toggle[1].first = e2.first;
		toggle[1].second = e1.second;
		toggle[2] = e1;
		toggle[3] = e2;
	}


	inline double logRatio(){return 0.0;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	//changed continuous variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,double> > >& contVarChanges(){
		return std::vector<std::pair<int,std::pair<int,double> > >();
	}

	//changed discrete variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,int> > >& disVarChanges(){
		return std::vector<std::pair<int,std::pair<int,int> > >();
	}

};
*/


/*!
 * Toggleing for RDS data.
 */
template<class  Engine >
class Rds {
protected:
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	boost::shared_ptr<std::vector< std::pair<int,int> > > edges;
	std::vector<bool> observed;
	std::vector<int> unobservedNodes;

	int nUnobsEdges;

	bool wasTetrad;
	int e1Index, e2Index;
	int lastIndex;

	double logProbRatio;
public:

	Rds() : nUnobsEdges(0), wasTetrad(false), e1Index(0), e2Index(0), lastIndex(0),
		logProbRatio(0.0){}

	Rds(Rcpp::List l) : nUnobsEdges(0), wasTetrad(false), e1Index(0), e2Index(0), lastIndex(0),
			logProbRatio(0.0){}

	Rds( BinaryNet<Engine> & network) : nUnobsEdges(0), wasTetrad(false), e1Index(0), e2Index(0), lastIndex(0),
			logProbRatio(0.0){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
	}

	virtual ~Rds(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		if(net->isDirected())
			Rf_error("RDS is only applicable for undirected networks");
		edges = net->edgelist();

		observed = std::vector<bool>(net->size(),false);
		for(int i=0;i<edges->size();i++){
			if(!net->isMissing(edges->at(i).first,edges->at(i).second)){
				observed.at(edges->at(i).first) = true;
				observed.at(edges->at(i).second) = true;
				edges->at(i) = edges->at(edges->size()-1);
				edges->pop_back();
				i--;
			}
		}
		for(int i=0;i<net->size();i++){
			for(int j=0;j<net->discreteVarNames().size();j++){
				if(net->discreteVariableObserved(j,i))
					observed.at(i) = true;
			}
		}
		nUnobsEdges = 0;
		for(int i=0;i<edges->size();i++){
			if(!observed.at(edges->at(i).first) &&
					!observed.at(edges->at(i).second)){
				nUnobsEdges++;
			}

		}
		unobservedNodes.clear();
		for(int i=0;i<observed.size();i++){
			if(!observed[i])
				unobservedNodes.push_back(i);
		}
		//cout<<unobservedNodes.size()<<"\n";
	}
	void togglesAccepted(bool apply){
		if(!apply)
			return;
		if(wasTetrad){
			if(!observed[edges->at(e1Index).first] && !observed[edges->at(e1Index).second])
				nUnobsEdges--;
			if(!observed[edges->at(e2Index).first] && !observed[edges->at(e2Index).second])
				nUnobsEdges--;
			if(!observed[toggle.at(0).first] && !observed[toggle.at(0).second])
				nUnobsEdges++;
			if(!observed[toggle.at(1).first] && !observed[toggle.at(1).second])
				nUnobsEdges++;
			edges->at(e1Index) = toggle.at(0);
			edges->at(e2Index) = toggle.at(1);
		}else{
			if(lastIndex<0){
				// edge has been added
				edges->push_back(toggle[0]);
				nUnobsEdges++;
			}else{
				//edge has been removed
				int size = edges->size();
				(*edges)[lastIndex] = (*edges)[size-1];
				edges->pop_back();
				nUnobsEdges--;
			}
		}
	}

	inline void generate(){
		if(Rf_runif(0.0,1.0)<.5){
			generateTieDyad();
			wasTetrad=false;
		}else{
			generateTetrad();
			wasTetrad=true;
		}
	}

	int pickEdge(){
		int i=0;
		int index = -1;
		while(i<100000){
			i++;
			index = floor(Rf_runif(0.0,(double)edges->size()));
			if(!observed[edges->at(index).first] && !observed[edges->at(index).second]){
				break;
			}
		}
		if(i>=100000)
			Rf_error("RDSToggle: could not find edge between unobserved nodes");
		return index;
	}


	inline void generateTieDyad(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		bool pickEdge = Rf_runif(0.0,1.0)>0.5;
		double nEdges = nUnobsEdges;



		double ndyads = unobservedNodes.size() * (unobservedNodes.size()-1.0) / 2.0;
		if(nEdges==0 && pickEdge){
			pickEdge=false;
		}

		if(pickEdge){

			int index = this->pickEdge();
			toggle[0] = edges->at(index);
			if(nEdges==1){
				this->logProbRatio = log(1.0/(ndyads + .5));
			}else{
				this->logProbRatio = log(nEdges / (ndyads + nEdges));
				// << nEdges<<" "<<ndyads<<" "<< this->logProbRatio<<" "<<log(nEdges / (ndyads + nEdges));
			}
			this->lastIndex=index;
		}else{
			int i1 = floor(Rf_runif(0.0,(double)unobservedNodes.size()));
			int i2 = floor(Rf_runif(0.0,(double)unobservedNodes.size()-1.0));
			if(i2>=i1)
				i2++;
			toggle[0].first = unobservedNodes[i1];
			toggle[0].second = unobservedNodes[i2];
			bool hasEdge = net->hasEdge(toggle[0].first,toggle[0].second);
			if(hasEdge){
				//just re-picking randomly is faster than find
				int index = this->pickEdge();
				toggle[0] = edges->at(index);
				this->lastIndex=index;
			}else{
				this->lastIndex=-1;
			}

			if(hasEdge){
				this->logProbRatio = log((nEdges==1 ? 2.0/(ndyads + 1.0) :
									nEdges / (ndyads + nEdges)));

			}else{
				this->logProbRatio = log((nEdges==0 ? .5*(ndyads + 1.0) :
				        1.0 + (ndyads)/(nEdges + 1.0)));
			}
		}

	}

	inline void generateTetrad(){
		std::pair<int,int>  e1, e2;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		int i = 0;
		double n = edges->size();
		while(i<100000000){
			i++;
			e1Index = floor(Rf_runif(0,n));
			e2Index = floor(Rf_runif(0,n-1.0));
			if(e2Index>=e1Index)
				e2Index++;
			e1 = edges->at(e1Index);
			e2 = edges->at(e2Index);
			if(e1.first!=e2.first && e1.first != e2.second &&
					e1.second!=e2.first && e1.second!=e2.second &&
					!net->hasEdge(e1.first,e2.second) &&
					!net->hasEdge(e2.first,e1.second))
				break;
		}
		if(i>=100000)
			Rf_error("TetradToggle: could not find tetrad");

		toggle[0].first = e1.first;
		toggle[0].second = e2.second;
		toggle[1].first = e2.first;
		toggle[1].second = e1.second;
		toggle[2] = e1;
		toggle[3] = e2;
		this->logProbRatio = 0.0;
	}


	inline double logRatio(){return logProbRatio;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "Rds";
	}

};

typedef DyadToggle<Directed, Rds<Directed> > DirectedRdsToggle;
typedef DyadToggle<Undirected, Rds<Undirected> > UndirectedRdsToggle;

/*
template<class  Engine >
class RDSTmpToggle : public AbstractDyadToggle<Engine> {
protected:
	typedef boost::shared_ptr<std::vector< std::pair<int,int> > > EdgeList;
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	EdgeList obsEdges;
	EdgeList unobsEdges;
	std::vector<bool> observed;
	std::vector<int> unobservedNodes;

	double probTetrad;

	double logProbRatio;

	bool wasTetrad;
	int e1Index, e2Index;
	bool e1obs, e2obs;
	int lastIndex;
public:
	DYAD_TOGGLE_BOILERPLATE(RDSTmpToggle)

	RDSTmpToggle(Rcpp::List l){}

	RDSTmpToggle( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		obsEdges = EdgeList(new std::vector< std::pair<int,int> >);
		unobsEdges = EdgeList(new std::vector< std::pair<int,int> >);
	}

	virtual ~RDSTmpToggle(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		obsEdges = EdgeList(new std::vector< std::pair<int,int> >);
		unobsEdges = EdgeList(new std::vector< std::pair<int,int> >);
		if(net->isDirected())
			Rf_error("RDS is only applicable for undirected networks");
		EdgeList edges = net->edgelist();

		observed = std::vector<bool>(net->size(),false);
		for(int i=0;i<edges->size();i++){
			if(!net->isMissing(edges->at(i).first,edges->at(i).second)){
				observed.at(edges->at(i).first) = true;
				observed.at(edges->at(i).second) = true;
				edges->at(i) = edges->at(edges->size()-1);
				edges->pop_back();
				i--;
			}
		}
		for(int i=0;i<net->size();i++){
			for(int j=0;j<net->discreteVarNames().size();j++){
				if(net->discreteVariableObserved(j,i))
					observed.at(i) = true;
			}
		}

		for(int i=0;i<edges->size();i++){
			if(!observed[edges->at(i).first] && !observed[edges->at(i).first]){
				unobsEdges->push_back(edges->at(i));
			}else{
				obsEdges->push_back(edges->at(i));
			}
		}

		unobservedNodes.clear();
		for(int i=0;i<observed.size();i++){
			if(!observed[i])
				unobservedNodes.push_back(i);
		}
		probTetrad = 1.0 - (unobservedNodes.size() / (double)net->size() );
		//cout<<unobservedNodes.size()<<"\n";
	}
	void togglesAccepted(bool apply){
		if(!apply)
			return;
		if(wasTetrad){
			bool o11 = observed[toggle.at(0).first];
			bool o12 = observed[toggle.at(0).second];
			bool o21 = observed[toggle.at(1).first];
			bool o22 = observed[toggle.at(1).second];
			if(o11 || o12)
				obsEdges->push_back(toggle[0]);
			else
				unobsEdges->push_back(toggle[0]);

			if(o21 || o22)
				obsEdges->push_back(toggle[1]);
			else
				unobsEdges->push_back(toggle[1]);

			if(e1obs){//o11 || o22){
				obsEdges->at(e1Index) = obsEdges->at(obsEdges->size()-1);
				obsEdges->pop_back();
			}else{
				unobsEdges->at(e1Index) = unobsEdges->at(unobsEdges->size()-1);
				unobsEdges->pop_back();
			}

			if(e2obs){//o21 || o12){
				if(e2Index == obsEdges->size())
					e2Index = e1Index;
				obsEdges->at(e2Index) = obsEdges->at(obsEdges->size()-1);
				obsEdges->pop_back();
			}else{
				if(e2Index == unobsEdges->size())
					e2Index = e1Index;
				unobsEdges->at(e2Index) = unobsEdges->at(unobsEdges->size()-1);
				unobsEdges->pop_back();
			}
		}else{
			if(lastIndex<0){
				// edge has been added
				unobsEdges->push_back(toggle[0]);
			}else{
				//edge has been removed
				int size = unobsEdges->size();
				unobsEdges->at(lastIndex) = unobsEdges->at(size-1);
				unobsEdges->pop_back();
			}
		}
	}

	inline void generate(){
		double rand = Rf_runif(0.0,1.0);
		if(rand< (1.0 - probTetrad)){
			generateTieDyad();
			return;
		}
		generateTetrad();
	}



	inline void generateTieDyad(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		bool pickEdge = Rf_runif(0.0,1.0)>0.5;
		double nEdges = unobsEdges->size();
		double nDyads = unobservedNodes.size() * (unobservedNodes.size()-1.0) / 2.0;
		if(nEdges==0 && pickEdge){
			pickEdge=false;
		}

		if(pickEdge){
			int index = floor(Rf_runif(0.0,nEdges));
			toggle[0].first = unobsEdges->at(index).first;
			toggle[0].second = unobsEdges->at(index).second;
			if(nEdges==1){
				this->logProbRatio = log(1.0/(nDyads + .5));
			}else{
				this->logProbRatio = log(nEdges / (nDyads + nEdges));
				// << nEdges<<" "<<ndyads<<" "<< this->logProbRatio<<" "<<log(nEdges / (ndyads + nEdges));
			}
			this->lastIndex=index;
		}else{
			//toggle[0] = net->randomNonEdge();
			int i1 = floor(Rf_runif(0.0,unobservedNodes.size()));
			int i2 = floor(Rf_runif(0.0,unobservedNodes.size()-1));
			if(i2>=i1)
				i2++;
			toggle[0].first = unobservedNodes.at(i1);
			toggle[0].second = unobservedNodes.at(i2);

			bool hasEdge = net->hasEdge(toggle[0].first,toggle[0].second);
			if(hasEdge){
				//just re-picking randomly is faster than find
				int index = floor(Rf_runif(0.0,nEdges));
				toggle[0].first = unobsEdges->at(index).first;
				toggle[0].second = unobsEdges->at(index).second;
				this->lastIndex=index;
			}else{
				this->lastIndex=-1;
			}

			if(hasEdge){
				this->logProbRatio = log((nEdges==1 ? 2.0/(nDyads + 1.0) :
									nEdges / (nDyads + nEdges)));

			}else{
				this->logProbRatio = log((nEdges==0 ? .5*(nDyads + 1.0) :
				        1.0 + (nDyads)/(nEdges + 1.0)));
			}
		}

		wasTetrad=false;
	}

	inline void generateTetrad(){
		std::pair<int,int>  e1, e2;
		toggle = std::vector< std::pair<int,int> >(4,std::pair<int,int>(-1,-1));
		int i = 0;
		int nobs = obsEdges->size();
		int nunobs = unobsEdges->size();
		//cout << nobs << " "<<nunobs<<"\n";
		while(i<1000000){
			i++;
			int i1 = floor(Rf_runif(0,nobs+nunobs));
			int i2 = floor(Rf_runif(0,nobs+nunobs-1.0));
			if(i2>=i1)
				i2++;
			e1obs = i1<nobs;
			if(!e1obs){
				e1Index = i1 - nobs;
				e1 = unobsEdges->at(e1Index);
			}else{
				e1Index = i1;
				e1 = obsEdges->at(e1Index);
			}
			e2obs = i2<nobs;
			if(!e2obs){
				e2Index = i2 - nobs;
				e2 = unobsEdges->at(e2Index);
			}else{
				e2Index = i2;
				e2 = obsEdges->at(e2Index);
			}
			if(e1.first!=e2.first && e1.first != e2.second &&
					e1.second!=e2.first && e1.second!=e2.second &&
					!net->hasEdge(e1.first,e2.second) &&
					!net->hasEdge(e2.first,e1.second))
				break;
		}
		if(i>=1000000){
			cout << nobs << " "<<nunobs<<"\n";
			cout<<e1.first<<" "<<e1.second<<"   "<<e2.first<<" "<<e2.second<<"\n";
			Rf_error("RDSToggle: could not find tetrad");
		}

		toggle[0].first = e1.first;
		toggle[0].second = e2.second;
		toggle[1].first = e2.first;
		toggle[1].second = e1.second;
		toggle[2] = e1;
		toggle[3] = e2;

		wasTetrad=true;
		logProbRatio=0.0;
	}

	inline double logRatio(){return logProbRatio;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	//changed continuous variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,double> > >& contVarChanges(){
		return std::vector<std::pair<int,std::pair<int,double> > >();
	}

	//changed discrete variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,int> > >& disVarChanges(){
		return std::vector<std::pair<int,std::pair<int,int> > >();
	}

};
*/

/*!
 * alternate between two toggles.
 *
 * note: toggles must be memoryless.
 */
template<class  ToggleType1,class ToggleType2, class Engine >
class CompoundToggle {
protected:
	ToggleType1 toggle1;
	ToggleType2 toggle2;
	bool useFirst;
public:

	CompoundToggle(){
		toggle1 = ToggleType1();
		toggle2 = ToggleType2();
		useFirst = true;
	}

	CompoundToggle(Rcpp::List l){
		toggle1 = ToggleType1(l);
		toggle2 = ToggleType2(l);
		useFirst = true;
	}



	CompoundToggle( BinaryNet<Engine> & network){
		toggle1 = ToggleType1(network);
		toggle2 = ToggleType2(network);
		useFirst = true;
	}

	virtual ~CompoundToggle(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		toggle1.setNetwork(n);
		toggle2.setNetwork(n);
	}

	void initialize(){
		toggle1.initialize();
		toggle2.initialize();
	}
	void togglesAccepted(bool apply){}

	inline void generate(){
		useFirst = !useFirst;
		if(useFirst)
			toggle1.generate();
		else
			toggle2.generate();

	}


	inline double logRatio(){
		if(useFirst)
			return toggle1.logRatio();
		else
			return toggle2.logRatio();
	}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		if(useFirst)
			return toggle1.dyadToggles();
		else
			return toggle2.dyadToggles();
	}

	inline std::string name(){
		return std::string("Compound_") + toggle1.name() + std::string("_") + toggle2.name();
	}

};

template <class Engine>
class CompoundNodeTieDyadNieghborhood :
		public CompoundToggle<NodeTieDyad<Engine>, Neighborhood<Engine>, Engine >{
protected:
public:
	CompoundNodeTieDyadNieghborhood() : CompoundToggle<NodeTieDyad<Engine>, Neighborhood<Engine>, Engine >(){}
	CompoundNodeTieDyadNieghborhood(Rcpp::List l) :CompoundToggle<NodeTieDyad<Engine>, Neighborhood<Engine>, Engine >(l){}
	CompoundNodeTieDyadNieghborhood( BinaryNet<Engine> & network) : CompoundToggle<NodeTieDyad<Engine>, Neighborhood<Engine>, Engine >(network){}
	virtual ~CompoundNodeTieDyadNieghborhood(){}
};

typedef DyadToggle<Directed, CompoundNodeTieDyadNieghborhood<Directed> > DirectedNtdNbrToggle;
typedef DyadToggle<Undirected, CompoundNodeTieDyadNieghborhood<Undirected> > UndirectedNtdNbrToggle;

/*!
 * random toggle of unobserved dyads. (ergm default. poor mixing)
 */
template<class  Engine >
class RandomDyadMissing {
protected:
	typedef  boost::shared_ptr< std::vector< std::pair<int,int> > > DyadVecPtr;
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	DyadVecPtr unobservedDyads;
public:

	RandomDyadMissing(){}

	RandomDyadMissing(Rcpp::List l){}

	RandomDyadMissing( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
	}

	virtual ~RandomDyadMissing(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
		unobservedDyads = net->missingDyads();
	}
	void togglesAccepted(bool apply){}

	inline void generate(){
		if(unobservedDyads->size()==0)
			Rf_error("Can not toggle unobserved dyads in fully observed network");
		int index = floor(Rf_runif(0.0,(double)unobservedDyads->size()));
		toggle[0].first = (*unobservedDyads)[index].first;
		toggle[0].second = (*unobservedDyads)[index].second;
		//cout << toggle[0].first << " " <<toggle[0].second<< " - ";
	}


	inline double logRatio(){return 0.0;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "RandomMissingDyad";
	}
};

typedef DyadToggle<Directed, RandomDyadMissing<Directed> > DirectedRandomDyadMissingToggle;
typedef DyadToggle<Undirected, RandomDyadMissing<Undirected> > UndirectedRandomDyadMissingToggle;


/*!
 * Nodal Tie-dyad toggling among missing dyads
 */
template<class  Engine >
class NodeTieDyadMissing {
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	boost::shared_ptr< BinaryNet<Engine> > net;
	std::vector< std::pair<int,int> > toggle;
	std::vector<int> nmissing;
	std::vector<int> verts;
	double logProbRatio;
public:

	NodeTieDyadMissing(Rcpp::List l) : logProbRatio(0.0){}

	NodeTieDyadMissing() : logProbRatio(0.0){}

	NodeTieDyadMissing( BinaryNet<Engine> & network) : logProbRatio(0.0){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;

	}

	virtual ~NodeTieDyadMissing(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net = n;
	}

	void initialize(){
		nmissing = std::vector<int>();
		verts = std::vector<int>();
		int ind = -1;
		for(int i=0;i<net->size();i++){
			int vertAdded = false;
			for(int j=0;j<net->size();j++){
				if(i!=j && net->isMissing(i,j)){
					if(!vertAdded){
						vertAdded=true;
						nmissing.push_back(0);
						verts.push_back(i);
						ind++;
					}
					nmissing.at(ind)++;
				}
			}
		}
		toggle.resize(1,std::make_pair(-1,-1));
	}
	void togglesAccepted(bool apply){}

	inline void generate(){
		if(verts.size()==0)
			::Rf_error("NTDNonObservedToggle: No missing dyads");
		int index = floor(Rf_runif(0.0,(double)nmissing.size()));
		int node = verts[index];

		int neighbor;
		std::vector<int> missingEdges;
		//const Set* nbs;
		//if(!net->isDirected())
		//	nbs = &net->neighbors(node);
		//else
		//	nbs = &net->outneighbors(node);
		//Set::iterator it = nbs->begin();
		NeighborIterator it;
		NeighborIterator end;
		int degree;
		if(!net->isDirected()){
			it = net->begin(node);
			end = net->end(node);
			degree = net->degree(node);
			//nbs = &net->neighbors(node);
		}else{
			it = net->outBegin(node);
			end = net->outEnd(node);
			degree = net->outdegree(node);
			//nbs = &net->outneighbors(node);
		}
		for(;it!=end;it++){
			if(net->isMissing(node,*it))
				missingEdges.push_back(*it);
		}
		double nedges = missingEdges.size();
		double ndyads = nmissing[index];
		bool pickEdge = Rf_runif(0.0,1.0)>0.5;
		if(nedges==0)
			pickEdge=false;
		bool hasEdge;
		if(pickEdge){
			hasEdge=true;
			int nbrIndex = floor(Rf_runif(0.0,(double)nedges));
			neighbor=missingEdges[nbrIndex];
			toggle[0].first = node;
			toggle[0].second = neighbor;
			assert(node != neighbor);

		}else{
			neighbor = net->randomDyad(node,true);
			toggle[0].first = node;
			toggle[0].second = neighbor;
			hasEdge = net->hasEdge(node,neighbor);
			assert(node != neighbor);
			//cout << neighbor;
		}

		double tForward, tReverse;
		if(hasEdge){
			if(nedges<1.5)
				tReverse = 1.0/ndyads;
			else
				tReverse = 0.5/ndyads;
			tForward = 0.5/nedges + 0.5/ndyads;

		}else{
			if(nedges<.5)
				tForward = 1.0/ndyads;
			else
				tForward = 0.5/ndyads;
			tReverse = 0.5/(nedges + 1.0) + 0.5/ndyads;

		}
		//cout<<ndyads;
		this->logProbRatio = log( tReverse/tForward );
		//cout<<this->logProbRatio << " ";

		//calculate log ratio for neighbor
		//do we need this?
/*		if(!net->isDirected()){
			it = net->begin(neighbor);
			end = net->end(neighbor);
			//nbs = &net->neighbors(neighbor);
			//std::vector<int> missingEdgesNbr;
			nedges = 0.0;
			//it = nbs->begin();
			for(;it!=end;it++){
				if(net->isMissing(neighbor,*it))
					nedges++;
			}
			std::vector<int>::iterator iter = lower_bound(verts.begin(),verts.end(),neighbor);
			index = iter - verts.begin();
			//nedges = missingEdgesNbr.size();
			ndyads = nmissing.at(index);
			//cout << verts[index] << neighbor <<" ";
			double tForwardNbr, tReverseNbr;
			if(hasEdge){
				if(nedges<1.5)
					tReverseNbr = 1.0/ndyads;
				else
					tReverseNbr = 0.5/ndyads;
				tForwardNbr = 0.5/nedges + 0.5/ndyads;
				//cout <<tReverse <<" "<< tReverseNbr << " | ";
			}else{
				if(nedges<.5)
					tForwardNbr = 1.0/ndyads;
				else
					tForwardNbr = 0.5/ndyads;
				tReverseNbr = 0.5/(nedges + 1.0) + 0.5/ndyads;

			}

			this->logProbRatio = log( (tReverse + tReverseNbr) / (tForward + tForwardNbr) );
			//cout<<this->logProbRatio << "\n";
		}
*/
	}


	inline double logRatio(){
		return logProbRatio;
	}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "NodeTieDyadMissing";
	}
};

typedef DyadToggle<Directed, NodeTieDyadMissing<Directed> > DirectedNodeTieDyadMissingToggle;
typedef DyadToggle<Undirected, NodeTieDyadMissing<Undirected> > UndirectedNodeTieDyadMissingToggle;

/*
 * toggle dyads
 *
 *   o
 *	|  \
 *	|    b
 *	|      \
 *	o---a---o
 *	 \      /
 *	  \    /
 *	   \  /
 *	    o   <--- a randomly chosen node
 *
 * 1. choose a random node
 * 2. choose two neighbors at random. the dyad between these two is denoted as "a"
 * 3. choose a random neighbor of the first node chosen in step 2. The connection between this node and the second node is denoted as "b"
 * 4. Toggle "a" with probability .5 and "b" with probability .5
 *
 * tries this up to 10 times, rejecting if the dyad is not missing. If at the end of
 * 10 tries, no valid dyad is found, a random node is selected, and then a random missing dyad
 * connected to that node is selected.
 *
 */
template<class  Engine >
class NeighborhoodMissing{
protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	std::vector< std::pair<int,int> > toggle;
	boost::shared_ptr< BinaryNet<Engine> > net;
	std::vector<int> verts;
	bool twoSteps;

/*	int ithElement(const Set& s,int& i){
		Set::iterator it = s.begin();
		std::advance(it,i);
		return *it;
	}
*/
public:
	NeighborhoodMissing(Rcpp::List l) : twoSteps(false){}

	NeighborhoodMissing() : twoSteps(false){}

	NeighborhoodMissing( BinaryNet<Engine> & network) : twoSteps(false){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		toggle = std::vector< std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
	}

	virtual ~NeighborhoodMissing(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		verts = std::vector<int>();
		for(int i=0;i<net->size();i++){
			int vertAdded = false;
			for(int j=0;j<net->size();j++){
				if(i!=j && net->isMissing(i,j)){
					if(!vertAdded){
						vertAdded=true;
						verts.push_back(i);
						continue;
					}
				}
			}
		}
		toggle.resize(1,std::make_pair(-1,-1));
	}
	void togglesAccepted(bool apply){

	}

	void generate(){
		for(int i=0;i<10;i++){
			if(generateToggle())
				return;
		}
		int node = verts.at(floor(Rf_runif(0.0,verts.size())));
		int neighbor = net->randomDyad(node,true);
		toggle[0].first = node;
		toggle[0].second = neighbor;
	}

	bool generateToggle(){
		int from, to;
		int node = floor(Rf_runif(0.0,net->size()));
		//const Set* nbs;
		//if(net->isDirected())
		//	nbs = &net->outneighbors(node);
		//else
		//	nbs = &net->neighbors(node);
		NeighborIterator it;
		NeighborIterator end;
		int degree;
		if(!net->isDirected()){
			it = net->begin(node);
			end = net->end(node);
			degree = net->degree(node);
			//nbs = &net->neighbors(node);
		}else{
			it = net->outBegin(node);
			end = net->outEnd(node);
			degree = net->outdegree(node);
			//nbs = &net->outneighbors(node);
		}
		if(degree<2){
			from = floor(Rf_runif(0.0,net->size()-1));
			to = floor(Rf_runif(0.0,net->size()-2));
			if(from>=node)
				node++;
			if(to>=std::min(node,from))
				to++;
			if(to>=std::max(node,from))
				to++;
		}else{
			int fromIndex = floor(Rf_runif(0.0,degree));
			std::advance(it,fromIndex);
			from = *it;
			std::advance(it,-fromIndex);
			//from = ithElement(*nbs,fromIndex);
			int toIndex = floor(Rf_runif(0.0,degree-1));
			if(toIndex>=fromIndex)
				toIndex++;
			std::advance(it,toIndex);
			to = *it;
			std::advance(it,-toIndex);
			//to = ithElement(*nbs,toIndex);
		}
		int first = from;
		int second = to;
		if(twoSteps){
			std::vector<int> v(3,0);
			v[0] = node;
			v[1] = to;
			v[2] = from;
			std::sort(v.begin(),v.end());
			//const Set* nbs1;
			//if(net->isDirected())
			//	nbs1 = &net->outneighbors(from);
			//else
			//	nbs1 = &net->neighbors(from);
			if(!net->isDirected()){
				it = net->begin(from);
				end = net->end(from);
				degree = net->degree(from);
			}else{
				it = net->outBegin(from);
				end = net->outEnd(from);
				degree = net->outdegree(from);
				//nbs = &net->outneighbors(node);
			}
			int reqSize = 1;
			if(net->hasEdge(from,node))
				reqSize++;
			if(net->hasEdge(from,to))
				reqSize++;
			if(degree<reqSize){
				first = floor(Rf_runif(0.0,net->size()-3));
				if(first>=v[0])
					first++;
				if(first>=v[1])
					first++;
				if(first>=v[2])
					first++;
			}else{
				int ind = floor(Rf_runif(0.0,degree+1-reqSize));
				//Set::iterator it = nbs1->begin();
				std::advance(it,ind);
				if(net->hasEdge(from,std::min(to,node)) && *it>=std::min(to,node))
					it++;
				if(net->hasEdge(from,std::max(to,node)) && *it>=std::max(to,node))
					it++;
				first = *it;
			}
			/*
			v.clear();
			v.push_back(node);
			v.push_back(from);
			v.push_back(first);
			v.push_back(to);
			std::sort(v.begin(),v.end());

			second = floor(Rf_runif(0.0,net->size()-4));
			if(second>=v[0])
				second++;
			if(second>=v[1])
				first++;
			if(second>=v[2])
				second++;
			if(second>=v[3])
				second++;
			*/
		}
		twoSteps = !twoSteps;
		if(net->isMissing(first,second)){
			toggle[0].first = first;
			toggle[0].second = second;
			return true;
		}else
			return false;

	}


	inline double logRatio(){return 0.0;}

	//changed ties
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return toggle;
	}

	inline std::string name(){
		return  "NeighborhoodMissing";
	}
};

typedef DyadToggle<Directed, NeighborhoodMissing<Directed> > DirectedNeighborhoodMissingToggle;
typedef DyadToggle<Undirected, NeighborhoodMissing<Undirected> > UndirectedNeighborhoodMissingToggle;

template <class Engine>
class CompoundNodeTieDyadNieghborhoodMissing :
		public CompoundToggle<NodeTieDyadMissing<Engine>, NeighborhoodMissing<Engine>, Engine >{
protected:
public:
	CompoundNodeTieDyadNieghborhoodMissing() : CompoundToggle<NodeTieDyadMissing<Engine>, NeighborhoodMissing<Engine>, Engine >(){}
	CompoundNodeTieDyadNieghborhoodMissing(Rcpp::List l) :CompoundToggle<NodeTieDyadMissing<Engine>, NeighborhoodMissing<Engine>, Engine >(l){}
	CompoundNodeTieDyadNieghborhoodMissing( BinaryNet<Engine> & network) : CompoundToggle<NodeTieDyadMissing<Engine>, NeighborhoodMissing<Engine>, Engine >(network){}
	virtual ~CompoundNodeTieDyadNieghborhoodMissing(){}
};

typedef DyadToggle<Directed, CompoundNodeTieDyadNieghborhoodMissing<Directed> > DirectedNtdNbrMissingToggle;
typedef DyadToggle<Undirected, CompoundNodeTieDyadNieghborhoodMissing<Undirected> > UndirectedNtdNbrMissingToggle;




}


#endif /* DYADTOGGLES_H_ */
