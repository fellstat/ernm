/*
 * VertexToggles.h
 *
 *  Created on: Jan 9, 2014
 *      Author: ianfellows
 */

#ifndef VERTEXTOGGLES_H_
#define VERTEXTOGGLES_H_

#include "VertexToggle.h"

namespace ernm{


/*!
 * Random toggling for a categorical vertex. Adaptive gaussian proposals for continuous vertex variates
 * with an acceptance rate of .234 as the goal. (see Roberts, G.O., Gelman, A., Gilks, W.R. (1997).
 * Weak Convergence and Optimal Scaling of Random Walk Metropolis Algorithms. Ann. Appl. Probab. 7, 110-20.)
 *
 */
template<class Engine>
class DefaultVertex{
protected:
	boost::shared_ptr< BinaryNet<Engine> > net;
	std::vector<int> contVars;
	std::vector<double> lowerLim;
	std::vector<double> upperLim;
	std::vector<int> disVars;
	std::vector<int> nlevels;
	std::vector<std::pair<int,std::pair<int,int> > > disToggle;
	std::vector<std::pair<int,std::pair<int,double> > > contToggle;

	std::vector<double> delta; //scale of continuous toggle change. adaptively chosen
	std::vector<int> nAccepted;
	std::vector<int> nRejected;
	int lastContIndex;
public:

	DefaultVertex() : lastContIndex(-1){
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
	}

	DefaultVertex(Rcpp::List l) : lastContIndex(-1){
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
	}

	DefaultVertex(BinaryNet<Engine>& network,
						std::vector<int>& cvars,std::vector<int>& dvars) : lastContIndex(-1){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		contVars = cvars;
		disVars = dvars;
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
	}

	DefaultVertex(BinaryNet<Engine> & network) : lastContIndex(-1){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
					std::make_pair(-1,std::make_pair(-1,-1)));
	}

	virtual ~DefaultVertex(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		lastContIndex = -1;
		nlevels = std::vector<int>();
		for(int i =0;i<disVars.size();i++){
			nlevels.push_back(net->discreteVariableAttributes(
					disVars.at(i)).labels().size());
		}
		bool initAdaptive = delta.size() != contVars.size();
		if(initAdaptive){
			delta = std::vector<double>(contVars.size(),1.0);
			nAccepted = std::vector<int>(contVars.size(),0);
			nRejected = std::vector<int>(contVars.size(),0);
		}

		lowerLim.clear();
		upperLim.clear();
		for(int i =0;i<contVars.size();i++){
			ContinAttrib attr = net->continVariableAttributes(contVars.at(i));
			double l = attr.hasLowerBound() ? attr.lowerBound() : -std::numeric_limits<double>::infinity();
			double u = attr.hasUpperBound() ? attr.upperBound() : std::numeric_limits<double>::infinity();
			//Rcpp::Rcout << "lower:" << l << " upper: " << u << "\n";
			lowerLim.push_back(l);
			upperLim.push_back(u);
			if(initAdaptive && attr.hasLowerBound() && attr.hasUpperBound())
				delta[i] = .1 * (u - l);
			else if(initAdaptive){
				double n = net->size();
				int varInd = contVars[i];
				double s = 0.0;
				double ssq = 0.0;
				for(int j=0;j<net->size();j++){
					double val = net->continVariableValue(varInd,j);
					s += val;
					ssq += pow(val,2.0);
				}
				delta[i] = sqrt(ssq/n - pow(s/n,2.0));
			}
			if(delta[i]<.00001)
				delta[i] = 1.0;
			//Rcpp::Rcout << "delta = " << delta[i]<< "\n";
		}
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
		contToggle = std::vector<std::pair<int,std::pair<int,double> > >(1,
						std::make_pair(-1,std::make_pair(-1,-1.0)));
	}

	void togglesAccepted(bool apply){
		if(lastContIndex>=0){
			if(apply){
				nAccepted[lastContIndex]++;
			}else{
				nRejected[lastContIndex]++;
			}
			int na = nAccepted[lastContIndex];
			int nr = nRejected[lastContIndex];
			if(na + nr > 100){
				double arate = (double)na / (double) (na + nr);
				double trate = .234;
				//if(contVars.size()==1)
				//	trate = .44;
				//Rcpp::Rcout <<"acc. rate: " << arate << " delta: " << delta[lastContIndex] <<"\n";
				if(arate > trate)
					delta[lastContIndex] *= 1.2;
				else
					delta[lastContIndex] *= .85;
				if(delta[lastContIndex] > upperLim[lastContIndex] - lowerLim[lastContIndex])
					delta[lastContIndex] = upperLim[lastContIndex] - lowerLim[lastContIndex];
				if(delta[lastContIndex] >= std::numeric_limits<double>::max()/100.0)
					delta[lastContIndex] = std::numeric_limits<double>::max()/100.0;
				if(delta[lastContIndex] < .00001)
					delta[lastContIndex] = .00001;
				nAccepted[lastContIndex] = nRejected[lastContIndex] = 0;
			}
		}
	}

	inline double logRatio(){return 0.0;}

	inline void generate(){
		lastContIndex = -1;
		if(disVars.size() + contVars.size() == 0 )
			Rf_error("DefaultVertexToggle: no vertex variables specified.");
		int vertex = floor(Rf_runif(0.0,net->size()));
		int index = floor(Rf_runif(0.0,contVars.size() + disVars.size()));
		//cout << "\ngenerate:"<<vertex<<" "<<index<<"\n";
		if(index >= contVars.size()){
			index = index - contVars.size();
			int varInd = disVars.at(index);
			int curVal = net->discreteVariableValue(varInd,vertex);
			int newVal = floor(Rf_runif(1.0,nlevels[index]));
			if(newVal>=curVal)
				newVal++;
			contToggle.clear();
			disToggle.clear();
			disToggle.push_back(std::make_pair(vertex,std::make_pair(varInd,newVal)));
		}else{
			lastContIndex = index;
			int varInd = contVars.at(index);
			double curVal = net->continVariableValue(varInd,vertex);
			double change = Rf_rnorm(0.0,delta[index]);
			double newVal = curVal + change;
			newVal = std::min(newVal, std::numeric_limits<double>::max());
			newVal = std::max(newVal, -std::numeric_limits<double>::max());
			while(newVal > upperLim[index])
				newVal = newVal - (upperLim[index] - lowerLim[index]);
			while(newVal < lowerLim[index])
				newVal = newVal + (upperLim[index] - lowerLim[index]);
			//Rcpp::Rcout<< "ul: " << upperLim[index] << " old value:" << curVal << " propose: " << newVal << "\n";
			disToggle.clear();
			contToggle.clear();
			contToggle.push_back( std::make_pair(vertex,std::make_pair(varInd,newVal)) );
		}
	}

	void setDiscreteVars(std::vector<int>& inds){
		disVars = inds;
	}

	void setContinuousVars(std::vector<int>& inds){
		contVars = inds;
	}

	//changed continuous variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,double> > >& contVarChanges(){
		return contToggle;
	}

	//changed discrete variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,int> > >& disVarChanges(){
		return disToggle;
	}

	inline std::string name(){
		return "DefaultVertex";
	}

};

typedef VertexToggle<Directed, DefaultVertex<Directed> > DirectedDefaultVertexToggle;
typedef VertexToggle<Undirected, DefaultVertex<Undirected> > UndirectedDefaultVertexToggle;


/*!
 * toggles randomly among missing vertex variables.
 */
template<class Engine>
class VertexMissing {
protected:
	boost::shared_ptr< BinaryNet<Engine> > net;
	std::vector<int> contVars;
	std::vector< std::pair<int,int> > contUnobserved; // std::vector< std::pair<variable id,vertex id>
	std::vector<int> disVars;
	std::vector< std::pair<int,int> > disUnobserved; // std::vector< std::pair<variable id,vertex id>
	std::vector<int> nlevels;
	std::vector<std::pair<int,std::pair<int,int> > > disToggle;
	std::vector<std::pair<int,std::pair<int,double> > > contToggle;
	std::vector<double> lowerLim;
	std::vector<double> upperLim;

	std::vector<double> delta; //scale of continuous toggle change. adaptively chosen
	std::vector<int> nAccepted;
	std::vector<int> nRejected;
	int lastContIndex;
public:

	VertexMissing() : lastContIndex(-1){
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
	}

	VertexMissing(Rcpp::List l) : lastContIndex(-1){
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
	}

	VertexMissing(BinaryNet<Engine>& network,
						std::vector<int> cvars,std::vector<int> dvars) : lastContIndex(-1){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		contVars = cvars;
		disVars = dvars;
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
	}

	VertexMissing(BinaryNet<Engine> & network) : lastContIndex(-1){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		net = n;
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
					std::make_pair(-1,std::make_pair(-1,-1)));
	}

	virtual ~VertexMissing(){}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net=n;
	}

	void initialize(){
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
		nlevels = std::vector<int>();
		for(int i =0;i<disVars.size();i++){
			nlevels.push_back(net->discreteVariableAttributes(
					disVars.at(i)).labels().size());
		}

		for(int i=0;i<net->size();i++){
			for(int j=0;j<disVars.size();j++){
				if(!net->discreteVariableObserved(disVars[j],i)){
					disUnobserved.push_back(std::make_pair(disVars[j],i));
				}
			}

			for(int j=0;j<contVars.size();j++){
				if(!net->continVariableObserved(contVars[j],i)){
					contUnobserved.push_back(std::make_pair(contVars[j],i));
				}
			}

		}

		bool initAdaptive = delta.size() != contVars.size();
		if(initAdaptive){
			delta = std::vector<double>(contVars.size(),1.0);
			nAccepted = std::vector<int>(contVars.size(),0);
			nRejected = std::vector<int>(contVars.size(),0);
		}

		lowerLim.clear();
		upperLim.clear();
		for(int i =0;i<contVars.size();i++){
			ContinAttrib attr = net->continVariableAttributes(contVars.at(i));
			double l = attr.hasLowerBound() ? attr.lowerBound() : -std::numeric_limits<double>::infinity();
			double u = attr.hasUpperBound() ? attr.upperBound() : std::numeric_limits<double>::infinity();
			//Rcpp::Rcout << "lower:" << l << " upper: " << u << "\n";
			lowerLim.push_back(l);
			upperLim.push_back(u);
			if(initAdaptive && attr.hasLowerBound() && attr.hasUpperBound())
				delta[i] = .1 * (u - l);
			else if(initAdaptive){
				double n = net->size();
				int varInd = contVars[i];
				double s = 0.0;
				double ssq = 0.0;
				for(int j=0;j<net->size();j++){
					double val = net->continVariableValue(varInd,j);
					s += val;
					ssq += pow(val,2.0);
				}
				delta[i] = sqrt(ssq/n - pow(s/n,2.0));
			}
			if(delta[i]<.00001)
				delta[i] = 1.0;
		}
		disToggle = std::vector<std::pair<int,std::pair<int,int> > >(1,
				std::make_pair(-1,std::make_pair(-1,-1)));
		contToggle = std::vector<std::pair<int,std::pair<int,double> > >(1,
						std::make_pair(-1,std::make_pair(-1,-1.0)));

	}

	void togglesAccepted(bool apply){
		if(lastContIndex>=0){
			if(apply){
				nAccepted[lastContIndex]++;
			}else{
				nRejected[lastContIndex]++;
			}
			int na = nAccepted[lastContIndex];
			int nr = nRejected[lastContIndex];
			if(na + nr > 100){
				double arate = (double)na / (double) (na + nr);
				double trate = .234;
				if(contVars.size()==1)
					trate = .44;
				//Rcpp::Rcout <<"acc. rate: " << arate << " delta: " << delta[lastContIndex] <<"\n";
				if(arate > trate)
					delta[lastContIndex] *= 1.2;
				else
					delta[lastContIndex] *= .85;
				if(delta[lastContIndex] > upperLim[lastContIndex] - lowerLim[lastContIndex])
					delta[lastContIndex] = upperLim[lastContIndex] - lowerLim[lastContIndex];
				if(delta[lastContIndex] >= std::numeric_limits<double>::max()/100.0)
					delta[lastContIndex] = std::numeric_limits<double>::max()/100.0;
				if(delta[lastContIndex] < .00001)
					delta[lastContIndex] = .00001;
				nAccepted[lastContIndex] = nRejected[lastContIndex] = 0;
			}
		}
	}

	inline double logRatio(){return 0.0;}

	inline void generate(){
		lastContIndex = -1;
		if(disVars.size() + contVars.size() == 0 )
			Rf_error("DefaultVertexToggle: no vertex variables specified.");
		bool toggleDiscrete = floor(Rf_runif(0.0,contVars.size() + disVars.size())) >= contVars.size();
		if(toggleDiscrete || contUnobserved.size()==0){
			if(disUnobserved.size()==0)
				::Rf_error("No unobserved variables");
			int tmp = floor(Rf_runif(0.0,disUnobserved.size()));
			int varInd = disUnobserved[tmp].first;
			int vertex = disUnobserved[tmp].second;
			int curVal = net->discreteVariableValue(varInd,vertex);
			int index=-1;
			//cout << varInd << " " << vertex << " " << disVars[0] << "\n";
			for(int i=0;i<disVars.size();i++){

				if(disVars[i]==varInd){
					index=i;
					break;
				}
			}
			assert(index > -1);
			int newVal = floor(Rf_runif(1.0,nlevels[index]));
			if(newVal>=curVal)
				newVal++;
			disToggle.clear();
			contToggle.clear();
			disToggle.push_back(std::make_pair(vertex,std::make_pair(varInd,newVal)));
		}else{
			int tmp = floor(Rf_runif(0.0,contUnobserved.size()));
			int varInd = contUnobserved[tmp].first;
			int index = lastContIndex = indexOf(varInd, contVars);
			int vertex = contUnobserved[tmp].second;
			double curVal = net->continVariableValue(varInd,vertex);
			double change = Rf_rnorm(0.0,delta[index]);
			double newVal = curVal + change;
			newVal = std::min(newVal, std::numeric_limits<double>::max());
			newVal = std::max(newVal, -std::numeric_limits<double>::max());
			while(newVal > upperLim[index])
				newVal = newVal - (upperLim[index] - lowerLim[index]);
			while(newVal < lowerLim[index])
				newVal = newVal + (upperLim[index] - lowerLim[index]);
			//Rcpp::Rcout<< "index: " << upperLim[index] << " old value:" << curVal << " propose: " << newVal << "\n";
			disToggle.clear();
			contToggle.clear();
			contToggle.push_back( std::make_pair(vertex,std::make_pair(varInd,newVal)) );
		}
	}

	void setDiscreteVars(std::vector<int>& inds){
		disVars = inds;
	}

	void setContinuousVars(std::vector<int>& inds){
		contVars = inds;
	}

	//changed continuous variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,double> > >& contVarChanges(){
		return contToggle;
	}

	//changed discrete variables <vertex<varId,new value>>
	inline std::vector<std::pair<int,std::pair<int,int> > >& disVarChanges(){
		return disToggle;
	}

	inline std::string name(){
		return "VertexMissing";
	}

};

typedef VertexToggle<Directed, VertexMissing<Directed> > DirectedVertexMissingToggle;
typedef VertexToggle<Undirected, VertexMissing<Undirected> > UndirectedVertexMissingToggle;


}


#endif /* VERTEXTOGGLES_H_ */
