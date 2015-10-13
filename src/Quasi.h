/*
 * Quasi.h
 *
 *  Created on: Jun 23, 2013
 *      Author: ianfellows
 */

#ifndef QUASI_H_
#define QUASI_H_
#include "BinaryNet.h"
#include "Model.h"
#include <vector>

namespace ernm {

/*!
 * Experimental partialy implement quasi-likelihood
 */
template<class Engine>
class Quasi {
protected:
	double logSamplingWeight;
	double logQuasiWeight;
	std::vector<std::pair<int,int> > dyads; 			//* std::pair<from,to>
	std::vector<std::pair<int,int> > catVerts;		//* std::pair<category, node id>
	std::vector< std::vector<double> > denomStats; 	//* the g(t') terms of sum_{t'}(exp(theta*(g(t')-g(t)) + o(t')-o(t)))
	std::vector<std::vector<double> > denomOff;		//* the o(t')
	std::vector<double> obsStats;				//*the g(t) terms
	std::vector<double> obsOff;					//*the o(t) terms
public:
	Quasi(){
		logSamplingWeight = 0.0;
		logQuasiWeight = 0.0;
	}
	virtual ~Quasi(){}

	double getLogSamplingWeight(){
		return this->logSamplingWeight;
	}

	double getLogQuasiWeight(){
		return this->logQuasiWeight;
	}

	void setLogQuasiWeight(double qw){
		this->logQuasiWeight = qw;
	}

	virtual void initialize(Model<Engine>& model) = 0;

	void calculate(Model<Engine>& model){
		obsStats = model.statistics();
		obsOff = model.offset();
		std::vector<int> ncat;
		for(int i=0;i<catVerts.size();i++){
			ncat.push_back(model.network()->
					discreteVariableAttributes(catVerts[i].first).labels().size());
		}
		int dn = dyads.size();
		int cn = catVerts.size();
		class Recur{
		public:
			void calc(int ind,int& dn,int& cn,
					std::vector<int>& ncat,
					std::vector< std::vector<double> >& denomStats,
					std::vector< std::vector<double> >& denomOff,
					std::vector<std::pair<int,int> >& dyads,
					std::vector<std::pair<int,int> >& catVerts,
					Model<Engine >& model){
				if(ind<=dn-1){
					if(ind==dn+cn-1){
						denomStats.push_back(model.statistics());
						denomOff.push_back(model.offset());
					}else
						calc(ind+1,dn,cn,ncat,denomStats,denomOff,dyads,catVerts,model);
					model.dyadUpdate(dyads[ind].first, dyads[ind].second);
					model.network()->toggle(dyads[ind].first, dyads[ind].second);
					if(ind==dn+cn-1){
						denomStats.push_back(model.statistics());
						denomOff.push_back(model.offset());
					}else
						calc(ind+1,dn,cn,ncat,denomStats,denomOff,dyads,catVerts,model);
					model.dyadUpdate(dyads[ind].first, dyads[ind].second);
					model.network()->toggle(dyads[ind].first, dyads[ind].second);
				}else if(ind<=dn+cn-1){
					int i = ind-dn;
					int n = ncat[catVerts[i].first];
					int tmp = model.network()->discreteVariableValue(
							catVerts[i].first,catVerts[i].second);
					int k=tmp;
					if(ind==dn+cn-1){
						denomStats.push_back(model.statistics());
						denomOff.push_back(model.offset());
					}else
						calc(ind+1,dn,cn,ncat,denomStats,denomOff,dyads,catVerts,model);

					for(int j=0;j<n;j++){
						k++;
						if(k>n)
							k=k-n;
						//Language("print",wrap(catVerts[i].first)).eval();
						//Language("print",wrap(k)).eval();
						model.discreteVertexUpdate(
							catVerts[i].second,
							catVerts[i].first ,
							k);
						model.network()->setDiscreteVariableValue(
								catVerts[i].first ,
								catVerts[i].second,
								k);
						if(k==tmp){

						}else if(ind==dn+cn-1){
							denomStats.push_back(model.statistics());
							denomOff.push_back(model.offset());
						}else
							calc(ind+1,dn,cn,ncat,denomStats,denomOff,dyads,catVerts,model);

					}
				}else
					Rf_error("Quasi: logic error");
			}
		};
		Recur a;
		a.calc(0,dn,cn,ncat,denomStats,denomOff,dyads,catVerts,model);

	}

	double logLik(std::vector<double> thetas){
		double ll = 0.0;
		for(int i=0;i<denomStats.size();i++){
			//cout<<".";
			double expon = 0.0;
			for(int j=0;j<obsStats.size();j++)
				expon += thetas[j] * (denomStats[i][j] - obsStats[j]);
			for(int j=0;j<obsOff.size();j++)
				expon += denomOff[i][j] - obsOff[j];
			ll += exp(expon);
		}
		return -log(ll);
	}

};

template<class Engine>
class QuasiRandomPatch : public Quasi< Engine > {
protected:
	int nSamp;
public:
	QuasiRandomPatch(int n){
		nSamp = n;
	}
	QuasiRandomPatch(){nSamp=1;}
	virtual ~QuasiRandomPatch(){}

	virtual void initialize(Model<Engine>& model){
		std::vector<int> randVars = model.randomVariables(true);
		bool randDyads = model.hasRandomGraph();
		int size = model.network()->size();
		for(int i=0;i<nSamp;i++){

			double tmp = Rf_runif(0.0,size*size + size*randVars.size());
			if(tmp<size*size){
				int from = floor(Rf_runif(0.0,size));
				int to = floor(Rf_runif(0.0,size-1));
				if(to>=from)
					to++;
				this->dyads.push_back(std::make_pair(from,to));
			}else{
				int varInd = floor(Rf_runif(0.0,randVars.size()));
				int vertInd = floor(Rf_runif(0.0,size));
				this->catVerts.push_back(std::make_pair(randVars[varInd],vertInd));
			}
		}
	}
};



} /* namespace ernm */
#endif /* QUASI_H_ */
