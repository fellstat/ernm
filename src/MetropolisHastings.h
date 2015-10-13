/*
 * MetropolisHastings.h
 *
 *  Created on: Jun 28, 2011
 *      Author: ianfellows
 */

#ifndef METROPOLISHASTINGSH_
#define METROPOLISHASTINGSH_

#include "Model.h"
#include "VertexToggles.h"
#include "ToggleController.h"
#include "DyadToggles.h"
#include <cmath>
#include <Rcpp.h>
#include <assert.h>

namespace ernm {

/*!
 * MCMC sampling.
 */
template<class Engine>
class MetropolisHastings{

protected:
	typedef boost::shared_ptr< Model<Engine> > ModelPtr;
	typedef boost::shared_ptr< AbstractDyadToggle<Engine> > DyadTogglePtr;
	typedef boost::shared_ptr< AbstractVertexToggle<Engine> > VertexTogglePtr;
	ModelPtr model;
	DyadTogglePtr dyadToggle;
	VertexTogglePtr vertToggle;
	double probDyad;

public:

	//Constructors.
	//Note that the defaults toggles are *CompoundNodeTieDyadNieghborhood and DefaultVertexTobble
	MetropolisHastings(){
		dyadToggle = DyadTogglePtr(new DyadToggle<Engine, CompoundNodeTieDyadNieghborhood<Engine> >());
		vertToggle = VertexTogglePtr(new VertexToggle<Engine, DefaultVertex<Engine> >());
		probDyad = .8;
	}

	MetropolisHastings(Model<Engine>& mod,DyadTogglePtr tp,
			VertexTogglePtr vtp){
		model = mod.clone();
		dyadToggle = tp;
		vertToggle = vtp;
		probDyad=.8;
	}

	MetropolisHastings(Model<Engine>& mod,AbstractDyadToggle<Engine>& tp,
			AbstractVertexToggle<Engine>& vtp){
		model = mod.clone();
		dyadToggle = DyadTogglePtr(tp.vCloneUnsafe());
		vertToggle = VertexTogglePtr(vtp.vCloneUnsafe());
		probDyad=.8;
	}

	MetropolisHastings(Model<Engine> mod){
		model = mod.clone();
		dyadToggle = DyadTogglePtr(new DyadToggle<Engine, CompoundNodeTieDyadNieghborhood<Engine> >(*mod.network()));
		vertToggle = VertexTogglePtr(new VertexToggle<Engine, DefaultVertex<Engine> >(*mod.network()));
		probDyad=.8;
	}

	MetropolisHastings(Model<Engine> mod, double pdyad){
		model = mod.clone();
		dyadToggle = DyadTogglePtr(new DyadToggle<Engine, CompoundNodeTieDyadNieghborhood<Engine> >(*mod.network()));
		vertToggle = VertexTogglePtr(new VertexToggle<Engine, DefaultVertex<Engine> >(*mod.network()));
		probDyad=pdyad;
	}



	~MetropolisHastings(){

	}

	/*!
	 * Initialize sampler
	 */
	void initialize(){

		//make sure model and toggles agree
		std::vector<int> tmp = model->randomVariables(true);
		vertToggle->vSetDiscreteVars(tmp);
		tmp = model->randomVariables(false);
		vertToggle->vSetContinuousVars(tmp);
		dyadToggle->vSetNetwork(model->network());
		vertToggle->vSetNetwork(model->network());

		//initialize toggles
		dyadToggle->vInitialize();
		vertToggle->vInitialize();
	}

	/*!
	 * Set exponential family model
	 */
	void setModel(const Model<Engine>& mod){
		model = mod.clone();
		dyadToggle->vSetNetwork(mod.network());
		vertToggle->vSetNetwork(mod.network());
		std::vector<int> tmp = mod.randomVariables(true);
		vertToggle->vSetDiscreteVars(tmp);
		tmp = mod.randomVariables(false);
		vertToggle->vSetContinuousVars(tmp);

	}

	ModelPtr getModel(){
		return model;
	}

	/*!
	 * Get model exposed to R
	 */
	SEXP getModelR(){
		return wrap(*model);
	}

	void setDyadToggle(AbstractDyadToggle<Engine>& tp){
		dyadToggle = tp;
	}

	void setDyadToggleType(std::string name, Rcpp::List params){
		dyadToggle = DyadTogglePtr(ToggleController<Engine>::getDyadToggle(name,params));
	}

	void setVertexToggle(AbstractVertexToggle<Engine>& vtp){
		vertToggle = vtp;
	}

	void setVertexToggleType(std::string name, Rcpp::List params){
		vertToggle = VertexTogglePtr(ToggleController<Engine>::getVertexToggle(name,params));
	}

	/*!
	 * Set the probability of a dyad toggle (as opposed to a vertex toggle)
	 */
	void setDyadProbability(double d){
		assert(d>=0.0 && d<=1.0);
		probDyad=d;
	}

	/*!
	 * Run MCMC for steps steps
	 */
	double run(int steps){
		//Rcpp::RNGScope scope;
		int accepted = 0;
		std::vector<std::pair<int,std::pair<int,double> > > contToggles;
		std::vector<std::pair<int,std::pair<int,int> > > disToggles;
		std::vector<std::pair<int,int> > tieToggles;

		double pd = probDyad;
		if(!model->hasAnyRandomVariables())
			pd = 1.0;
		if(!model->hasRandomGraph())
			pd = 0.0;

		bool isDyadToggle = false;
		for(int i=0;i<steps;i++){
			isDyadToggle = Rf_runif(0.0,1.0) < pd;
			if(isDyadToggle){
				dyadToggle->vGenerate();

				tieToggles = dyadToggle->vDyadToggles();
				contToggles.clear();
				disToggles.clear();
				//cout <<"toggles: " << tieToggles[0].first<<" "<<tieToggles[0].second<<"\n";
				//cout << "is tie: " << model->network()->hasEdge(tieToggles[0].first,tieToggles[0].second)<<"\n";
				double lastLik = model->logLik();
				for(int j=0;j<tieToggles.size();j++){
					model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
					model->network()->toggle(tieToggles[j].first, tieToggles[j].second);
				}

				double lr = model->logLik() - lastLik;
				//cout << "lr:"<<lr<<"\n";
				lr += dyadToggle->vLogRatio();
				//cout << "v toggle lr:"<<dyadToggle.logRatio()<<"\n";
				if(lr>log(Rf_runif(0.0,1.0))){
					accepted++;
					dyadToggle->vTogglesAccepted(true);
				}else{
					for(int j=0;j<tieToggles.size();j++){
						model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
						model->network()->toggle(tieToggles[j].first,tieToggles[j].second);
					}
					dyadToggle->vTogglesAccepted(false);
				}
			}else{

				vertToggle->vGenerate();
				tieToggles.clear();
				contToggles = vertToggle->vContVarChanges();
				disToggles = vertToggle->vDisVarChanges();
				double lastLik = model->logLik();

				std::vector<int> oldDisValues = std::vector<int>(disToggles.size(),-1);
				for(int j=0;j<disToggles.size();j++){
					oldDisValues[j] =model->network()->discreteVariableValue(
											disToggles[j].second.first,
											disToggles[j].first);
					model->discreteVertexUpdate(
							disToggles[j].first,
							disToggles[j].second.first ,
							disToggles[j].second.second);
					model->network()->setDiscreteVariableValue(
							disToggles[j].second.first,
							disToggles[j].first,
							disToggles[j].second.second);
				}


				std::vector<double> oldContValues = std::vector<double>(contToggles.size(),-1);
				for(int j=0;j<contToggles.size();j++){
					//std::cout << contToggles[j].second.first << " " << contToggles[j].first << "\n";
					oldContValues[j] =model->network()->continVariableValue(
											contToggles[j].second.first,
											contToggles[j].first);
					//std::cout << oldContValues[j] << " " << contToggles[j].second.second << "\n";
					model->continVertexUpdate(
							contToggles[j].first,
							contToggles[j].second.first ,
							contToggles[j].second.second);
					model->network()->setContinVariableValue(
							contToggles[j].second.first,
							contToggles[j].first,
							contToggles[j].second.second);
				}


				double lr = model->logLik() - lastLik;
				lr += vertToggle->vLogRatio();
				if(lr>log(Rf_runif(0.0,1.0))){
					accepted++;
					vertToggle->vTogglesAccepted(true);
				}else{
					for(int j=((int)disToggles.size())-1;j>=0;j--){
						model->discreteVertexUpdate(
							disToggles[j].first,
							disToggles[j].second.first ,
							oldDisValues.at(j));
						model->network()->setDiscreteVariableValue(
							disToggles[j].second.first,
							disToggles[j].first,
							oldDisValues.at(j));
					}

					for(int j=((int)contToggles.size())-1;j>=0;j--){
						model->continVertexUpdate(
							contToggles[j].first,
							contToggles[j].second.first ,
							oldContValues.at(j));
						model->network()->setContinVariableValue(
							contToggles[j].second.first,
							contToggles[j].first,
							oldContValues.at(j));
					}

					vertToggle->vTogglesAccepted(false);
				}

			}
		}
		return ((double)accepted)/((double)steps);
	}

	/*!
	 * Generates a list of BinaryNet objects using MCMC.
	 *
	 * Exposed to R
	 */
	Rcpp::List generateSample(int burnIn,int interval,int sampleSize){
		model->calculate();
		GetRNGstate();
		initialize();
		this->run(burnIn);
		Rcpp::List lis;
		double ar = 0.0;
		for(int i=0;i<sampleSize-1;i++){
			R_CheckUserInterrupt();
			lis.push_back(model->network()->cloneR());
			ar += this->run(interval) / (sampleSize-1.0);
		}
		lis.push_back(model->network()->cloneR());
		lis.attr("acceptRatio") = ar;
		PutRNGstate();
		return lis;
	}

	/*!
	 * Generates a matrix of model statistics. Offset values are output as
	 * an attribute.
	 *
	 * Exposed to R
	 */
	NumericMatrix generateSampleStatistics(int burnIn,int interval,int sampleSize){
		std::vector<double> offs;
		std::vector<double> stats;
		model->calculate();
		NumericMatrix m(sampleSize,model->statistics().size());
		NumericMatrix off(sampleSize,model->offset().size());
		GetRNGstate();
		initialize();
		this->run(burnIn);
		double ar = 0.0;
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			if(i!=0)
				ar += this->run(interval) / (sampleSize-1.0);
			stats = model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			offs = model->offset();
			for(int j=0;j<offs.size();j++)
				off(i,j) = offs[j];
		}
		PutRNGstate();
	    List lis;
	    lis.push_back(R_NilValue);
		lis.push_back(model->names());
		m.attr("dimnames") = lis;
		if(offs.size()>0)
			m.attr("offset") = off;
		m.attr("acceptRatio") = ar;
		return m;
	}

	/*!
	 * Generates a matrix of sample statistics. Additionally it applies
	 * the R Function supplimentalFunction to each sample and returns the
	 * result.
	 *
	 * Exported to R. Really shows the power of Rcpp
	 */
	List generateSampleStatisticsSupplimental(int burnIn, int interval,
			int sampleSize, Function supplimentalFunction){
		model->calculate();
		NumericMatrix m(sampleSize,model->statistics().size());
		List sup;
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			if(i!=0)
				this->run(interval);
			std::vector<double> stats = model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			sup.push_back(supplimentalFunction(*model->network()));
		}
		PutRNGstate();
		List lis;
		lis.push_back(m);
		lis.push_back(sup);
		return lis;
	}


};


}

#endif /* METROPOLISHASTINGSH_ */
