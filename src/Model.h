/*
 * Model.h
 *
 *  Created on: Jun 27, 2011
 *      Author: ianfellows
 */

#ifndef MODELH_
#define MODELH_

#include "Stat.h"
#include "Offset.h"
#include "StatController.h"
#include <vector>
#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <Rcpp.h>
#include <RcppCommon.h>

namespace ernm{


/*!
 * a representation of an ernm model
 */
template<class Engine>
class Model{
protected:
	typedef boost::shared_ptr< AbstractStat<Engine > > StatPtr;
	typedef std::vector< StatPtr >  StatVector;

	typedef boost::shared_ptr< AbstractOffset<Engine > > OffsetPtr;
	typedef std::vector< OffsetPtr >  OffsetVector;

	StatVector stats;						//!statistics
	OffsetVector offsets;
	boost::shared_ptr< BinaryNet<Engine> > net;			//!the relevant network

	//The domain
	boost::shared_ptr< bool > randomGraph;
	boost::shared_ptr< std::vector<int> > randomDiscreteVariables;
	boost::shared_ptr< std::vector<int> > randomContinVariables;
public:
	Model(){
		boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>());
		net=n;
		randomGraph = boost::shared_ptr< bool >(new bool);
		randomDiscreteVariables = boost::shared_ptr< std::vector<int> >(new std::vector<int>) ;
		randomContinVariables = boost::shared_ptr< std::vector<int> >(new std::vector<int>) ;
		*randomGraph = true;
	}

	Model(BinaryNet<Engine>& network){
		boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>(network));
		net = n;
		randomGraph = boost::shared_ptr< bool >(new bool);
		randomDiscreteVariables = boost::shared_ptr< std::vector<int> >(new std::vector<int>) ;
		randomContinVariables = boost::shared_ptr< std::vector<int> >(new std::vector<int>) ;
		*randomGraph = true;
	}

	Model(const Model& mod){
		stats = mod.stats;
		offsets = mod.offsets;
		net = mod.net;
		randomGraph = mod.randomGraph;
		randomDiscreteVariables = mod.randomDiscreteVariables;
		randomContinVariables = mod.randomContinVariables;
	}

	/*!
	 * if deep, then the model statistics are de-aliased
	 */
	Model(const Model& mod, bool deep){
		stats = mod.stats;
		offsets = mod.offsets;
		net = mod.net;
		randomGraph = mod.randomGraph;
		randomDiscreteVariables = mod.randomDiscreteVariables;
		randomContinVariables = mod.randomContinVariables;
		if(deep){
			for(int i=0;i<stats.size();i++)
				stats[i] = stats[i]->vClone();
			for(int i=0;i<offsets.size();i++)
				offsets[i] = offsets[i]->vClone();
			randomGraph = boost::shared_ptr< bool >(new bool);
			randomDiscreteVariables = boost::shared_ptr< std::vector<int> >(new std::vector<int>) ;
			randomContinVariables = boost::shared_ptr< std::vector<int> >(new std::vector<int>) ;
			*randomGraph = *mod.randomGraph;
			*randomDiscreteVariables = *mod.randomDiscreteVariables;
			*randomContinVariables = *mod.randomContinVariables;
		}
	}

	virtual ~Model(){}

	/*!
	 * R constructor for RCPP
	 *
	 */
	Model(SEXP sexp){
		boost::shared_ptr<Model> xp = unwrapRobject< Model<Engine> >(sexp);
		stats = xp->stats;
		offsets = xp->offsets;
		net = xp->net;
		randomGraph = xp->randomGraph;
		randomDiscreteVariables = xp->randomDiscreteVariables;
		randomContinVariables = xp->randomContinVariables;
	}

	/*!
	 * coerce to R object. for RCPP
	 */
	operator SEXP() const{
		return wrapInReferenceClass(*this,Engine::engineName() + "Model");
	}

	/*!
	 * clone the model (but not the network assosiated with it)
	 */
	boost::shared_ptr< Model<Engine> > clone() const{
		return boost::shared_ptr< Model<Engine> >(new Model<Engine>(*this, true));
	}

	void copy(Model<Engine>& mod){
		stats = mod.stats;
		offsets = mod.offsets;
		net = mod.net;
		randomGraph = mod.randomGraph;
		randomDiscreteVariables = mod.randomDiscreteVariables;
		randomContinVariables = mod.randomContinVariables;
	}

	void copy(Model<Engine>& mod,bool deep){
		net = mod.net;
		if(deep){
			stats.resize(mod.stats.size());
			offsets.resize(mod.offsets->vSize());
			for(int i=0;i<stats.size();i++)
				stats[i] = mod.stats[i]->vClone();
			for(int i=0;i<offsets->vSize();i++)
				offsets[i] = mod.offsets[i]->vClone();
			*randomGraph = *mod.randomGraph;
			*randomDiscreteVariables = *mod.randomDiscreteVariables;
			*randomContinVariables = *mod.randomContinVariables;
		}else{
			stats = mod.stats;
			offsets = mod.offsets;
			randomGraph = mod.randomGraph;
			randomDiscreteVariables = mod.randomDiscreteVariables;
			randomContinVariables = mod.randomContinVariables;
		}
	}

	/*!
	 * Does the domain of the model include the graph
	 */
	bool hasRandomGraph() const{
		return *randomGraph;
	}

	void setRandomGraph(bool random){
		*randomGraph = random;
	}

	/*!
	 * Which nodal variables are included in the domain of the model
	 */
	std::vector<int> randomVariables(bool discrete) const{
		if(discrete)
			return *randomDiscreteVariables;
		else
			return *randomContinVariables;
	}

	void setRandomVariables(std::vector<int> variables, bool discrete){
		if(discrete)
			*randomDiscreteVariables = variables;
		else
			*randomContinVariables = variables;
	}

	bool hasAnyRandomVariables() const{
		return randomDiscreteVariables->size()>0 || randomContinVariables->size()>0;
	}

	void setRandomVariablesR(std::vector<std::string> vars){
		std::vector<std::string> dv = net->discreteVarNames();
		std::vector<std::string> cv = net->continVarNames();
		std::vector<int> di,ci;
		int ind;
		for(int i=0;i<vars.size();i++){
			ind = indexOf(vars[i],dv);
			if(ind>=0){
				di.push_back(ind);
			}else{
				ind = indexOf(vars[i],cv);
				if(ind>=0)
					ci.push_back(ind);
				else
					Rf_error("Model::setRandomVariables : Unknown variable");
			}
		}
		*randomDiscreteVariables = di;
		*randomContinVariables = ci;
	}

	std::vector<std::string> getRandomVariablesR() const{
		std::vector<std::string> vars;
		std::vector<std::string> dv = net->discreteVarNames();
		std::vector<std::string> cv = net->continVarNames();
		for(int i=0;i<randomDiscreteVariables->size();i++){
			vars.push_back(dv.at(randomDiscreteVariables->at(i)));
		}
		for(int i=0;i<randomContinVariables->size();i++){
			vars.push_back(cv.at(randomContinVariables->at(i)));
		}
		return vars;
	}



	/*!
	 * the model terms
	 */
	std::vector<double> terms(){
		int n=0;
		for(int i=0;i<stats.size();i++){
			n += stats.at(i)->vSize();
		}
		for(int i=0;i<offsets.size();i++){
			n += offsets.at(i)->vSize();
		}
		std::vector<double> v(n,0.0);
		int c=0;
		for(int i=0;i<stats.size();i++){
			std::vector<double> vals = stats.at(i)->vValues();
			for(int j=0;j<vals.size();j++){
				v[c] = vals[j];
				c++;
			}
		}
		for(int i=0;i<offsets.size();i++){
			std::vector<double> vals = offsets.at(i)->vValues();
			for(int j=0;j<vals.size();j++){
				v[c] = vals[j];
				c++;
			}
		}
		return v;
	}

	/*!
	 * the model parameters
	 */
	std::vector<double> thetas(){
		int n=0;
		for(int i=0;i<stats.size();i++){
			n += stats.at(i)->vTheta().size();
		}
		std::vector<double> v(n,0.0);
		int c=0;
		for(int i=0;i<stats.size();i++){
			std::vector<double> vals = stats.at(i)->vTheta();
			for(int j=0;j<vals.size();j++){
				v[c] = vals[j];
				//cout << stats.at(i)->theta()[j];
				c++;
			}
		}
		return v;
	}

	/*!
	 * set the model paramters
	 */
	void  setThetas(std::vector<double> newThetas){
		int n=0;
		for(int i=0;i<stats.size();i++){
			n += stats.at(i)->vTheta().size();
		}
		if(newThetas.size()!= n){
			//Rcpp::Rcout << n  << " " << newThetas.size() << " ";
			::Rf_error("Model.setThetas: size mismatch:");
		}
		int c=0;
		for(int i=0;i<stats.size();i++){
			std::vector<double>* vals = &stats.at(i)->vTheta();
			for(int j=0;j<vals->size();j++){
				(*vals)[j] = newThetas[c];
				//cout << stats.at(i)->theta()[j];
				c++;
			}
		}
	}


	/*!
	 * the model statistics
	 */
	std::vector<double> statistics(){
		int n=0;
		for(int i=0;i<stats.size();i++){
			n += stats.at(i)->vSize();
		}
		std::vector<double> v(n,0.0);
		int c=0;
		for(int i=0;i<stats.size();i++){
			//std::vector<double> vals = stats.at(i)->vStatistics();
			for(int j=0;j<stats.at(i)->vStatistics().size();j++){
				v[c] = stats.at(i)->vStatistics()[j];
				c++;
			}
		}
		return v;
	}
	/*!
	 * Copy model statistics into v
	 */
	void statistics(std::vector<double>& v){
			int c=0;
			for(int i=0;i<stats.size();i++){
				//std::vector<double> vals = stats.at(i)->vStatistics();
				for(int j=0;j<stats.at(i)->vStatistics().size();j++){
					v[c] = stats.at(i)->vStatistics()[j];
					c++;
				}
			}
		}

	/*!
	 * returns statistics with names for R
	 */
	NumericVector statisticsR(){
		NumericVector res = wrap(statistics());
		res.attr("names") = wrap(names());
		return res;
	}

	/*!
	 * returns thetas with names for R
	 */
	NumericVector thetasR(){
		NumericVector res = wrap(thetas());
		res.attr("names") = wrap(names());
		return res;
	}

	/*!
	 * the model statistic names
	 */
	std::vector<std::string> names(){
		int n=0;
		for(int i=0;i<stats.size();i++){
			n += stats.at(i)->vSize();
		}
		std::vector<std::string> v(n,"??");
		int c=0;
		for(int i=0;i<stats.size();i++){
			std::vector<std::string> vals = stats.at(i)->vStatNames();
			for(int j=0;j<vals.size();j++){
				v[c] = vals[j];
				c++;
			}
		}
		return v;
	}

	/*!
	 * the model offsets
	 */
	std::vector<double> offset(){
		int n=0;
		for(int i=0;i<offsets.size();i++){
			n += offsets.at(i)->vSize();
		}
		std::vector<double> v(n,0.0);
		int c=0;
		for(int i=0;i<offsets.size();i++){
			std::vector<double> vals = offsets.at(i)->vValues();
			for(int j=0;j<vals.size();j++){
				v[c] = vals[j];
				c++;
			}
		}
		return v;
	}

	/*!
	 * the log likelihood of the model
	 */
	double logLik(){
		double ll = 0.0;
		for(int i=0;i<stats.size();i++){
			ll += stats[i]->vLogLik();
		}
		for(int i=0;i<offsets.size();i++){
			ll += offsets[i]->vLogLik();
		}
		return ll;
		//std::vector<double> vals = terms();
		//double d = 0.0;
		//for(int i=0;i<vals.size();i++)
		//	d += vals[i];
		//return d;
	}

	/*!
	 * add a statistic
	 */
	void addStatPtr(StatPtr  s){
		stats.push_back(s);
		s->vCalculate(*net);
	}

	/*!
	 * add a statistic to the model
	 */
	void addStat(const AbstractStat<Engine>&  s){
		StatPtr ps((&s)->clone());
		ps->vCalculate(*net);
		stats.push_back(ps);
	}

	/*!
	 * add a offset
	 */
	void addOffsetPtr(OffsetPtr  o){
		offsets.push_back(o);
		o->vCalculate(*net);
	}

	/*!
	 * add a offset to the model
	 */
	void addOff(const AbstractOffset<Engine>&  o){
		OffsetPtr ps((&o)->vClone());
		ps->vCalculate(*net);
		offsets.push_back(ps);
	}

	/*!
	 * add a statistic by name. uses StatController
	 */
	void addStatistic(const std::string name, Rcpp::List params){
		AbstractStat<Engine>* ps = StatController<Engine>::getStat(name, params);
		if(ps==NULL){
			::Rf_error("Invalid stat");
			return;
		}
		ps->vCalculate(*net);
		stats.push_back(StatPtr(ps));
	}

	/*!
	 * add a offset by name. uses StatController
	 */
	void addOffset(const std::string name, Rcpp::List params){
		AbstractOffset<Engine>* ps = StatController<Engine>::getOffset(name, params);
		if(ps==NULL){
			::Rf_error("Invalid offset");
			return;
		}
		ps->vCalculate(*net);
		offsets.push_back(OffsetPtr(ps));
	}

	/*!
	 * calculates statistics and offsets
	 */
	void calculate(){
		calculateStatistics();
		calculateOffsets();
	}

	/*!
	 * calculate the statistics
	 */
	void calculateStatistics(){
		for(int i=0;i<stats.size();i++){
			stats[i]->vCalculate(*net);
		}
	}

	/*!
	 * calculate the statistics
	 */
	void calculateOffsets(){
		for(int i=0;i<offsets.size();i++){
			offsets[i]->vCalculate(*net);
		}
	}

	void dyadUpdate(int from, int to){
		for(int k=0;k<stats.size();k++){
			stats[k]->vDyadUpdate(*net, from, to);
		}
		for(int k=0;k<offsets.size();k++){
			offsets[k]->vDyadUpdate(*net, from, to);
		}
	}


	void discreteVertexUpdate(int vertex, int variable, int newValue){
		for(int k=0;k<stats.size();k++)
			stats[k]->vDiscreteVertexUpdate(*net,vertex, variable, newValue);
		for(int k=0;k<offsets.size();k++)
			offsets[k]->vDiscreteVertexUpdate(*net,vertex, variable, newValue);
	}

	void continVertexUpdate(int vertex, int variable, double newValue){
		for(int k=0;k<stats.size();k++)
			stats[k]->vContinVertexUpdate(*net,vertex, variable, newValue);
		for(int k=0;k<offsets.size();k++)
			offsets[k]->vContinVertexUpdate(*net,vertex, variable, newValue);
	}

	/*!
	 * get the network
	 */
	boost::shared_ptr< BinaryNet<Engine> > network() const{
		return net;
	}

	SEXP getNetworkR() const{
		return wrap(*net);
	}

	/*!
	 * set the network
	 */
	void setNetwork(const BinaryNet<Engine>& network){
		boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>(network));
		net = n;
	}
	void setNetworkR(const BinaryNet<Engine>& network){
		boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>(network));
		net = n;
	}
	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > network){
		net = network;
	}

};

#include <Rcpp.h>

}

#endif /* MODELH_ */
