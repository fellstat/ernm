/*
 * SimpleToggle.h
 *
 *  Created on: Jun 28, 2011
 *      Author: ianfellows
 */

#ifndef SIMPLETOGGLEH_
#define SIMPLETOGGLEH_

#include "Rcpp.h"
#include <cmath>
#include <vector>
#include <assert.h>
#include <memory>
#include <boost/shared_ptr.hpp>
#include "util.h"
#include "BinaryNet.h"

namespace ernm{


template<class Engine>
class AbstractDyadToggle{
public:
	virtual void vSetNetwork(const boost::shared_ptr< BinaryNet<Engine> > n) = 0;
	virtual void vInitialize() = 0;
	virtual void vTogglesAccepted(bool apply) = 0;
	virtual void vGenerate() = 0;
	virtual double vLogRatio() = 0;
	virtual std::vector< std::pair<int,int> >& vDyadToggles() = 0;
	virtual std::string vName() = 0;
	virtual AbstractDyadToggle<Engine>* vCreateUnsafe(Rcpp::List params) = 0;
	virtual AbstractDyadToggle<Engine>* vCloneUnsafe() = 0;
	virtual ~AbstractDyadToggle(){};
};



/*!
 * The templated interface for dyad toggles.
 * Engine is the BinaryNetwork engine, and ToggleEngine is the toggle implementation class.
 */
template<class NetworkEngine, class ToggleEngine >
class DyadToggle : public AbstractDyadToggle<NetworkEngine>{
protected:
	ToggleEngine tog;

public:
	DyadToggle() : tog(){}

	DyadToggle(Rcpp::List l) : tog(l){}

	DyadToggle( BinaryNet<NetworkEngine> & network) : tog(network){}

	virtual ~DyadToggle(){}

	/*!
	 * Changes the network to generate toggles for.
	 */
	inline void setNetwork(const boost::shared_ptr< BinaryNet<NetworkEngine> > n){
		tog.setNetwork(n);
	}

	/*!
	 * Initializes the toggle. This should be called before generate.
	 */
	inline void initialize(){
		tog.initialize();
	}

	/*!
	 * Tells the toggle whether the proposed changes to network were accepted or not
	 */
	inline void togglesAccepted(bool apply){
		tog.togglesAccepted(apply);
	}

	/*!
	 * generates a new proposed change to the network
	 */
	inline void generate(){
		tog.generate();
	}


	/*!
	 * The MCMC log probability ratio for the last proposed change.
	 */
	inline double logRatio(){
		return tog.logRatio();
	}

	inline std::string name(){
		return tog.name();
	}

	/*!
	 * Gets a reference to the proposed dyad change.
	 */
	inline std::vector< std::pair<int,int> >& dyadToggles(){
		return tog.dyadToggles();
	}

	//virtual functions
	virtual void vSetNetwork(const boost::shared_ptr< BinaryNet<NetworkEngine> > n){
		setNetwork(n);
	}
	virtual void vInitialize(){
		initialize();
	}
	virtual void vTogglesAccepted(bool apply){
		togglesAccepted(apply);
	}
	virtual void vGenerate(){
		generate();
	}
	virtual double vLogRatio(){
		return logRatio();
	}
	virtual std::vector< std::pair<int,int> >& vDyadToggles(){
		return dyadToggles();
	}
	virtual std::string vName(){
		return name() ;
	}
	virtual AbstractDyadToggle<NetworkEngine>* vCreateUnsafe(Rcpp::List params){
		return new DyadToggle<NetworkEngine,ToggleEngine>(params);
	}
	virtual AbstractDyadToggle<NetworkEngine>* vCloneUnsafe(){
		return new DyadToggle<NetworkEngine,ToggleEngine>(*this);
	}

};






}

#endif /* SIMPLETOGGLEH_ */
