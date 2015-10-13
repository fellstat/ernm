/*
 * VertexToggles.h
 *
 *  Created on: Jul 20, 2011
 *      Author: ianfellows
 */

#ifndef VERTEXTOGGLESH_
#define VERTEXTOGGLESH_

#include "BinaryNet.h"
#include <vector>


namespace ernm{

template<class Engine>
class AbstractVertexToggle{
public:
	virtual void vSetNetwork(const boost::shared_ptr< BinaryNet<Engine> > n) = 0;
	virtual void vInitialize() = 0;
	virtual void vTogglesAccepted(bool apply) = 0;
	virtual void vGenerate() = 0;
	virtual double vLogRatio() = 0;
	virtual std::vector<std::pair<int,std::pair<int,double> > >& vContVarChanges() = 0;
	virtual std::vector<std::pair<int,std::pair<int,int> > >& vDisVarChanges() = 0;
	virtual void vSetDiscreteVars(std::vector<int>& inds) = 0;
	virtual void vSetContinuousVars(std::vector<int>& inds) = 0;
	virtual std::string vName() = 0;
	virtual AbstractVertexToggle<Engine>* vCreateUnsafe(Rcpp::List params) = 0;
	virtual AbstractVertexToggle<Engine>* vCloneUnsafe() = 0;
	virtual ~AbstractVertexToggle(){};
};


/*!
 * The templated interface for vertex toggles.
 * Engine is the BinaryNetwork engine, and ToggleEngine is the toggle implementation class.
 */
template<class NetworkEngine, class ToggleEngine >
class VertexToggle : public AbstractVertexToggle<NetworkEngine>{
protected:
	ToggleEngine tog;

public:
	VertexToggle() : tog(){}

	VertexToggle(Rcpp::List l) : tog(l){}

	VertexToggle( BinaryNet<NetworkEngine> & network) : tog(network){}

	VertexToggle(BinaryNet<NetworkEngine>& network,
							std::vector<int>& cvars,std::vector<int>& dvars) : tog(network,cvars,dvars){}

	virtual ~VertexToggle(){}

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

	/*!
	 * The name of the toggle to be exposed in R
	 */
	inline std::string name(){
		return tog.name();
	}

	/*!
	 * sets which discrete indices should be toggled
	 */
	void setDiscreteVars(std::vector<int>& inds){
		tog.setDiscreteVars(inds);
	}

	/*!
	 * Sets the continuous variables to be toggled
	 */
	void setContinuousVars(std::vector<int>& inds){
		tog.setContinuousVars(inds);
	}

	/*!
	 * changed continuous variables <vertex, <varId, new value> >
	 */
	inline std::vector<std::pair<int,std::pair<int,double> > >& contVarChanges(){
		return tog.contVarChanges();
	}

	/*!
	 * changed discrete variables <vertex, <varId, new value> >
	 */
	inline std::vector<std::pair<int,std::pair<int,int> > >& disVarChanges(){
		return tog.disVarChanges();
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

	virtual std::vector<std::pair<int,std::pair<int,double> > >& vContVarChanges(){
		return contVarChanges();
	}

	virtual std::vector<std::pair<int,std::pair<int,int> > >& vDisVarChanges(){
		return disVarChanges();
	}

	virtual void vSetDiscreteVars(std::vector<int>& inds){
		setDiscreteVars(inds);
	}

	virtual void vSetContinuousVars(std::vector<int>& inds){
		setContinuousVars(inds);
	}

	virtual std::string vName(){
		return name() ;
	}

	virtual AbstractVertexToggle<NetworkEngine>* vCreateUnsafe(Rcpp::List params){
		return new VertexToggle<NetworkEngine,ToggleEngine>(params);
	}

	virtual AbstractVertexToggle<NetworkEngine>* vCloneUnsafe(){
		return new VertexToggle<NetworkEngine,ToggleEngine>(*this);
	}


};




}

#endif /* VERTEXTOGGLESH_ */
