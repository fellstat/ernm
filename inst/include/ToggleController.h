/*
 * ToggleController.h
 *
 *  Created on: Jan 3, 2014
 *      Author: ianfellows
 */

#ifndef TOGGLECONTROLLER_H_
#define TOGGLECONTROLLER_H_

#include <string>
#include <map>
#include "DyadToggle.h"
#include "VertexToggle.h"

namespace ernm {

template<class Engine>
class ToggleController {
private:
	typedef typename boost::shared_ptr< AbstractDyadToggle<Engine> > DtPtr;
	typedef typename boost::shared_ptr< std::map< std::string, DtPtr > > DtMapPtr;
	typedef typename boost::shared_ptr< AbstractVertexToggle<Engine> > VtPtr;
	typedef typename boost::shared_ptr< std::map< std::string, VtPtr > > VtMapPtr;
protected:
	static DtMapPtr dyadMapPtr;
	static VtMapPtr vertexMapPtr;

public:

	ToggleController(){};

	virtual ~ToggleController(){};
	static void init(){
		if(!dyadMapPtr)
			dyadMapPtr = DtMapPtr(new std::map< std::string, DtPtr >);
		if(!vertexMapPtr)
			vertexMapPtr = VtMapPtr(new std::map< std::string, VtPtr >);
	}
	static void addToggle(DtPtr pS){
		init();
		dyadMapPtr->insert(std::make_pair(pS->vName(),pS));
	}

	static void addToggle(VtPtr pS){
		init();
		vertexMapPtr->insert(std::make_pair(pS->vName(),pS));
	}

	static std::vector<std::string> getAvailableDyadToggles(){
		std::vector<std::string> r;
		typename std::map< std::string, DtPtr >::iterator it;
		for(it = dyadMapPtr->begin();
				it != dyadMapPtr->end();++it){
			r.push_back(it->first);
		}
		return r;
	}

	static std::vector<std::string> getAvailableVertexToggles(){
		std::vector<std::string> r;
		typename std::map< std::string, VtPtr >::iterator it;
		for(it = vertexMapPtr->begin();
				it != vertexMapPtr->end();++it){
			r.push_back(it->first);
		}
		return r;
	}

	static AbstractDyadToggle<Engine>* getDyadToggle(std::string name, Rcpp::List params){
		DtPtr pS;
		try{
			pS = dyadMapPtr->at(name);

		}catch(...){
			std::string warning = std::string("Could not find dyad toggler: ") + name + ". Available togglers are: ";
			std::vector<std::string> dt = getAvailableDyadToggles();
			for(int i=0;i< dt.size();i++)
				warning += std::string(i==0 ? "" : ", ") + dt[i];
			Rf_error("%s", warning.c_str());
			return NULL;
		}
		if(pS==NULL){
			Rf_error("%s",(std::string("Could not find dyad toggler: ") + name).c_str());
			return NULL;
		}
		return pS->vCreateUnsafe(params);
	}

	static AbstractVertexToggle<Engine>* getVertexToggle(std::string name, Rcpp::List params){
		VtPtr pS;
		try{
			pS = vertexMapPtr->at(name);
		}catch(...){
			Rf_error("%s",(std::string("Could not find vertex toggler: ") + name).c_str());
			return NULL;
		}
		if(pS==NULL){
			Rf_error("%s",(std::string("Could not find vertex toggler: ") + name).c_str());
			return NULL;
		}
		return pS->vCreateUnsafe(params);
	}
};

template<class Engine>
void registerToggle(boost::shared_ptr< AbstractDyadToggle<Engine> > ps){
	ToggleController<Engine>::addToggle(ps);
}

template<class Engine>
void registerToggle(boost::shared_ptr< AbstractVertexToggle<Engine> > ps){
	ToggleController<Engine>::addToggle(ps);
}

} /* namespace ernm */


/*!
 * Called upon loading ERNM, registering built-in togglers
 */
void initToggles();

void registerDirectedDyadToggle(Rcpp::XPtr< ernm::AbstractDyadToggle<ernm::Directed> > ps);

void registerUndirectedDyadToggle(Rcpp::XPtr< ernm::AbstractDyadToggle<ernm::Undirected> > ps);

void registerDirectedVertexToggle(Rcpp::XPtr< ernm::AbstractVertexToggle<ernm::Directed> > ps);

void registerUndirectedVertexToggle(Rcpp::XPtr< ernm::AbstractVertexToggle<ernm::Undirected> > ps);

#endif /* TOGGLECONTROLLER_H_ */
