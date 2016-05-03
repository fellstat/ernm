/*
 * LatentOrderLikelihood.h
 *
 *  Created on: Jun 22, 2015
 *      Author: goodfellow
 */

#ifndef LATENTORDERLIKELIHOOD_H_
#define LATENTORDERLIKELIHOOD_H_

#include "Model.h"
#include "VertexToggles.h"
#include "ToggleController.h"
#include "DyadToggles.h"
#include "ShallowCopyable.h"

#include <cmath>
#include <Rcpp.h>
#include <assert.h>
#include <vector>

namespace ernm{

template<class Engine>
class LatentOrderLikelihood : public ShallowCopyable{
protected:
	typedef boost::shared_ptr< Model<Engine> > ModelPtr;
	ModelPtr model;
	ModelPtr noTieModel;

	template<class T>
	void shuffle(std::vector<T>& vec, long offset){
		for( int i=0; i < offset; i++){
			long ind = floor(Rf_runif(0.0,1.0)*offset);
			T tmp = vec[i];
			vec[i] = vec[ind];
			vec[ind] = tmp;
		}
	}

	void removeEdges(ModelPtr mod){
		long n = mod->network()->size();
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				mod->network()->removeEdge(i,j);
			}
		}
		mod->calculate();
	}
public:

	LatentOrderLikelihood(){}

	LatentOrderLikelihood(Model<Engine> mod){
		model = mod.clone();
		noTieModel = mod.clone();
		noTieModel->setNetwork(mod.network()->clone());
		removeEdges(noTieModel);
	}

	/*!
	 * R constructor for RCPP
	 *
	 */
	LatentOrderLikelihood(SEXP sexp){
		boost::shared_ptr<LatentOrderLikelihood> xp = unwrapRobject< LatentOrderLikelihood<Engine> >(sexp);
		model = xp->model;
		noTieModel = xp->noTieModel;
	}

	/*!
	 * coerce to R object. for RCPP
	 */
	operator SEXP() const{
		return wrapInReferenceClass(*this,Engine::engineName() + "LatentOrderLikelihood");
	}

	virtual ShallowCopyable* vShallowCopyUnsafe() const{
		return new LatentOrderLikelihood(*this);
	}

	~LatentOrderLikelihood(){}

	void setModel(const Model<Engine>& mod){
		model = mod.clone();
		noTieModel = mod.clone();
		noTieModel->setNetwork(mod.network()->clone());
		removeEdges(noTieModel);
	}

	void setThetas(std::vector<double> newThetas){
		model->setThetas(newThetas);
		noTieModel->setThetas(newThetas);
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

	List fullLogLik(int nOrders, double downsampleRate){
		List result;
		long n = model->network()->size();
		std::vector<long> vertices(n);
		for(int i=0; i<n;i++){
			//long ind = floor(Rf_runif(0.0,1.0)*n);
			vertices[i] = i;
		}
		for( int i=0; i<nOrders; i++){
			GetRNGstate();
			this->shuffle(vertices, n);
			PutRNGstate();
			result.push_back(this->conditionalLogLik(downsampleRate, vertices));
			//result.push_back(this->conditionalLogLikRandomDyad(downsampleRate));
		}
		return result;
	}

	List fullLogLikWithFunc(int nOrders, double downsampleRate, Function vertexOrderingFunction){
		List result;
		for( int i=0; i<nOrders; i++){
			GetRNGstate();
			std::vector<long> vertices = as< std::vector<long> >(vertexOrderingFunction());
			PutRNGstate();
			result.push_back(this->conditionalLogLik(downsampleRate, vertices));
		}
		return result;
	}

	List conditionalLogLikRandomDyad(double downsampleRate){
		std::cout << "enter conditionalLogLikRandomDyad\n";
		GetRNGstate();
		long n = model->network()->size();
		long nStats = model->thetas().size();
		std::vector<long> vertices(n);
		std::vector<std::pair<long, long> > dyad_order;
		for(int i=0; i<n;i++){
			//long ind = floor(Rf_runif(0.0,1.0)*n);
			vertices[i] = i;
			for(int j=i+1;j<n;j++){
				dyad_order.push_back(std::make_pair(i,j));
			}
		}
		this->shuffle(dyad_order,dyad_order.size());
		//this->shuffle(vertices, n);

		ModelPtr runningModel = noTieModel->clone();
		runningModel->setNetwork(noTieModel->network()->clone());
		runningModel->calculate();
		//std::cout << "n2 edges: " << noTieModel->network()->nEdges();
		List samples;
		double llik = runningModel->logLik();
		double llikInit = runningModel->logLik();
		double llikChange, ldenom, probTie;
		std::vector<double> terms = runningModel->statistics();
		std::vector<double>  newTerms = runningModel->statistics();
		std::vector<double> deriv(nStats, 0.0);
		NumericMatrix hessian(nStats, nStats);
		bool sample;
		double lpartition = 0.0;
		for(int i=0; i < dyad_order.size(); i++){
			int vertex = dyad_order[i].first;
			int alter = dyad_order[i].second;
			sample = Rf_runif(0.0,1.0) < downsampleRate;
			runningModel->statistics(terms);
			//std::cout <<"(" << vertex << ", " << alter << ")\n";
			if(runningModel->network()->hasEdge(vertex,alter)){
				Rf_error("Logic error: edge found where there should be none.\n");
			}
			llik = runningModel->logLik();
			//Calculate likelihood for vertex --> alter
			runningModel->dyadUpdate(vertex, alter);
			runningModel->network()->toggle(vertex, alter);
			llikChange = runningModel->logLik() - llik;
			//std::cout << llikChange;
			ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
			lpartition += ldenom;
			probTie = exp(llikChange - ldenom);
			runningModel->statistics(newTerms);
			bool hasEdge = model->network()->hasEdge(vertex, alter);
			for(int k=0; k < nStats; k++){
				double changeK = newTerms[k] - terms[k];
				deriv[k] += (hasEdge ? changeK : 0.0) - probTie * changeK;
				for(int l=0; l < nStats; l++){
					double changeL = newTerms[l] - terms[l];
					hessian(k,l) -= changeK * changeL * probTie * (1.0 - probTie);
				}
			}

			if(hasEdge){
				llik += llikChange;
			}else{
				runningModel->dyadUpdate(vertex, alter);
				runningModel->network()->toggle(vertex, alter);
			}

			if(sample){
				for(int k=0; k<terms.size(); k++){
					terms[k] = newTerms[k] - terms[k];
				}
				List s;
				s["hasEdge"] = wrap(model->network()->hasEdge(vertex, alter));
				s["changeStats"] = wrap(terms);
				samples.push_back(s);
			}

			if(runningModel->network()->isDirected()){
				//Calculate likelihood for alter --> vertex
				sample = Rf_runif(0.0,1.0) < downsampleRate;
				runningModel->statistics(terms);

				runningModel->dyadUpdate(alter, vertex);
				runningModel->network()->toggle(alter, vertex);
				llikChange = runningModel->logLik() - llik;
				ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
				lpartition += ldenom;
				probTie = exp(llikChange - ldenom);
				runningModel->statistics(newTerms);
				bool hasEdge = model->network()->hasEdge(alter, vertex);
				for(int k=0; k < nStats; k++){
					double changeK = newTerms[k] - terms[k];
					deriv[k] += (hasEdge ? changeK : 0.0) - probTie * changeK;
					for(int l=0; l < nStats; l++){
						double changeL = newTerms[l] - terms[l];
						hessian(k,l) -= changeK * changeL * probTie * (1.0 - probTie);
					}
				}
				if(hasEdge){
					llik += llikChange;
				}else{
					runningModel->dyadUpdate(alter, vertex);
					runningModel->network()->toggle(alter, vertex);
				}

				if(sample){
					for(int k=0; k<terms.size(); k++){
						terms[k] = newTerms[k] - terms[k];
					}
					List s;
					s["hasEdge"] = wrap(model->network()->hasEdge(alter, vertex));
					s["changeStats"] = wrap(terms);
					samples.push_back(s);
				}


			}
		}
		//runningModel->calculate();
		//newTerms = runningModel->statistics();
		//for(int k=0; k < nStats; k++){
		//	deriv[k] = newTerms[k] - deriv[k];
		//}
		//std::cout << "db:" << newTerms[0];
		PutRNGstate();
		List result;
		result["logLik"] = llik - llikInit;
		result["logPartition"] = lpartition;
		result["derivative"] = wrap(deriv);
		result["hessian"] = hessian;
		result["samples"] = samples;
		return result;
	}

	List conditionalLogLik(double downsampleRate, std::vector<long> vert_order){
		//std::cout << "enter conditionalLogLik\n";
		GetRNGstate();
		long n = model->network()->size();
		long nStats = model->thetas().size();
		std::vector<long> vertices(n);
		/*std::vector<std::pair<long, long> > dyad_order;
		for(int i=0; i<n;i++){
			//long ind = floor(Rf_runif(0.0,1.0)*n);
			vertices[i] = i;
			for(int j=i+1;j<n;j++){
				dyad_order.push_back(std::make_pair(i,j));
			}
		}
		this->shuffle(dyad_order,dyad_order.size());
		//this->shuffle(vertices, n);
		*/
		ModelPtr runningModel = noTieModel->clone();
		runningModel->setNetwork(noTieModel->network()->clone());
		runningModel->calculate();
		//std::cout << "n2 edges: " << noTieModel->network()->nEdges();
		List samples;
		double llik = runningModel->logLik();
		double llikChange, ldenom, probTie;
		std::vector<double> terms = runningModel->statistics();
		std::vector<double>  newTerms = runningModel->statistics();
		std::vector<double> deriv(nStats, 0.0);
		NumericMatrix hessian(nStats, nStats);
		for(int k=0; k < nStats; k++){
			for(int l=0; l < nStats; l++){
				hessian(k,l) = 0.0;
			}
		}
		bool sample;
		double lpartition = 0.0;
		for(int i=0; i < n; i++){
			int vertex = vert_order[i];
			this->shuffle(vert_order,i);

			for(int j=0; j < i; j++){
				int alter = vert_order[j];
				sample = Rf_runif(0.0,1.0) < downsampleRate;
				runningModel->statistics(terms);
				//std::cout <<"(" << vertex << ", " << alter << ")\n";
				if(runningModel->network()->hasEdge(vertex,alter)){
					Rf_error("Logic error: edge found where there should be none.\n");
				}
				llik = runningModel->logLik();
				//Calculate likelihood for vertex --> alter
				runningModel->dyadUpdate(vertex, alter);
				runningModel->network()->toggle(vertex, alter);
				llikChange = runningModel->logLik() - llik;
				//std::cout << llikChange;
				ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
				lpartition += ldenom;
				probTie = exp(llikChange - ldenom);
				runningModel->statistics(newTerms);
				bool hasEdge = model->network()->hasEdge(vertex, alter);
				for(int k=0; k < nStats; k++){
					double changeK = newTerms[k] - terms[k];
					deriv[k] += changeK * (hasEdge - probTie);  //(hasEdge ? changeK : 0.0) - probTie * changeK;
					for(int l=0; l < nStats; l++){
						double changeL = newTerms[l] - terms[l];
						hessian(k,l) -= changeK * changeL * probTie * (1.0 - probTie);
					}
				}

				if(hasEdge){
					llik += llikChange;
				}else{
					runningModel->dyadUpdate(vertex, alter);
					runningModel->network()->toggle(vertex, alter);
				}

				if(sample){
					for(int k=0; k<terms.size(); k++){
						terms[k] = newTerms[k] - terms[k];
					}
					List s;
					s["hasEdge"] = wrap(model->network()->hasEdge(vertex, alter));
					s["changeStats"] = wrap(terms);
					samples.push_back(s);
				}

				if(runningModel->network()->isDirected()){
					//Calculate likelihood for alter --> vertex
					sample = Rf_runif(0.0,1.0) < downsampleRate;
					runningModel->statistics(terms);

					runningModel->dyadUpdate(alter, vertex);
					runningModel->network()->toggle(alter, vertex);
					llikChange = runningModel->logLik() - llik;
					ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
					lpartition += ldenom;
					probTie = exp(llikChange - ldenom);
					runningModel->statistics(newTerms);
					bool hasEdge = model->network()->hasEdge(alter, vertex);
					for(int k=0; k < nStats; k++){
						double changeK = newTerms[k] - terms[k];
						deriv[k] += (hasEdge ? changeK : 0.0) - probTie * changeK;
						for(int l=0; l < nStats; l++){
							double changeL = newTerms[l] - terms[l];
							hessian(k,l) -= changeK * changeL * probTie * (1.0 - probTie);
						}
					}
					if(hasEdge){
						llik += llikChange;
					}else{
						runningModel->dyadUpdate(alter, vertex);
						runningModel->network()->toggle(alter, vertex);
					}

					if(sample){
						for(int k=0; k<terms.size(); k++){
							terms[k] = newTerms[k] - terms[k];
						}
						List s;
						s["hasEdge"] = wrap(model->network()->hasEdge(alter, vertex));
						s["changeStats"] = wrap(terms);
						samples.push_back(s);
					}
				}


			}
		}
		//runningModel->calculate();
		//newTerms = runningModel->statistics();
		//for(int k=0; k < nStats; k++){
		//	deriv[k] = newTerms[k] - deriv[k];
		//}
		//std::cout << "db:" << newTerms[0];
		PutRNGstate();
		List result;
		result["logLik"] = llik;
		result["logPartition"] = lpartition;
		result["derivative"] = wrap(deriv);
		result["hessian"] = hessian;
		result["samples"] = samples;
		return result;
	}

	SEXP generateNetwork(){
		GetRNGstate();
		long n = model->network()->size();
		std::vector<long> vertices(n);
		for(int i=0; i<n;i++){
			//long ind = floor(Rf_runif(0.0,1.0)*n);
			vertices[i] = i;
		}
		this->shuffle(vertices, n);
		PutRNGstate();
		return this->generateNetworkWithOrder(vertices);
		//return this->generateNetworkRandomDyad();
	}

	SEXP generateNetworkRandomDyad(){
		//std::cout << "enter conditionalLogLik\n";
		GetRNGstate();
		long n = model->network()->size();
		std::vector<std::pair<long, long> > dyad_order;
		for(int i=0; i<n;i++){
			//long ind = floor(Rf_runif(0.0,1.0)*n);
			for(int j=i+1;j<n;j++){
				dyad_order.push_back(std::make_pair(i,j));
			}
		}
		this->shuffle(dyad_order,dyad_order.size());

		long nStats = model->thetas().size();
		ModelPtr runningModel = noTieModel->clone();
		runningModel->setNetwork(noTieModel->network()->clone());
		runningModel->calculate();
		std::vector<double> eStats = runningModel->statistics();
		std::vector<double> terms = runningModel->statistics();
		std::vector<double>  newTerms = runningModel->statistics();
		NumericMatrix grad(nStats, nStats);
		//std::cout << "n2 edges: " << noTieModel->network()->nEdges();
		double llik = runningModel->logLik();
		double llikChange, ldenom, probTie;
		bool hasEdge = false;
		for(int i=0; i < dyad_order.size(); i++){
			int vertex = dyad_order[i].first;
			int alter = dyad_order[i].second;
			//std::cout <<"(" << vertex << ", " << alter << ")\n";
			if(runningModel->network()->hasEdge(vertex,alter)){
				Rf_error("Logic error: edge found where there should be none.\n");
			}
			llik = runningModel->logLik();
			runningModel->statistics(terms);
			runningModel->dyadUpdate(vertex, alter);
			runningModel->network()->toggle(vertex, alter);
			runningModel->statistics(newTerms);
			llikChange = runningModel->logLik() - llik;
			//std::cout << llikChange;
			ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
			probTie = exp(llikChange - ldenom);
			//std::cout << probTie << "\n";
			hasEdge = true;
			if(probTie < Rf_runif(0.0, 1.0)){
				runningModel->dyadUpdate(vertex, alter);
				runningModel->network()->toggle(vertex, alter);
				hasEdge = false;
			}

			for(int m=0; m<terms.size(); m++){
				eStats[m] += (newTerms[m] - terms[m]) * probTie;
			}

			for(int k=0; k < nStats; k++){
				double changeK = newTerms[k] - terms[k];
				for(int l=0; l < nStats; l++){
					double changeL = newTerms[l] - terms[l];
					double actStat = hasEdge ? changeL : 0.0;
					grad(k,l) += changeK * changeL * probTie * (1.0 - probTie) + changeK * probTie
							;
				}
			}
			if(runningModel->network()->isDirected()){
				if(runningModel->network()->hasEdge(alter, vertex)){
					Rf_error("Logic error: edge found where there should be none.\n");
				}
				llik = runningModel->logLik();
				runningModel->statistics(terms);
				runningModel->dyadUpdate(alter, vertex);
				runningModel->network()->toggle(alter, vertex);
				runningModel->statistics(newTerms);
				llikChange = runningModel->logLik() - llik;
				//std::cout << llikChange;
				ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
				probTie = exp(llikChange - ldenom);
				//std::cout << probTie << "\n";
				if(probTie < Rf_runif(0.0, 1.0)){
					runningModel->dyadUpdate(alter, vertex);
					runningModel->network()->toggle(alter, vertex);
				}
				for(int m=0; m<terms.size(); m++){
					eStats[m] += (newTerms[m] - terms[m]) * probTie;
				}
			}

		}
		PutRNGstate();
		List result;
		result["network"] = runningModel->network()->cloneR();
		result["stats"] = wrap(runningModel->statistics());
		result["expectedStats"] = wrap(eStats);
		return result;
	}

	SEXP generateNetworkWithOrder(std::vector<long> vert_order){
		//std::cout << "enter conditionalLogLik\n";
		GetRNGstate();
		long n = model->network()->size();
		/*std::vector<long> vertices(n);
		for(int i=0; i<n;i++){
			//long ind = floor(Rf_runif(0.0,1.0)*n);
			vertices[i] = i;
		}
		this->shuffle(vertices, n);
		std::vector<long> vert_order = vertices;*/
		long nStats = model->thetas().size();
		ModelPtr runningModel = noTieModel->clone();
		runningModel->setNetwork(noTieModel->network()->clone());
		runningModel->calculate();
		std::vector<double> eStats = std::vector<double>(nStats, 0.0);//runningModel->statistics();
		std::vector<double> stats = std::vector<double>(nStats, 0.0);
		std::vector<double> terms = runningModel->statistics();
		std::vector<double>  newTerms = runningModel->statistics();
		std::vector<double>  emptyStats = runningModel->statistics();
		NumericMatrix grad(nStats, nStats);
		NumericMatrix grad2(nStats, nStats);
		for(int k=0; k < nStats; k++){
			for(int l=0; l < nStats; l++){
				grad(k,l) = 0.0;
				grad2(k,l) = 0.0;
			}
		}
		//std::cout << "n2 edges: " << noTieModel->network()->nEdges();
		double llik = runningModel->logLik();
		double llikChange, ldenom, probTie;
		bool hasEdge = false;
		for(int i=0; i < n; i++){
			int vertex = vert_order[i];
			this->shuffle(vert_order,i);

			for(int j=0; j < i; j++){
				int alter = vert_order[j];
				//std::cout <<"(" << vertex << ", " << alter << ")\n";
				if(runningModel->network()->hasEdge(vertex,alter)){
					Rf_error("Logic error: edge found where there should be none.\n");
				}
				llik = runningModel->logLik();
				runningModel->statistics(terms);
				runningModel->dyadUpdate(vertex, alter);
				runningModel->network()->toggle(vertex, alter);
				runningModel->statistics(newTerms);
				llikChange = runningModel->logLik() - llik;
				//std::cout << llikChange;
				ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
				probTie = exp(llikChange - ldenom);
				//std::cout << probTie << "\n";
				hasEdge = true;
				if(probTie < Rf_runif(0.0, 1.0)){
					runningModel->dyadUpdate(vertex, alter);
					runningModel->network()->toggle(vertex, alter);
					hasEdge = false;
				}

				for(int k=0; k < nStats; k++){
					double changeK = newTerms[k] - terms[k];
					for(int l=0; l < nStats; l++){
						double changeL = newTerms[l] - terms[l];
						double actStat = hasEdge ? changeL : 0.0;
						grad(k,l) -= changeK * changeL * probTie * (1.0 - probTie) + changeK * probTie * (stats[l] - eStats[l]);
						grad2(k,l) -= changeK * probTie * (stats[l] - eStats[l]);
					}
				}

				for(int m=0; m<terms.size(); m++){
					eStats[m] += (newTerms[m] - terms[m]) * probTie;
					if(hasEdge)
						stats[m] += newTerms[m] - terms[m];
				}
				if(runningModel->network()->isDirected()){
					if(runningModel->network()->hasEdge(alter, vertex)){
						Rf_error("Logic error: edge found where there should be none.\n");
					}
					llik = runningModel->logLik();
					runningModel->statistics(terms);
					runningModel->dyadUpdate(alter, vertex);
					runningModel->network()->toggle(alter, vertex);
					runningModel->statistics(newTerms);
					llikChange = runningModel->logLik() - llik;
					//std::cout << llikChange;
					ldenom = R::log1pexp(llikChange);//log(1.0 + exp(llikChange));
					probTie = exp(llikChange - ldenom);
					//std::cout << probTie << "\n";
					if(probTie < Rf_runif(0.0, 1.0)){
						runningModel->dyadUpdate(alter, vertex);
						runningModel->network()->toggle(alter, vertex);
					}
					for(int k=0; k < nStats; k++){
						double changeK = newTerms[k] - terms[k];
						for(int l=0; l < nStats; l++){
							double changeL = newTerms[l] - terms[l];
							double actStat = hasEdge ? changeL : 0.0;
							grad(k,l) -= changeK * changeL * probTie * (1.0 - probTie) +
									changeK * probTie * (stats[l] - eStats[l]);
						}
					}

					for(int m=0; m<terms.size(); m++){
						eStats[m] += (newTerms[m] - terms[m]) * probTie;
						if(hasEdge)
							stats[m] += newTerms[m] - terms[m];
					}
				}
			}
		}
		PutRNGstate();
		List result;
		result["network"] = runningModel->network()->cloneR();
		result["emptyNetworkStats"] = wrap(emptyStats);
		result["stats"] = wrap(stats);
		result["expectedStats"] = wrap(eStats);
		result["gradient"] = grad;
		result["gradient1"] = grad;
		result["gradient2"] = grad2;
		return result;
	}
};


}
#endif /* LATENTORDERLIKELIHOOD_H_ */
