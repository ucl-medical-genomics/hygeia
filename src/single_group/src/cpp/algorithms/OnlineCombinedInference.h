/// \file
/// \brief Combines both marginal smoothing and online parameter estimates
///
/// This file contains the functions for implementing joint online parameter estimation and
/// online fixed-lag approximation of some smoothed functionals.


#ifndef __ONLINECOMBINEDINFERENCE_H
#define __ONLINECOMBINEDINFERENCE_H

#include "algorithms/OnlineMarginalSmoothing.h"
#include "algorithms/OnlineParameterEstimation.h"

/// Class for 
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle> class OnlineCombinedInference
{
public:
  
  /// Initialises the class.
  OnlineCombinedInference
  (
    Rng& rng,
    Model<ModelParameters, LatentVariable, Covariate, Observation>& model, 
    Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc,
    OnlineMarginalSmoothing<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& onlineMarginalSmoothing,
    OnlineParameterEstimation<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& onlineParameterEstimation
  ) : 
    rng_(rng), 
    model_(model),
    smc_(smc),
    onlineMarginalSmoothing_(onlineMarginalSmoothing),
    onlineParameterEstimation_(onlineParameterEstimation)
  {
    // Empty
  }
  
  /// Specifies total number of SMC steps.
  void setNSteps(const unsigned int nSteps) {nSteps_ = nSteps;}
  /// Specifies if we should estimate expectations of certain 
  /// test functions under the joint smoothing distribution.
  void setUseOnlineMarginalSmoothing(const bool useOnlineMarginalSmoothing) {useOnlineMarginalSmoothing_ = useOnlineMarginalSmoothing;}
  /// Specifies if we should estimate the model parameters via 
  /// an online stochastic-gradient ascent scheme.
  void setUseOnlineParameterEstimation(const bool useOnlineParameterEstimation) {useOnlineParameterEstimation_ = useOnlineParameterEstimation;}
  /// Runs the algorithm for joint online estimation of the model parameters
  /// and expectations of certain test functions under marginals of the joint
  /// smoothing distribution.
  void run
  (
    std::vector<arma::colvec>& functionalEstimates, 
    std::vector<arma::colvec>& thetaEstimates
  )
  {
    // stores the parameter estimates
    std::vector<arma::colvec> functionalEstimatesAux; // stores the smoothed test function estimates associated with different timeIndex which can be read from the following vector:
    std::vector<unsigned int> timeIndicesAux; // stores the time indices associated with the smoothed test function estimates
    
    unsigned int T = getNSteps();
    // unsigned int R = getNTestFunctions();
    unsigned int S = T;
    smc_.initialise();
    if (useOnlineMarginalSmoothing_)
    {
      functionalEstimatesAux.reserve(S);
      timeIndicesAux.reserve(S);
      onlineMarginalSmoothing_.initialise(functionalEstimatesAux, timeIndicesAux);
    }
    if (useOnlineParameterEstimation_)
    {
      thetaEstimates.reserve(T);
      onlineParameterEstimation_.initialise(thetaEstimates); 
    }
    
    for (unsigned int t = 1; t < T; t++)
    {
      
      if (t % 1000 == 0)
      {
        std::cout << "Iteration/time: " << t << "; nParticles: " << smc_.getNParticlesCurr()  << std::endl;
        if (useOnlineParameterEstimation_)
       {
          std::cout << "P: " << model_.getModelParameters().getP() << std::endl;
          std::cout << "omega: " << model_.getModelParameters().getOmega().t() << std::endl;
        } 
      }
      smc_.iterate();
      smc_.evaluateBackwardKernels();
      
      if (useOnlineMarginalSmoothing_)
      {
        if (t == T-1)
        {
          onlineMarginalSmoothing_.setIsFinalStep(true);
        }
        onlineMarginalSmoothing_.update(functionalEstimatesAux, timeIndicesAux);
      }
      if (useOnlineParameterEstimation_)
      {
        onlineParameterEstimation_.update(thetaEstimates);
      }
    }
    
    // Rearranging the output:
    // NOTE: If we are really storing all the estimates in memory, then this can be simplified by
    // simply storing the estimates for time step $t$ in the $t$th position of the vector.
    if (useOnlineMarginalSmoothing_)
    {
      functionalEstimates.resize(T);
      unsigned int t;
      for (unsigned int s = 0; s < S; s++)
      {
        t = timeIndicesAux[s];
        functionalEstimates[t] = functionalEstimatesAux[s];
      }
    }
  }
  
  
private:
   
  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  
  /// Returns the total number of SMC steps.
  unsigned int getNSteps() const {return nSteps_;}
  /// Returns the total number of test functions whose
  /// expectations we want to approximate at each time step.
  unsigned int getNTestFunctions() const {return onlineMarginalSmoothing_.getNTestFunctions();}

  Rng& rng_; // random number generation.
  Model<ModelParameters, LatentVariable, Covariate, Observation>& model_; // the targeted model.
  Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc_; // the SMC algorithm
  OnlineMarginalSmoothing<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& onlineMarginalSmoothing_; // online approximation of marginal smoothing functionals
  OnlineParameterEstimation<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& onlineParameterEstimation_; // online stochastic-gradient ascent estimation of the model parameters
  
  unsigned int nSteps_; // maximum number SMC steps
  bool useOnlineMarginalSmoothing_; // should we estimate expectations of certain integrals under marginals of the joint smoothing distribution?
  bool useOnlineParameterEstimation_; // should the model parameters be updated via an online stochastic-gradient ascent scheme?
  
};
#endif
