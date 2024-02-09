/// \file
/// \brief Implements an the adaptive online smoother from 
///
/// This file contains the functions for implementing 


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
    std::vector<double> functionalEstimatesAux; // stores the smoothed test function estimates associated with different (timeIndex, testFunctionIndex)-pairs which can be read from the following two vectors:
    std::vector<unsigned int> timeIndicesAux; // stores the time indices associated with the smoothed test function estimates
    std::vector<unsigned int> testFunctionIndicesAux; // stores the test function indices associated with the smoothed test function estimates
    
    unsigned int T = getNSteps();
    unsigned int R = getNTestFunctions();
    unsigned int S = T * R;
    smc_.initialise();
    if (useOnlineMarginalSmoothing_)
    {
//       std::cout << "useOnlineMarginalSmoothing at time /*0*/" << std::endl;
      functionalEstimatesAux.reserve(S);
      timeIndicesAux.reserve(S);
      testFunctionIndicesAux.reserve(S);
//       std::cout << "finished resizing marginal smoothing output vectors" << std::endl;
      
      onlineMarginalSmoothing_.initialise(functionalEstimatesAux, timeIndicesAux, testFunctionIndicesAux);
    }
    if (useOnlineParameterEstimation_)
    {
//       std::cout << "useOnlineParameterEstimation at time 0" << std::endl;
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
//         std::cout << "resample type: " << smc_.getSmcResampleType() << std::endl;
//         std::cout << "proposal type: " << smc_.getSmcProposalType() << std::endl;
      }

      smc_.iterate();
//       std::cout << "smc_.iterate() complete at time " << t << std::endl;
      smc_.evaluateBackwardKernels();
//       std::cout << "smc_.evaluateBackwardKernels() complete at time " << t << std::endl;
      
      if (useOnlineMarginalSmoothing_)
      {
//         std::cout << "useOnlineMarginalSmoothing at time " << t << std::endl;
        if (t == T-1)
        {
          onlineMarginalSmoothing_.setIsFinalStep(true);
        }
//         std::cout << "START: onlineMarginalSmoothing_.update() at time " << t << std::endl;
        onlineMarginalSmoothing_.update(functionalEstimatesAux, timeIndicesAux, testFunctionIndicesAux);
//         std::cout << "END: onlineMarginalSmoothing_.update() at time " << t << std::endl;
      }
      if (useOnlineParameterEstimation_)
      {
//         std::cout << "START: onlineParameterEstimation_.update() at time " << t << std::endl;
        onlineParameterEstimation_.update(thetaEstimates);
//         std::cout << "END: onlineParameterEstimation_.update() at time " << t << std::endl;
      }
    }
    
//     std::cout << "OnlineCombinedInference iterations complete!" << std::endl;
    
    // Rearranging the output:
    if (useOnlineMarginalSmoothing_)
    {
      functionalEstimates.resize(T);
      for (unsigned int t = 0; t < T; t++)
      {
        functionalEstimates[t].set_size(R);
      }
      unsigned int t, r;
      for (unsigned int s = 0; s < S; s++)
      {
        t = timeIndicesAux[s];
        r = testFunctionIndicesAux[s];
        functionalEstimates[t](r) = functionalEstimatesAux[s];
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
