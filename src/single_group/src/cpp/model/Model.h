/// \file
/// \brief An abstract model class for time-series models.
///
/// This file contains a class template for an abstract class that 
/// implements a generic time-series model

#ifndef __MODEL_H
#define __MODEL_H

#include "misc/misc.h"
#include "rng/Rng.h"

/// Class for the model specification.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation> class Model
{
public:
  
  /// Initialises the class.
  Model
  (
    Rng& rng
  ) : 
    rng_(rng)
  {
    nCores_ = 1;
  }
  
  /// Returns the number of observations.
  unsigned int getNObservations() {return observations_.size();}
  /// Returns the member holding the model parameters.
  const ModelParameters& getModelParameters() const {return modelParameters_;}
  /// Returns the latent variables.
  const std::vector<LatentVariable>& getLatentVariables() const {return latentVariables_;}
  /// Returns the $t$th latent variable.
  const LatentVariable& getLatentVariables(const unsigned int t) const {return latentVariables_[t];}
  /// Returns the covariates.
  const std::vector<Covariate>& getCovariates() const {return covariates_;}
  /// Returns the $t$th covariate.
  const Covariate& getCovariates(const unsigned int t) const {return covariates_[t];}
  /// Returns the observations.
  const std::vector<Observation>& getObservations() const {return observations_;}
  /// Returns the $t$th observation.
  const Observation& getObservations(const unsigned int t) const {return observations_[t];}
  /// Returns the length of the parameter vector.
  unsigned int getDimTheta() const {return dimTheta_;}
  /// Specifies the length of the parameter vector.
  void setDimTheta(const unsigned int dimTheta) 
  {
    dimTheta_ = dimTheta;
  }
  /// Specifies the covariates.
  void setCovariates(const std::vector<Covariate>& covariates)
  {
    covariates_ = covariates;
  }
  /// Specifies the observations.
  void setObservations(const std::vector<Observation>& observations)
  {
    observations_ = observations;
  }
  /// Generates synthetic data.
  void simulateData(const unsigned int nObservations, std::vector<LatentVariable>& latentVariables, std::vector<Covariate>& covariates, std::vector<Observation>& observations)
  {
    latentVariables.resize(nObservations);
    observations.resize(nObservations);
    latentVariables[0] = sampleFromInitialDistribution();
    for (unsigned int t=1; t<nObservations; t++)
    {
      latentVariables[t] = sampleFromTransitionEquation(t, latentVariables[t-1]);
    }
    for (unsigned int t=0; t<nObservations; t++)
    {
      observations[t] = sampleFromObservationEquation(t, latentVariables[t]);
    }
  }
  /// Generates synthetic data.
  void simulateData(const unsigned int nObservations)
  {
    simulateData(nObservations, latentVariables_, covariates_, observations_);
  }
  /// Samples a latent variable from the prior distribution at time $t=0$.
  LatentVariable sampleFromInitialDistribution()
  {
    LatentVariable latentVariableCurr;
    sampleFromInitialDistribution(latentVariableCurr);
    return latentVariableCurr;
  }
  /// Samples a latent variable from the transition equation at time $t>0$. 
  LatentVariable sampleFromTransitionEquation(const unsigned int t, const LatentVariable& latentVariablePrev)
  {
    LatentVariable latentVariableCurr;
    sampleFromTransitionEquation(t, latentVariableCurr, latentVariablePrev);
    return latentVariableCurr;
  }
  /// Samples a single observation from the observation equation at time $t>=0$.
  Observation sampleFromObservationEquation(const unsigned int t, const LatentVariable& latentVariableCurr)
  {
    Observation observation;
    sampleFromObservationEquation(t, observation, latentVariableCurr);
    return(observation);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user
  
  /// Specifies the model parameters which are to be inferred.
  void setUnknownParameters(const arma::colvec& theta);
  /// Specifies the known model parameters (i.e. those which we do not infer).
  void setKnownParameters(const arma::colvec& vartheta);
  /// Samples a single latent value from the prior on the transformed (generative-)model parameters.
  arma::colvec sampleFromParameterPrior();
  /// Samples a latent variable from the prior distribution at time $t=0$.
  void sampleFromInitialDistribution(LatentVariable& latentVariableCurr);
  /// Samples a latent variable from the transition equation at time $t>0$. 
  void sampleFromTransitionEquation(const unsigned int t, LatentVariable& latentVariableCurr, const LatentVariable& latentVariablePrev);
  /// Samples a single observation from the observation equation at time $t>=0$.
  void sampleFromObservationEquation(const unsigned int t, Observation& observation, const LatentVariable& latentVariableCurr);
  /// Evaluates the log-density of the initial state.
  double evaluateLogInitialDensity(const LatentVariable& latentVariable);
  /// Evaluates the log-transition density at some time $t > 0$.
  double evaluateLogTransitionDensity(const unsigned int t, const LatentVariable& latentVariableCurr, const LatentVariable& latentVariablePrev);
  /// Evaluates the log-observation density at some time $t$.
  double evaluateLogObservationDensity(const unsigned int t, const LatentVariable& latentVariableCurr);
  /// Evaluates the gradient of the logarithm of the prior density for
  /// the initial state.
  arma::colvec evaluateGradThetaLogInitialDensity(const LatentVariable& latentVariable);
  /// Evaluates the gradient of the logarithm of the 
  /// transition density
  arma::colvec evaluateGradThetaLogTransitionDensity(const unsigned int t, const LatentVariable& latentVariableCurr, const LatentVariable& latentVariablePrev);
  /// Evaluates the gradient of the logarithm of the observation density.
  arma::colvec evaluateGradThetaLogObservationDensity(const unsigned int t, const LatentVariable& latentVariableCurr);
  ////////////////////////////////////////////////////////////////////////////////////////////////////

private:
  
  Rng& rng_; // random number generation.
  std::vector<LatentVariable> latentVariables_; // stores the latent variables
  std::vector<Covariate> covariates_; // stores the covariates
  std::vector<Observation> observations_; // stores the observations
  ModelParameters modelParameters_; // holds the parameters of the (generative) model
  unsigned int dimTheta_; // dimension of the vector of (transformed) (generative-)model parameters
  unsigned int nCores_; // number of cores (this parameter is currently not used)
  
};

#endif 
