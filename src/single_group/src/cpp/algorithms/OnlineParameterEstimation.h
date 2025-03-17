/// \file
/// \brief Implements a SMC-based online gradient-ascent algorithm for parameter estimation.
///
/// This file contains the functions for implementing the online, sequential Monte Carlo based
/// stochastic-gradient-ascent type static-parameter estimation method from
/// Poyiadjis, G., Doucet, A., & Singh, S. S. (2011). Particle approximations of the score and observed information matrix in state space models with application to parameter estimation. Biometrika, 98(1), 65-80.


#ifndef __ONLINEPARAMETERESTIMATION_H
#define __ONLINEPARAMETERESTIMATION_H

#include "algorithms/Smc.h"
#include "misc/GradientAscent.h"

/// Class for the online SMC-based stochastic gradient-ascent procedure.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle> class OnlineParameterEstimation
{
public:
  
  /// Initialises the class.
  OnlineParameterEstimation
  (
    Rng& rng,
    Model<ModelParameters, LatentVariable, Covariate, Observation>& model, 
    Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc,
    GradientAscent& gradientAscent
  ) : 
    rng_(rng), 
    model_(model),
    smc_(smc),
    gradientAscent_(gradientAscent)
    // NOTE: we need to specify theta_ before using this class
  {
    // Empty
  }
  
  /// Specifies the parameter vector.
  void setTheta(const arma::colvec& theta) {theta_ = theta;}
  /// Specifies the number of SMC steps between parameter estimates.
  void setNStepsWithoutParameterUpdate(const unsigned int nStepsWithoutParameterUpdate) {nStepsWithoutParameterUpdate_ = nStepsWithoutParameterUpdate;}
  /// Computes the marginal smoothing estimates at the initial time step.
  void initialise(std::vector<arma::colvec>& thetaEstimates)
  {
    gradientAscent_.initialise(getDimTheta());
    initialisePhi();
    initialiseGradients();
    thetaEstimates.push_back(theta_);
    // NOTE: for simplicity, we do not update $\theta$ at the initial step
  }
  /// Computes and updates all smoothing estimates at the current time step.
  void update(std::vector<arma::colvec>& thetaEstimates)
  {
    updatePhi();
    if (getStep() % nStepsWithoutParameterUpdate_ == 0)
    {
      updateGradients();
      gradientAscent_.iterate(theta_, gradientCurr_ - gradientPrev_);
      model_.setUnknownParameters(theta_);
    }
    thetaEstimates.push_back(theta_);
  }
  
  
private:
   
  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user

  //////////////////////////////////////////////////////////////////////////////////////////////////// 

  /// Returns the number of particles from the current time step.
  unsigned int getNParticlesCurr() const {return smc_.getNParticlesCurr();}
  /// Returns the number of particles from the previous time step.
  unsigned int getNParticlesPrev() const {return smc_.getNParticlesPrev();}
  /// Returns the maximum number of particles.
  unsigned int getNParticlesMax() const {return smc_.getNParticlesMax();}
  /// Returns the length of the parameter vector $\theta$.
  unsigned int getDimTheta() const {return model_.getDimTheta();}
  /// Returns SMC step counter.
  unsigned int getStep() const {return smc_.getStep();}
  /// Returns the type of proposal kernel used by the SMC algorithm.
  const SmcProposalType& getSmcProposalType() const {return smc_.getSmcProposalType();}
  /// Returns the ancestor index of the $n$th particle from the previous time step.
  const Particle& getParticlesCurr(const unsigned int n) const {return smc_.getParticlesCurr()[n];}
  /// Returns the $n$th particle from the previous time step.
  const Particle& getParticlesPrev(const unsigned int n) const {return smc_.getParticlesPrev()[n];}
  /// Returns the $n$th particle from the previous time step.
  unsigned int getAncestorIndicesCurr(const unsigned int n) const {return smc_.getAncestorIndicesCurr()(n);}
  /// Returns the ancestor index of the $n$th particle from the previous time step.
  unsigned int getAncestorIndicesPrev(const unsigned int n) const {return smc_.getAncestorIndicesPrev()(n);}
  /// Returns the collection of backward kernels associated with the $n$th current particle.
  const std::vector<arma::colvec>& getBackwardKernels() const {return smc_.getBackwardKernels();}
  /// Returns the latent variable associated with the 
  /// $n$th particle from the current time step.
  const LatentVariable& getLatentVariableCurr(const unsigned int n) const {return smc_.getParticlesCurr()[n].getLatentVariable();}
  /// Returns the latent variable associated with
  /// the $n$th particle from the previous time step.
  const LatentVariable& getLatentVariablePrev(const unsigned int n) const {return smc_.getParticlesPrev()[n].getLatentVariable();}
  
  /// Returns the gradient of the unnormalised target density at the initial time step.
  /// evaluated at the latent variable associated with the $n$th particle.
  arma::colvec evaluateGradThetaLogInitialVarGamma(const unsigned int n)
  {
    return model_.evaluateGradThetaLogInitialDensity(getLatentVariableCurr(n)) +  model_.evaluateGradThetaLogObservationDensity(getStep(), getLatentVariableCurr(n));
  }
  /// Returns the gradient of the unnormalised target density at the current time step.
  /// evaluated at the latent variables associated with the $n$th particle at the current 
  /// time step and $m$th particle at the previous time step.
  arma::colvec evaluateGradThetaLogVarGamma(const unsigned int n, const unsigned int m)
  {
    return model_.evaluateGradThetaLogTransitionDensity(getStep(), getLatentVariableCurr(n), getLatentVariablePrev(m)) + 
    model_.evaluateGradThetaLogObservationDensity(getStep(), getLatentVariableCurr(n));
  }
  /// Initialises the auxiliary parameters from which
  /// the gradient estimates are computed.
  void initialisePhi()
  {
    unsigned int N = getNParticlesCurr();
    phiCurr_.resize(getNParticlesMax());
    phiPrev_.resize(getNParticlesMax());
    for (unsigned int n = 0; n < N; n++)
    {
      phiCurr_[n] = evaluateGradThetaLogInitialVarGamma(n);
      phiPrev_[n].set_size(getDimTheta());
    }
  }
  /// Initialises the gradient estimates.
  void initialiseGradients()
  {
    gradientPrev_.set_size(getDimTheta());
    gradientPrev_.zeros(); 
    gradientCurr_ = smc_.computeFilteredMean(phiCurr_);
  }
  /// Updates the auxiliary parameters from which the
  /// gradient estimates are computed.
  void updatePhi()
  {
    phiPrev_ = phiCurr_; // TODO: make this more efficient
    unsigned int N = getNParticlesCurr();
    phiCurr_.resize(N); // NOTE: this means every element of phiCurr_ needs to be resized too!
    
    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      unsigned int R = smc_.getNRegimes();
      unsigned int M = N - R;
      
      for (unsigned int n = 0; n < M; n++)
      {
        phiCurr_[n] = phiPrev_[getAncestorIndicesCurr(n)] + evaluateGradThetaLogVarGamma(n, getAncestorIndicesCurr(n));
      }
      for (unsigned int r = 0; r < R; r++)
      {
        phiCurr_[M + r].zeros(getDimTheta());
        for (unsigned int n = 0; n < getNParticlesPrev(); n++)
        {
          phiCurr_[M + r] = phiCurr_[M + r] + getBackwardKernels()[r](n) * (phiPrev_[n] + evaluateGradThetaLogVarGamma(M + r, n));
        }
      }
    }
    else // i.e. for non-change point model SMC algorithms
    {
      for (unsigned int n = 0; n < N; n++)
      {
        phiCurr_[n].zeros();
        for (unsigned int m = 0; m < getNParticlesPrev(); m++)
        {
          phiCurr_[n] = phiCurr_[n] + getBackwardKernels()[n](m) * (phiPrev_[m] + evaluateGradThetaLogVarGamma(n, m));
        }
      }
    }
  }
  /// Updates the gradient estimates.
  void updateGradients()
  {
    gradientPrev_ = gradientCurr_;
    gradientCurr_ = smc_.computeFilteredMean(phiCurr_);
  }

  Rng& rng_; // random number generation.
  Model<ModelParameters, LatentVariable, Covariate, Observation>& model_; // the targeted model.
  Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc_; // the SMC algorithm
  GradientAscent& gradientAscent_; // a class for dealing with gradient-ascent and ADAM optimisation
  arma::colvec theta_; // the current parameter vector
  bool normaliseGradients_; // should the gradients be normalised according to their $L_1$ norm?
  bool useAdam_; // should the ADAM optimiser be used (otherwise it is plain stochastic gradient ascent)
  unsigned int nStepsWithoutParameterUpdate_; // the number of SMC steps between parameter estimates
  std::vector<arma::colvec> phiCurr_, phiPrev_; // auxiliary quantities used for computing the estimates the gradient
  arma::colvec gradientCurr_, gradientPrev_; // the current and previous $\theta$-gradient estimates
  
};
#endif
