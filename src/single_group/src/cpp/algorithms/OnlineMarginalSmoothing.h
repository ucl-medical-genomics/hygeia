/// \file
/// \brief Implements an adaptive fixed-lag smoothing algorithm
///
/// This file contains the functions for implementing (a variant of) the adaptive
///  fixed-lag smoothing procedure from:
///
/// Alenlöv, J., & Olsson, J. (2019). Particle-based adaptive-lag online marginal smoothing in general state-space models. IEEE Transactions on Signal Processing, 67(21), 5571-5582.


#ifndef __ONLINEMARGINALSMOOTHING_H
#define __ONLINEMARGINALSMOOTHING_H

#include "algorithms/Smc.h"

/// Class for the online adaptive fixed-lag smoothing procedure.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation,  class SmcParameters, class Particle> class OnlineMarginalSmoothing
{
public:
  
  /// Initialises the class.
  OnlineMarginalSmoothing
  (
    Rng& rng,
    Model<ModelParameters, LatentVariable, Covariate, Observation>& model, 
    Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc
  ) : 
    rng_(rng), 
    model_(model),
    smc_(smc)
  {
    isFinalStep_ = false;
  }
  
  /// Specifies the threshold-parameter epsilon which governs when 
  /// we stop updating the smoothing estimates.
  void setEpsilon(const double epsilon) {epsilon_ = epsilon;}
  /// Specifies whether the final SMC step has been reached
  void setIsFinalStep(const bool isFinalStep) {isFinalStep_ = isFinalStep;}
  /// Computes the marginal smoothing estimates at the initial time step.
  void initialise
  (
    std::vector<arma::colvec>& functionalEstimatesAux,
    std::vector<unsigned int>& timeIndicesAux
  )
  {
    psiCurr_.resize(0);
    psiTimeIndices_.resize(0);
    initialisePsi();
    storeEstimates(functionalEstimatesAux, timeIndicesAux);
  }
  /// Computes and updates all smoothing estimates at the current time step.
  void update
  (
    std::vector<arma::colvec>& functionalEstimatesAux,
    std::vector<unsigned int>& timeIndicesAux
  )
  {
    psiPrev_ = psiCurr_; // can't do .swap() here because both psiCurr_ needs to have the same size as psiPrev_
    updatePsi();
    initialisePsi();
    storeEstimates(functionalEstimatesAux, timeIndicesAux);
    
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user

  /// Returns the number of test functions associated with each time step.
  unsigned int getNTestFunctions();
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  
private:
   
  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user

  /// Evaluates the test function at some time $t$.
  double evaluateTestFunction(const unsigned int timeIndex, const unsigned int testFunctionIndex, const LatentVariable& latentVariableCurr);
  //////////////////////////////////////////////////////////////////////////////////////////////////// 

  /// Returns the number of particles from the current time step.
  unsigned int getNParticlesCurr() const {return smc_.getNParticlesCurr();}
  /// Returns the number of particles from the previous time step.
  unsigned int getNParticlesPrev() const {return smc_.getNParticlesPrev();}
  /// Returns the maximum number of particles.
  unsigned int getNParticlesMax() const {return smc_.getNParticlesMax();}
  /// Returns SMC step counter.
  unsigned int getStep() const {return smc_.getStep();}
  /// Returns the threshold-parameter epsilon which governs when 
  /// we stop updating the smoothing estimates.
  double getEpsilon() const {return epsilon_;}
  /// Returns whether the final SMC step has been reached
  unsigned int getIsFinalStep() const {return isFinalStep_;}
  /// Returns the type of proposal kernel used by the SMC algorithm.
  const SmcProposalType& getSmcProposalType() const {return smc_.getSmcProposalType();}
  /// Returns the ancestor index of the $n$th particle from the previous time step.
  Particle& getParticlesCurr(const unsigned int n) const {return smc_.getParticlesCurr()[n];}
  /// Returns the $n$th particle from the previous time step.
  Particle& getParticlesPrev(const unsigned int n) const {return smc_.getParticlesPrev()[n];}
  /// Returns the $n$th particle from the previous time step.
  unsigned int getAncestorIndicesCurr(const unsigned int n) const {return smc_.getAncestorIndicesCurr()(n);}
  /// Returns the ancestor index of the $n$th particle from the previous time step.
  unsigned int getAncestorIndicesPrev(const unsigned int n) const {return smc_.getAncestorIndicesPrev()(n);}
  /// Returns the latent variable associated with the 
  /// $n$th particle from the current time step.
  const LatentVariable& getLatentVariableCurr(const unsigned int n) const {return smc_.getParticlesCurr()[n].getLatentVariable();}
  /// Returns the latent variable associated with
  /// the $n$th particle from the previous time step.
  const LatentVariable& getLatentVariablePrev(const unsigned int n) const {return smc_.getParticlesPrev()[n].getLatentVariable();}
  /// Returns the collection of backward kernels associated with the $n$th current particle.
  const std::vector<arma::colvec>& getBackwardKernels() const {return smc_.getBackwardKernels();}
  /*
  /// Writes the final smoothing estimates to a file along with 
  /// the associated time and test-function indices.
  void writeToFile(const unsigned int timeIndex, const unsigned int testFunctionIndex, const double estimate)
  {
    // TODO
  }
  */
  /// Evaluates the test function at some time 
  void initialisePsi()
  {
    unsigned int t = getStep();
    unsigned int N = getNParticlesCurr();
    unsigned int R = getNTestFunctions();
    
    unsigned int S = psiTimeIndices_.size();
    if (psiCurr_.capacity() <= psiCurr_.size() + 1)
    {
      psiCurr_.reserve(S + 2);
    }
    
    if (psiTimeIndices_.capacity() <= psiTimeIndices_.size() + 1)
    {
      psiTimeIndices_.reserve(S + 2);
    }
    std::vector<std::vector<double>> psiCurrAux(R);
    for (unsigned int r = 0; r < R; r++)
    {
      psiCurrAux[r].resize(getNParticlesMax());
      for (unsigned int n = 0; n < N; n++)
      {
        psiCurrAux[r][n] = evaluateTestFunction(t, r, getLatentVariableCurr(n));
      }
    }
    psiTimeIndices_.push_back(t);
    psiCurr_.push_back(psiCurrAux);
  }
  /// Evaluates the test function at some time 
  void updatePsi()
  {
    unsigned int N = getNParticlesCurr();
    unsigned int S = psiCurr_.size();
    unsigned int R = getNTestFunctions();

    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      unsigned int R = getNTestFunctions();
      unsigned int M = N - R;
      for (unsigned int s = 0; s < S; s++)
      {
        for (unsigned int r = 0; r < R; r++)
        {
          for (unsigned int n = 0; n < M; n++)
          {
            psiCurr_[s][r][n] = psiPrev_[s][r][getAncestorIndicesCurr(n)];
          }
          for (unsigned int q = 0; q < R; q++)
          {
            psiCurr_[s][r][M + q] = 0;

            for (unsigned int n = 0; n < getNParticlesPrev(); n++)
            {
              psiCurr_[s][r][M + q] = psiCurr_[s][r][M + q] + getBackwardKernels()[q](n) * psiPrev_[s][r][n];
            }
          }
        }
      }
    }
    else // i.e. for non-change point model SMC algorithms
    {
      for (unsigned int s = 0; s < S; s++)
      {
        for (unsigned int r = 0; r < R; r++)
        {
          for (unsigned int n = 0; n < N; n++)
          {
            psiCurr_[s][r][n] = 0;
            for (unsigned int m = 0; m < getNParticlesPrev(); m++)
            {
              psiCurr_[s][r][n] = psiCurr_[s][r][n] + getBackwardKernels()[n](m) * psiPrev_[s][r][m];
            }
          }
        }
      }
    }
  }
  /// Computes the estimates of the smoothed test function.
  void storeEstimates
  (
    std::vector<arma::colvec>& functionalEstimatesAux,
    std::vector<unsigned int>& timeIndicesAux
  )
  {

    unsigned int S = psiCurr_.size();
    unsigned int R = getNTestFunctions();
    
    std::vector<unsigned int> psiTimeIndicesAux(0);
    std::vector<std::vector<std::vector<double>>> psiCurrAux(0);
    
    psiTimeIndicesAux.reserve(S);
    psiCurrAux.reserve(S);
    
    for (unsigned int s = 0; s < S; s++)
    {

      bool storeEstimateAux = true;
      if (!getIsFinalStep())
      {
        /// NOTE: In slight contrast to Alenlöv, J., & Olsson, J. (2019), we only store the estimate
        /// of the $r$th functional from time $t$ once the variances of /all/ functionals
        /// from time $t$ are below the threshold $\varepsilon$. In the change-point model ,
        /// this ensures that the regime-probability estimates are guaranteed to sum to $1$.
        unsigned int r = 0;
        while ((r < R) && (smc_.computeFilteredVariance(psiCurr_[s][r]) < getEpsilon()))
        {
          r++;
        }
        if (r < R) {
          storeEstimateAux = false;
        }
      }


      if (storeEstimateAux)
      {
        arma::colvec functionalEstimateAuxSingle(R);
        for (unsigned int r = 0; r < R; ++r)
        {
          functionalEstimateAuxSingle(r) = smc_.computeFilteredMean(psiCurr_[s][r]);
        }
        functionalEstimatesAux.push_back(functionalEstimateAuxSingle);
        timeIndicesAux.push_back(psiTimeIndices_[s]);
      }
      else // Otherwise, keep the following estimates.
      {
        psiTimeIndicesAux.push_back(psiTimeIndices_[s]);
        psiCurrAux.push_back(psiCurr_[s]);
      }
    }
    

    psiTimeIndices_.swap(psiTimeIndicesAux);
    psiCurr_.swap(psiCurrAux);
    
  }

  Rng& rng_; // random number generation.
  Model<ModelParameters, LatentVariable, Covariate, Observation>& model_; // the targeted model.
  Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc_; // the SMC algorithm
  double epsilon_; // the variance threshold 
  std::vector<std::vector<std::vector<double>>> psiCurr_, psiPrev_; // auxiliary quantities used for computing the estimates of smoothed functionals
  std::vector<unsigned int> psiTimeIndices_; // stores the auxiliary time indices of the estimated test functions
  bool isFinalStep_; // has the final SMC step been reached?
};
#endif
