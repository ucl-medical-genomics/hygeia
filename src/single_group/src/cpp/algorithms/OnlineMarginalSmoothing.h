/// \file
/// \brief Implements an the adaptive online smoother from 
///
/// This file contains the functions for implementing 


#ifndef __ONLINEMARGINALSMOOTHING_H
#define __ONLINEMARGINALSMOOTHING_H

#include "algorithms/Smc.h"

/// Class for 
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
    std::vector<double>& functionalEstimatesAux, 
    std::vector<unsigned int>& timeIndicesAux, 
    std::vector<unsigned int>& testFunctionIndicesAux
  )
  {
    psiCurr_.resize(0);
    psiTimeIndices_.resize(0);
    psiTestFunctionIndices_.resize(0);
    initialisePsi();
    storeEstimates(functionalEstimatesAux, timeIndicesAux, testFunctionIndicesAux);
  }
  /// Computes and updates all smoothing estimates at the current time step.
  void update
  (
    std::vector<double>& functionalEstimatesAux, 
    std::vector<unsigned int>& timeIndicesAux, 
    std::vector<unsigned int>& testFunctionIndicesAux
  )
  {
//     std::cout << "psiCurr_.swap(psiPrev_):" << std::endl;
    psiPrev_ = psiCurr_; // can't do .swap() here because both psiCurr_ needs to have the same size as psiPrev_
//       std::cout << "updatePsi():" << std::endl;
    updatePsi();
//       std::cout << "initialisePsi():" << std::endl;
    initialisePsi();
//       std::cout << "storeEstimates():" << std::endl;
    storeEstimates(functionalEstimatesAux, timeIndicesAux, testFunctionIndicesAux);
//     std::cout << "end storeEstimates()" << std::endl;
    
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
//     std::cout << "start initialisePsi()" << std::endl;
    unsigned int t = getStep();
    unsigned int N = getNParticlesCurr();
    unsigned int R = getNTestFunctions();
    
    unsigned int S = psiTimeIndices_.size();
    if (psiCurr_.capacity() <= psiCurr_.size() + R)
    {
      psiCurr_.reserve(S + 2 * R);
    }
    
    if (psiTimeIndices_.capacity() <= psiTimeIndices_.size() + R)
    {
      psiTimeIndices_.reserve(S + 2 * R);
    }
    if (psiTestFunctionIndices_.capacity() <= psiTestFunctionIndices_.size() + R)
    {
      psiTestFunctionIndices_.reserve(S + 2 * R);
    }
    std::vector<double> psiCurrAux(getNParticlesMax());
    for (unsigned int r = 0; r < R; r++)
    {
      psiTimeIndices_.push_back(t);
      psiTestFunctionIndices_.push_back(r);
      for (unsigned int n = 0; n < N; n++)
      {
        psiCurrAux[n] = evaluateTestFunction(t, r, getLatentVariableCurr(n));
      }
      psiCurr_.push_back(psiCurrAux);
    }
//     std::cout << "end initialisePsi()" << std::endl;
  }
  /// Evaluates the test function at some time 
  void updatePsi()
  {
//     std::cout << "start updatePsi()" << std::endl;
    unsigned int N = getNParticlesCurr();
    unsigned int S = psiCurr_.size();
    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      unsigned int R = getNTestFunctions();
      unsigned int M = N - R;
      for (unsigned int s = 0; s < S; s++)
      {
        for (unsigned int n = 0; n < M; n++)
        {
          psiCurr_[s][n] = psiPrev_[s][getAncestorIndicesCurr(n)];
        }
        for (unsigned int r = 0; r < R; r++)
        {
//           std::cout << "psiCurr_[s].size(): " << psiCurr_[s].size() << "; nParticlesCurr_: " << N << std::endl;
//           std::cout << "psiPrev_[s].size(): " << psiPrev_[s].size() << "; getBackwardKernels()[r].size(): " << getBackwardKernels()[r].size() << std::endl;
     
          psiCurr_[s][M + r] = 0;
          
          for (unsigned int n = 0; n < getNParticlesPrev(); n++)
          {
            psiCurr_[s][M + r] = psiCurr_[s][M + r] + getBackwardKernels()[r](n) * psiPrev_[s][n]; 
          }
        }
      }
    }
    else // i.e. for non-change point model SMC algorithms
    {
      for (unsigned int s = 0; s < S; s++)
      {
        for (unsigned int n = 0; n < N; n++)
        {
          psiCurr_[s][n] = 0;
          for (unsigned int m = 0; m < getNParticlesPrev(); m++)
          {
            psiCurr_[s][n] = psiCurr_[s][n] + getBackwardKernels()[n](m) * psiPrev_[s][m];
          }
        }
      }
    }
//     std::cout <<  "end updatePsi()" << std::endl;
  }
  /// Computes the estimates of the smoothed test function.
  void storeEstimates
  (
    std::vector<double>& functionalEstimatesAux, 
    std::vector<unsigned int>& timeIndicesAux, 
    std::vector<unsigned int>& testFunctionIndicesAux
  )
  {
//     std::cout << "start storeEstimates()" << std::endl;
    
    unsigned int S = psiCurr_.size();
    
//     std::cout << "S: " << S << std::endl;
//     std::cout << "psiCurr_.capacity(): " << psiCurr_.capacity() << std::endl;
    
    std::vector<unsigned int> psiTimeIndicesAux(0);
    std::vector<unsigned int> psiTestFunctionIndicesAux(0);
    std::vector<std::vector<double>> psiCurrAux(0);
    
//     std::cout << "storeEstimates(), reserve" << std::endl;
    psiTimeIndicesAux.reserve(S);
    psiTestFunctionIndicesAux.reserve(S);
    psiCurrAux.reserve(S);
    
//     std::cout << "storeEstimates(), loop" << std::endl;
    for (unsigned int s = 0; s < S; s++)
    {
//       std::cout << "s: " << s << std::endl;
      // Store the those estimates whose variance is below the threshold.
      if (smc_.computeFilteredVariance(psiCurr_[s]) < getEpsilon() || getIsFinalStep())
      {
        functionalEstimatesAux.push_back(smc_.computeFilteredMean(psiCurr_[s]));
        timeIndicesAux.push_back(psiTimeIndices_[s]);
        testFunctionIndicesAux.push_back(psiTestFunctionIndices_[s]);
      }
      else // Otherwise, keep the following estimates.
      {
        psiTimeIndicesAux.push_back(psiTimeIndices_[s]);
        psiTestFunctionIndicesAux.push_back(psiTestFunctionIndices_[s]);
        psiCurrAux.push_back(psiCurr_[s]);
      }
    }
    
//     std::cout << "soreEstimates(), swap" << std::endl;
    // Put the "kept" estimates back into the vectors.
    
//     std::cout << "psiTimeIndices_.capacity(): " << psiTimeIndices_.capacity()  << "; psiTimeIndicesAux.size(): " << psiTimeIndicesAux.size() << std::endl;
//     std::cout << "psiTestFunctionIndices_.capacity(): " << psiTestFunctionIndices_.capacity()  << "; psiTestFunctionIndicesAux.size(): " << psiTestFunctionIndicesAux.size() << std::endl;
//     std::cout << "psiCurr_.capacity(): " << psiCurr_.capacity()  << "; psiCurrAux.size(): " << psiCurrAux.size() << std::endl;
    

    psiTimeIndices_.swap(psiTimeIndicesAux);
    psiTestFunctionIndices_.swap(psiTestFunctionIndicesAux);
    psiCurr_.swap(psiCurrAux);
    
//     psiTimeIndices_ = psiTimeIndicesAux;
//     psiTestFunctionIndices_ = psiTestFunctionIndicesAux;
//     psiCurr_ = psiCurrAux;
    
//     std::cout << "end storeEstimates()" << std::endl;
  }

  Rng& rng_; // random number generation.
  Model<ModelParameters, LatentVariable, Covariate, Observation>& model_; // the targeted model.
  Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>& smc_; // the SMC algorithm
  double epsilon_; // the variance threshold 
  std::vector<std::vector<double>> psiCurr_, psiPrev_; // auxiliary quantities used for computing the estimates of smoothed functionals
  std::vector<unsigned int> psiTimeIndices_; // stores the auxiliary time indices of the estimated test functions
  std::vector<unsigned int> psiTestFunctionIndices_; // stores the auxiliary indices of the type of estimated test function  
  bool isFinalStep_; // has the final SMC step been reached?
};
#endif
