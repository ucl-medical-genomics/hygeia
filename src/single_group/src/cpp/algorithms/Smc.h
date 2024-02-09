/// \file
/// \brief Implements a sequential Monte Carlo algorithm for change-point models.
///
/// This file implements a sequential Monte Carlo algorithm for 
/// models for change-point models.


#ifndef __SMC_H
#define __SMC_H

#include "model/Model.h"
#include "misc/resample.h"

/// Specifiers for various resampling algorithms used by the 
/// sequential Monte Carlo algorithm for change-point models.
enum SmcResampleType
{ 
  SMC_RESAMPLE_MULTINOMIAL = 0, 
  SMC_RESAMPLE_SYSTEMATIC,
  SMC_RESAMPLE_OPTIMAL_FINITE_STATE
};
/// Specifiers for various proposal kernels used in 
/// the simple (i.e. non-change-point model)
/// sequential Monte Carlo algorithm.
enum SmcProposalType 
{ 
  SMC_PROPOSAL_PRIOR = 0,
  SMC_PROPOSAL_CHANGE_POINT_MODEL
};

/// Class for running sequential Monte Carlo algorithms.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle> class Smc
{
public:
  
  /// Initialises the class.
  Smc
  (
    Rng& rng,
    Model<ModelParameters, LatentVariable, Covariate, Observation>& model
  ) : 
    rng_(rng), 
    model_(model)
  {
    // Empty
  }
  
  /// Returns the SMC step counter.
  unsigned int getStep() const {return step_;}
  /// Returns the SMC parameters.
  const SmcParameters& getSmcParameters() const {return smcParameters_;}
  /// Returns the SMC parameters (non-constant refernece allows modification).
  SmcParameters& getReferenceToSmcParameters() {return smcParameters_;}
  /// Returns the current number of particles.
  unsigned int getNParticlesCurr() const {return nParticlesCurr_;}
  /// Returns the previous number of particles.
  unsigned int getNParticlesPrev() const {return nParticlesPrev_;}
  /// Returns the maximum number of particles.
  unsigned int getNParticlesMax() const {return nParticlesMax_;}
  /// Returns the number of regimes.
  unsigned int getNRegimes() const {return nRegimes_;}
  /// Returns the type of resample kernel to use.
  const SmcResampleType& getSmcResampleType() const {return smcResampleType_;}
  /// Returns the type of proposal kernel to use.
  const SmcProposalType& getSmcProposalType() const {return smcProposalType_;}
  /// Returns the particles associated with the previous time step.
  const std::vector<Particle>& getParticlesPrev() const {return particlesPrev_;}
  /// Returns the particles associated with the current time step.
  const std::vector<Particle>& getParticlesCurr() const {return particlesCurr_;}
  /// Returns a single particle associated with the previous time step
  const Particle& getParticlesPrev(const unsigned int n) const {return particlesPrev_[n];}
  /// Returns a single particle associated with the current time step
  const Particle& getParticlesCurr(const unsigned int n) const {return particlesCurr_[n];}
  /// Returns the ancestor indices associated with the previous time step.
  const arma::uvec& getAncestorIndicesPrev() const {return ancestorIndicesPrev_;}
  /// Returns the ancestor indices associated with the current time step.
  const arma::uvec& getAncestorIndicesCurr() const {return ancestorIndicesCurr_;}
  /// Returns the backward kernels for the current time step.
  const std::vector<arma::colvec>& getBackwardKernels() const {return backwardKernels_;}
  /// Returns the log-unnormalised particle weights 
  /// associated with the previous time step.
  const arma::colvec& getLogUnnormalisedWeightsPrev() const {return logUnnormalisedWeightsPrev_;}
  /// Returns the log-unnormalised particle weights 
  /// associated with the current time step.
  const arma::colvec& getLogUnnormalisedWeightsCurr() const {return logUnnormalisedWeightsCurr_;}
  /// Returns the self-normalised particle weights 
  /// associated with the previous time step.
  const arma::colvec& getSelfNormalisedWeightsPrev() const {return selfNormalisedWeightsPrev_;}
  /// Returns the self-normalised particle weights 
  /// associated with the current time step.
  const arma::colvec& getSelfNormalisedWeightsCurr() const {return selfNormalisedWeightsCurr_;}
  /// Returns the log of the normalising-constant estimate
  /// associated with the previous time step.
  double getLogSumOfUnnormalisedWeightsPrev() const {return logSumOfUnnormalisedWeightsPrev_;}
  /// Returns the log of the normalising-constant estimate
  /// associated with the current time step.
  double getLogSumOfUnnormalisedWeightsCurr() const {return logSumOfUnnormalisedWeightsCurr_;}
  
  /// Specifies the type of resampling scheme to use.
  void setSmcResampleType(const SmcResampleType smcResampleType) {smcResampleType_ = smcResampleType;}
  /// Specifies the type of proposal kernel to use.
  void setSmcProposalType(const SmcProposalType smcProposalType) {smcProposalType_ = smcProposalType;}
  /// Specifies the total number of particles to use.
  void setNParticlesCurr(const unsigned int nParticlesCurr) {nParticlesCurr_ = nParticlesCurr;}
  /// Specifies the total number of particles to use.
  void setNParticlesPrev(const unsigned int nParticlesPrev) {nParticlesPrev_ = nParticlesPrev;}
  /// Specifies the maximum number of particles to use (same as nParticles for a standard particle filter).
  void setNParticlesMax(const unsigned int nParticlesMax) {nParticlesMax_ = nParticlesMax;}
  /// Specifies the number of regimes.
  void setNRegimes(const unsigned int nRegimes) {nRegimes_ = nRegimes;}
  
  /// Initialises various containers used by the SMC algorithm.
  /// and performs the first step of the algorithm.
  void initialise()
  {
    step_ = 0;
//     std::cout << "SMC step: " << step_ << std::endl;
    unsigned int N;
    unsigned int R;
    unsigned int M;
    
    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      setNParticlesCurr(getNRegimes());
      N = getNParticlesCurr();
      R = getNRegimes();
      M = N - R;
    }
    else
    {
      setNParticlesCurr(getNParticlesMax());
      N = getNParticlesCurr();
      R = N;
      M = N;
    }
    
//     std::cout << "N: " << N << ", " << getNParticlesCurr() <<  "; R: " << R << ", " << getNRegimes() << "; M: " << M << "; NMax: " << getNParticlesMax() << std::endl;
    
    particlesPrev_.resize(N);
    particlesCurr_.resize(N);
    logUnnormalisedWeightsPrev_.set_size(N);
    logUnnormalisedWeightsCurr_.set_size(N);
    selfNormalisedWeightsPrev_.set_size(N);
    selfNormalisedWeightsCurr_.set_size(N);
    ancestorIndicesPrev_.set_size(M);
    ancestorIndicesCurr_.set_size(M);
    backwardKernels_.resize(R);
    for (unsigned int r = 0; r < R; r++)
    {
      backwardKernels_[r].set_size(N);
    }
    
//     std::cout << "start: sampling initial particles" << std::endl; 
    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      sampleInitialParticlesCp();
      computeInitialWeightsCp();
    }
    else
    {
      sampleInitialParticles();
      computeInitialWeights();
    }
    selfNormaliseWeights();  
    
//     std::cout << "finished: sampling initial particles" << std::endl;
    
    if (logUnnormalisedWeightsCurr_.has_nan()) 
    {
      std::cout << "log-unnormalised weights at time " << step_ << ": " << logUnnormalisedWeightsCurr_.t() << std::endl;
      
      for (unsigned int n = 0; n < getNParticlesCurr(); n++)
      {
        std::cout << "n: " << n << "; dCurr=" << particlesCurr_[n].getDistance() << "; rCurr=" << particlesCurr_[n].getRegime() << "; a: " << ancestorIndicesCurr_[n] << " ";
      }
      std::cout << " " << std::endl;
      for (unsigned int n = 0; n < getNParticlesPrev(); n++)
      {
        std::cout << "n: " << n << "; dPrev=" << particlesPrev_[n].getDistance() << "; rPrev=" << particlesPrev_[n].getRegime() << " ";
      }
      std::cout << " " << std::endl;
    }
    if (selfNormalisedWeightsCurr_.has_nan()) 
    {
      std::cout << "self-normalised weights at time " << step_ << ": " << selfNormalisedWeightsCurr_.t() << std::endl;
    }
//     std::cout << "self-normalised weights at time " << step_ << ": " << selfNormalisedWeightsCurr_.t() << std::endl;
  }
  /// Performs the $t+1$th iteration of the SMC algorithm.
  void iterate()
  {
    step_++;
//     std::cout << "SMC step: " << step_ << std::endl;
    
    setNParticlesPrev(getNParticlesCurr());
    
    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      if (getNParticlesPrev() + getNRegimes() > getNParticlesMax()) {
        setNParticlesCurr(getNParticlesMax());
      }
      else 
      {
        setNParticlesCurr(getNParticlesPrev() + getNRegimes());
      }
    }
        
        // TODO: put the swap()s back in to improve computational efficiency
    
//     if (particlesPrev_.size() == getNParticlesCurr()) 
//     {
//       particlesPrev_.swap(particlesCurr_);
//       logUnnormalisedWeightsPrev_.swap(logUnnormalisedWeightsCurr_);
//       selfNormalisedWeightsPrev_.swap(selfNormalisedWeightsCurr_);
//       ancestorIndicesPrev_.swap(ancestorIndicesCurr_); 
//     }
//     else
//     {
      particlesPrev_ = particlesCurr_;
      logUnnormalisedWeightsPrev_ = logUnnormalisedWeightsCurr_; 
      selfNormalisedWeightsPrev_ = selfNormalisedWeightsCurr_; 
      ancestorIndicesPrev_ = ancestorIndicesCurr_; 
      
      particlesCurr_.resize(getNParticlesCurr());
      logUnnormalisedWeightsCurr_.resize(getNParticlesCurr());
      selfNormalisedWeightsCurr_.resize(getNParticlesCurr());
      ancestorIndicesCurr_.resize(getNParticlesCurr() - getNRegimes());
//     }
    logSumOfUnnormalisedWeightsPrev_ = logSumOfUnnormalisedWeightsCurr_;
    
    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {

      
      // Generates the ancestor indices, i.e. resamples the particles
      
//       std::cout << "selfNnormalised weights at time " << step_ << ": " << selfNormalisedWeightsCurr_.t() << std::endl;
      
//       std::cout << "log-unnormalised weights at time " << step_ << ": " << logUnnormalisedWeightsCurr_.t() << std::endl;
     
//       std::cout << "resampleCp()" << std::endl;
//       std::cout << "getNParticlesCurr(): " << getNParticlesCurr() << "; getNParticlesPrev(): " << getNParticlesPrev() << std::endl;
      resampleCp(); 
      // Generates the particles at the current time step
//         std::cout << "sampleParticlesCp()" << std::endl;
      sampleParticlesCp();
      // Computes the particle weights at the current time step
//         std::cout << "computeWeightsCp()" << std::endl;
      computeWeightsCp();
//         std::cout << "finished: computeWeightsCp()" << std::endl;
    }
    else
    {

      // Generates the ancestor indices, i.e. resamples the particles
      resample();
      // Generates the particles at the current time step
      sampleParticles();
      // Computes the particle weights at the current time step
      computeWeights();
    }
    selfNormaliseWeights();
    
    ///////////////////////
  if (logUnnormalisedWeightsCurr_.has_nan()) 
    {
      std::cout << "log-unnormalised weights at time " << step_ << ": " << logUnnormalisedWeightsCurr_.t() << std::endl;
      
      for (unsigned int n = 0; n < nParticlesCurr_; n++)
      {
        std::cout << "n: " << n << "; dCurr=" << particlesCurr_[n].getDistance() << "; rCurr=" << particlesCurr_[n].getRegime() << " ";
        if (n < ancestorIndicesCurr_.size())
        {
          std::cout << "; a: " << ancestorIndicesCurr_[n] << " ";
        }
        
      }
      std::cout << " " << std::endl;
      for (unsigned int n = 0; n < nParticlesPrev_; n++)
      {
        std::cout << "n: " << n << "; dPrev=" << particlesPrev_[n].getDistance() << "; rPrev=" << particlesPrev_[n].getRegime() << " ";
      }
      std::cout << " " << std::endl;
    }
    /////////////////////////
  }
  /// Evaluates the backward kernel at the most recent time step.
  void evaluateBackwardKernels()
  {
    unsigned int N = particlesCurr_.size();

    if (getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL)
    {
      unsigned int R = getNRegimes();
      unsigned int M = N - R;
      for (unsigned int r = 0; r < R; r++)
      {
        backwardKernels_[r].set_size(particlesPrev_.size());
        for (unsigned int n = 0; n < particlesPrev_.size(); n++)
        {
          backwardKernels_[r](n) = logUnnormalisedWeightsPrev_(n) + model_.evaluateLogTransitionDensity(step_, particlesCurr_[M + r].getLatentVariable(), particlesPrev_[n].getLatentVariable());
        }
        
        backwardKernels_[r] = arma::exp(normaliseExp(backwardKernels_[r]));
        if (backwardKernels_[r].has_nan())
        {
          backwardKernels_[r].zeros();
        }
      }
    }
    else
    {
      for (unsigned int n = 0; n < N; n++)
      {
        for (unsigned int m = 0; m < N; m++)
        {
          backwardKernels_[n](m) = logUnnormalisedWeightsPrev_(m) + model_.evaluateLogTransitionDensity(step_, particlesCurr_[n].getLatentVariable(), particlesPrev_[m].getLatentVariable());
        }
        backwardKernels_[n] = arma::exp(normaliseExp(backwardKernels_[n]));
        if (backwardKernels_[n].has_nan())
        {
          backwardKernels_[n].zeros();
        }
      }
    }
  }
  /// Computes a weighted average of a collection of column vectors (each of the same size);
  /// the weights are the current filter weights.
  double computeFilteredMean(const std::vector<double>& x)
  {
    double est = 0.0;
//     for (unsigned int n = 0; n < x.size(); n++)
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      est = est + selfNormalisedWeightsCurr_(n) * x[n];
    }
    return est;
  }
  /// Overload for arma::colvec.
  arma::colvec computeFilteredMean(const std::vector<arma::colvec>& x)
  {
    arma::colvec est(x[0].size(), arma::fill::zeros);
//     for (unsigned int n = 0; n < x.size(); n++)
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      est = est + selfNormalisedWeightsCurr_(n) * x[n];
    }
    return est;
  }
  /// Computes an empirical variance;
  /// the weights are the current filter weights.
  double computeFilteredVariance(const std::vector<double>& x)
  {
    double est = 0.0;
    double meanEst = computeFilteredMean(x);
//     for (unsigned int n = 0; n < x.size(); n++)
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      est += selfNormalisedWeightsCurr_(n) * std::pow(x[n] - meanEst, 2.0);
    }
    return est;
  }

    
  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user
  
  /// Specifies the known proposal parameters (i.e. those which we do not infer).
  void setKnownParameters(const arma::colvec& varphi);
  /// Specifies the (proposal) parameters which are to be inferred.
  void setUnknownParameters(const arma::colvec& phi);
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
private:

  //////////////////////////////////////////////////////////////////////////////////////////////////// NOTE: these are model-specific and need to be implemented by the user
  
  /// Evaluates the logarithm of the proposal density at the first time step.
  double evaluateLogInitialProposalDensity(const Particle& particleCurr);
  /// Evaluates the logarithm of the proposal density at the current time step.
  double evaluateLogProposalDensity(const unsigned int step, const Particle& particleCurr, const Particle& particlePrev);
  /// Generates a particle from the proposal kernel at time $t=0$.
  void sampleFromInitialProposalDistribution(Particle& particleCurr);
  /// Generates a particle from the proposal kernel at some time $t>0$.
  void sampleFromProposalDistribution(const unsigned int step, Particle& particleCurr, const Particle& particlePrev);
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Resamples the particles.
  /// Only used for non-change-point models.
  void resample() // TODO: add option for ESS-based adaptive resampling
  {
    if (smcResampleType_ == SMC_RESAMPLE_SYSTEMATIC)
    {
      resample::systematic(ancestorIndicesCurr_, logUnnormalisedWeightsCurr_, selfNormalisedWeightsPrev_, nParticlesCurr_, logUnnormalisedWeightsPrev_, logSumOfUnnormalisedWeightsPrev_);
    }
    else if (smcResampleType_ == SMC_RESAMPLE_MULTINOMIAL) // i.e. multinomial resampling
    {
      resample::multinomial(ancestorIndicesCurr_, logUnnormalisedWeightsCurr_, selfNormalisedWeightsPrev_, nParticlesCurr_, logUnnormalisedWeightsPrev_, logSumOfUnnormalisedWeightsPrev_);
    }
    else
    {
      std::cout << "ERROR: only non-adaptive systematic and multinomial resampling are currently supported for non-change point models!" << std::endl;
    }
  }
  /// Resamples the particles.
  /// Only used for change-point models (otherwise, resample() is used).
  void resampleCp()
  {
    unsigned int M = getNParticlesCurr() - getNRegimes();
    
    if (getNParticlesCurr() < getNParticlesPrev() + getNRegimes())
    {
      // Compute indices of particles with non-zero weights.
      arma::uvec indices = arma::find_finite(logUnnormalisedWeightsPrev_);
      if (indices.size() > M)
      {
        if (smcResampleType_ == SMC_RESAMPLE_OPTIMAL_FINITE_STATE)
        {
          resample::optimalFiniteState(ancestorIndicesCurr_, logUnnormalisedWeightsCurr_, selfNormalisedWeightsPrev_, M, logUnnormalisedWeightsPrev_, logSumOfUnnormalisedWeightsPrev_);
        }
        else if (smcResampleType_ == SMC_RESAMPLE_SYSTEMATIC || smcResampleType_ == SMC_RESAMPLE_MULTINOMIAL)
        {
          if (smcResampleType_ == SMC_RESAMPLE_SYSTEMATIC)
          {
            resample::systematic(ancestorIndicesCurr_, logUnnormalisedWeightsCurr_, selfNormalisedWeightsPrev_, M, logUnnormalisedWeightsPrev_, logSumOfUnnormalisedWeightsPrev_);
          }
          else if (smcResampleType_ == SMC_RESAMPLE_MULTINOMIAL) // i.e. multinomial resampling
          {
            resample::multinomial(ancestorIndicesCurr_, logUnnormalisedWeightsCurr_, selfNormalisedWeightsPrev_, M, logUnnormalisedWeightsPrev_, logSumOfUnnormalisedWeightsPrev_);
          }
        }
      }
      else 
      {
        // TODO: make this more efficient by avoiding the sorting operation.
        
        // Simply keep the particles with the largest weights.
        // By construction, this includes all particles with non-zero weights.
        indices = arma::sort_index(logUnnormalisedWeightsPrev_, "descend");
        ancestorIndicesCurr_ = indices(arma::span(0, M - 1));
        logUnnormalisedWeightsCurr_(arma::span(0, M - 1)) = logUnnormalisedWeightsPrev_.elem(indices(arma::span(0, M - 1)));
      }
    }
    else
    {
      ancestorIndicesCurr_ = arma::linspace<arma::uvec>(0, nParticlesPrev_ - 1, nParticlesPrev_);
      logUnnormalisedWeightsCurr_(arma::span(0, nParticlesPrev_ - 1)) = logUnnormalisedWeightsPrev_;
    }
    

  }
  /// Samples particles at the initial time step.
  /// Used for any non-change point model.
  void sampleInitialParticles()
  {
//     std::cout << "no of particles at the initial time step in the standard PF: " << getNParticlesCurr() << std::endl;
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      sampleFromInitialProposalDistribution(n);
    }
  }
  /// Samples particles at the initial time step.
  /// Used for the change point model.
  void sampleInitialParticlesCp()
  {
//     std::cout << "no of particles at the initial time step in the change-point PF algorithm: " << getNParticlesCurr() << std::endl;
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      particlesCurr_[n].setDistance(1);
      particlesCurr_[n].setRegime(n);
    }
  }
  /// Evaluates the particle weights at the initial time step.
  void computeInitialWeights()
  {
    logUnnormalisedWeightsCurr_.fill(- std::log(getNParticlesCurr()));
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      logUnnormalisedWeightsCurr_(n) += 
        evaluateLogInitialVarGamma(n) -
        evaluateLogInitialProposalDensity(n);
    } 
  }
  /// Evaluate the particle weights at the initial time step.
  /// Used for the change point model.
  void computeInitialWeightsCp()
  {
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      logUnnormalisedWeightsCurr_(n) = evaluateLogInitialVarGamma(n);
    }
  }
  /// Samples particles at the current time step.
  /// Used for any non-change point model.
  void sampleParticles()
  {
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      sampleFromProposalDistribution(n, ancestorIndicesCurr_[n]);
    }
  }
  /// Samples particles at the current time step.
  /// Only used for change-point models (otherwise, this 
  /// sample() needs to be implemented by the user.
  void sampleParticlesCp() 
  {
    unsigned int N = getNParticlesCurr();
    unsigned int R = getNRegimes();
    unsigned int M = N - R;

    // Generating the first $N$ particles
    for (unsigned int n = 0; n < M; n++)
    {
      particlesCurr_[n].setDistance(particlesPrev_[ancestorIndicesCurr_[n]].getDistance() + 1);
      particlesCurr_[n].setRegime(particlesPrev_[ancestorIndicesCurr_[n]].getRegime());
    }
    // Determining the remaining $R$ particles
    for (unsigned int r = 0; r < R; r++)
    {
      particlesCurr_[M + r].setDistance(1);
      particlesCurr_[M + r].setRegime(r); // NOTE: regimes are counted from $0$ to $R-1$
    }
  }
  /// Computes the particle weights at the current time step.
  /// Needs to be implemented by the user (for non-change point models);
  void computeWeights()
  {
    for (unsigned int n = 0; n < getNParticlesCurr(); n++)
    {
      logUnnormalisedWeightsCurr_(n) +=
        evaluateLogVarGamma(n, ancestorIndicesCurr_[n]) -
        evaluateLogProposalDensity(n, ancestorIndicesCurr_[n]);
    }
  }
  /// Computes the particle weights at the current time step.
  /// Only used for change-point models (otherwise, computeWeights() is used)
  void computeWeightsCp() 
  {
//     std::cout << "NParticlesCurr: " << getNParticlesCurr() << std::endl;
//     std::cout << "NParticlesPrev: " << getNParticlesPrev() << std::endl;
//     std::cout << "logUnnormalisedWeightsCurr_.size(): " << logUnnormalisedWeightsCurr_.size() << std::endl;
//     std::cout << "logUnnormalisedWeightsPrev_.size(): " << logUnnormalisedWeightsPrev_.size() << std::endl;
//     std::cout << "particlesCurr_.size(): " << particlesCurr_.size() << std::endl;
//     std::cout << "particlesPrev_.size(): " << particlesPrev_.size() << std::endl;
    
//     std::cout << "logUnnormalisedWeightsPrev_.t(): " <<  logUnnormalisedWeightsPrev_.t() << std::endl;
    unsigned int N = getNParticlesCurr();
    unsigned int R = getNRegimes();
    unsigned int M = N - R;
    
    // Computing the first $M$ weights:
//     logUnnormalisedWeightsCurr_.set_size(N);
    for (unsigned int n = 0; n < M; n++)
    {
      logUnnormalisedWeightsCurr_(n) += 
        evaluateLogVarGamma(n, ancestorIndicesCurr_[n]);
    }
//     std::cout << " " << std::endl;
//     std::cout << "arma::trans(logUnnormalisedWeightsCurr_(arma::span(0, M-1))): " <<  arma::trans(logUnnormalisedWeightsCurr_(arma::span(0, M - 1))) << std::endl;
        
        
    // Computing the remaining $R$ weights:
    arma::colvec logAux(logUnnormalisedWeightsPrev_.size());
    for (unsigned int r = 0; r < R; r++)
    {
      for (unsigned int n = 0; n < logUnnormalisedWeightsPrev_.size(); n++)
      {
        logAux(n) = evaluateLogVarGamma(M + r, n) + logUnnormalisedWeightsPrev_[n];
//         std::cout << "1st component of logAux(n) for n = " << n << ": " << evaluateLogVarGamma(M + r, n) << std::endl;
      }
//       std::cout << " " << std::endl;
//       std::cout << "logAux: " << logAux.t() << std::endl;
      logUnnormalisedWeightsCurr_(M + r) = sumExp(logAux); // summation in log-space for numerical stability
    }
  }
  /// Self-normalises the weights.
  void selfNormaliseWeights()
  {
    selfNormalisedWeightsCurr_ = arma::exp(normaliseExp(logUnnormalisedWeightsCurr_, logSumOfUnnormalisedWeightsCurr_));
  }
  /// Evaluates the logarithm of the unnormalised target
  /// density at the initial time step.
  double evaluateLogInitialVarGamma(const unsigned int n)
  {
    return model_.evaluateLogInitialDensity(particlesCurr_[n].getLatentVariable()) + 
    model_.evaluateLogObservationDensity(0, particlesCurr_[n].getLatentVariable());
  }
  /// Evaluates the logarithm of the ratio between the unnormalised 
  /// unnormalised target density between the current and previous time step.
  double evaluateLogVarGamma(const unsigned int n, const unsigned int m)
  {
//     std::cout << "evaluateLogTransitionDensity: " << model_.evaluateLogTransitionDensity(step_, particlesCurr_[n].getLatentVariable(), particlesPrev_[m].getLatentVariable()) << "; evaluateLogObservationDensity: " << model_.evaluateLogObservationDensity(step_, particlesCurr_[n].getLatentVariable()) << " ";
    
    return 
    model_.evaluateLogTransitionDensity(step_, particlesCurr_[n].getLatentVariable(), particlesPrev_[m].getLatentVariable()) + 
    model_.evaluateLogObservationDensity(step_, particlesCurr_[n].getLatentVariable());
  }
  /// Evaluates the logarithm of the proposal
  /// density at the initial time step.
  double evaluateLogInitialProposalDensity(const unsigned int n)
  {
    return evaluateLogInitialProposalDensity(particlesCurr_[n]);
  }
  /// Evaluates the logarithm of the proposal density at some time $t > 0$.
  double evaluateLogProposalDensity(const unsigned int n, const unsigned int m)
  {
    return evaluateLogProposalDensity(step_, particlesCurr_[n], particlesPrev_[m]);
  }
  /// Generates a particle from the proposal kernel at time $t=0$.
  void sampleFromInitialProposalDistribution(const unsigned int n)
  {
    sampleFromInitialProposalDistribution(particlesCurr_[n]);
  }
  /// Generates a particle from the proposal kernel at some time $t>0$.
  void sampleFromProposalDistribution(const unsigned int n, const unsigned int m)
  {
    sampleFromProposalDistribution(step_, particlesCurr_[n], particlesPrev_[m]);
  }
  
  Rng& rng_; // random number generation.
  Model<ModelParameters, LatentVariable, Covariate, Observation>& model_; // the targeted model.
  SmcParameters smcParameters_; // parameters associated with the importance-sampling scheme
  SmcResampleType smcResampleType_; // type of resampling scheme used
  SmcProposalType smcProposalType_; // type of proposal kernel used
  std::vector<Particle> particlesCurr_, particlesPrev_; // stores the particles
  arma::uvec ancestorIndicesCurr_, ancestorIndicesPrev_; // stores the ancestor indices (if resampling/subsampling is used)
  arma::colvec logUnnormalisedWeightsCurr_, logUnnormalisedWeightsPrev_; // stores the log-unnormalised particle weights
  arma::colvec selfNormalisedWeightsCurr_, selfNormalisedWeightsPrev_; // stores the self-normalised particle weights
  double logSumOfUnnormalisedWeightsCurr_, logSumOfUnnormalisedWeightsPrev_; // stores the log-estimate of the normalising constant
  std::vector<arma::colvec> backwardKernels_; // for the change-point model algorithm the $r$th element of this vector contains the backward kernel associated with the (nParticles_-nRegimes+r)th current particle; otherwise the $n$th element contains the $n$th backward kernel 
//   unsigned int nParticles_; // total number of particles 
  unsigned int nParticlesCurr_, nParticlesPrev_; // total number of particles at current and previous time steps
  unsigned int nParticlesMax_; // maximum number of particles to use (same as nParticlesCurr_ = nParticlesPrev_ for a standard particle filter)
  unsigned int step_; // the current step time index
  
  
  /// NOTE: the following members are only used if getSmcProposalType() == SMC_PROPOSAL_CHANGE_POINT_MODEL
  unsigned int nRegimes_; // number of regimes in the case of a change point model
  
};

#endif
