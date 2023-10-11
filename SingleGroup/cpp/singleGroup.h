/// \file
/// \brief A change-point model for DNA methylation measurements.
///
/// Implements a change-point model for DNA methylation 
/// measurements for data from a single group (i.e. without
/// a control group).


#ifndef __SINGLEGROUP_H
#define __SINGLEGROUP_H

#include <cstdlib> // provides abs()
#include "algorithms/OnlineCombinedInference.h"


/// Type of methylation regime.
enum MethylationRegimeType
{ 
  METHYLATION_REGIME_LARGE_MEAN_SMALL_SD = 0,
  METHYLATION_REGIME_SMALL_MEAN_SMALL_SD,
  METHYLATION_REGIME_LARGE_MEAN_LARGE_SD,
  METHYLATION_REGIME_SMALL_MEAN_LARGE_SD,
  METHYLATION_REGIME_MEDIUM_MEAN_LARGE_SD,
  METHYLATION_REGIME_CHAOTIC
};

///////////////////////////////////////////////////////////////////////////////
// Implementation of the Model class template
///////////////////////////////////////////////////////////////////////////////

/// Holds the latent variables associated with a single "time" step.
class LatentVariable
{
public:
  
  /// Returns the index of the current regime.
  MethylationRegimeType getMethylationRegimeType() const {return methylationRegimeType_;}
  /// Returns the sojourn time for the current regime.
  unsigned int getSojournTime() const {return sojournTime_;}
  /// Specifies the index of the current regime.
  void setMethylationRegimeType(const MethylationRegimeType methylationRegimeType) {methylationRegimeType_ = methylationRegimeType;}
  /// Specifies the sojourn time for the current regime.
  void setSojournTime(const unsigned int sojournTime) {sojournTime_ = sojournTime;}
  
private:
  
  MethylationRegimeType methylationRegimeType_; // the index of the current regime
  unsigned int sojournTime_; // the sojourn time in the current regime

};

/// Holds the covariates associated with a single "time" step.
class Covariate
{
public:
  
  /// Returns the genomic position of the CpG site.
  unsigned int getGenomicPosition() const {return genomicPosition_;}
  /// Returns the number of reads taken at this CpG site.
  const arma::uvec& getNTotalReads() const {return nTotalReads_;}
  /// Specifies the genomic position of the CpG site.
  void setGenomicPosition(const unsigned int genomicPosition) {genomicPosition_ = genomicPosition;}
  /// Specifies the number of reads taken at this CpG site.
  void setNTotalReads(const arma::uvec& nTotalReads) {nTotalReads_ = nTotalReads;}
  
private:
  
  unsigned int genomicPosition_; // the genomic position of the CpG site
  arma::uvec nTotalReads_; // the number of reads (methylated + unmethylated) taken at this CpG site

};

/// Holds the observations associated with a single "time" step.
/// Here: each component of the vector represents the number of 
/// methylated reads in a single sample; the length of this vector 
/// thus constitutes the number of samples.
typedef arma::uvec Observation;

/// Holds all the model parameters.
class ModelParameters
{
public:
  
  /// Returns the minimum distance between change points.
  const double getU() const {return u_;}
  /// Returns the transition matrix for the regimes.
  const arma::mat& getP() const {return P_;}
  /// Returns the transpose of the transition matrix for the regimes.
  const arma::mat& getPTrans() const {return PTrans_;}
  /// Returns the first set of beta-binomial distribution parameters.
  const arma::colvec& getAlpha() const {return alpha_;}
  /// Returns the second set of beta-binomial distribution parameters.
  const arma::colvec& getBeta() const {return beta_;}
  /// Returns the $r$th entry of the first set of 
  /// beta-binomial distribution parameters.
  double getAlpha(const unsigned int r) const {return alpha_(r);}
  /// Returns the $r$th entry of the second set of
  /// beta-binomial distribution parameters.
  double getBeta(const unsigned int r) const {return beta_(r);}
  /// Returns the first set of negative-binomial distribution parameters.
  const arma::colvec& getOmega() const {return omega_;}
  /// Returns the second set of negative-binomial distribution parameters.
  const arma::colvec& getKappa() const {return kappa_;}
  /// Returns the $r$th entry of the first set of 
  /// negative-binomial distribution parameters.
  double getOmega(const unsigned int r) const {return omega_(r);}
  /// Returns the $r$th entry of the second set of
  /// negative-binomial distribution parameters.
  double getKappa(const unsigned int r) const {return kappa_(r);}
  /// Returns whether the regime-specific negative-binomial parameters
  /// $\kappa_r$ are fixed (otherwise, they are estimated)
  bool getIsKappaFixed() const {return isKappaFixed_;}
  /// Returns the number of different methylation regimes.
  unsigned int getNMethylationRegimes() const {return nMethylationRegimes_;}
  /// Returns indices of components of $\theta$.
  /// NOTE: as usual, we always start counting from zero,
  /// i.e. $r=0$ represents the first "regime".
  unsigned int getIndexPMin(const unsigned int r) const {return r * (nMethylationRegimes_ - 1);}
  unsigned int getIndexPMax(const unsigned int r) const {return (r + 1) * (nMethylationRegimes_ - 1) - 1;}
  unsigned int getIndexPMin() const {return getIndexPMin(0);}
  unsigned int getIndexPMax() const {return getIndexPMax(nMethylationRegimes_ - 1);}
  unsigned int getIndexOmega(const unsigned int r) const {return nMethylationRegimes_ * (nMethylationRegimes_ - 1) + r;}
  unsigned int getIndexOmegaMin() const {return nMethylationRegimes_ * (nMethylationRegimes_ - 1);}
  unsigned int getIndexOmegaMax() const {return nMethylationRegimes_ * nMethylationRegimes_ - 1;}
  unsigned int getIndexKappa(const unsigned int r) const {return nMethylationRegimes_ * nMethylationRegimes_ + r;}
  unsigned int getIndexKappaMin() const {return nMethylationRegimes_ * nMethylationRegimes_;}
  unsigned int getIndexKappaMax() const {return nMethylationRegimes_ * (nMethylationRegimes_ + 1) - 1;}
  /// Returns the pre-calculated values of rho.
  double getRho(const unsigned int dPrev, const unsigned int rPrev)
  {
    if (dPrev > dMax_[rPrev])
    {
      extendAuxiliaryQuantities(dPrev, rPrev);
    }
    return rho_[rPrev][dPrev-1];
  }
  /// Returns the pre-calculated values of the derivative of $\log \rho$ w.r.t. the component
  /// of $\theta$ which corresponds to the $\omega$-parameter for the $r$th regime.
  double getGradThetaOmegaLogRho(const unsigned int dPrev, const unsigned int rPrev) 
  {
    if (dPrev > dMax_[rPrev])
    {
      extendAuxiliaryQuantities(dPrev, rPrev);
    }
    return gradThetaOmegaLogRho_[rPrev][dPrev-1];
  }
  /// Returns the pre-calculated values of the derivative of $\log \rho$ w.r.t. the component
  /// of $\theta$ which corresponds to the $\kappa$-parameter for the $r$th regime.
  double getGradThetaKappaLogRho(const unsigned int dPrev, const unsigned int rPrev)
  {
    if (dPrev > dMax_[rPrev])
    {
      extendAuxiliaryQuantities(dPrev, rPrev);
    }
    return gradThetaKappaLogRho_[rPrev][dPrev-1];
  }
  /// Returns the exit status which indicates numerical problems in the calculation of rho.
  double getExitStatus(const unsigned int dPrev, const unsigned int rPrev) const
  {
    return exitStatus_[rPrev][dPrev-1];
  }
  
  
//   /// Returns the number of samples at each time point.
//   /// NOTE: this is only used for generating data from the model.
//   unsigned int getNSamples() const {return nSamples_;}
//   /// Returns the average number of reads at each time point.
//   /// NOTE: this is only used for generating data from the model.
//   unsigned int getMeanNTotalReads() const {return meanNTotalReads_;}
  
  /// Specifies whether the regime-specific negative-binomial parameters
  /// $\kappa_r$ are fixed (otherwise, they are estimated)
  void setIsKappaFixed(const bool isKappaFixed) {isKappaFixed_ = isKappaFixed;}
  
  /*
  /// Specifies the number of samples at each time point.
  /// NOTE: this is only used for generating data from the model.
  void setNSamples(const unsigned int nSamples) {nSamples_ = nSamples;}
  /// Specifies the average number of reads at each time point.
  /// NOTE: this is only used for generating data from the model.
  void setMeanNTotalReads(const unsigned int meanNTotalReads) {meanNTotalReads_ = meanNTotalReads;}
  */
  /// Specifies the known/fixed model parameters.
  void setKnownParameters(const arma::colvec& vartheta)
  {
    /// vartheta = (u_, nMethylationRegimes_ alpha_, beta_, isKappaFixed_, kappa_) if isKappaFixed_ == false
    /// vartheta = (u_, nMethylationRegimes_ alpha_, beta_, isKappaFixed_) if isKappaFixed_ == true
    u_ = static_cast<unsigned int>(vartheta(0));
    nMethylationRegimes_ = static_cast<unsigned int>(vartheta(1));
    unsigned int R = nMethylationRegimes_;
    alpha_ = vartheta(arma::span(2,R+1));
    beta_ = vartheta(arma::span(R+2,2*R+1));
//     setNSamples(static_cast<unsigned int>(vartheta(2*R+2)));
//     setMeanNTotalReads(vartheta(2*R+3));
    setIsKappaFixed(static_cast<bool>(vartheta(2*R+2)));
    if (isKappaFixed_)
    {
      kappa_ = vartheta(arma::span(2*R+3,3*R+2));
      dimTheta_ = R * R;
    }
    else
    {
      dimTheta_ = R * (R + 1);
    }
    P_.set_size(R,R);
    PTrans_.set_size(R,R);
    omega_.set_size(R);
  }
  /// Specifies the model parameters which are inferred by the algorithm.
  void setUnknownParameters(const arma::colvec& theta)
  {
    unsigned int R = getNMethylationRegimes();
    
    PTrans_.zeros();
    unsigned int idxMin, idxMax;
    arma::colvec aux(R-1);
    for (unsigned int r=0; r<R; r++) 
    {
      idxMin = getIndexPMin(r);
      idxMax = getIndexPMax(r);
      aux = arma::exp(normaliseExp(theta(arma::span(idxMin,idxMax))));
      aux.insert_rows(r,1);
      PTrans_.col(r) = aux;
    }
    P_ = PTrans_.t();
    
    idxMin = getIndexOmegaMin();
    idxMax = getIndexOmegaMax();
    omega_ = inverseLogit(theta(arma::span(idxMin,idxMax)));
    
    if (!isKappaFixed_)
    {
      idxMin = getIndexKappaMin();
      idxMax = getIndexKappaMax();
      kappa_ = arma::exp(theta(arma::span(idxMin,idxMax)));
    }
    
//     std::cout << "P: " << P_ << std::endl;
//     std::cout << "omega: " << omega_.t() << std::endl;


//     std::cout << "initialising auxiliary quantities"<< std::endl;
 
    /// Calculates some quantities to avoid re-calculating them for each particle.
    dMax_.resize(R);
    littleH_.resize(R);
    bigH_.resize(R);
    rho_.resize(R);
    exitStatus_.resize(R);
    gradThetaOmegaLogLittleH_.resize(R);
    gradThetaOmegaBigH_.resize(R);
    gradThetaOmegaLogRho_.resize(R);
    if (!isKappaFixed_)
    {
      gradThetaKappaLogLittleH_.resize(R);
      gradThetaKappaBigH_.resize(R);
      gradThetaKappaLogRho_.resize(R);
    }
    
    for (unsigned int r=0; r<R; r++)
    {
      dMax_[r] = 10; 
      littleH_[r].reserve(0);
      bigH_[r].reserve(0);
      rho_[r].reserve(0);
      exitStatus_[r].reserve(0);
      gradThetaOmegaLogLittleH_[r].reserve(0);
      gradThetaOmegaBigH_[r].reserve(0);
      gradThetaOmegaLogRho_[r].reserve(0);
      if (!isKappaFixed_)
      {
        gradThetaKappaLogLittleH_[r].reserve(0);
        gradThetaKappaBigH_[r].reserve(0);
        gradThetaKappaLogRho_[r].reserve(0);
      }
    }
    
//         std::cout << "first extension of auxiliary quantities"<< std::endl;
    for (unsigned int r=0; r<R; r++)
    {
      extendAuxiliaryQuantities(dMax_[r], r);
    }
  }
  /// Calculates additional elements for the auxiliary quantities.
  void extendAuxiliaryQuantities(const unsigned int dMaxNew, const unsigned int r)
  {
//             std::cout << "START: extension of auxiliary quantities"<< std::endl;
    unsigned int dMaxOld = littleH_[r].size();
    
    littleH_[r].reserve(dMaxNew);
    bigH_[r].reserve(dMaxNew);
    rho_[r].reserve(dMaxNew);
    exitStatus_[r].reserve(dMaxNew);
    gradThetaOmegaLogLittleH_[r].reserve(dMaxNew);
    gradThetaOmegaBigH_[r].reserve(dMaxNew);
    gradThetaOmegaLogRho_[r].reserve(dMaxNew);
    if (!isKappaFixed_)
    {
      gradThetaKappaLogLittleH_[r].reserve(dMaxNew);
      gradThetaKappaBigH_[r].reserve(dMaxNew);
      gradThetaKappaLogRho_[r].reserve(dMaxNew);
    }
    
    for (unsigned int d=dMaxOld; d<std::min(u_-1, dMaxNew); d++)
    {
      littleH_[r][d] = 0.0;
      bigH_[r][d] = 0.0;
      rho_[r][d] = 0.0;
      exitStatus_[r][d] = false;
      gradThetaOmegaLogLittleH_[r][d] = 0.0;
      gradThetaOmegaBigH_[r][d] = 0.0;
      gradThetaOmegaLogRho_[r][d] = 0.0;
      if (!isKappaFixed_)
      {
        gradThetaKappaLogLittleH_[r][d] = 0.0;
        gradThetaKappaBigH_[r][d] = 0.0;
        gradThetaKappaLogRho_[r][d] = 0.0;
      }
    }
    for (unsigned int d=u_-1; d<dMaxNew; d++)
    {
      littleH_[r][d] = std::exp(evaluateLogNegativeBinomialDensity(d+1-u_, kappa_(r), omega_(r)));
      
      if (exitStatus_[r][d-1] == true || bigH_[r][d-1] >= 1.0)
      {
        bigH_[r][d-1] = 0.99999;
        rho_[r][d] = 1.0;
        exitStatus_[r][d] = true;
      }
      else
      {
        bigH_[r][d] = bigH_[r][d-1] + littleH_[r][d];
        rho_[r][d] = littleH_[r][d] / (1.0 - bigH_[r][d-1]);
        exitStatus_[r][d] = false;
      }
      
      gradThetaOmegaLogLittleH_[r][d] = (static_cast<double>(d+1-u_)/omega_(r) - kappa_(r) /(1.0 - omega_(r))) * gradLogitEvaluatedAtInverseLogit(omega_(r));
      gradThetaOmegaBigH_[r][d] = gradThetaOmegaBigH_[r][d-1] + littleH_[r][d] *  gradThetaOmegaLogLittleH_[r][d];
      gradThetaOmegaLogRho_[r][d] = gradThetaOmegaLogLittleH_[r][d] + gradThetaOmegaBigH_[r][d-1] / (1.0 - bigH_[r][d-1]);
      
      if (!isKappaFixed_)
      {
        gradThetaKappaLogLittleH_[r][d] = kappa_(r) * (digamma(static_cast<double>(d+1-u_) + kappa_(r)) - digamma(kappa_(r)) - std::log(1.0 - omega_(r)));
        gradThetaKappaBigH_[r][d] = gradThetaOmegaBigH_[r][d-1] + littleH_[r][d] *  gradThetaKappaLogLittleH_[r][d];
        gradThetaKappaLogRho_[r][d] = gradThetaKappaLogLittleH_[r][d] + gradThetaKappaBigH_[r][d-1] / (1.0 - bigH_[r][d-1]);
      }
    }
    dMax_[r] = dMaxNew;
    
//          std::cout << "END: extension of auxiliary quantities"<< std::endl;
  }

  /// Returns the probability of a change point 
  /// occuring at the current time step, i.e. this
  /// is evaluates the probability $\rho$.
  double evaluateChangePointProbability
  (
    const unsigned int t, // time index
    const unsigned int dPrev, // the sojourn time at the previous time step
    const unsigned int rPrev, // the regime at the previous time step
    bool exitStatus // is set to TRUE if the parameters are not in a (numerically) valid range 
  ) const
  {

    if (dPrev >= u_)
    {
      double bigH = 0.0;
      if (dPrev > u_)
      {
        for (unsigned int i=u_; i<dPrev; i++)
        {
          bigH += evaluateLittleH(i,rPrev);
        }
      }
      if (bigH < 1.0)
      {     
        exitStatus = false;
        return evaluateLittleH(dPrev,rPrev) / (1.0 - bigH);
      }
      else
      { 
        std::cout << "WARNING: bigH at time t: " << t << ": " << bigH << "; rPrev: " << rPrev << "; dPrev: " << dPrev <<std::endl;
        
        std::cout << "evaluateLittleH(dPrev-1,rPrev): " << evaluateLittleH(dPrev-1,rPrev) << "; omega_(rPrev): " << omega_(rPrev) << std::endl;
        exitStatus = true;
        return 1.0; /// TODO
      }
    }
    else
    {
      exitStatus = false;
      return 0.0;
    }
  }
  /// Evaluates the shifted negative-binomial density,
  /// i.e. the function $h$ used in the construction of the 
  /// change-point probabilities.
  double evaluateLittleH
  (
    const unsigned int d, // the sojourn time
    const unsigned int r // the regime
  ) const
  {
    if (d >= u_)
    { 
      return std::exp(evaluateLogNegativeBinomialDensity(d-u_, kappa_(r), omega_(r)));
    }
    else
    {
      return 0.0;
    }
  }
  /// Evaluates gradient of the logarithm of the shifted negative-binomial density,
  /// i.e. of the logarithm of the function $h$ used in the construction of the 
  /// change-point probabilities, w.r.t. to the parameter governing the
  /// "success probability" parameter.
  double evaluateGradThetaOmegaLogLittleH
  (
    const unsigned int d, // the sojourn time
    const unsigned int r
  ) const
  {
    if (d >= u_)
    { 
      return (static_cast<double>(d-u_)*(1.0 - omega_(r)) - kappa_(r)*omega_(r));
//       return (static_cast<double>(d-u_)/omega_(r) - kappa_(r) /(1.0 - omega_(r))) * gradLogitEvaluatedAtInverseLogit(omega_(r));
    }
    else
    {
      return 0.0;
    }
  }
  /// Evaluates gradient of the logarithm of the shifted negative-binomial density,
  /// i.e. of the logarithm of the function $h$ used in the construction of the 
  /// change-point probabilities, w.r.t. to the parameter governing the
  /// "number of failures until success" parameter.
  double evaluateGradThetaKappaLogLittleH
  (
    const unsigned int d, // the sojourn time
    const unsigned int r
  ) const
  {
    if (d >= u_)
    { 
      /// TODO: check this!
      return kappa_(r) * (digamma(static_cast<double>(d-u_) + kappa_(r)) - digamma(kappa_(r)) - std::log(1.0 - omega_(r)));
    }
    else
    {
      return 0.0;
    }
  }
  
  
private:
  
  unsigned int nMethylationRegimes_; // the number of regimes
  arma::colvec alpha_; // regime-specific values for the first parameters of the beta laws governing the observation densities (assumed known/fixed)
  arma::colvec beta_; // regime-specific values for the second parameters of the beta laws governing the observation densities (assumed known/fixed)
  arma::mat P_, PTrans_; // the transition matrix for the methylation regimes (has zeros on its diagonal)
  unsigned int u_; // shift-parameter for the negative-binomial distribution (must be $\geq 1$)
  arma::colvec omega_; // regime-specific values for the second parameters of the negative-binomial laws governing the distribution of the sojourn times
  bool isKappaFixed_; // are the regime-specific negative-binomial parameters kappa_ fixed (otherwise they are estimated)?
  arma::colvec kappa_; // regime-specific values for the first parameters of the negative-binomial laws governing the distribution of the sojourn times
  unsigned int dimTheta_; // length of the unknown-parameter vector theta
  
  // The following quantities are only needed in order to avoid repeated calculations:
  std::vector<std::vector<double>> littleH_, bigH_, rho_, gradThetaOmegaLogLittleH_, gradThetaOmegaBigH_, gradThetaOmegaLogRho_, gradThetaKappaLogLittleH_, gradThetaKappaBigH_, gradThetaKappaLogRho_;
  std::vector<std::vector<bool>> exitStatus_;
  std::vector<unsigned int> dMax_; // maximum (regime-specific) distances for which the above quantities are calculated
  
  
  // The following constants are only used for generating
  // synthetic data from the model:
//   unsigned int nSamples_; // the number of samples at each time point
//   double meanNTotalReads_; // the average number of reads at each time point (the numer of reads is then distributed according to a Poisson distribution with mean meanNTotalReads_).
};

/// Specifies the known hyperparameters.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
void Model<ModelParameters, LatentVariable, Covariate, Observation>::setKnownParameters(const arma::colvec& vartheta)
{
//   std::cout << "START: setKnownModelParameters()" << std::endl;
    
  modelParameters_.setKnownParameters(vartheta);
  unsigned int R = modelParameters_.getNMethylationRegimes();
  if (modelParameters_.getIsKappaFixed())
  {
    setDimTheta(R * R);
  }
  else
  {
    setDimTheta(R * (R + 1));
  }
//   std::cout << "END: setKnownModelParameters()" << std::endl;
}
/// Determines the model parameters from arma::colvec theta.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
void Model<ModelParameters, LatentVariable, Covariate, Observation>::setUnknownParameters(const arma::colvec& theta)
{
//   std::cout << "START: setUnknownModelParameters()" << std::endl;
  modelParameters_.setUnknownParameters(theta);
//   std::cout << "END: setUnknownModelParameters()" << std::endl;
}
/// Samples a single latent value from the prior on the transformed (generative-)model parameters.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
arma::colvec Model<ModelParameters, LatentVariable, Covariate, Observation>::sampleFromParameterPrior()
{
  return arma::randn<arma::colvec>(getDimTheta());
}
/// Samples a latent variable from the prior distribution at time $t=0$.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
void Model<ModelParameters, LatentVariable, Covariate, Observation>::sampleFromInitialDistribution
(
  LatentVariable& latentVariableCurr
)
{
  unsigned int R = modelParameters_.getNMethylationRegimes();
//   unsigned int dCurr = 1;
  unsigned int rCurr = rng_.randomUniformInt(0,R-1);
  
  latentVariableCurr.setSojournTime(1);
  latentVariableCurr.setMethylationRegimeType(static_cast<MethylationRegimeType>(rCurr));
//   std::cout << "proposed regime at time 0: " << static_cast<MethylationRegimeType>(rCurr) << std::endl;
}
/// Samples a latent variable from the transition equation at time $t>0$. 
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
void Model<ModelParameters, LatentVariable, Covariate, Observation>::sampleFromTransitionEquation
(
  const unsigned int t, 
  LatentVariable& latentVariableCurr, 
  const LatentVariable& latentVariablePrev
)
{
  unsigned int dPrev = latentVariablePrev.getSojournTime();
  unsigned int rPrev = static_cast<unsigned int>(latentVariablePrev.getMethylationRegimeType());
  unsigned int dCurr;
  unsigned int rCurr;
  arma::colvec prob;
  
//   std::cout << "START: sampleFromTransitionEquation" << std::endl;
//     std::cout << "rPrev: " << rPrev << "; dPrev: " << dPrev << std::endl;
  
//   std::cout << "START: computing rho" << std::endl;
    
//   bool exitStatus = false;
  double rho = modelParameters_.getRho(dPrev, rPrev);
//   double rho = modelParameters_.evaluateChangePointProbability(t,dPrev,rPrev,exitStatus); // probability of a change point
  
  /*
  if (std::abs(rho - modelParameters_.getRho(dPrev, rPrev)) > 0.000001)
  {
    std::cout << "rho difference in sampleFromTransitionEquation(): " << rho - modelParameters_.getRho(dPrev, rPrev) << std::endl;
  }
  */
  
  if (modelParameters_.getExitStatus(dPrev, rPrev))
  {
    std::cout << "WARNING: non-zero exit status of evaluateChangePointProbability() called by sampleFromTransitionEquation()" << std::endl; 
  }
  
//   std::cout << "rho: " << rho << std::endl;
  
  if (arma::randu() <= rho)
  {
    dCurr = 1;
    prob = modelParameters_.getPTrans().col(rPrev);
//     std::cout << "prob: " << prob.t() << std::endl;
    rCurr = sampleInt(prob);
  }
  else
  {
    dCurr = dPrev + 1;
    rCurr = rPrev;
  }
  latentVariableCurr.setSojournTime(dCurr);
  latentVariableCurr.setMethylationRegimeType(static_cast<MethylationRegimeType>(rCurr));
}
/// Samples a single observation from the observation equation at time $t>=0$.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
void Model<ModelParameters, LatentVariable, Covariate, Observation>::sampleFromObservationEquation
(
  const unsigned int t, 
  Observation& observation, 
  const LatentVariable& latentVariableCurr
)
{
//   std::cout << "latentVariableCurr.getSojournTime(): " << latentVariableCurr.getSojournTime() << std::endl;
//   std::cout << "latentVariableCurr.getMethylationRegimeType(): " << latentVariableCurr.getMethylationRegimeType() << std::endl;
  

  /*
  unsigned int S = 1; // make these accessible from the outside modelParameters_.getNSamples();
  double meanNTotalReads = 100.0; // TODO: make these accessible from the outside modelParameters_.getMeanNTotalReads();
  */
  
  unsigned int r = static_cast<unsigned int>(latentVariableCurr.getMethylationRegimeType());
  
//   std::cout << "r at time " << t << " :" << r << std::endl;
  double alpha = modelParameters_.getAlpha(r);
  double beta  = modelParameters_.getBeta(r);
  
//     std::cout << "alpha, beta: " << alpha << "; "  << beta << std::endl;
  
//   arma::uvec nTotalReads(S); 
  
  observation.set_size(covariates_[t].getNTotalReads().size()); 
//   arma::uvec nMethylatedReads(S);
  for (unsigned s=0; s<observation.size(); s++)
  {
//     nTotalReads(s)      = rng_.randomPoisson(meanNTotalReads);
    observation(s) = rng_.randomBetaBinomial(covariates_[t].getNTotalReads()(s), alpha, beta);
  }
//   observation.setNTotalReads(nTotalReads);
//   observation.setNMethylatedReads(nMethylatedReads);
}
/*
/// Samples latent variables and observations from the model.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
void Model<ModelParameters, LatentVariable, Covariate, Observation>::simulateData
(
  const unsigned int nObservations
)
{
  
  latentVariables_.resize(nObservations);
  observations_.resize(nObservations);
  
  latentVariables_[1] = sampleFromInitialDistribution();
  for (unsigned int t=1; t<nObservations; t++)
  {
    latentVariables_[t] = sampleFromObservationEquation(t, latentVariables_[t], latentVariables_[t-1]);
  }
  for (unsigned int t=0; t<nObservations; t++)
  {
    observations_[t] = sampleFromObservationEquation(t, latentVariables_[t]);
  }
  
  double meanNTotalReads = 30.0; // TODO: make this available from the outside
  unsigned int nSamples = 4; // TODO: make this available from the outside
  

  unsigned int R = modelParameters_.getNMethylationRegimes();

  unsigned int d = 1;
  MethylationRegimeType r = rng_.randomUniformInt(0,R-1);
  latentVariables_[0].setSojournTime(d);
  latentVariables_[0].setMethylationRegimeType(r);
  
  if (nObservations > 1)
  {
    for (unsigned int t=1; t<nObservations; t++)
    {
      if (arma::randu() <= modelParameters_.evaluateChangePointProbability(t,d,r))
      {
        d = 1;
        r = sampleInt(modelParameters_.getPTrans().col(r));
      }
      else
      {
        d++;
      }
      latentVariables_[t].setSojournTime(d);
      latentVariables_[t].setMethylationRegimeType(r);
    }
  }
  
  arma::uvec nTotalReads(nSamples);
  arma::uvec nMethylatedReads(nSamples);
  double alpha, beta;
  for (unsigned int t=0; t<nObservations; t++)
  {
    for (unsigned s=0; s<nSamples; s++)
    {
      r     = latentVariables_[t].getMethylationRegimeType();
      alpha = modelParameters_.getAlpha(r);
      beta  = modelParameters_.getBeta(r);
      nTotalReads(s)      = rng_.randomPoisson(meanNTotalReads);
      nMethylatedReads(s) = rng_.randomBetaBinomial(nTotalReads(s), alpha, beta);
    }
    observations_[t].setNTotalReads(nTotalReads);
    observations_[t].setNMethylatedReads(nMethylatedReads);
  }
}
*/
/// Evaluates the log-prior density associated with the initial state.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
double Model<ModelParameters, LatentVariable, Covariate, Observation>::evaluateLogInitialDensity
(
  const LatentVariable& latentVariableCurr
)
{
  return - std::log(modelParameters_.getNMethylationRegimes());
}
/// Evaluates the log-prior density associated with the initial state.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
double Model<ModelParameters, LatentVariable, Covariate, Observation>::evaluateLogTransitionDensity
(
  const unsigned int t,
  const LatentVariable& latentVariableCurr,
  const LatentVariable& latentVariablePrev
)
{
  unsigned int dCurr = latentVariableCurr.getSojournTime();
  unsigned int dPrev = latentVariablePrev.getSojournTime();
  MethylationRegimeType rCurr = latentVariableCurr.getMethylationRegimeType();
  MethylationRegimeType rPrev = latentVariablePrev.getMethylationRegimeType();
  double rho; 
  
  double logDensity = - std::numeric_limits<double>::infinity();
//   bool exitStatus = false;
  
  if (dCurr == 1 && rCurr != rPrev && dPrev >= modelParameters_.getU())
  {
//     rho = modelParameters_.evaluateChangePointProbability(t, dPrev, rPrev, exitStatus);
    rho = modelParameters_.getRho(dPrev, rPrev);
    
    /*
      if (std::abs(rho - modelParameters_.getRho(dPrev, rPrev)) > 0.000001)
      {
        std::cout << "rho difference in evaluateLogTransitionDensity() with dCurr == 1: " << rho - modelParameters_.getRho(dPrev, rPrev) << "; exitStatus: " << exitStatus << "; rho: " << rho << "; getRho(): " << modelParameters_.getRho(dPrev, rPrev) << std::endl;
      }
    */
    
    if (modelParameters_.getExitStatus(dPrev, rPrev))
    {
      logDensity = std::log(modelParameters_.getP()(rPrev,rCurr));
    }
    else 
    {
      logDensity = std::log(rho) + std::log(modelParameters_.getP()(rPrev,rCurr));
    }
  } 
  else if (dCurr > 1 && rCurr == rPrev)
  {
//     rho = modelParameters_.evaluateChangePointProbability(t, dPrev, rPrev, exitStatus);
    rho = modelParameters_.getRho(dPrev, rPrev);
    
    /*
      if (std::abs(rho - modelParameters_.getRho(dPrev, rPrev)) > 0.000001)
      {
        std::cout << "rho difference in evaluateLogTransitionDensity() with dCurr > 1: " << rho - modelParameters_.getRho(dPrev, rPrev) << "; exitStatus: " << exitStatus << "; rho: " << rho << "; getRho(): " << modelParameters_.getRho(dPrev, rPrev) << std::endl;
      }
    */
      
    if (!modelParameters_.getExitStatus(dPrev, rPrev) && rho <= 1) 
    {
      logDensity = std::log(1.0 - rho);
    }
  }
  return logDensity;
  
      ////////////////////////////////
      /*
    if (std::isnan(logDensity))
      {
        std::cout << "WARNING: " << 
        << "rho: " 
        << rho
        << "; omega: " 
        << modelParameters_.getOmega(rPrev)
        << "; kappa: " 
        << modelParameters_.getKappa(rPrev)
        << "; dCurr: " << dCurr
        << "; rCurr: " << rCurr
        << "; dPrev: " << dPrev
        << "; rPrev: " << rPrev
        << "; logDensity: "
        << logDensity
        << std::endl;
      }
      */
  ///////////////////////
  
}
/// Evaluates the log-observation density
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
double Model<ModelParameters, LatentVariable, Covariate, Observation>::evaluateLogObservationDensity
(
  const unsigned int t, 
  const LatentVariable& latentVariableCurr
)
{
  unsigned int rCurr = static_cast<unsigned int>(latentVariableCurr.getMethylationRegimeType());
  double alpha = modelParameters_.getAlpha(rCurr);
  double beta  = modelParameters_.getBeta(rCurr);
  double logDensity = 0.0;
  for (unsigned int s=0; s<getObservations(t).size(); s++)
  {
    logDensity += 
    evaluateLogBetaBinomialDensity(getObservations(t)(s), getCovariates(t).getNTotalReads()(s), alpha, beta);
  }
  return logDensity;
}
/// Evaluates the gradient of the log-density of the initial state.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
arma::colvec Model<ModelParameters, LatentVariable, Covariate, Observation>::evaluateGradThetaLogInitialDensity
(
  const LatentVariable& latentVariable
)
{ 
//   std::cout << "evaluateGradThetaLogInitialVarGamma" << std::endl;
  arma::colvec grad(getDimTheta(), arma::fill::zeros);
  return grad;
  
}
/// Evaluates the gradient of the log-transition density.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
arma::colvec Model<ModelParameters, LatentVariable, Covariate, Observation>::evaluateGradThetaLogTransitionDensity
(
  const unsigned int t,
  const LatentVariable& latentVariableCurr, 
  const LatentVariable& latentVariablePrev
)
{
  /// TODO: clean up this function and check correctness!
  
//   unsigned int R = modelParameters_.getNMethylationRegimes();
  unsigned int dCurr = latentVariableCurr.getSojournTime();
  unsigned int dPrev = latentVariablePrev.getSojournTime();
  unsigned int rCurr = static_cast<unsigned int>(latentVariableCurr.getMethylationRegimeType());
  unsigned int rPrev = static_cast<unsigned int>(latentVariablePrev.getMethylationRegimeType());
  unsigned int idxMin, idxMax, idx;
  
  arma::colvec grad(getDimTheta(), arma::fill::zeros);

  /*
  // Auxiliary quantities needed for the rest of the gradient.
  double bigH = 0.0;
  arma::colvec littleH(dPrev, arma::fill::zeros); // NOTE: the first element of this vector corresponds to distance $1$
  for (unsigned int d=modelParameters_.getU()-1; d<dPrev; d++)
  {
    littleH(d) = modelParameters_.evaluateLittleH(d+1,rPrev);
  }
  if (dPrev > 1)
  {
    bigH = arma::sum(littleH(arma::span(0,dPrev-2)));
  }
  if (bigH >= 1.0) // TODO: this is not correct
  { 
    std::cout << "WARNING: bigH >= 1!" << std::endl;
    bigH = 0.99999; // TODO
  }
  
  // Gradient for the parameters governing the 
  // "success probability" parameters of the negative-binomial
  // distribution.
  arma::colvec gradThetaLogLittleH(dPrev);
  for (unsigned int d=0; d<dPrev; d++)
  {
    gradThetaLogLittleH(d) = modelParameters_.evaluateGradThetaOmegaLogLittleH(d+1,rPrev);
  }
  idx = modelParameters_.getIndexOmega(rPrev);
  
  grad(idx) = gradThetaLogLittleH(dPrev-1);
  if (dPrev > 1)
  {
    for (unsigned int d=0; d<dPrev-1; d++)
    {
      grad(idx) = grad(idx) + (littleH(d) * gradThetaLogLittleH(d)) / (1.0 - bigH);
    }
//     grad(idx) = grad(idx) + arma::sum(littleH(arma::span(0,dPrev-2)) % gradThetaLogLittleH(arma::span(0,dPrev-2)),1)/ (1.0 - bigH);
  }
  
  if (std::abs(grad(idx) - modelParameters_.getGradThetaOmegaLogRho(dPrev, rPrev)) > 0.000001)
  {
    std::cout << "rho difference in sampleFromTransitionEquation(): " << grad(idx) - modelParameters_.getGradThetaOmegaLogRho(dPrev, rPrev) << std::endl;
  }
  */
  
  idx = modelParameters_.getIndexOmega(rPrev);
  grad(idx) = modelParameters_.getGradThetaOmegaLogRho(dPrev, rPrev);
  
  
  // Gradient for the parameters governing the 
  // "number of failures until success" parameters 
  // of the negative-binomial distribution.
  if (!modelParameters_.getIsKappaFixed())
  {
    /*
    for (unsigned int d=0; d<dPrev; d++)
    {
      gradThetaLogLittleH(d) = modelParameters_.evaluateGradThetaKappaLogLittleH(d+1,rPrev);
    }
    idx = modelParameters_.getIndexKappa(rPrev);
    grad(idx) = gradThetaLogLittleH(dPrev-1);
    if (dPrev > 1)
    {
      for (unsigned int d=0; d<dPrev-1; d++)
      {
        grad(idx) = grad(idx) + (littleH(d) * gradThetaLogLittleH(d)) / (1.0 - bigH);
      }
//       grad(idx) = grad(idx) + arma::sum(littleH(arma::span(0,dPrev-2)) % gradThetaLogLittleH(arma::span(0,dPrev-2)),1)/ (1.0 - bigH);
    }
    */
    idx = modelParameters_.getIndexOmega(rPrev);
    grad(idx) = modelParameters_.getGradThetaKappaLogRho(dPrev, rPrev);
  }
  
  if (dCurr == 1 && rCurr != rPrev && dPrev >= modelParameters_.getU())
  {
    arma::colvec aux = (-1.0) * modelParameters_.getPTrans().col(rPrev);
    aux(rCurr) = aux(rCurr) + 1;
    idxMin = modelParameters_.getIndexPMin(rPrev);
    idxMax = modelParameters_.getIndexPMax(rPrev);
    aux.shed_row(rPrev);
    grad(arma::span(idxMin,idxMax)) = aux;
  } 
  else if (dCurr > 1 && rCurr == rPrev)
  {
//     bool exitStatus = false;
//     double rho = modelParameters_.evaluateChangePointProbability(t, dPrev, rPrev, exitStatus); // = littleH(dPrev-1) / (1.0 - bigH);
    
    double rho = modelParameters_.getRho(dPrev, rPrev);
    
    if (!modelParameters_.getExitStatus(dPrev, rPrev) && rho < 1.0) // TODO: check if this will cause problems!
    {
      idx = modelParameters_.getIndexOmega(rPrev);
      grad(idx) = - grad(idx) * rho / (1.0 - rho);
      if (!modelParameters_.getIsKappaFixed())
      {
        idx = modelParameters_.getIndexKappa(rPrev);
        grad(idx) = - grad(idx) * rho / (1.0 - rho);
      }
    }
    else
    {
      grad.zeros(); // TODO: check if this is valid
    }
  }
  else // in this case, the transition density is zero!
  {
    grad.zeros();
  }
  
  /*
  if (grad.has_nan()) 
  {
    
    std::cout << "bigH calc manually: " << bigH << std::endl;
  }
  */
  
  return grad;
}
/// Evaluates the gradient of the log-observation density.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation>
arma::colvec Model<ModelParameters, LatentVariable, Covariate, Observation>::evaluateGradThetaLogObservationDensity
(
  const unsigned int t,
  const LatentVariable& latentVariableCurr
)
{
  arma::colvec grad(getDimTheta(), arma::fill::zeros);
  return grad;
}


///////////////////////////////////////////////////////////////////////////////
// Implementation of the Smc class template
///////////////////////////////////////////////////////////////////////////////

/// Holds a single particle used by the 
/// sequential Monte Carlo algorithm for change-point models.
template<class LatentVariable> class Particle
{
public:
  
  /// Returns the latent variable.
  const LatentVariable& getLatentVariable() const {return latentVariable_;}
  /// Returns the index of the current regime.
  unsigned int getRegime() const {return static_cast<unsigned int>(latentVariable_.getMethylationRegimeType());}
  /// Returns the sojourn time for the current regime.
  unsigned int getDistance() const {return latentVariable_.getSojournTime();}
  /// Specifies the latent variable.
  void setLatentVariable(const LatentVariable& latentVariable) {latentVariable_ = latentVariable;}
  /// Specifies the index of the current regime.
  void setRegime(const unsigned int regime) {latentVariable_.setMethylationRegimeType(static_cast<MethylationRegimeType>(regime));}
  /// Specifies the sojourn time for the current regime.
  void setDistance(const unsigned int distance) {latentVariable_.setSojournTime(distance);}
  
private:
  
  LatentVariable latentVariable_; // a value of the latent variable in the original parametrisation
  
};

/// Holds additional parameters used by the change-point SMC scheme.
class SmcParameters
{
public:
  // Empty
private:
  // Empty
};
/// Specifies the proposal parameters which do not change over the course of the algorithm.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
void Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::setKnownParameters
(
  const arma::colvec& varphi
)
{
  setNRegimes(model_.getModelParameters().getNMethylationRegimes());
}
/// Determines the proposal parameters from arma::colvec phi.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
void Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::setUnknownParameters
(
  const arma::colvec& phi
)
{
  // Empty
}
  /// Evaluates the log-initial proposal density.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
double Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::evaluateLogInitialProposalDensity
( 
  const Particle & particleCurr
)
{
  return model_.evaluateLogInitialDensity(particleCurr.getLatentVariable());
}
/// Evaluates the log-proposal density at time $t>1$
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
double Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::evaluateLogProposalDensity
(
  const unsigned int t,
  const Particle& particleCurr,
  const Particle& particlePrev
)
{
  return model_.evaluateLogTransitionDensity(t, particleCurr.getLatentVariable(), particlePrev.getLatentVariable());
}
/// Samples a particle at time $1$.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
void Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::sampleFromInitialProposalDistribution
(
  Particle& particleCurr
)
{
  particleCurr.setLatentVariable(model_.sampleFromInitialDistribution());
}
/// Samples a particle at time $t>1$.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
void Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::sampleFromProposalDistribution
(
  const unsigned int t, 
  Particle& particleCurr,
  const Particle& particlePrev
)
{
  particleCurr.setLatentVariable(model_.sampleFromTransitionEquation(t, particlePrev.getLatentVariable()));
}

///////////////////////////////////////////////////////////////////////////////
// Implementation of the OnlineMarginalSmoothing class template
///////////////////////////////////////////////////////////////////////////////

/// Evaluates the test function at some time $t$.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
double OnlineMarginalSmoothing<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::evaluateTestFunction
(
  const unsigned int timeIndex, 
  const unsigned int testFunctionIndex, 
  const LatentVariable& latentVariableCurr
)
{
//   std::cout << "storing " << testFunctionIndex << "th test function -- regime: " << latentVariableCurr.getMethylationRegimeType() << std::endl;
  if (testFunctionIndex == static_cast<unsigned int>(latentVariableCurr.getMethylationRegimeType()))
  {
    return 1.0;
  }
  else
  {
    return 0.0;
  }
}
/// Returns the number of test functions associated with each time step.
template <class ModelParameters, class LatentVariable, class Covariate, class Observation, class SmcParameters, class Particle>
unsigned int OnlineMarginalSmoothing<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle>::getNTestFunctions()
{
  return smc_.getNRegimes();
}

///////////////////////////////////////////////////////////////////////////////
// Other auxiliary functions.
///////////////////////////////////////////////////////////////////////////////

/// Converts LatentVariable to arma::uvec.
arma::uvec convertLatentVariableToArmaUvec(const LatentVariable& latentVariable)
{
  arma::uvec latVec(2);
  latVec(0) = latentVariable.getSojournTime();
  latVec(1) = static_cast<unsigned int>(latentVariable.getMethylationRegimeType());
  return latVec;
}
/// Converts std::vector<LatentVariable> to arma::umat.
arma::umat convertLatentVariablesToArmaUmat(const std::vector<LatentVariable>& latentVariables)
{
//     std::cout << "START: convertLatentVariablesToArmaUmat()" << std::endl; 
  unsigned int T = latentVariables.size();
  arma::umat latMat(2, T);
  for (unsigned int t=0; t<T; t++)
  {
    latMat.col(t) = convertLatentVariableToArmaUvec(latentVariables[t]);
  }
//       std::cout << "END: convertLatentVariablesToArmaUmat()" << std::endl;
  return latMat;
}
/// Converts arma::uvec to LatentVariable.
LatentVariable convertArmaUvecToLatentVariable(const arma::uvec& latVec)
{
  LatentVariable latentVariable;
  latentVariable.setSojournTime(latVec(0));
  latentVariable.setMethylationRegimeType(static_cast<MethylationRegimeType>(latVec(1)));
  return latentVariable;
}
/// Converts arma::umat to std::vector<LatentVariable>.
std::vector<LatentVariable> convertArmaUmatToLatentVariables(const arma::umat& latMat)
{
  unsigned int T = latMat.n_cols;
  std::vector<LatentVariable> latentVariables(T);
  for (unsigned int t=0; t<T; t++)
  {
    latentVariables[t] = convertArmaUvecToLatentVariable(latMat.col(t));
  }
  return latentVariables;
}
/// Converts arma::umat to std::vector<Covariate>.
std::vector<Covariate> convertArmaUmatToCovariates
(
  const arma::uvec& genomicPositions, 
  const arma::umat& nTotalReads
)
{
  unsigned int T = nTotalReads.n_cols;
//   unsigned int S = nTotalReads.n_rows;
  std::vector<Covariate> covariates(T);
  for (unsigned int t=0; t<T; t++)
  {
    covariates[t].setGenomicPosition(genomicPositions(t));
    covariates[t].setNTotalReads(nTotalReads.col(t));
  }
  return covariates;
}
/// Converts std::vector<Observation> to arma::umat.
arma::umat convertObservationsToArmaUmat(const std::vector<Observation>& observations)
{
  return convertStdVecToArmaMat(observations);
}
/// Converts arma::umat to std::vector<Observation>.
std::vector<Observation> convertArmaUmatToObservations(const arma::umat& obsMat)
{
  return convertArmaMatToStdVec(obsMat);
}

#endif
