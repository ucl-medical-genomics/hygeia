/// \file
/// \brief Gradient-ascent methods
///
/// A class for gradient-ascent type optimisation

#ifndef __GRADIENTASCENT_H
#define __GRADIENTASCENT_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <random>
#include <vector>
#include "rng/Rng.h"
#include "misc/misc.h"

class GradientAscent
{
  
public:
  
  /// Initialises the class.
  GradientAscent
  (
    // Empty
  ) 
  {
    setUseAdam(true);
    setNormaliseGradients(false);
    setLearningRateFactor(1.0);
  }
  
  /// Returns the iteration number.
  unsigned int getIter() const {return iter_;}
  /// Returns the size of the parameter vector.
  unsigned int getDimParam() const {return dimParam_;}
  /// Returns the exponent which determines the learning-rate decay.
  double getLearningRateExponent() const {return learningRateExponent_;}
  /// Returns the exponent which determines the learning-rate decay.
  double getLearningRateFactor() const {return learningRateFactor_;}
  /// Specifies the size of the parameter vector.
  void setDimParam(const unsigned int dimParam) {dimParam_ = dimParam;}
  /// Specifies whether the gradients should be normalised 
  /// according to their $L_1$ norm.
  void setNormaliseGradients(const bool normaliseGradients) {normaliseGradients_ = normaliseGradients;}
  /// Specifies if the ADAM optimiser should be used 
  /// (otherwise it is plain stochastic gradient ascent)
  void setUseAdam(const bool useAdam) {useAdam_ = useAdam;}
  /// Specifies the exponent which determines the learning-rate decay.
  void setLearningRateExponent(const double learningRateExponent) {learningRateExponent_ = learningRateExponent;}
  /// Specifies the factor which determines the learning-rate decay.
  void setLearningRateFactor(const double learningRateFactor) {learningRateFactor_ = learningRateFactor;}
  /// Initialise the class
  void initialise(const unsigned int dimParam)
  {
//      std::cout << " BEGIN: initialise(): learningRateExponent_: " << learningRateExponent_ << std::endl;
    setDimParam(dimParam);
    iter_ = 0;
//     learningRateExponent_ = 0.5; // the learning rate will be std::pow(1.0 / static_cast<double>(iter_), learningRateExponent_); // TODO: make this accessible/modifiable from the outside
    if (useAdam_)
    {
      // Use default parameters for ADAM: // TODO: make these accessible/modifiable from the outside
      beta1_ = 0.9;
      beta2_ = 0.999;
      epsilon_ = std::exp(-8.0*std::log(10));
      
      // Initialise the mean and vARIANCE AUXILIARY PARAMETERS:
      adamMean_.zeros(dimParam_);
      adamVar_.zeros(dimParam_);
    }
//         std::cout << " END: initialise(): learningRateExponent_: " << learningRateExponent_ << std::endl;
  }
  /// Updates the estimate using a single iteration of gradient-ascent/ADAM.
  void iterate(arma::colvec& param, const arma::colvec gradient)
  {
//         std::cout << " in iterate(): learningRateExponent_: " << learningRateExponent_ << std::endl;
    if (useAdam_)
    {
      /*
      if (param.has_nan())
      {
        std::cout << "param pre ADAM: " << param.t() << std::endl;
      }
      */
      adam(evaluateLearningRate(), gradient, adamMean_, adamVar_, param);
      
      /*
      if (param.has_nan())
      {
        std::cout << "param post ADAM: " << param.t() << std::endl;
      }
      */
    }
    else
    {
      if (normaliseGradients_)
      {
        param = param + evaluateLearningRate() * arma::normalise(gradient, 1); 
      }
      else
      {
        param = param + evaluateLearningRate() * gradient; 
      }
    }
    iter_++;
  }
    
private:
  
  /// Returns the learning rate.
  double evaluateLearningRate()
  { 
    return learningRateFactor_ / std::pow(static_cast<double>(iter_+1.0), learningRateExponent_);
  }
  /// Performs the ADAM update for a single stochastic gradient-ascent iteration.
  void adam
  (
    const double learningRate,
    const arma::colvec& gradient, 
    arma::colvec& mean, 
    arma::colvec& var,
    arma::colvec& param
  )
  {
         
    
//     if (gradient.has_nan() || param.has_nan())
//     {
//       std::cout << "iter_: " << iter_ << "; learningRate: " << learningRate << "; param: " << param.t()  << "gradient: " << gradient.t() << std::endl;
//         std::cout << "omega elements of grad: " << arma::trans(gradient(arma::span(30,35))) << " ";
//     }
    
    
    mean = beta1_ * mean + (1.0-beta1_) * gradient;
    var  = beta2_ * var  + (1.0-beta2_) * gradient % gradient;
    
        /*
    if (mean.has_nan()) 
    {
      std::cout << "updated ADAM mean: " << mean.t() << std::endl;
    }
    *
    if (var.has_nan()) 
    {
      std::cout << "updated ADAM var: " << var.t() << std::endl;
    }
    */
        
  
//       std::cout << "updated ADAM mean: " << mean.t() << std::endl;
//       std::cout << "updated ADAM var: " << var.t() << std::endl;
    
    
//     std::cout << "ADAM update difference: " << arma::trans(learningRate * mean % arma::pow(arma::sqrt(var/(1.0 - std::pow(beta2_, iter_+1))) + epsilon_, -1.0) / (1.0 - std::pow(beta1_, iter_+1))) << std::endl;
    
    param = param + learningRate * mean % arma::pow(arma::sqrt(var/(1.0 - std::pow(beta2_, iter_+1))) + epsilon_, -1.0) / (1.0 - std::pow(beta1_, iter_+1));
  }
  
  unsigned int iter_; // iteration number
  unsigned int dimParam_; // dimension of the parameter vector
  bool normaliseGradients_; // should the gradients be normalised according to their $L_1$ norm?
  double learningRateExponent_, learningRateFactor_; // the learning rate will be learningRateFactor_ / std::pow(static_cast<double>(iter_+1.0), learningRateExponent_) 
  bool useAdam_; // should the ADAM optimiser be used (otherwise it is plain stochastic gradient ascent)
  arma::colvec adamMean_, adamVar_; // auxiliary variables needed for the ADAM optimiser
  double beta1_, beta2_, epsilon_; // tuning parameters used by ADAM
  
};
#endif
