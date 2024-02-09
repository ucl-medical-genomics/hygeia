/// \file
/// \brief Defines number of small helper functions
///
/// This file contains several helper functions used by the 
/// model or by the algorithms.

#ifndef __MISC_H
#define __MISC_H

#include <RcppArmadillo.h>
#include <random>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// Some global constants
////////////////////////////////////////////////////////////////////////////////
/// Constant $\log(2\pi)$ needed e.g. for log-Gaussian densities.
const double log2Pi = std::log(2.0 * M_PI);
////////////////////////////////////////////////////////////////////////////////
// Logit functions and related transforms 
////////////////////////////////////////////////////////////////////////////////

/// Digamma function.
double digamma(const double x)
{
  return R::digamma(x); 
}
/// Logit transform:
double logit(const double x)
{
  return std::log(x) - std::log(1.0-x);
}
/// Inverse logit ("logistic", "sigmoidal") transform:
double inverseLogit(const double x)
{
  return 1.0 / (1.0 + std::exp((-1.0) * x));
}
/// Modified logit transform which maps
/// from $[-1,1]$ to the real line.
double biLogit(const double x)
{
  return std::log(1.0+x) - std::log(1.0-x);
}
/// Inverse of the modified logit transform
/// i.e. this maps from the real line to $[-1,1]$.
double inverseBiLogit(const double x)
{
  return 2.0 / (1.0 + std::exp((-1.0) * x)) - 1.0;
}

/// Digamma function.
arma::mat digamma(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return digamma(val);});
}
/// Logit transform:
arma::mat logit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return logit(val);});
}
/// Inverse logit ("logistic", "sigmoidal") transform:
arma::mat inverseLogit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return inverseLogit(val);});
}
/// Modified logit transform which maps
/// from $[-1,1]$ to the real line.
arma::mat biLogit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return biLogit(val);});
}
/// Inverse of the modified logit transform
/// i.e. this maps from the real line to $[-1,1]$
arma::mat inverseBiLogit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return inverseBiLogit(val);});
}


/// Derivative of logit transform:
double gradLogit(const double x)
{
  return 1.0/x + 1.0/(1.0-x);
}
/// Derivative of logit transform, evaluated at 
/// inverseLogit(x);
double gradLogitEvaluatedAtInverseLogit(const double x)
{
  return (2.0 + std::exp(-x) + std::exp(x));
}
/// Derivative of inverse logit ("logistic", "sigmoidal") transform:
double gradInverseLogit(const double x)
{
  return 1.0 / ((1.0 + std::exp(-x)) * (1.0 + std::exp(x)));
}
/// Derivative of inverse logit ("logistic", "sigmoidal") transform,
/// evaluated at logit(x);
double gradInverseLogitEvaluatedAtLogit(const double x)
{
  return x * (1.0 - x);
}
/// Derivative of inverse bilogit transform,
/// evaluated at bilogit(x);
double gradInverseBiLogitEvaluatedAtBiLogit(const double x)
{
  return (1.0 - x * x) / 2.0;
}
/// Derivative of logit transform:
arma::mat gradLogit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return gradLogit(val);});
}
/// Derivative of logit transform, evaluated at 
/// inverseLogit(x);
arma::mat gradLogitEvaluatedAtInverseLogit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return gradLogitEvaluatedAtInverseLogit(val);});
}
/// Derivative of nverse logit ("logistic", "sigmoidal") transform:
arma::mat gradInverseLogit(const arma::mat& x)
{
  arma::mat y = x;
  return y.transform([](double val) { return gradInverseLogit(val);});
}
/// Derivative of inverse logit ("logistic", "sigmoidal") transform,
/// evaluated at logit(x).
arma::mat gradInverseLogitEvaluatedAtLogit(const arma::mat x)
{
  arma::mat y = x;
  return y.transform([](double val) { return gradInverseLogitEvaluatedAtLogit(val);});
}

/// Declaration of normaliseExp needed here.
arma::vec normaliseExp(const arma::vec&);
/// Evaluates the gradient of the softmax function.
arma::mat gradSoftmax(const arma::colvec& x)
{
  arma::colvec y = arma::exp(normaliseExp(x));
  unsigned int n = y.size();
  return arma::diagmat(y) - arma::repmat(y, 1, n) * arma::repmat(y.t(), n, 1); 
}

////////////////////////////////////////////////////////////////////////////////
// Sampling from some well-known distributions
////////////////////////////////////////////////////////////////////////////////
/// Samples a single value from a categorical distribution (with size $1$)
/// (i.e. outputs a value in {0, ..., length(W)-1}).
unsigned int sampleInt(const arma::colvec& W) 
{ 
  unsigned int x = arma::conv_to<unsigned int>::from(arma::find(arma::cumsum(W) >= arma::randu(), 1, "first"));
  return x;
}
/// Samples a single value from a categoridal distribution (with size $1$)
/// (i.e. outputs a value in {0, ..., length(W)-1}).
unsigned int sampleInt(const arma::rowvec& W) 
{ 
  unsigned int x = arma::conv_to<unsigned int>::from(arma::find(arma::cumsum(W) >= arma::randu(), 1, "first"));
  return x;
}
/// Samples multiple values from a categorical distribution (with size $1$)
/// (i.e. outputs N values in {0, ..., length(W)-1})
arma::uvec sampleInt(unsigned int N, const arma::colvec& W) 
{ 
  arma::colvec cumW = arma::cumsum(W);
  arma::uvec x(N);
  for (unsigned int n=0; n<N; n++)
  {
    x(n) =  arma::conv_to<unsigned int>::from(arma::find(cumW >= arma::randu(), 1, "first"));
  }
  return x;
}
////////////////////////////////////////////////////////////////////////////////
// Evaluating log-densities of some well-known distributions
// as well as derivatives of these w.r.t. some of their parameters.
////////////////////////////////////////////////////////////////////////////////

/// Evaluates the logarithm of the difference between two evaluations of the 
/// cumulative distribution function of a univariate Gaussian distribution
/// at two different points.
double evaluateLogDifferenceOfDiscreteNormalCdfs(const double x1, const double x2, const double mean, const double sd)
{
  double y = std::log(R::pnorm(x1, mean, sd, true, false) - R::pnorm(x2, mean, sd, true, false));
  if (!std::isfinite(y))
  {
    y = std::log(R::pnorm(x2, mean, sd, false, false) - R::pnorm(x1, mean, sd, false, false));
  }
  return y;
}
/// Evaluates the logarithm of the difference between 1.0 and an evaluation of the 
/// cumulative distribution function of a univariate Gaussian distribution.
double evaluateLogOfOneMinusDiscreteNormalCdf(const double x, const double mean, const double sd)
{
  double y = std::log(1.0 - R::pnorm(x, mean, sd, true, false));
  if (!std::isfinite(y))
  {
    y = R::pnorm(x, mean, sd, false, true);
  }
  return y;
}

/// Evaluates the log-density of discrete normal distribution.
double evaluateLogDiscreteNormalDensity(const int x, const double mean, const double sd)
{
  return evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd);
//   return std::log(R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false));
}
/// Evaluates the derivative of the log-density of a
/// discrete normal distribution w.r.t. 
/// the mean parameter.
double evaluateGradMeanLogDiscreteNormalDensity(const int x, const double mean, const double sd)
{
  
  return (
  (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
  ) / sd;
//   return (
//   (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//   ) / sd;
}
/// Evaluates the derivative of the log-density of a
/// truncated discrete normal distribution w.r.t. 
/// the standard-deviation parameter.
double evaluateGradSdLogDiscreteNormalDensity(const int x, const double mean, const double sd)
{
//   double grad = (
//   (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//   ) / (sd * sd);
  
    double grad = (
  (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
  ) / (sd * sd);
  
  if (std::isfinite(grad))
  {
    return grad;
  }
  else
  {
    return 0.0;
  }
}
/// Evaluates the log-density of left-truncated discrete normal distribution.
double evaluateLogLeftTruncatedDiscreteNormalDensity(const int x, const int lower, const double mean, const double sd)
{
  if (x < lower)
  { 
    std::cout << "WARNING: x < lower in evaluateLogLeftTruncatedDiscreteNormalDensity" << std::endl;
    return - std::numeric_limits<double>::infinity();
  }
  else
  { 
    return evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd) - evaluateLogOfOneMinusDiscreteNormalCdf(lower, mean, sd);
//     return std::log(R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false)) - std::log(1.0 - R::pnorm(lower, mean, sd, true, false));
  }
  
  
}
/// Evaluates the derivative of the log-density of a
/// right-truncated discrete normal distribution w.r.t. 
/// the mean parameter.
double evaluateGradMeanLogLeftTruncatedDiscreteNormalDensity(const int x, const int lower, const double mean, const double sd)
{

  if (x < lower)
  {
    return 0.0;
  }
  else
  {
//     double grad = (
//     (0.0 - R::dnorm(lower, mean, sd, false))/ (1.0 - R::pnorm(lower, mean, sd, true, false)) - 
//     (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//     ) / sd;
    
    double grad = (
    (0.0 - R::dnorm(lower, mean, sd, false))/ std::exp(evaluateLogOfOneMinusDiscreteNormalCdf(lower, mean, sd)) - 
    (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
    ) / sd;
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the derivative of the log-density of a
/// right-truncated discrete normal distribution w.r.t. 
/// the variance parameter.
double evaluateGradVarLogLeftTruncatedDiscreteNormalDensity(const int x, const int lower, const double mean, const double sd)
{
  if (x < lower)
  {
    return 0.0;
  }
  else
  {
    
//     double grad = (
//     (0.0 - R::dnorm(lower, mean, sd, false)*(lower-mean))/ (1.0 - R::pnorm(lower, mean, sd, true, false)) - 
//     (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//     ) / (2.0 * sd * sd * sd);
    
        double grad = (
    (0.0 - R::dnorm(lower, mean, sd, false)*(lower-mean))/ std::exp(evaluateLogOfOneMinusDiscreteNormalCdf(lower, mean, sd)) - 
    (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
    ) / (2.0 * sd * sd * sd);
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the derivative of the log-density of a
/// right-truncated discrete normal distribution w.r.t. 
/// the standard-deviation parameter.
double evaluateGradSdLogLeftTruncatedDiscreteNormalDensity(const int x, const int lower, const double mean, const double sd)
{
  if (x < lower)
  {
    return 0.0;
  }
  else
  {
    double grad = (
    (0.0 - R::dnorm(lower, mean, sd, false)*(lower-mean))/ std::exp(evaluateLogOfOneMinusDiscreteNormalCdf(lower, mean, sd)) - 
    (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
    ) / (sd * sd);
    
//     double grad = (
//     (0.0 - R::dnorm(lower, mean, sd, false)*(lower-mean))/ (1.0 - R::pnorm(lower, mean, sd, true, false)) - 
//     (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//     ) / (sd * sd);
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the log-density of right-truncated discrete normal distribution.
double evaluateLogRightTruncatedDiscreteNormalDensity(const int x, const int upper, const double mean, const double sd)
{
  if (x > upper)
  {
    return - std::numeric_limits<double>::infinity();
  }
  else
  {
    return std::log(R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false)) - R::pnorm(upper+1, mean, sd, true, true);
//     return std::log(R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false)) - std::log(R::pnorm(upper+1, mean, sd, true, false));
  }
}
/// Evaluates the derivative of the log-density of a
/// right-truncated discrete normal distribution w.r.t. 
/// the mean parameter.
double evaluateGradMeanLogRightTruncatedDiscreteNormalDensity(const int x, const int upper, const double mean, const double sd)
{
  if (x > upper)
  {
    return 0.0;
  }
  else
  {
    double grad = (
    (R::dnorm(upper+1, mean, sd, false) - 0.0) / (R::pnorm(upper+1, mean, sd, true, false) - 0.0) - 
    (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
    ) / sd;
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the derivative of the log-density of a
/// right-truncated discrete normal distribution w.r.t. 
/// the standard-deviation parameter.
double evaluateGradSdLogRightTruncatedDiscreteNormalDensity(const int x, const int upper, const double mean, const double sd)
{
  if (x > upper)
  {
    return 0.0;
  }
  else
  {
    double grad = (
    (R::dnorm(upper+1, mean, sd, false)*(upper+1-mean) - 0.0)/ (R::pnorm(upper+1, mean, sd, true, false) - 0.0) - 
    (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
    ) / (sd * sd);
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the log-density of right-truncated discrete normal distribution.
double evaluateLogTruncatedDiscreteNormalDensity(const int x, const int lower, const int upper, const double mean, const double sd)
{
  if (x < lower || x > upper)
  {
    return - std::numeric_limits<double>::infinity();
  }
  else
  {
//     return std::log(R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false)) - std::log(R::pnorm(upper+1, mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false));
    return evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd) - evaluateLogDifferenceOfDiscreteNormalCdfs(upper+1, lower, mean, sd);
  }
}
/// Evaluates the derivative of the log-density of a
/// truncated discrete normal distribution w.r.t. 
/// the mean parameter.
double evaluateGradMeanLogTruncatedDiscreteNormalDensity(const int x, const int lower, const int upper, const double mean, const double sd)
{
  if (x < lower || x > upper)
  {
    return 0.0;
  }
  else
  {
    double grad = (
    (R::dnorm(upper+1, mean, sd, false) - R::dnorm(lower, mean, sd, false))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(upper+1, lower, mean, sd)) - 
    (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
    ) / sd;
    
//     double grad = (
//     (R::dnorm(upper+1, mean, sd, false) - R::dnorm(lower, mean, sd, false))/ (R::pnorm(upper+1, mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false)) - 
//     (R::dnorm(x+1, mean, sd, false) - R::dnorm(x, mean, sd, false))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//     ) / sd;
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the derivative of the log-density of a
/// truncated discrete normal distribution w.r.t. 
/// the standard-deviation parameter.
double evaluateGradSdLogTruncatedDiscreteNormalDensity(const int x, const int lower, const int upper, const double mean, const double sd)
{
  if (x < lower || x > upper)
  {
    return 0.0;
  }
  else
  {
    
    
    double grad = (
    (R::dnorm(upper+1, mean, sd, false)*(upper+1-mean) - R::dnorm(lower, mean, sd, false)*(lower-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(upper+1, lower, mean, sd)) - 
    (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
    ) / (sd * sd);
//     double grad = (
//     (R::dnorm(upper+1, mean, sd, false)*(upper+1-mean) - R::dnorm(lower, mean, sd, false)*(lower-mean))/ (R::pnorm(upper+1, mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false)) - 
//     (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//     ) / (sd * sd);
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the derivative of the log-density of a
/// truncated discrete normal distribution w.r.t. 
/// the standard-deviation parameter.
double evaluateGradVarLogTruncatedDiscreteNormalDensity(const int x, const int lower, const int upper, const double mean, const double sd)
{
  if (x < lower || x > upper)
  {
    return 0.0;
  }
  else
  {
    
    double grad = (
    (R::dnorm(upper+1, mean, sd, false)*(upper+1-mean) - R::dnorm(lower, mean, sd, false)*(lower-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(upper+1, lower, mean, sd)) - 
    (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ std::exp(evaluateLogDifferenceOfDiscreteNormalCdfs(x+1, x, mean, sd))
    ) / (2.0 * sd * sd * sd);
    
//     double grad = (
//     (R::dnorm(upper+1, mean, sd, false)*(upper+1-mean) - R::dnorm(lower, mean, sd, false)*(lower-mean))/ (R::pnorm(upper+1, mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false)) - 
//     (R::dnorm(x+1, mean, sd, false)*(x+1-mean) - R::dnorm(x, mean, sd, false)*(x-mean))/ (R::pnorm(x+1, mean, sd, true, false) - R::pnorm(x, mean, sd, true, false))
//     ) / (2.0 * sd * sd * sd);
    
    if (std::isfinite(grad))
    {
      return grad;
    }
    else
    {
      return 0.0;
    }
  }
}
/// Evaluates the log-density of a Poisson distribution for some
/// specified mean parameter
double evaluateLogPoissonDensity(const unsigned int x, const double mean)
{
  if (mean > 0.0) 
  {
    return x * std::log(mean) - mean - std::lgamma(x+1);
  }
  else
  {
    return - std::numeric_limits<double>::infinity();
  }
}
/// Evaluates the derivative of the log-density of a Poisson distribution 
/// w.r.t. to the mean parameter
double evaluateGradMeanLogPoissonDensity(const unsigned int x, const double mean)
{
  return static_cast<double>(x) / mean - 1.0;
}
/// Evaluates the log-density of a binomial distribution for some
/// specified range and probability parameters.
double evaluateLogBinomialDensity
(
  const unsigned int x, 
  const unsigned int range, 
  const double prob
)
{
  if ((x == 0) && (prob == 0))
  {
    return 0.0;
  }
  else if (x <= range)
  {
    return std::lgamma(range+1) - std::lgamma(x+1) - std::lgamma(range-x+1) + x * std::log(prob) + (range-x) * std::log(1.0 - prob);
  }
  else
  {
    return - std::numeric_limits<double>::infinity();
  }
}
/// Evaluates the derivative of the log-density of a binomial distribution
/// w.r.t. the probability parameter.
double evaluateGradProbLogBinomialDensity(const unsigned int x, const unsigned int range, const double prob)
{
  if (x <= range)
  {
    return static_cast<double>(x) / prob - static_cast<double>(range - x) / (1.0 - prob);
  }
  else
  {
    return 0.0;
  }
}
/// Evaluates the log-density of a multinomial distribution for some
/// specified size and probability parameters.
double evaluateLogMultinomialDensity
(
  const arma::uvec& x, 
  const unsigned int size, 
  const arma::colvec& prob
)
{
  if (arma::sum(x) == size)
  {
    return std::lgamma(size+1) - arma::accu(arma::lgamma(x + arma::ones(arma::size(x)))) + arma::accu(x % arma::log(prob));
  }
  else
  {
    return - std::numeric_limits<double>::infinity();
  }
}
/// Evaluates the gradient of the log-density of a 
/// multinomial distribution w.r.t. a parameter vector which 
/// parametrises the probabilities through a softmax fuction.
arma::colvec evaluateGradSoftmaxLogMultinomialDensity
(
  const arma::uvec& x, 
  const unsigned int size, 
  const arma::colvec& prob
)
{
  if (arma::sum(x) == size)
  {
    return prob * arma::accu(x) + x % (arma::ones<arma::colvec>(x.size()) - 2.0 * prob);
  }
  else
  {
    return arma::zeros<arma::colvec>(x.size());
  }
}
/// Evaluates the log-density of a beta-binomial distribution
/// specified range parameter (from the binomial distribution) 
/// and the two shape parameters from the beta distribution.
/// NOTE: this is the parametrisation used by Wikipedia; 
/// in R, the negative-binomial distribution is parametrised with
/// failure probability 1-prob.
double evaluateLogBetaBinomialDensity(const unsigned int x, const unsigned int range, const double shape1, const double shape2)
{
  if (x <= range)
  {
    return std::lgamma(range+1) - std::lgamma(x+1) - std::lgamma(range-x+1) + std::lgamma(x+shape1) + std::lgamma(range - x + shape2) - std::lgamma(range + shape1 + shape2) + std::lgamma(shape1 + shape2) - std::lgamma(shape1) - std::lgamma(shape2);
  }
  else
  {
    return - std::numeric_limits<double>::infinity();
  }
}
/// Evaluates the derivative of the log-density of a 
/// beta-binomial distribution w.r.t. the first shape parameter
double evaluateGradShape1LogBetaBinomialDensity(const unsigned int x, const unsigned int range, const double shape1, const double shape2)
{
  if (x <= range)
  {
    return digamma(x+shape1) - digamma(range + shape1 + shape2) + digamma(shape1 + shape2) - digamma(shape1);
  }
  else
  {
    return 0.0;
  }
}
/// Evaluates the derivative of the log-density of a 
/// beta-binomial distribution w.r.t. the second shape parameter
double evaluateGradShape2LogBetaBinomialDensity(const unsigned int x, const unsigned int range, const double shape1, const double shape2)
{
  if (x <= range)
  {
    return digamma(range - x + shape2) - digamma(range + shape1 + shape2) + digamma(shape1 + shape2) - digamma(shape2);
  }
  else
  {
    return 0.0;
  }
}
/// Evaluates the log-density of a negative-binomial distribution
/// with specified number of failures (definition extended to allow for
/// this parameter to be real-valued) and specified success probability.
/// NOTE: this is the parametrisation used by tWikipedia; 
/// in R, the negative-binomial distribution is parametrised with
/// failure probability 1-prob.
double evaluateLogNegativeBinomialDensity
(
  const unsigned int x, 
  const double size,
  const double prob
)
{
  
  if ((x == 0) && (prob == 0))
  {
    return 0.0;
  }
  else if (prob == 0)
  {
    return - std::numeric_limits<double>::infinity();
  }
  else
  {
    return std::lgamma(x + size) - std::lgamma(size) - std::lgamma(x+1) + size * std::log(1-prob) + x * std::log(prob);
  }
}
/// Evaluates the derivative of the log-density of a negative-binomial 
/// distribution w.r.t. the number-of-failures parameter.
/// NOTE: this is the parametrisation used by tWikipedia; 
/// in R, the negative-binomial distribution is parametrised with
/// failure probability 1-prob.
double evaluateGradSizeLogNegativeBinomialDensity
(
  const unsigned int x, 
  const double size,
  const double prob
)
{
  return digamma(x + size) - digamma(size) + log(1.0-prob);
}
/// Evaluates the derivative of the log-density of a negative-binomial 
/// distribution w.r.t. the success-probability parameter.
double evaluateGradProbLogNegativeBinomialDensity
(
  const unsigned int x, 
  const double size,
  const double prob
)
{
  return static_cast<double>(x) / prob - size / (1.0 - prob);
}
/// Evaluates the density of a multivariate Gaussian distribution.
double evaluateLogMultivariateNormalDensity
(
  const arma::colvec& x,  
  const arma::colvec& mean,  
  const arma::mat& sigma,
  bool isChol = false // are we passing the Cholesky decomposition of the covariance matrix (assumed to be upper(!) triangular)
) 
{ 
  arma::mat rooti; // Inverse root of the covariance matrix
  if (isChol == false) 
  {
    rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  } 
  else 
  {
    rooti = arma::trans(arma::inv(arma::trimatu(sigma)));
  }
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(x.n_rows)/2.0) * log2Pi;
  arma::colvec z = rooti * (x - mean);    
  return constants - 0.5 * arma::sum(z%z) + rootisum;     
}
////////////////////////////////////////////////////////////////////////////////
// Normalising vectors in log-space
////////////////////////////////////////////////////////////////////////////////

/// Computes a sum of a vector in log-space
/// (returns logarithm of the sum).
double sumExp(const arma::vec& logW) 
{ 
  double logWMax = arma::max(logW);
  if (logWMax > - std::numeric_limits<double>::infinity())
  {
    return logWMax + log(arma::sum(arma::exp(logW - logWMax)));
  }
  else
  {
//     std::cout << "NOTE: sumExp() returned -std::numeric_limits<double>::infinity()" << std::endl;
    return - std::numeric_limits<double>::infinity();
  }
}
/// Computes a sum of each column of a matrix in log-space.
/// (returns a column vector containing the logarithms of each sum).
arma::colvec sumExpMult(const arma::mat& logW) 
{ 
  unsigned int n = logW.n_cols; // number of columns to normalise
  arma::colvec y(n);
  double logWMax;
  
  for (unsigned int i = 0; i < n; i++)
  {
    logWMax = arma::max(logW.col(i));
    if (logWMax > - std::numeric_limits<double>::infinity())
    {
      y(i) = logWMax + log(arma::sum(arma::exp(logW.col(i) - logWMax)));
    }
    else
    {
      y(i) = - std::numeric_limits<double>::infinity();
    }
  }
  return y;
}
/// Normalises a single distribution in log-space (returns 
/// normalised weights in log space).
arma::vec normaliseExp(const arma::vec& logW) 
{ 
  double logWMax = arma::max(logW);
  double logZ = logWMax + log(arma::sum(arma::exp(logW - logWMax)));
  return(logW - logZ);
}
/// Normalises a single distribution in log-space (returns 
/// normalised weights and normalising constant in log space).
arma::vec normaliseExp(const arma::vec& logW, double& logZ) 
{ 
  double logWMax = arma::max(logW);
  logZ = logWMax + log(arma::sum(arma::exp(logW - logWMax)));
  return(logW - logZ);
}
/// Normalises a single distribution in log-space and overwrites the vector of weights
void normaliseExpInplace(arma::colvec& logW) 
{ 
  double logWMax = arma::max(logW);
  double logZ = logWMax + log(arma::sum(arma::exp(logW - logWMax)));
  // return the 1-unit norm (to make sure the elements of the vector sum to 1)
  logW = arma::normalise(arma::exp(logW - logZ), 1); 
}
////////////////////////////////////////////////////////////////////////////////
// Conversions between std::vector<arma::colvec> and arma::mat
////////////////////////////////////////////////////////////////////////////////
/// Converts a std::vector<arma::colvec> of length T 
/// (in which each element is an N-dimensional arma::colvec)
/// to an (N,T)-dimensional arma::mat.
void convertStdVecToArmaMat(const std::vector<arma::colvec>& x, arma::mat& y)
{
  unsigned int N = x[0].n_rows;
  unsigned int T = x.size();
  
  y.set_size(N, T);
  for (unsigned int t=0; t<T; t++)
  {
    y.col(t) = x[t];
  }
}
/// Converts a std::vector<arma::colvec> of length T 
/// (in which each element is an N-dimensional arma::colvec)
/// to an (N,T)-dimensional arma::mat.
arma::mat convertStdVecToArmaMat(const std::vector<arma::colvec>& x)
{
  unsigned int N = x[0].n_rows;
  unsigned int T = x.size();
  
  arma::mat y(N, T);
  for (unsigned int t=0; t<T; t++)
  {
    y.col(t) = x[t];
  }
  return y;
}
/// Converts a std::vector<arma::uvec> of length T 
/// (in which each element is an N-dimensional arma::uvec)
/// to an (N,T)-dimensional arma::umat.
void convertStdVecToArmaMat(const std::vector<arma::uvec>& x, arma::umat& y)
{
  unsigned int N = x[0].n_rows;
  unsigned int T = x.size();
  
  y.set_size(N, T);
  for (unsigned int t=0; t<T; t++)
  {
    y.col(t) = x[t];
  }
}
/// Converts a std::vector<arma::uvec> of length T 
/// (in which each element is an N-dimensional arma::uvec)
/// to an (N,T)-dimensional arma::umat.
arma::umat convertStdVecToArmaMat(const std::vector<arma::uvec>& x)
{
  unsigned int N = x[0].n_rows;
  unsigned int T = x.size();
  
  arma::umat y(N, T);
  for (unsigned int t=0; t<T; t++)
  {
    y.col(t) = x[t];
  }
  return y;
}
/// Converts an (N,T)-dimensional arma::mat to 
/// a std::vector<arma::colvec> of length T.
void convertArmaMatToStdVec(const arma::mat& y, std::vector<arma::colvec>& x)
{
  unsigned int T = y.n_cols;
  x.resize(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y.col(t);
  }
}
/// Converts an (N,T)-dimensional arma::mat to 
/// a std::vector<arma::colvec> of length T.
std::vector<arma::colvec> convertArmaMatToStdVec(const arma::mat& y)
{
  unsigned int T = y.n_cols;
  std::vector<arma::colvec> x(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y.col(t);
  }
  return x;
}
/// Converts an (N,T)-dimensional arma::umat to 
/// a std::vector<arma::uvec> of length T.
void convertArmaMatToStdVec(const arma::umat& y, std::vector<arma::uvec>& x)
{
  unsigned int T = y.n_cols;
  x.resize(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y.col(t);
  }
}
/// Converts an (N,T)-dimensional arma::umat to 
/// a std::vector<arma::uvec> of length T.
std::vector<arma::uvec> convertArmaMatToStdVec(const arma::umat& y)
{
  unsigned int T = y.n_cols;
  std::vector<arma::uvec> x(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y.col(t);
  }
  return x;
}
/// Converts an (N,T)-dimensional arma::colvec to 
/// a std::vector<double> of length T.
std::vector<double> convertArmaVecToStdVec(const arma::colvec& y)
{
  unsigned int T = y.size();
  std::vector<double> x(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y(t);
  }
  return x;
}
/// Converts an (N,T)-dimensional arma::rowevec to 
/// a std::vector<double> of length T.
std::vector<double> convertArmaVecToStdVec(const arma::rowvec& y)
{
  unsigned int T = y.size();
  std::vector<double> x(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y(t);
  }
  return x;
}
/// Converts an (N,T)-dimensional arma::colvec to 
/// a std::vector<unsigned int> of length T.
std::vector<unsigned int> convertArmaVecToStdVec(const arma::ucolvec& y)
{
  unsigned int T = y.size();
  std::vector<unsigned int> x(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y(t);
  }
  return x;
}
/// Converts an (N,T)-dimensional arma::rowevec to 
/// a std::vector<unsigned int> of length T.
std::vector<unsigned int> convertArmaVecToStdVec(const arma::urowvec& y)
{
  unsigned int T = y.size();
  std::vector<unsigned int> x(T);
  for (unsigned int t=0; t<T; t++)
  {
    x[t] = y(t);
  }
  return x;
}
/// Converts a vector (of $d * (d+1)/2$ elements) into a
/// lower-triangular $d \times d$ matrix by filling the matrix
/// column-wise.
arma::mat convertVecToTrimatL(const arma::colvec& x)
{
  unsigned int d = - 0.5 + std::sqrt(0.25 + 2 * x.size());
  arma::mat A(d, d, arma::fill::zeros);
  unsigned int k=0;
  for (unsigned int j=0; j<d; j++)
  {
    for (unsigned int i=j; i<d; i++)
    {
      A(i,j) = x(k);
      k++;
    }
  }
  return A;
}
/// Converts a lower-triangular $d \times d$ matrix into a 
/// vector (of $d * (d+1)/2$ elements) by taking elements
/// column-wise from the matrix.
arma::colvec convertTrimatLToVec(const arma::mat& A)
{
  unsigned int d = A.n_rows;
  arma::colvec x(d*(d+1)/2);
  unsigned int k=0;
  for (unsigned int j=0; j<d; j++)
  {
    for (unsigned int i=j; i<d; i++)
    {
      x(k) = A(i,j);
      k++;
    }
  }
  return x;
}
/// Removes the elements/rows of the matrix X which
/// correspond to the upper triangular elements of the matrices
/// arma::reshape(X.col(i), d, d), where d = std::sqrt(X.n_rows).
arma::mat shedTrimatURows(const arma::mat& X)
{
  unsigned int d = std::sqrt(X.n_rows);
  arma::umat sub(2, d*(d-1)/2);
  unsigned int k=0;
  
  for (unsigned int j=0; j<d; j++)
  {
    for (unsigned int i=0; i<j; i++)
    {
      sub(0,k) = i;
      sub(1,k) = j;
      k++;
    }
  }
  arma::mat Y = X;
  arma::uvec rowsToShed = arma::sub2ind(arma::size(d,d), sub);
  Y.shed_rows(rowsToShed);
  return Y;
}
////////////////////////////////////////////////////////////////////////////////
// Printing useful information to the command line
////////////////////////////////////////////////////////////////////////////////
/// Prints a matrix to the command line in such a way that it can easily be 
/// copied-and-pasted into R.
void printMatrixToR(const arma::mat& A)
{
  std::cout << "A <- matrix(c(";
  for (unsigned int j=0; j<A.n_cols; j++)
  {
    for (unsigned int i=0; i<A.n_rows; i++)
    {
      std::cout << A(i,j);
      if (i<A.n_rows-1 || j<A.n_cols)
      {
        std::cout << ",";
      }
      else
      {
        std::cout << ")"; 
      }
    }
  }
  std::cout << A.n_rows << "," << A.n_cols << ")" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
// Rootfinding
////////////////////////////////////////////////////////////////////////////////
/// The bisection method.
double bisection(
  bool& isBracketing,
  const std::function<double(const double)>& f,
  const double lb, 
  const double ub,
  const double tolX, // tolerance: interval boundaries
  const double tolF, // tolerance: values of f
  const unsigned int nIterations
)
{
  double a = lb;
  double b = ub;
  double fa = f(lb);
  double fb = f(ub);
  unsigned int i = 0;
  double x = a;
  double fx = fa;

  if (fa * fb > 0) 
  { // CASE I: interval not bracketing
    isBracketing = false;
    std::cout << "Error: interval not bracketing!" << std::endl;
  }
  else
  { // CASE II: interval is bracketing
    isBracketing = true;

    while ((i == 0) || ((std::abs(a - b) > tolX) && (std::abs(fx) > tolF) && (i < nIterations)))
    {
      x = (a + b) / 2.0; // new midpoint
      fx = f(x);
      i++;
      
      if (fa * fx <= 0) // i.e. the root must be in [a,x]
      {
        b  = x;
//           fb = fx;
      }
      else // i.e. the root must be in (x,b]
      {
        a  = x;
        fa = fx;
      }
    }
  }
  if (i == nIterations)
  {
    std::cout << "Error: maximum number of iterations exceeded!" << std::endl;
  }
  return x;
}


#endif
