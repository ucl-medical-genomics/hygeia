/// \file
/// \brief A class for random number generation.
///
/// This file contains the Rng class which generates random numbers from 
/// various distributions.

#ifndef __RNG_H
#define __RNG_H

#include <RcppArmadillo.h>
#include <vector>
#include <iostream>
#include <random>

class Rng
{
  
public:
  
  //////////////////////////////////////////////////////////////////////////////
  // Constructors and Destructor
  //////////////////////////////////////////////////////////////////////////////
  
  /// Initialises the random number generator with engine specified by the user
  /// and with a pseudo-random system-generated seed.
  Rng()
  {
    randomiseSeed(); 
  }
  /// Initialises the random number generator with a seed
  /// specified by the user.
  Rng(const unsigned int seed) : seed_(seed)
  {
    engine_.seed(seed_); 
  }

  //////////////////////////////////////////////////////////////////////////////
  // Member functions for setting or obtaining information about engine/seed
  //////////////////////////////////////////////////////////////////////////////
  
  /// Changes the seed.
  void setSeed(const unsigned int seed)
  {
    seed_ = seed;
    engine_.seed(seed_); 
  }
  /// Randomises the seed.
  void randomiseSeed()
  {
    seed_ = std::random_device{}(); 
    engine_.seed(seed_);  
  }
  /// Returns the seed.
  unsigned int & getSeed() {return seed_;}
  
  //////////////////////////////////////////////////////////////////////////////
  // Member functions for sampling from specific parametrised distributions
  //////////////////////////////////////////////////////////////////////////////
  
  /// Returns a random number from a Bernoulli distribution with specified
  /// success probability.
  bool randomBernoulli(const double prob);
  /// Returns a random number from a binomial distribution with specified
  /// range and success probability parameters.
  unsigned int randomBinomial(const unsigned int range, const double prob);
  /// Returns a random number from a discrete distribution on 
  /// \f$\{0, 1, 2, ..., \mathit{weights.size()}-1\}\f$ for a specified weight vector.
  unsigned int randomDiscrete(std::vector<double>& weights);
  /// Returns a random number from an exponential distribution with specified
  /// rate parameter.
  double randomExponential(const double rate);
  /// Returns a random number from Fisher's F-distribution with specified
  /// degrees-of-freedom parameters.
  double randomFisher(const double df1, const double df2);
  /// Returns a random number from a geometric distribution on 
  /// \f$\{0, 1, 2, \dotsc \}\f$ with specified success probability.
  double randomGeometric(const double prob);
  /// Returns a random number from a gamma distribution with specified shape
  /// and scale parameters.
  double randomGamma(const double shape, const double scale);
  /// Returns a random number from a lognormal distribution with specified 
  /// location and scale parameters.
  double randomLognormal(const double location, const double scale);
  /// Returns a random number from a normal distribution with specified mean 
  /// and standard deviation.
  double randomNormal(const double mean, const double stdDev);
  /// Returns a random number from a discrete normal distribution 
  /// with specified mean and standard deviation parameters.
  int randomDiscreteNormal(const double mean, const double stdDev);
  /// Returns a random number from a left-truncated discrete normal distribution 
  /// with specified mean and standard deviation parameters.
  int randomLeftTruncatedDiscreteNormal(const int lower, const double mean, const double stdDev);
  /// Returns a random number from a righ-truncated discrete normal distribution 
  /// with specified mean and standard deviation parameters.
  int randomRightTruncatedDiscreteNormal(const int upper, const double mean, const double stdDev);
  /// Returns a random number from a (left- and right-)
  /// truncated discrete normal distribution 
  /// with specified mean and standard deviation parameters.
  int randomTruncatedDiscreteNormal(const int lower, const int upper, const double mean, const double stdDev);
  /// Returns a random number from a Poisson distribution with specified mean.
  unsigned int randomPoisson(const double mean);
  /// Returns a random number from a Student-t distribution with specified 
  /// degrees-of-freedom parameter.
  double randomStudent(const double df);
  /// Returns a random number from a Uniform distribution on 
  /// \f$\{\mathit{from}, \mathit{from}+1, \dotsc, \mathit{thru}-1, \mathit{thru}\}\f$.
  int randomUniformInt(const int from, const int thru);
  /// Returns a random number from a Uniform distribution on 
  /// \f$[\mathit{from}, \mathit{to}]\f$.
  double randomUniformReal(const double from, const double to);
  /// Returns a random number from a negative-binomial
  /// distribution with specified number of failures
  /// (extended to include any positive real number)
  /// and success probability.
  unsigned int randomNegativeBinomial(const double num, const double prob);
  /// Returns a random number from a beta distribution with two shape parameters.
  double randomBeta(const double shape1, const double shape2);
  /// Returns a random number from a beta-binomial
  /// distribution with specified range parameter (from the binomial distribution) 
  /// and the two shape parameters from the beta distribution.
  unsigned int randomBetaBinomial(const unsigned int range, const double shape1, const double shape2);
  
private:

  unsigned int seed_;
  std::mt19937 engine_;
  
};


////////////////////////////////////////////////////////////////////////////////
// Member functions for sampling from specific parametrised distributions
////////////////////////////////////////////////////////////////////////////////

/// Returns a random number from a Bernoulli distribution with specified
/// success probability.
bool Rng::randomBernoulli(const double prob)
{
  static std::bernoulli_distribution d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(prob));
}
/// Returns a random number from a binomial distribution with specified
/// range and success probability parameters.
unsigned int Rng::randomBinomial(const unsigned int range, const double prob)
{
  static std::binomial_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(range, prob));
}
/// Returns a random number from a discrete distribution on 
/// \f$\{0, 1, 2, ..., \mathit{weights.size()}-1\}\f$ for a specified weight vector.
unsigned int Rng::randomDiscrete(std::vector<double>& weights)
{
  //static std::discrete_distribution<> d{};
  //using parameterType = decltype(d)::param_type;
  //return d(engine_, parameterType(weights) );
  static std::discrete_distribution<unsigned int> d(weights.begin(), weights.end());
  return d(engine_);
}
/// Returns a random number from an exponential distribution with specified
/// rate parameter.
double Rng::randomExponential(const double rate)
{
  static std::exponential_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(rate));
}
/// Returns a random number from Fisher's F-distribution with specified
/// degrees-of-freedom parameters.
double Rng::randomFisher(const double df1, const double df2)
{
  static std::fisher_f_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(df1, df2));
}
/// Returns a random number from a geometric distribution on 
/// \f$\{0, 1, 2, \dotsc \}\f$ with specified success probability.
double Rng::randomGeometric(const double prob)
{
  static std::geometric_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(prob));
}
/// Returns a random number from a gamma distribution with specified shape
/// and scale parameters.
double Rng::randomGamma(const double shape, const double scale)
{
  static std::gamma_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(shape, scale));
}
/// Returns a random number from a lognormal distribution with specified 
/// location and scale parameters.
double Rng::randomLognormal(const double location, const double scale)
{
  static std::lognormal_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(location, scale));
}
/// Returns a random number from a normal distribution with specified mean 
/// and standard deviation.
double Rng::randomNormal(const double mean, const double stdDev)
{
  static std::normal_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(mean, stdDev));
}
/// Returns a random number from a discrete normal distribution 
/// with specified mean and standard deviation parameters.
int Rng::randomDiscreteNormal(const double mean, const double stdDev)
{
  static std::normal_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return std::floor(d(engine_, parameterType(mean, stdDev)));
}
/// Returns a random number from a left-truncated discrete normal distribution 
/// with specified mean and standard deviation parameters.
int Rng::randomLeftTruncatedDiscreteNormal(const int lower, const double mean, const double stdDev)
{
  double minU = R::pnorm(lower, mean, stdDev, true, false);
  double maxU = 1.0;
  double u = randomUniformReal(minU, maxU);
  
  if (u < 1.0)
  {
    int x = std::floor(R::qnorm(u, mean, stdDev, true, false));
    if (x < lower)
    {
      return lower;
    }
    else
    {
      return x;
    }
  }
  else
  {
    std::cout << "WARNING: u = " << "; minU: " << minU << "; maxU: " << maxU << "! Output of Rng::randomLeftTruncatedDiscreteNormal(lower, upper, mean, stdDev) set to lower!" << std::endl;
    return lower;
  }
}
/// Returns a random number from a righ-truncated discrete normal distribution 
/// with specified mean and standard deviation parameters.
int Rng::randomRightTruncatedDiscreteNormal(const int upper, const double mean, const double stdDev)
{
  double minU = 0.0;
  double maxU = R::pnorm(upper+1, mean, stdDev, true, false);
  double u = randomUniformReal(minU, maxU);
    
  if (u > 0.0)
  {
    int x = std::floor(R::qnorm(u, mean, stdDev, true, false));
    if (x > upper)
    {
      return upper;
    }
    else
    {
      return x;
    }
  }
  else
  {
    std::cout << "WARNING: u = " << u << "; minU: " << minU << "; maxU: " << maxU << "! Output of Rng::randomRightTruncatedDiscreteNormal(lower, upper, mean, stdDev) set to upper!" << std::endl;
    return upper;
  }
}
/// Returns a random number from a (left- and right-)
/// truncated discrete normal distribution 
/// with specified mean and standard deviation parameters.
int Rng::randomTruncatedDiscreteNormal(const int lower, const int upper, const double mean, const double stdDev)
{
  if (std::isfinite(mean) && std::isfinite(stdDev))
  {
    if (lower == upper)
    {
      return lower;
    }
    else if (lower > upper)
    {
      std::cout << "ERROR: lower > upper in randomTruncatedDiscreteNormal()!" << std::endl;
      return -1;
    }
    else
    {
      
      double minU = R::pnorm(lower, mean, stdDev, true, false);
      double maxU = R::pnorm(upper+1, mean, stdDev, true, false);
      double u = randomUniformReal(minU, maxU);
      
      if (u < 1.0 && u > 0.0)
      {
        int x = std::floor(R::qnorm(u, mean, stdDev, true, false));
        if (x < lower)
        {
          return lower;
        }
        else if (x > upper)
        {
          return upper;
        }
        else
        {
          return x;
        }
      }
      else if (u > 0.0) // i.e. if u is numerically 1.0
      {
        std::cout << "WARNING: u = " << u << "; minU: " << minU << "; maxU: " << maxU << "! Output of Rng::randomTruncatedDiscreteNormal(lower, upper, mean, stdDev) set to lower!" << std::endl;
        return lower;
      }
      else // i.e. if u is numerically 0
      {
        std::cout << "WARNING: u = " << u << "; minU: " << minU << "; maxU: " << maxU << "! Output of Rng::randomTruncatedDiscreteNormal(lower, upper, mean, stdDev) set to upper!" << std::endl;
        return upper;
      }
    }
  }
  else
  {
    return - std::numeric_limits<int>::infinity();
  }
}
/// Returns a random number from a Poisson distribution with specified mean.
unsigned int Rng::randomPoisson(const double mean)
{
  static std::poisson_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(mean));
}
/// Returns a random number from a Student-t distribution with specified 
/// degrees-of-freedom parameter.
double Rng::randomStudent(const double df)
{
  static std::student_t_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(df));
}
/// Returns a random number from a Uniform distribution on 
/// \f$\{\mathit{from}, \mathit{from}+1, \dotsc, \mathit{thru}-1, \mathit{thru}\}\f$.
int Rng::randomUniformInt(const int from, const int thru)
{
  static std::uniform_int_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(from, thru));
}
/// Returns a random number from a Uniform distribution on 
/// \f$[\mathit{from}, \mathit{to}]\f$.
double Rng::randomUniformReal(const double from, const double to)
{
  static std::uniform_real_distribution<> d{};
  using parameterType = decltype(d)::param_type;
  return d(engine_, parameterType(from, to));
}
/// Returns a random number from a negative-binomial
/// distribution with specified number of failures
/// (extended to include any positive real number)
/// and success probability.
/// NOTE: this is the parametrisation used by tWikipedia; 
/// in R, the negative-binomial distribution is parametrised with
/// failure probability 1-prob.
unsigned int Rng::randomNegativeBinomial(const double size, const double prob)
{
  double mean = randomGamma(size, prob/(1.0-prob));
  return randomPoisson(mean);
}
/// Returns a random number from a beta distribution with two shape parameters.
double Rng::randomBeta(const double shape1, const double shape2)
{
  double x = randomGamma(shape1, 1.0);
  double y = randomGamma(shape2, 1.0);
  return x / (x + y);
}
/// Returns a random number from a beta-binomial
/// distribution with specified range parameter (from the binomial distribution) 
/// and the two shape parameters from the beta distribution.
unsigned int Rng::randomBetaBinomial(const unsigned int range, const double shape1, const double shape2)
{
  double prob = randomBeta(shape1, shape2);
  return randomBinomial(range, prob);
}
#endif
