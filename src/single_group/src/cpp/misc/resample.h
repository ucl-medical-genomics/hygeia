/// \file
/// \brief Some helper functions for resampling and dealing particle weights.
///
/// This file contains the functions for performing resampling in SMC algorithms
/// as well as some auxiliary functions for normalising the weights in log-space.

#ifndef __RESAMPLE_H
#define __RESAMPLE_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <random>
#include <iostream>
#include <vector>
#include "misc/misc.h"

namespace resample
{
  
  ////////////////////////////////////////////////////////////////////////////////
  // Standard multinomial resampling
  ////////////////////////////////////////////////////////////////////////////////
  
  /// \brief Performs multinomial resampling (using already normalised weights)
  void multinomialBase
  (
    arma::uvec& parentIndices,   // stores the post-resampling particle labels
    const arma::colvec& selfNormalisedWeightsPrev, // self-normalised particle weights
    unsigned int N               // desired number of offspring
  )        
  {
    parentIndices = sampleInt(N, selfNormalisedWeightsPrev);
  }
  /// \brief Performs multinomial resampling (using already normalised weights)
  void multinomial
  (
    arma::uvec& parentIndices,   // stores the post-resampling particle labels
    const arma::colvec& selfNormalisedWeightsPrev, // self-normalised particle weights
    unsigned int N               // desired number of offspring
  )        
  {
    parentIndices = sampleInt(N, selfNormalisedWeightsPrev);
  }
  /// \brief Performs multinomial resampling (and self-normalises the weights)
  void multinomial
  (
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    arma::colvec& logUnnormalisedWeightsCurr, // the log-unnormalised post-resampling weights
    const arma::colvec& selfNormalisedWeightsPrev, // self-normalised particle weights
    const unsigned int N, // total desired number of offspring
    const arma::colvec& logUnnormalisedWeightsPrev, // the the log-unnormalised particle weights
    const double logSumOfUnnormalisedWeightsPrev // the logarithm of the sum of the unnormalised weights from the previous time step
  )        
  {
    multinomialBase(parentIndices, selfNormalisedWeightsPrev, N);
    logUnnormalisedWeightsCurr.fill(logSumOfUnnormalisedWeightsPrev - std::log(N));
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Conditional multinomial resampling
  ////////////////////////////////////////////////////////////////////////////////

  /// \brief Performs conditional systematic resampling given a particular
  /// uniform random number.
  void conditionalMultinomialBase
  (
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    unsigned int& b,       // particle index of the distinguished particle
    const arma::colvec& w, // particle weights
    const unsigned int N,  // total number of offspring
    const unsigned int a   // parent index
  )
  { 
    parentIndices = sampleInt(N, w);
    b = 0; // NOTE: for simplicity, we set the index of the reference input path to 0 when using multinomial resampling!
    parentIndices(b) = a;
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  // Standard systematic resampling
  ////////////////////////////////////////////////////////////////////////////////
  
  /// \brief Performs systematic resampling given a particular
  /// uniform random number  (using already normalised weights)
  void systematicBase
  (
    const double u,              // unform random variable supplied by the user
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    const arma::colvec& selfNormalisedWeightsPrev, // self-normalised particle weights
    unsigned int N               // total number of offspring
  )        
  {
    
//     std::cout << "systematic resampling 1" << std::endl;
    arma::colvec T = (arma::linspace(0, N - 1, N) + u) / N;
//     std::cout << "systematic resampling 2" << std::endl;
    arma::colvec Q = arma::cumsum(selfNormalisedWeightsPrev);

//         std::cout << "systematic resampling 3" << std::endl;
//         std::cout << "N: " << N << " " << parentIndices.size() << " " << w.t() << " " << Q.t() << " " << T.t() << std::endl;
    unsigned int i = 0;
    unsigned int j = 0;
      
    while (j < N) 
    {
      if (T(j) <= Q(i)) 
      {
        parentIndices(j) = i;
        ++j;
      } 
      else 
      {
        ++i;
      }
    }
//         std::cout << "systematic resampling 4" << std::endl;
  }
  /// \brief Performs systematic resampling (using already normalised weights)
  void systematic
  (
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    const arma::colvec& selfNormalisedWeightsPrev,       // self-normalised particle weights
    const unsigned int N               // total number of offspring
  )
  {
    systematicBase(arma::randu(), parentIndices, selfNormalisedWeightsPrev, N);
  }
  /// \brief Performs systematic resampling (and self-normalises the weights)
  void systematic
  (
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    arma::colvec& logUnnormalisedWeightsCurr, // the log-unnormalised post-resampling weights
    const arma::colvec& selfNormalisedWeightsPrev, // self-normalised particle weights
    const unsigned int N, // total desired number of offspring
    const arma::colvec& logUnnormalisedWeightsPrev, // the the log-unnormalised particle weights
    const double logSumOfUnnormalisedWeightsPrev // the logarithm of the sum of the unnormalised weights from the previous time step
  )        
  {
    systematicBase(arma::randu(), parentIndices, selfNormalisedWeightsPrev, N);
    logUnnormalisedWeightsCurr.fill(logSumOfUnnormalisedWeightsPrev - std::log(N));
  }

  
  
  ////////////////////////////////////////////////////////////////////////////////
  // Conditional systematic resampling
  ////////////////////////////////////////////////////////////////////////////////
  
  /// \brief Performs conditional systematic resampling given a particular
  /// uniform random number.
  void conditionalSystematicBase
  (
    double u,              // uniform random variable
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    unsigned int& b,       // particle index of the distinguished particle
    const arma::colvec& w, // particle weights
    const unsigned int N,  // total number of offspring
    const unsigned int a   // parent index
  )
  { 
    arma::colvec Q = arma::cumsum(N * w);
    
    // determines the highest posstible stratum for each cumulated particle weight
    arma::uvec bins = arma::conv_to<arma::uvec>::from(arma::ceil(Q) - 1);
    bins.elem(arma::find(bins == N)).fill(N - 1);
    arma::colvec wAux; // stores the probabilities of sampling b from a particular stratum
      
    // First step: obtain the index of the distinguished particle at 
    // the current step given the entire history of the particle system
    // and in particular, given that the parent of this particle has index a.
    
    if (a == 0 || bins(a) == bins(a - 1)) 
    {
      b = bins(a);
    }
    else
    {
      wAux.zeros(N); 
      if (bins(a) > bins(a - 1) + 1)
      {
        // Assigning equal weights to all possible strata associated with b:
        wAux(arma::span(bins(a - 1) + 1, bins(a) - 1)).fill(1);
      }
        
      // Dealing with the last possibly stratum:
      wAux(bins(a)) = Q(a) - bins(a); 
      
      // Dealing with the first possible stratum:
      wAux(bins(a-1)) = bins(a - 1) - Q(a - 1) + 1;
      
      // Self-normalising the weights:
      wAux = arma::normalise(wAux, 1);

      // Determining the index of the distinguished particle:
      b = sampleInt(wAux);
    }
    
    // Second step: obtain a single uniform random variable compatible
    // with the indices a and b.
    
    double lb = 0; // lower bound for u
    double ub = 1; // upper bound for u
    if (a > 0 && b == bins(a - 1)) // i.e. if b falls into the first possible stratum
    {
      lb = Q(a-1) - bins(a - 1);
    } 
    if (b == bins(a)) // i.e. if b falls into the last possible stratum
    {
      ub = Q(a) - bins(a);
    } 
    
    u = lb + (ub - lb) * u;
    
    // Third step: perform standard systematic resampling given u.
    arma::colvec T = (arma::linspace(0, N - 1, N) + u);

    unsigned int i = 0;
    unsigned int j = 0;  
      
    while (j <= b) 
    {
      if (T(j) <= Q(i)) 
      {
        parentIndices(j) = i;
        ++j;
      } 
      else 
      {
        ++i;
      }
    }
    
    if (parentIndices(b) != a) 
    {
      std::cout << "Warning: conditional systematic resampling did not set the parent index of the conditioning path correctly!" << std::endl;
      std::cout << "Most likely cause: weight of the conditioning path is numerically zero!" << std::endl;
      //std::cout << "lb, ub, u, b, a, W(a), parentIndices(b): " << lb << ", " << ub << ", " << u << ", " << b << ", " << a << ", " << w(a) << ", " << parentIndices(b) << std::endl;          
      parentIndices(b) = a;
    }
    
    i = a;
    j = b+1;  
      
    while (j < N) 
    {
      /////////////////////////////////////////////////////////////////////////
      // This is only used to catch errors resulting from w(a) = 0 (which shouldn't normally happen)
      if (i == N) 
      {
        std::cout << "Warning: i = N!" << std::endl;
        std::cout << "########## Skipped resampling step! ##########" << std::endl;
        b = a;
        parentIndices = arma::linspace<arma::uvec>(0, N - 1, N);
        break;
      }
      /////////////////////////////////////////////////////////////////////////
      
      if (T(j) <= Q(i)) 
      {
        parentIndices(j) = i;
        ++j;
      } 
      else 
      {
        ++i;
      }
    }
    
  }
  /// \brief Performs conditional systematic resampling.
  void conditionalSystematic
  ( 
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    unsigned int& b,       // particle index of the distinguished particle
    const arma::colvec& w, // particle weights
    const unsigned int N,  // total number of offspring
    const unsigned int a   // parent index
  )
  { 
    conditionalSystematicBase(arma::randu(), parentIndices, b, w, N, a);
  }
  
    
  ////////////////////////////////////////////////////////////////////////////////
  // Optimal finite-state resampling from P. Fearnhead's PhD thesis (1998)
  ////////////////////////////////////////////////////////////////////////////////
  
  /// \brief Performs optimal finite-state resampling.
  void optimalFiniteState
  (
    arma::uvec& parentIndices, // stores the post-resampling particle labels
    arma::colvec& logUnnormalisedWeightsCurr, // the log-unnormalised post-resampling weights
    const arma::colvec& selfNormalisedWeightsPrev, // self-normalised particle weights
    const unsigned int M, // total desired number of offspring
    const arma::colvec& logUnnormalisedWeightsPrev, // the the log-unnormalised particle weights
    const double logSumOfUnnormalisedWeightsPrev // the logarithm of the sum of the unnormalised weights from the previous time step
  )        
  {

    parentIndices.set_size(M);
    logUnnormalisedWeightsCurr.zeros();
    
    unsigned int N = selfNormalisedWeightsPrev.size(); // number of from among which to choose the ancestors
    arma::uvec sortedIndices = arma::sort_index(selfNormalisedWeightsPrev, "descend"); // sorted indices, i.e. the $n$th entry of arma::sort_index(w, "descent") returns the position of the $n$th largest element in w, so that, e.g. the last element of arma::sort_index(w, "descent") returns the position of the smallest element in w
    arma::colvec q = selfNormalisedWeightsPrev.elem(sortedIndices); // self-normalised particle weights sorted in descending order
    arma::colvec logQ = arma::log(q);
    arma::colvec Q = arma::reverse(arma::cumsum(arma::reverse(q))); // "reverse" cumulative sum of sorted weights
   
    //////////////////////////////////////////
    //////////////////////////////////////////
    // alternative attempt to improve numerical stability
// // //     arma::colvec logQRev = arma::log(arma::reverse(q));
// // //     arma::colvec logQAlt(logQRev.size());
// // //     for (unsigned i = 0; i < logQAlt.size(); i++)
// // //     {
// // //       logQAlt(i) = sumExp(logQRev(arma::span(0, i)));
// // //     }
// // //     logQAlt = arma::reverse(logQAlt);
// // //     unsigned int kOldAlt = 1;
// // //     unsigned int kNewAlt = 0;
// // //     double logCAlt;
// // // 
// // //     while (kNewAlt != kOldAlt)
// // //     {
// // //       kOldAlt = kNewAlt;
// // //       logCAlt = std::log(M - kOldAlt) - logQAlt(kOldAlt);
// // //       kNewAlt = kOldAlt + arma::sum(logQ(arma::span(kOldAlt, N - 1)) > -logCAlt * arma::ones<arma::colvec>(N - kOldAlt));
// // // //       std::cout << "kNewAlt: " << kNewAlt << "; kOldAlt: " << kOldAlt << "; logCAlt: " << logCAlt << std::endl;
// // //     }
    //////////////////////////////////////////
    //////////////////////////////////////////
    
    unsigned int kOld = 1;
    unsigned int kNew = 0;
    double logC;
    while (kNew != kOld)
    {
      kOld = kNew;
      logC = std::log(M - kOld) - std::log(Q(kOld));
      kNew = kOld + arma::sum(logQ(arma::span(kOld, N - 1)) > -logC * arma::ones<arma::colvec>(N - kOld));
//       std::cout << "kNew: " << kNew << "; kOld: " << kOld << "; logC: " << logC << std::endl;
    }

    
    if (std::isfinite(logC)) 
    {
      unsigned int K = kNew;
      unsigned int L = M - K; // residual number of particles (which are not necessarily "kept")

      if (K > 0) {
        parentIndices(arma::span(0, K - 1)) = sortedIndices(arma::span(0, K - 1));
        logUnnormalisedWeightsCurr(arma::span(0, K - 1)) = logUnnormalisedWeightsPrev.elem(parentIndices(arma::span(0, K - 1)));
      }
      if (K < M) {
        arma::uvec indRes(L);
        arma::colvec residualWeights = arma::exp(normaliseExp(logQ(arma::span(K, N - 1))));
        resample::systematic(indRes, residualWeights, L);
        parentIndices(arma::span(K, M - 1)) = sortedIndices.elem(indRes + K * arma::ones<arma::uvec>(L)); 
      }
      
      for (unsigned int n = K; n < M; n++)
      {
        logUnnormalisedWeightsCurr(n) = logSumOfUnnormalisedWeightsPrev - logC;
      }
    }
    else // simply keep the particles with the largest weights 
    { 
      // Note that logC is infinity if there are fewer particles with non-zero weights at the
      // previous time step than ancestors that we wish to select. Of course, in this case, one should
      // simply keep all the particles with non-zero weights. This is what the following 
      // lines accomplish.
      
      arma::uvec indices = arma::sort_index(logUnnormalisedWeightsPrev, "descend");
      parentIndices = indices(arma::span(0, M - 1));
      logUnnormalisedWeightsCurr(arma::span(0, M - 1)) = logUnnormalisedWeightsPrev.elem(indices(arma::span(0, M - 1)));
      
      
//       std::cout << "arma::log(Q): " << arma::trans(arma::log(Q)) << std::endl;
//       std::cout << "logQAlt: " << logQAlt.t() << std::endl;
//       std::cout << "kNew: " << kNew << "; kOld: " << kOld << "; logC: " << logC << std::endl;
//       std::cout << "kNewAlt: " << kNewAlt << "; kOldAlt: " << kOldAlt << "; logCAlt: " << logCAlt << std::endl;
//       
//       std::cout << "WARNING: resorting to standard systematic resampling due to numerical issues!" << std::endl;
//       std::cout << "WARNING: numerical issues in optimal finite-state resampling! Keeping the particles with the largest weights instead!" << std::endl;
//       std::cout << "logC: " << logC << std::endl;
//       std::cout << "logUnnormalisedWeightsPrev.t(): " << logUnnormalisedWeightsPrev.t() << std::endl;
      // TODO: just keep the particles with the largest weights
      

      
//       for (unsigned m = 0; m < M; m++)
//       {
//         logUnnormalisedWeightsCurr(m) = logUnnormalisedWeightsPrev(indices(m));
//       }
      
      
//       std::cout << "logUnnormalisedWeightsCurr.t(): " <<  logUnnormalisedWeightsCurr.t() << std::endl;
      
//       resample::systematic
//       (
//         parentIndices, // stores the post-resampling particle labels
//         logUnnormalisedWeightsCurr, // the log-unnormalised post-resampling weights
//         selfNormalisedWeightsPrev, // self-normalised particle weights
//         M, // total desired number of offspring
//         logUnnormalisedWeightsPrev, // the the log-unnormalised particle weights
//         logSumOfUnnormalisedWeightsPrev
//       );
    }
  }
  
}
#endif
