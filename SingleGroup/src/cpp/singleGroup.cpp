#include <RcppArmadillo.h>
#include "time.h"
#include "singleGroup.h"

// TODO: disable range checks (by using at() for indexing elements of cubes/matrices/vectors)
// once the code is tested; 
// To that end, compile the code with ARMA_NO_DEBUG defined.

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


////////////////////////////////////////////////////////////////////////////////
// Samples from the prior of the parameters.
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::colvec sampleFromParameterPriorCpp
(
  const arma::colvec& vartheta // hyperparameters
)
{ 
  Rng rng;
  Model<ModelParameters, LatentVariable, Covariate, Observation> model(rng);
  model.setKnownParameters(vartheta);;
  return model.sampleFromParameterPrior();
}
////////////////////////////////////////////////////////////////////////////////
// Simulates latent variables and observations from the model.
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List simulateDataCpp
(
  const unsigned int nObservations, // number of observations
  const arma::colvec& vartheta, // hyperparameters 
  const arma::colvec& theta, // parameters
  const arma::uvec& genomicPositions,        // the genomic positions of each CpG site
  const arma::umat& nTotalReads              // the no of (methylated + unmethylated) reads for each sample at each CpG site
)
{
//   std::cout << "START: simulateDataCpp()" << std::endl;
  Rng rng;
  Model<ModelParameters, LatentVariable, Covariate, Observation> model(rng);
  model.setKnownParameters(vartheta);
  model.setUnknownParameters(theta);
  model.setCovariates(convertArmaUmatToCovariates(genomicPositions, nTotalReads));
  model.simulateData(nObservations);
  
//   std::cout << "END: simulateDataCpp()" << std::endl;
  return Rcpp::List::create
  (
    Rcpp::Named("latentVariables")  = convertLatentVariablesToArmaUmat(model.getLatentVariables()),
    Rcpp::Named("nMethylatedReads") = convertObservationsToArmaUmat(model.getObservations())
  ); 
}

////////////////////////////////////////////////////////////////////////////////
// Runs the online smoothing and parameter-estimation algorithm
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List runOnlineCombinedInferenceCpp
(
  const arma::colvec& vartheta,              // hyperparameters and other auxiliary model parameters
  const arma::colvec& thetaInit,             // initial value for the model-parameter vector theta
  const arma::uvec& genomicPositions,        // the genomic positions of each CpG site
  const arma::umat& nTotalReads,             // the no of (methylated + unmethylated) reads for each sample at each CpG site
  const arma::umat& nMethylatedReads,        // the no of methylated reads for each sample at each CpG site
  const unsigned int nParticlesMax,          // maximum number particles used by the SMC algorithm
  const unsigned int smcProposalType,        // type of particle proposals to use
  const unsigned int smcResampleType,        // type of resampling scheme to use
  const bool useOnlineMarginalSmoothing,     // determines if we should estimate the expectations of the regimes under the marginal smoothing distributions
  const double epsilon,                      // tuning parameter for the online-smoothing algorithm
  const bool useOnlineParameterEstimation,   // determines if the model parameters should be updated throughout the algorithm (otherwise they are fixed at their initial values)
  const bool normaliseGradients,             // should the gradients be normalised according to their $L_1$ norm?
  const bool useAdam,                        // should the ADAM optimiser be used (otherwise it is plain stochastic gradient ascent)
  const unsigned int nStepsWithoutParameterUpdate, // number of SMC steps without parameter updates
  const double learningRateExponent,         // the exponent governing the learning-rate decay
  const double learningRateFactor            // the factor governing the learning-rate decay
)
{

  // Starting the timer:
  clock_t t1,t2; // timing
  t1 = clock(); // start timer

  Rng rng;

  /////////////////////////////////////////////////////////////////////////////
  // Model class.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up Model class" << std::endl;
  
  Model<ModelParameters, LatentVariable, Covariate, Observation> model(rng);
  model.setKnownParameters(vartheta);
  model.setCovariates(convertArmaUmatToCovariates(genomicPositions, nTotalReads));
  model.setObservations(convertArmaUmatToObservations(nMethylatedReads));
  model.setUnknownParameters(thetaInit);

  /////////////////////////////////////////////////////////////////////////////
  // Class for running the SMC algorithm.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up Smc class" << std::endl;
  
  Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>> smc(rng, model);
  smc.setSmcResampleType(static_cast<SmcResampleType>(smcResampleType));
  smc.setSmcProposalType(static_cast<SmcProposalType>(smcProposalType));
  smc.setNParticlesMax(nParticlesMax);
  smc.setNRegimes(model.getModelParameters().getNMethylationRegimes()); 
  
  /////////////////////////////////////////////////////////////////////////////
  // Class for approximating expectations w.r.t. marginal smoothing distributions.
  /////////////////////////////////////////////////////////////////////////////

//   std::cout << "set up OnlineMarginalSmoothing class" << std::endl;
  
  OnlineMarginalSmoothing<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>>
  onlineMarginalSmoothing(rng, model, smc);
  
  if (useOnlineMarginalSmoothing)
  {
    onlineMarginalSmoothing.setEpsilon(epsilon);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Class for online parameter estimation.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up OnlineParameterEstimation class" << std::endl;
  
  GradientAscent gradientAscent;
  OnlineParameterEstimation<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>>
  onlineParameterEstimation(rng, model, smc, gradientAscent);
  
  if (useOnlineParameterEstimation)
  {
//     std::cout << "in CPP file: learningRateExponent: " << learningRateExponent << std::endl;
    gradientAscent.setNormaliseGradients(normaliseGradients);
    gradientAscent.setUseAdam(useAdam);
    gradientAscent.setLearningRateExponent(learningRateExponent);
    gradientAscent.setLearningRateFactor(learningRateFactor);
    onlineParameterEstimation.setTheta(thetaInit);
    onlineParameterEstimation.setNStepsWithoutParameterUpdate(nStepsWithoutParameterUpdate);
  }
   
  /////////////////////////////////////////////////////////////////////////////
  // Class for running the overall combined online inference algorithm.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up OnlineCombinedInference class" << std::endl;
  
  OnlineCombinedInference<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>>
  onlineCombinedInference(rng, model, smc, onlineMarginalSmoothing, onlineParameterEstimation);
  onlineCombinedInference.setNSteps(model.getNObservations());
  onlineCombinedInference.setUseOnlineMarginalSmoothing(useOnlineMarginalSmoothing);
  onlineCombinedInference.setUseOnlineParameterEstimation(useOnlineParameterEstimation);
  
  /////////////////////////////////////////////////////////////////////////////
  // Running the combined online inference algorithm.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "running the online inference algorithm" << std::endl;
  
  std::vector<arma::colvec> regimeProbabilityEstimates;
  std::vector<arma::colvec> thetaEstimates;
  onlineCombinedInference.run(regimeProbabilityEstimates, thetaEstimates);
  
  t2 = clock(); // stop timer 
  double cpuTime = (static_cast<double>(t2)-static_cast<double>(t1)) / CLOCKS_PER_SEC; // elapsed time in seconds
  
  return Rcpp::List::create(
    Rcpp::Named("regimeProbabilityEstimates") = regimeProbabilityEstimates,
    Rcpp::Named("thetaEstimates")  = thetaEstimates,
    Rcpp::Named("cpuTime") = cpuTime
  );

}
////////////////////////////////////////////////////////////////////////////////
// Empirically approximates (log-)normalising constant
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List approximateLogNormalisingConstantCpp
(
  const unsigned int nReplicates,            // number of independent estimates of the normalising constant
  const arma::colvec& vartheta,              // hyperparameters and other auxiliary model parameters
  const arma::colvec& theta,                 // value for the model-parameter vector theta
  const arma::uvec& genomicPositions,        // the genomic positions of each CpG site
  const arma::umat& nTotalReads,             // the no of (methylated + unmethylated) reads for each sample at each CpG site
  const arma::umat& nMethylatedReads,        // the no of methylated reads for each sample at each CpG site
  const unsigned int nParticlesMax,          // maximum number particles used by the SMC algorithm
  const unsigned int smcProposalType,        // type of proposal (kernel) for the latent states
  const unsigned int smcResampleType         // type of resampling scheme to use
)
{
  
  Rng rng;

  /////////////////////////////////////////////////////////////////////////////
  // Model class.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up Model class" << std::endl;
  
  Model<ModelParameters, LatentVariable, Covariate, Observation> model(rng);
  model.setKnownParameters(vartheta);
  model.setCovariates(convertArmaUmatToCovariates(genomicPositions, nTotalReads));
  model.setObservations(convertArmaUmatToObservations(nMethylatedReads));
  model.setUnknownParameters(theta);

  /////////////////////////////////////////////////////////////////////////////
  // Class for running the SMC algorithm.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up Smc class" << std::endl;
  
  Smc<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>> smc(rng, model);
  smc.setSmcResampleType(static_cast<SmcResampleType>(smcResampleType));
  smc.setSmcProposalType(static_cast<SmcProposalType>(smcProposalType));
  smc.setNParticlesMax(nParticlesMax);
  smc.setNRegimes(model.getModelParameters().getNMethylationRegimes()); 

  /////////////////////////////////////////////////////////////////////////////
  // Class for approximating expectations w.r.t. marginal smoothing distributions.
  /////////////////////////////////////////////////////////////////////////////

//   std::cout << "set up OnlineMarginalSmoothing class" << std::endl;
  
  OnlineMarginalSmoothing<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>>
  onlineMarginalSmoothing(rng, model, smc);

  /////////////////////////////////////////////////////////////////////////////
  // Class for online parameter estimation.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up OnlineParameterEstimation class" << std::endl;
  
  GradientAscent gradientAscent;
  OnlineParameterEstimation<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>>
  onlineParameterEstimation(rng, model, smc, gradientAscent);

  /////////////////////////////////////////////////////////////////////////////
  // Class for running the overall combined online inference algorithm.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "set up OnlineCombinedInference class" << std::endl;
  
  OnlineCombinedInference<ModelParameters, LatentVariable, Covariate, Observation, SmcParameters, Particle<LatentVariable>>
  onlineCombinedInference(rng, model, smc, onlineMarginalSmoothing, onlineParameterEstimation);
  onlineCombinedInference.setNSteps(model.getNObservations());
  onlineCombinedInference.setUseOnlineMarginalSmoothing(false);
  onlineCombinedInference.setUseOnlineParameterEstimation(false);
  
  /////////////////////////////////////////////////////////////////////////////
  // Running the combined online inference algorithm.
  /////////////////////////////////////////////////////////////////////////////
  
//   std::cout << "approximating the normalising constants" << std::endl;

  arma::colvec logNormalisingConstantEst(nReplicates);
  std::vector<arma::colvec> regimeProbabilityEstimates; // NOTE: not really used here
  std::vector<arma::colvec> thetaEstimates; // NOTE: not really used here
  for (unsigned int m = 0; m < nReplicates; m++)
  {
    onlineCombinedInference.run(regimeProbabilityEstimates, thetaEstimates);
    logNormalisingConstantEst(m) = smc.getLogSumOfUnnormalisedWeightsCurr();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("logNormalisingConstantEst") = logNormalisingConstantEst
  );
  
}
