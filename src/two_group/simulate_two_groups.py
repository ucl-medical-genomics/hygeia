########
#Script for generating data for the simulation study using the generative model
#with different methylation regimes
########


import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
import os
from absl import flags
from absl import app
import sys
import time
sys.path.insert(0, '..')
from hygeia.case_control_regime_model import CaseControlRegimeModel
import hygeia.filter_and_smoother_algorithm as filter_and_smoother_algorithm
from hygeia.case_control_proposal_mappings import CaseControlProposal
tfd = tfp.distributions
tfb = tfp.bijectors

dtype = tf.float32
flags.DEFINE_integer("regimes_config",
                     default=8,
                     help="regimes_configuration for generative model")
flags.DEFINE_integer("minimum_duration_min",
                     default=5,
                     help="minimum duration between change points minimum possible value")
flags.DEFINE_integer("minimum_duration_max",
                     default=5,
                     help="minimum duration between change points maximum possible value")
flags.DEFINE_float("omega_case_min",
                     default=.8,
                     help="omega parameter for duration of case and control group same for all regimes, must be in (0,1)")
flags.DEFINE_float("omega_case_max",
                     default=.8,
                     help="omega parameter for duration of case and control group same for all regimes, must be in (0,1)")
flags.DEFINE_integer("num_observations",
                     default=10000,
                     help="number of observations")
flags.DEFINE_integer("num_samples_control",
                     default=1,
                     help="number of samples in the control group")
flags.DEFINE_integer("num_samples_case",
                     default=1,
                     help="number of samples in the case group")
flags.DEFINE_float("average_reads",
                     default=100,
                     help="average read deapth for simulating data")
flags.DEFINE_float("prob_no_data",
                     default=.0,
                     help="probability of no data at a given site")
flags.DEFINE_float("merge_log_prob_min",
                     default=np.log(.1),
                     help="minimum value for merge probability")
flags.DEFINE_float("merge_log_prob_max",
                     default=np.log(.1),
                     help="maximum value for merge probability")
flags.DEFINE_float("split_prob_min",
                     default=.01,
                     help="minimum value for split probability")
flags.DEFINE_float("split_prob_max",
                     default=.01,
                     help="maximum value for split probability")
flags.DEFINE_float("kappa_min",
                     default=2.,
                     help="kappa value of regime duration")
flags.DEFINE_float("kappa_max",
                     default=2.,
                     help="kappa value of regime duration")
flags.DEFINE_bool("uniform_regime_transition",
                     default=True,
                     help="if regime transition matrix should be uniform")
flags.DEFINE_string(
    'model_dir',
    default=os.path.join(os.getcwd(),'simulation_study_known_model'),
    help="Directory to put the results.")
flags.DEFINE_integer("seed",
                     default=0,
                     help="seed used for sampling random variables")
flags.DEFINE_multi_integer("num_resampled_particles",
                     default=[10, 20],
                     help="number M of particles that are resampled in the marginal filter")
flags.DEFINE_integer("num_samples_backward",
                      default=25,
                      help="number of particles used for backward sampling algorithm")
flags.DEFINE_bool("multinomial",
                      default=False,
                      help="multinomial or residual resampling")


FLAGS = flags.FLAGS


def main(argv):
  del argv  # unused
  tf.random.set_seed(FLAGS.seed)
  np.random.seed(FLAGS.seed)


  ## save flags to file
  fv = flags._flagvalues.FlagValues()
  key_flags = FLAGS.get_key_flags_for_module(sys.argv[0])
  s = '\n'.join(f.serialize() for f in key_flags)
  print('specified flags:\n{}'.format(s))
  path = os.path.join(FLAGS.model_dir)
  if not os.path.exists(path):
      os.makedirs(path)
  flag_file = open(os.path.join(path, 'flags.txt'), "w")
  flag_file.write(s)
  flag_file.close()
  #####
  #sample parameters for dataset
  #####

  ###############################################################################
  ## PARAMETERS
  ###############################################################################
  ## The following parameters are fixed/assumed to be known, i.e.
  ## they are not estimated by the parameter-estimation scheme.
  ## They are also the same for both the case and control groups:
  regimes_config = FLAGS.regimes_config
  if regimes_config == 1:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.85, 0.15, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.10, 0.10, 0.08, 0.08, 0.15, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 2:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.85, 0.15, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.10, 0.10, 0.04, 0.04, 0.15, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 3:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.85, 0.15, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.15, 0.15, 0.08, 0.08, 0.15, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 4:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.85, 0.15, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.15, 0.15, 0.04, 0.04, 0.15, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 5:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.99, 0.01, 0.80, 0.20, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.05, 0.05, 0.05, 0.05, 0.1, 1 / np.sqrt(6)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 6:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.99, 0.01, 0.80, 0.20, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.05, 0.05, 0.1, 0.1, 0.1, 1 / np.sqrt(6)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 7:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.85, 0.15, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.05, 0.05, 0.05, 0.05, 0.05, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 8:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.8, 0.2, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.05, 0.05, 0.1, 0.1, 0.1, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 9:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.75, 0.25, 0.50, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.1, 0.1, 0.1, 0.1, 0.1, 1 / np.sqrt(12)]),
                              dtype = dtype, name = 'sigma_true')
  elif regimes_config == 10:
    # regime specific mean parameters beta laws
    mu_true = tf.Variable([0.95, 0.05, 0.8, 0.2, 0.50],
                          dtype = dtype, name = 'mu_true')
    # regime specific standard-deviation parameters beta laws
    sigma_true = tf.Variable(np.array([0.05, 0.05, 0.1, 0.1, 0.1]),
                              dtype = dtype, name = 'sigma_true')

  n_methylation_regimes = mu_true.shape[0]

  omega_case = np.random.uniform(low=FLAGS.omega_case_min, high=FLAGS.omega_case_max,
                                  size = n_methylation_regimes)
  minimum_duration = np.random.choice(range(FLAGS.minimum_duration_min,
                                            FLAGS.minimum_duration_max+1))
  log_q_merge = np.random.uniform(low=FLAGS.merge_log_prob_min, high=FLAGS.merge_log_prob_max)
  log_q_split = np.random.uniform(low=np.log(FLAGS.split_prob_min), high=np.log(FLAGS.split_prob_max))

  kappa = np.random.uniform(low=np.log(FLAGS.kappa_min), high=np.log(FLAGS.kappa_max))


  P_softmax_merged_ = np.array([[np.log(1.-np.exp(log_q_merge)), log_q_merge],
                                [log_q_split, np.log(1.-np.exp(log_q_split))]])
  P_softmax_merged = tf.Variable(P_softmax_merged_, dtype = dtype,
                                  name = 'P_softmax_merged')

  # transition matrix for the regimes of the control and case group
  # we parameterise the regime transitions with a softmax transform
  # we ignore the diagonal entries of the param matrix by setting it to -inf
  # before constructing the transition matrix by taking the softmax transform
  # we don't explicitly parameterise the regime transition matrix for the case group,
  # they are fixed to uniform
  if FLAGS.uniform_regime_transition:
    p = 1./(n_methylation_regimes-1)*np.ones([n_methylation_regimes,n_methylation_regimes])
    np.fill_diagonal(p, np.zeros([n_methylation_regimes]))
    p_softmax = np.log(p)
  else:
    alpha = 2/3 * np.ones([n_methylation_regimes])
    p = tfp.distributions.Dirichlet(concentration = alpha).sample(n_methylation_regimes).numpy()
    np.fill_diagonal(p, np.zeros([n_methylation_regimes]))
    p_softmax = np.log(p)
    print(p_softmax)
  omega_control = np.random.uniform(low=FLAGS.omega_case_min, high=FLAGS.omega_case_max,
                                  size = n_methylation_regimes)
  P_softmax_control = tf.Variable(p_softmax, dtype = dtype, name = 'P_softmax_control')
  ## success-probability parameters of the regime-specific negative-binomial distributions
  # governing the function h() which determines the change-point probabilities
  omega_inv_logit_control = tf.Variable(np.exp(omega_control) / (1 + np.exp(omega_control)),
                                          name = 'omega_logit_control', dtype = dtype)
  omega_inv_logit_case = tf.Variable(np.exp(omega_case) / (1 + np.exp(omega_case)),
                                      name = 'omega_logit_case', dtype = dtype)


  ## The parameters kappa are fixed/assumed to be known
  # number-of-failures parameters of the regime-specific negative-binomial distributions governing the function h()
  # which determines the change-point probabilities
  kappa_true_control = tf.Variable(kappa * np.ones([n_methylation_regimes]),
                                    dtype = dtype, name = 'kappa_true_control')
  kappa_true_case = tf.Variable(kappa * np.ones([n_methylation_regimes]),
                                dtype = dtype, name = 'kappa_true_case')

  ## number of total reads for the case and control groups
  n_samples_control = FLAGS.num_samples_control
  n_samples_case = FLAGS.num_samples_case
  n_observations = FLAGS.num_observations
  n_total_reads_control = tf.random.poisson(lam=FLAGS.average_reads, shape=(n_observations,n_samples_control))
  n_total_reads_case = tf.random.poisson(lam=FLAGS.average_reads, shape=(n_observations,n_samples_case))
  #adjust for missing data
  n_total_reads_control = n_total_reads_control * np.random.binomial(
      n = 1, p = 1- FLAGS.prob_no_data, size = (n_observations,n_samples_control))
  n_total_reads_case= n_total_reads_case * np.random.binomial(
      n = 1, p = 1- FLAGS.prob_no_data, size = (n_observations,n_samples_case))
  n_observations = n_total_reads_case.shape[0]

  #########
  ## Build generative model
  #########
  generative_model = CaseControlRegimeModel(n_methylation_regimes, mu_true, sigma_true, P_softmax_control,
                                            P_softmax_merged,
                                            omega_inv_logit_control, omega_inv_logit_case, minimum_duration,
                                            kappa_true_control, kappa_true_case,
                                            n_total_reads_control, n_total_reads_case)


  #########
  ## Simulate from generative model
  #########

  @tf.function
  def simulate_from_model():
    return generative_model.simulate(n_observations, FLAGS.seed)

  simulated_results = simulate_from_model()
  true_latent_states = simulated_results.latent_states
  true_merged_states = true_latent_states['merged_state'].numpy()
  true_control_states = true_latent_states['control_state'].numpy()[:,0, :]
  true_case_states = true_latent_states['case_state'].numpy()[:, 0, :]
  observations_control = simulated_results.observation['control'].numpy()[:,0, :]
  observations_case = simulated_results.observation['case'].numpy()[:,0, :]

  observations = {'control': observations_control,
              'case': observations_case}

  #save simulated data
  np.savetxt(os.path.join(path, 'true_merged_states.csv'),
      true_merged_states, delimiter=',')
  np.savetxt(os.path.join(path, 'true_control_states.csv'),
              true_control_states, delimiter = ',')
  np.savetxt(os.path.join(path, 'true_case_states.csv'),
              true_case_states, delimiter = ',')

  np.savetxt(os.path.join(path, 'observations_control.csv'),
        observations_control, delimiter=',')
  np.savetxt(os.path.join(path, 'observations_case.csv'),
        observations_case, delimiter=',')
  np.savetxt(os.path.join(path, 'n_total_reads_control.csv'),
        n_total_reads_control, delimiter=',')
  np.savetxt(os.path.join(path, 'n_total_reads_case.csv'),
        n_total_reads_case, delimiter=',')


  # define test functions for smoothing
  def test_function(state):
    merged_state = state['merged_state']#tf.squeeze(state['merged_state'], -1)
    regime_control = tf.gather(state['control_state'], 1, axis = len(state['control_state'].shape) - 1)
    regime_case = tf.gather(state['case_state'], 1, axis = len(state['control_state'].shape) - 1)
    return tf.cast(tf.concat([
        tf.stack([merged_state ==0],-1),
        tf.stack([regime_control == i for i in range(n_methylation_regimes)],-1),
        tf.stack([regime_case == i for i in range(n_methylation_regimes)],-1)], -1), tf.float32)

  #and evaluate test function on the true latent states
  true_functionals = tf.squeeze(test_function(true_latent_states), 1)
  np.savetxt(os.path.join(path, 'true_functionals.csv'),
              true_functionals, delimiter = ',')


  #proposal mapping for marginal filters
  case_control_proposal_mapping = CaseControlProposal(n_methylation_regimes)
  #####
  # run non-marginal filter with optimal resampling scheme
  #####
  @tf.function
  def run_non_marginal_pf_optimal(num_resampled_ancestors = 100):
    num_particles = num_resampled_ancestors * (2 * n_methylation_regimes + n_methylation_regimes ** 2)

    return filter_and_smoother_algorithm.run(
      observations = observations,
      initial_state_prior = generative_model.intitial_state_dist(batch_size = 1),
      transition_fn = generative_model.transition_fn,
      observation_fn = generative_model.observation_fn,
      proposal_fn = case_control_proposal_mapping.proposal_fn_standard_filter,
      initial_proposal = case_control_proposal_mapping.initial_proposal_fn_standard_filter,
      num_particles = num_particles,
      num_resampled_ancestors = num_resampled_ancestors,
      optimal_resampling = True,
      multinomial_resampling = FLAGS.multinomial,
      num_simulations = FLAGS.num_samples_backward)


  log_normalizing_constants_optimal = {}
  optimal_time = {}
  optimal_time_backward = {}

  for M in FLAGS.num_resampled_particles:
    print(M)
    N = M * (2 * n_methylation_regimes + n_methylation_regimes ** 2)
    start_time = time.time()
    optimal_filter_backward_results, unnormalized_optimal_log_weight = run_non_marginal_pf_optimal(M)
    optimal_time[N] = time.time() - start_time
    log_normalizing_constants_optimal[N] = tf.reduce_logsumexp(
      unnormalized_optimal_log_weight).numpy()

    backward_simulation_particles = optimal_filter_backward_results.particle

    marginal_functionals = tf.reduce_mean(test_function(backward_simulation_particles), axis = [1]).numpy()
    split_probs = marginal_functionals[:, 0]
    regime_probs = marginal_functionals[:, 1:]

    # save backward trajctories
    backward_particles_merged_state = backward_simulation_particles['merged_state'].numpy()
    backward_particles_control_state = backward_simulation_particles['control_state'].numpy()
    backward_particles_case_state = backward_simulation_particles['case_state'].numpy()

    np.save(os.path.join(path, 'optimal_backward_particles_merged_state_' + str(N) + '_' + str(FLAGS.seed)),
            backward_particles_merged_state.astype(np.int16))
    np.save(os.path.join(path, 'optimal_backward_particles_control_state_' + str(N) + '_' + str(FLAGS.seed)),
            backward_particles_control_state.astype(np.int16))
    np.save(os.path.join(path, 'optimal_backward_particles_case_state_' + str(N) + '_' + str(FLAGS.seed)),
            backward_particles_case_state.astype(np.int16))

    np.save(os.path.join(path, 'optimal_split_probs_' + str(N) + '_' + str(FLAGS.seed)),
            split_probs)
    np.save(os.path.join(path, 'optimal_regime_probs_' + str(N) + '_' + str(FLAGS.seed)),
            regime_probs)

  # save data
  with open(os.path.join(path, 'log_normalizing_constants_optimal_' + str(FLAGS.seed) + '.txt'), "w") as f:
    print(log_normalizing_constants_optimal, file = f)
  with open(os.path.join(path, 'optimal_time_' + str(FLAGS.seed) + '.txt'), "w") as f:
    print(optimal_time, file = f)
  with open(os.path.join(path, 'optimal_time_backward_' + str(FLAGS.seed) + '.txt'), "w") as f:
    print(optimal_time_backward, file = f)


if __name__ == '__main__':
  app.run(main)
