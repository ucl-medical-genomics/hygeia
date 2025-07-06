import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
import os
from absl import flags, app
import sys
import pandas as pd
import hygeia.filter_and_smoother_algorithm as filter_and_smoother_algorithm
from hygeia.case_control_regime_model import CaseControlRegimeModel
from hygeia.case_control_proposal_mappings import CaseControlProposal
from pathlib import Path
import time

tfd = tfp.distributions
tfb = tfp.bijectors

dtype = tf.float32

flags.DEFINE_list("mu",
                     default='0.95,0.05,0.80,0.20,0.50,0.50',
                     help="mu of the beta distribution")
flags.DEFINE_list("sigma",
                     default='0.05,0.05,0.1,0.1,0.1,0.2886751',
                     help="sigma of the beta distribution")
flags.DEFINE_integer("minimum_duration",
                     default=3,
                     help="minimum duration between change points")
flags.DEFINE_float("omega_case",
                     default=.8,
                     help="omega parameter for duration of case group same for all regimes, must be in (0,1)")
flags.DEFINE_float("merge_log_prob",
                     default=np.log(.1),
                     help=" value for merge probability")
flags.DEFINE_float("split_prob",
                     default=.01,
                     help=" value for split probability")
flags.DEFINE_multi_integer("num_resampled_particles",
                     default=[50],
                     help="number M of particles that are resampled in the marginal filter")
flags.DEFINE_integer("num_samples_backward",
                      default=25,
                      help="number of particles used for backward sampling algorithm")
flags.DEFINE_bool("multinomial",
                      default=False,
                      help="multinomial or residual resampling")
flags.DEFINE_string("chrom",
                      default="22",
                      help="The chromosome to analyze (chr22, or 22, as per input file)")
flags.DEFINE_string(
    'results_dir',
    default=os.path.join(Path(os.getcwd()).parents[0],'test'),
    help="Directory to put the results.")
flags.DEFINE_string(
    'data_dir',
    default=os.path.join('data'),
    help="Directory of the read data.")
flags.DEFINE_string(
    'single_group_dir',
    default=os.path.join('test_data', 'single_group_results'),
    help="Directory of the single group estimation results.")
flags.DEFINE_integer("seed",
                     default=0,
                     help="seed used for sampling random variables")
flags.DEFINE_integer("batch",
                      default=0,
                      help="index of the selected chromosome segment")
flags.DEFINE_integer("segment_size",
                      default=100000,
                      help="size of the selected chromosome segment (in CpG sites)")
flags.DEFINE_integer("buffer_size",
                      default=5000,
                      help="size of the buffer segment (in CpG sites)")
FLAGS = flags.FLAGS


def get_estimated_control_group_param(chromosome, n_methylation_regimes):
  theta_data = pd.read_table(os.path.join(FLAGS.single_group_dir, 'theta_' + str(chromosome) + '.csv.gz'), sep = ',')
  #estimated_params = theta_data.to_numpy()[0]
  estimated_params = pd.to_numeric(theta_data['data']).to_numpy()
  p_softmax = np.zeros([n_methylation_regimes, n_methylation_regimes])
  i=0
  for r in range(n_methylation_regimes):
    for r1 in range(n_methylation_regimes):
      if r !=r1:
        p_softmax[r,r1] = np.math.exp(estimated_params[i])
        i+=1
    p_softmax[r, :] = p_softmax[r, :] / np.sum(p_softmax[r,:])
  omega_logit_control = estimated_params[-n_methylation_regimes:]
  return np.log(p_softmax), omega_logit_control


def main(argv):
  del argv  # unused
  tf.random.set_seed(FLAGS.seed)
  np.random.seed(FLAGS.seed)
  ## save flags to file
  fv = flags._flagvalues.FlagValues()
  key_flags = FLAGS.get_key_flags_for_module(sys.argv[0])
  s = '\n'.join(f.serialize() for f in key_flags)
  print('specified flags:\n{}'.format(s))
  path = os.path.join(FLAGS.results_dir,
                      'chrom_{}_{}'.format(
                        flags.FLAGS.chrom, flags.FLAGS.batch))
  if not os.path.exists(path):
    os.makedirs(path)
  flag_file = open(os.path.join(path, 'flags'+ str(FLAGS.seed)+'.txt'), "w")
  flag_file.write(s)
  flag_file.close()

  split_log_prob = np.log(FLAGS.split_prob)

  ###############################################################################
  ## PARAMETERS
  ###############################################################################
  ## The following parameters are fixed/assumed to be known, i.e.
  ## they are not estimated by the parameter-estimation scheme.
  ## They are also the same for both the case and control groups:

  # regime specific mean parameters beta laws
  mu_true = tf.Variable(np.array(FLAGS.mu, dtype = np.float32),
                        dtype = dtype, name = 'mu_true')
  # regime specific standard-deviation parameters beta laws
  sigma_true = tf.Variable(np.array(FLAGS.sigma, dtype = np.float32),
                            dtype = dtype, name = 'sigma_true')

  n_methylation_regimes = mu_true.shape[0]

  # transition matrix for the regimes of the control and case group
  # we parameterise the regime transitions with a softmax transform
  # we ignore the diagonal entries of the param matrix by setting it to -inf
  # before constructing the transition matrix by taking the softmax transform
  # we don't explicitly parameterise the regime transition matrix for the case group,
  # they are fixed to uniform
  # p_softmax = np.zeros([n_methylation_regimes,n_methylation_regimes])
  p_softmax, omega_logit_control = get_estimated_control_group_param(FLAGS.chrom,
                                                               n_methylation_regimes)
  P_softmax_control = tf.Variable(p_softmax, dtype = dtype, name = 'P_softmax_control')

  print(P_softmax_control)
  print(tf.math.exp(P_softmax_control))

  ## success-probability parameters of the regime-specific negative-binomial distributions
  # governing the function h() which determines the change-point probabilities
  omega_case = tf.Variable(FLAGS.omega_case * np.ones([n_methylation_regimes]), dtype = dtype)
  omega_logit_control = tf.Variable(omega_logit_control, dtype = dtype)
  def inv_logit(x):
   return tf.math.exp(x)/(1+tf.math.exp(x))
  omega_inv_logit_case = inv_logit(omega_case)
  omega_control = inv_logit(omega_logit_control)
  omega_inv_logit_control = inv_logit(omega_control)

  # minimum duration between change points
  minimum_duration = FLAGS.minimum_duration

  ## The parameters kappa are fixed/assumed to be known
  # number-of-failures parameters of the regime-specific negative-binomial distributions governing the function h()
  # which determines the change-point probabilities
  kappa_true_control = tf.Variable(2 * np.ones([n_methylation_regimes]),
                                   dtype = dtype, name = 'kappa_true_control')
  kappa_true_case = tf.Variable(2 * np.ones([n_methylation_regimes]),
                                dtype = dtype, name = 'kappa_true_case')

  ## transition matrix for the merged_regime indicator
  P_softmax_merged_ = np.array([[np.log(1. - np.exp(FLAGS.merge_log_prob)), FLAGS.merge_log_prob],
                                [split_log_prob, np.log(1. - np.exp(split_log_prob))]])
  P_softmax_merged = tf.Variable(P_softmax_merged_, dtype = dtype,
                                 name = 'P_softmax_merged')



  ###############################################################################
  ## Get Data
  ###############################################################################
  data_dir = FLAGS.data_dir

  #full chromosome data
  positions = pd.read_table(
    os.path.join(data_dir,
                 'positions_{}.txt.gz'.format(FLAGS.chrom)), sep = ',', header = None)
  n_total_reads_control = pd.read_table(
    os.path.join(data_dir,
                 'n_total_reads_control_{}.txt.gz'.format(FLAGS.chrom)), sep = ',', header = None)
  n_methylated_reads_control = pd.read_table(
    os.path.join(data_dir,
                 'n_methylated_reads_control_{}.txt.gz'.format(FLAGS.chrom)), sep = ',', header = None)
  n_total_reads_case = pd.read_table(
    os.path.join(data_dir,
                 'n_total_reads_case_{}.txt.gz'.format(FLAGS.chrom)), sep = ',', header = None)
  n_methylated_reads_case = pd.read_table(
    os.path.join(data_dir,
                 'n_methylated_reads_case_{}.txt.gz'.format(FLAGS.chrom)), sep = ',', header = None)


  #Select the segment corresponding to the batch index (if valid, otherwise exit)
  if FLAGS.batch * FLAGS.segment_size > positions.shape[0]:
    print("Batch index is too large for the chromosome")
    sys.exit(0)

  index= range(max(0, FLAGS.batch * FLAGS.segment_size - FLAGS.buffer_size),
                min((FLAGS.batch + 1) * FLAGS.segment_size + FLAGS.buffer_size, positions.shape[0]))

  observations_control = tf.convert_to_tensor(n_methylated_reads_control.iloc[index, :], dtype)
  observations_case = tf.convert_to_tensor(n_methylated_reads_case.iloc[index, :], dtype)
  n_total_reads_control = tf.convert_to_tensor(n_total_reads_control.iloc[index, :], dtype)
  n_total_reads_case = tf.convert_to_tensor(n_total_reads_case.iloc[index, :], dtype)
  positions = tf.convert_to_tensor(positions.iloc[index, :], tf.int64)

  observations = {'control': observations_control,
              'case': observations_case}
  assert (np.sum(n_total_reads_case < observations_case) == 0)
  assert (np.sum(n_total_reads_control < n_total_reads_control) == 0)

  if FLAGS.batch == 0:
    #for incomplete segment: positions.shape may be smaller than segment size
    return_index = range(0, min(positions.shape[0], FLAGS.segment_size))
  else:
    #for incomplete segment: positions.shape may be smaller than segment size
    return_index = range(FLAGS.buffer_size, min(positions.shape[0], FLAGS.buffer_size + FLAGS.segment_size))


  #########
  ## Build generative model
  #########

  generative_model = CaseControlRegimeModel(n_methylation_regimes, mu_true, sigma_true, P_softmax_control, P_softmax_merged,
                                                 omega_inv_logit_control, omega_inv_logit_case, minimum_duration,
                                                 kappa_true_control, kappa_true_case,
                                                 n_total_reads_control, n_total_reads_case)



  # define test functions for smoothing
  def test_function(state):
    merged_state = state['merged_state']#tf.squeeze(state['merged_state'], -1)
    regime_control = tf.gather(state['control_state'], 1, axis = len(state['control_state'].shape) - 1)
    regime_case = tf.gather(state['case_state'], 1, axis = len(state['control_state'].shape) - 1)
    return tf.cast(tf.concat([
        tf.stack([merged_state ==0],-1),
        tf.stack([regime_control == i for i in range(n_methylation_regimes)],-1),
        tf.stack([regime_case == i for i in range(n_methylation_regimes)],-1)], -1), tf.float32)

  #proposal mapping for marginal filters
  case_control_proposal_mapping = CaseControlProposal(n_methylation_regimes)

  #removed buffered states and save
  np.savetxt(os.path.join(path, 'observations_control.csv.gz'),
        observations_control.numpy().astype(np.int16)[return_index], delimiter=',')
  np.savetxt(os.path.join(path, 'observations_case.csv.gz'),
        observations_case.numpy().astype(np.int16)[return_index], delimiter=',')
  np.savetxt(os.path.join(path, 'n_total_reads_control.csv.gz'),
        n_total_reads_control.numpy().astype(np.int16)[return_index], delimiter=',')
  np.savetxt(os.path.join(path, 'n_total_reads_case.csv.gz'),
        n_total_reads_case.numpy().astype(np.int16)[return_index], delimiter=',')
  np.savetxt(os.path.join(path, 'positions.csv.gz'),
        positions.numpy()[return_index], delimiter=',')


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

    #removed buffered states and save
    np.save(os.path.join(path, 'optimal_backward_particles_merged_state_' + str(N) + '_' + str(FLAGS.seed)),
            backward_particles_merged_state.astype(np.int16)[return_index])
    np.save(os.path.join(path, 'optimal_backward_particles_control_state_' + str(N) + '_' + str(FLAGS.seed)),
            backward_particles_control_state.astype(np.int16)[return_index])
    np.save(os.path.join(path, 'optimal_backward_particles_case_state_' + str(N) + '_' + str(FLAGS.seed)),
            backward_particles_case_state.astype(np.int16)[return_index])

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
