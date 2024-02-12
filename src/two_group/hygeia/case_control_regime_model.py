import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
from tensorflow_probability.python.util import SeedStream
import collections
from tensorflow_probability.python.internal import prefer_static
from hygeia import case_control_distributions

tfd = tfp.distributions
tfb = tfp.bijectors

###
#helper functions
###
def get_logits_from_probs(probs):
  return tf.math.log(probs) - tf.math.log1p(-probs)


def get_beta_distribution_params(mu, sigma):
  nu = (mu * (1 - mu) / (sigma ** 2)) - 1.
  alpha = mu * nu
  beta = (1 - mu) * nu
  return alpha, beta

## For saving simulated data
GenerativeModelStepResults = collections.namedtuple(
    'GenerativeModelStepResults',
    ['latent_states',
     'observation',
     'transition_log_prob',
     'observation_log_prob'
     ])

GenerativeModelLoopResults = collections.namedtuple(
    'GenerativeModelLoopResults',
    ['step',
     'previous_step_results',
     'accumulated_step_results'
     ])

class CaseControlRegimeModel:
  """
  Case Control Model with Three Regimes
  They expose a transition_fn and an observation_fn
  """

  def __init__(self,
               n_methylation_regimes,
               mu_true,
               sigma_true,
               P_softmax_control,
               P_softmax_merged,
               omega_inv_logit_control,
               omega_inv_logit_case,
               minimum_duration,
               kappa_control,
               kappa_case,
               n_total_reads_control,
               n_total_reads_case,
               dtype = tf.float32):
    """Initializes the model."""
    self._n_methylation_regimes = n_methylation_regimes
    self._mu_true = mu_true
    self._sigma_true = sigma_true
    self._P_softmax_control = P_softmax_control
    self._P_softmax_merged = P_softmax_merged
    self._omega_inv_logit_control = omega_inv_logit_control
    self._omega_inv_logit_case = omega_inv_logit_case
    self._minimum_duration = minimum_duration
    self._kappa_control = kappa_control
    self._kappa_case = kappa_case
    self._n_total_reads_control = n_total_reads_control
    self._n_total_reads_case = n_total_reads_case
    self._dtype = dtype





  def _next_merged_state_probs(self, step, previous_merged_state, previous_duration_control,
                              previous_duration_case):
    #tf.print('P merged ', tf.nn.softmax(self._P_softmax_merged, axis=1))
    probs = tf.nn.softmax(self._P_softmax_merged)
    probs = prefer_static.cond(step == 0, lambda: tf.constant([[0.,1.],[0.,1.]]), lambda: probs)
    return tf.where(tf.expand_dims(tf.minimum(previous_duration_case, previous_duration_control)>=self._minimum_duration, -1),
                    tf.gather(probs, previous_merged_state),
                    tf.gather(tf.eye(2), previous_merged_state))


  def _next_regime_control_probs(self, previous_regime_control):
    control_regime_transition_probs = tf.gather(
        tf.nn.softmax(tf.linalg.set_diag(self._P_softmax_control, -np.math.inf* tf.ones_like(
          self._P_softmax_control[-1])), axis = -1), previous_regime_control)
    return control_regime_transition_probs


  def transition_fn(self, step, state):
    """Define the transition distribution of the states"""

    previous_merged_state = prefer_static.cond(state['merged_state'].shape[-1]==1,
                                               lambda: tf.squeeze(state['merged_state'], -1),
                                               lambda: state['merged_state'])
    previous_duration_control = tf.gather(state['control_state'], 0, axis=len(state['control_state'].shape)-1)
    previous_regime_control = tf.gather(state['control_state'], 1, axis=len(state['control_state'].shape)-1)
    previous_duration_case = tf.gather(state['case_state'], 0, axis=len(state['case_state'].shape)-1)
    previous_regime_case =  tf.gather(state['case_state'], 1, axis=len(state['case_state'].shape)-1)




    def next_duration_control_rho(self, previous_duration_control):
      previous_duration_control = tf.cast(previous_duration_control, self._dtype)
      neg_bin_logits_control = tf.gather((self._omega_inv_logit_control),
                                     previous_regime_control)
      neg_bin_counts_control = tf.gather(self._kappa_control,
                                       previous_regime_control)
      h_control_dist = tfd.NegativeBinomial(total_count=neg_bin_counts_control,
                                          probs= tf.math.log(
                                              neg_bin_logits_control/(1-neg_bin_logits_control)))
      log_h_control = tf.where(previous_duration_control >= self._minimum_duration,
                           h_control_dist.log_prob(previous_duration_control - self._minimum_duration),
                           -np.Inf)
      log_h_survival_control = tf.where(previous_duration_control > self._minimum_duration,
                               h_control_dist.log_survival_function(previous_duration_control - self._minimum_duration - 1.),
                               0.)
      rho_control = tf.where(log_h_control == -np.Inf,
                             0.,
                             tf.math.exp(log_h_control - log_h_survival_control)
                             )
      #if log_h_survival_control or log_h_control become inf, set rho to fixed value
      fixed_value_inf = .1
      rho_control = tf.where(tf.math.is_finite(rho_control),
                             rho_control,
                             fixed_value_inf)

      rho_control = prefer_static.cond(step == 0, lambda: 1., lambda: rho_control)

      return rho_control


    def next_duration_case_rho(self, previous_duration_case):
      previous_duration_case = tf.cast(previous_duration_case, self._dtype)
      neg_bin_logits_case = tf.gather(self._omega_inv_logit_case,
                                    previous_regime_case)
      neg_bin_counts_case = tf.gather(self._kappa_case,
                                       previous_regime_case)
      h_case_dist = tfd.NegativeBinomial(total_count=neg_bin_counts_case,
                                         probs = tf.math.log(
                                             neg_bin_logits_case / (1 - neg_bin_logits_case)))
      log_h_case = tf.where(previous_duration_case >= self._minimum_duration,
                         h_case_dist.log_prob(previous_duration_case - self._minimum_duration),
                         -np.Inf)
      log_h_survival_case = tf.where(previous_duration_case > self._minimum_duration,
                             h_case_dist.log_survival_function(previous_duration_case - self._minimum_duration - 1.),
                             0.)
      rho_case = tf.where(log_h_case == -np.Inf,
                           0.,
                           tf.math.exp(log_h_case - log_h_survival_case)
                           )
      #if log_h_survival_case or log_h_case become inf, set rho to fixed value
      fixed_value_inf = .1
      rho_case = tf.where(tf.math.is_finite(rho_case),
                             rho_case,
                             fixed_value_inf)

      rho_case = prefer_static.cond(step == 0, lambda: 1., lambda: rho_case)

      return rho_case


    transition_dist = tfd.JointDistributionNamed(dict(
     merged_state = tfd.Independent(tfd.Categorical(probs = self._next_merged_state_probs(step, previous_merged_state,
                                                                    previous_duration_control,
                                                                    previous_duration_case)),
                                    len(previous_duration_control.shape)-1),
     control_state = lambda merged_state: case_control_distributions.ControlStateTransition(
       current_merged_state = merged_state,
       previous_duration_control = previous_duration_control,
       previous_regime_control = previous_regime_control,
       rho = next_duration_control_rho(self, previous_duration_control),
       regime_probs = self._next_regime_control_probs(previous_regime_control)),
     case_state = lambda merged_state, control_state: case_control_distributions.CaseStateTransition(
       previous_duration_case = previous_duration_case,
       previous_regime_case = previous_regime_case,
       previous_merged_state = previous_merged_state,
       current_merged_state = merged_state,
       current_control_states = control_state,
       rho = next_duration_case_rho(self, previous_duration_case),
       n_regimes = self._n_methylation_regimes)
      ))


    return transition_dist



  def observation_fn(self, step, state):
    """Defines the observation distribution of the model"""
    #get the number of reads
    n_reads_control = self._n_total_reads_control[tf.cast(step, tf.int32), :]
    n_reads_case = self._n_total_reads_case[tf.cast(step, tf.int32), :]
    r_control = tf.gather(state['control_state'], 1, axis = len(state['control_state'].shape) - 1)
    alpha_control, beta_control = get_beta_distribution_params(
      tf.gather(self._mu_true, r_control),
      tf.gather(self._sigma_true, r_control)
    )
    r_case = tf.gather(state['case_state'], 1, axis = len(state['case_state'].shape) - 1)
    alpha_case, beta_case = get_beta_distribution_params(
      tf.gather(self._mu_true, r_case),
      tf.gather(self._sigma_true, r_case)
    )
    #expand dimension
    n_reads_control = n_reads_control*tf.ones_like(tf.expand_dims(alpha_control,1))
    alpha_control = tf.tile(tf.expand_dims(alpha_control,-1),[1, tf.shape(n_reads_control)[-1]])
    beta_control = tf.tile(tf.expand_dims(beta_control,-1),[1, tf.shape(n_reads_control)[-1]])
    n_reads_case = n_reads_case*tf.ones_like(tf.expand_dims(alpha_case,1))
    alpha_case = tf.tile(tf.expand_dims(alpha_case,-1),[1, tf.shape(n_reads_case)[-1]])
    beta_case = tf.tile(tf.expand_dims(beta_case, -1), [1, tf.shape(n_reads_case)[-1]])


    control_observation_dist = tfd.BetaBinomial(total_count = n_reads_control,
                                            concentration1 = alpha_control,
                                            concentration0 = beta_control)
    case_observation_dist = tfd.BetaBinomial(total_count = n_reads_case,
                                            concentration1 = alpha_case,
                                            concentration0 = beta_case)
    observation_dist = tfd.JointDistributionNamed(dict(
      control = tfd.Independent(control_observation_dist,1),
      case = tfd.Independent(case_observation_dist,1)))

    return observation_dist


  def intitial_state_dist(self, batch_size = 1):
    """Defines the initial density of the latent states"""
    phantom_state_batch = tfd.JointDistributionNamed(dict(
      merged_state = tfd.Independent(tfd.Categorical(probs = [[0., 1.]]), 1),
      control_state = case_control_distributions.InitialControlState(self._n_methylation_regimes,
                                                                     init_duration = 0),
      case_state = lambda control_state: tfd.Deterministic(control_state)
    )).sample(batch_size)

    initial_state_dist = self.transition_fn(0, phantom_state_batch)
    return initial_state_dist

  #@tf.function
  def simulate(self, n_observations, seed):
    """Simulate from the model (given the number of reads)"""
    seed = SeedStream(seed, 'sample_generative_model')
    with tf.name_scope('sample_generative_model'):

      def sample_model_transitions(step,
                                   previous_step_results,
                                   accumulated_step_results,
                                   seed = seed):


        with tf.name_scope('generative_model'):
          seed = SeedStream(seed, 'generative_model')
          transition_kernel = self.transition_fn(
            step, previous_step_results.latent_states)
          next_latent_states = transition_kernel.sample(seed = seed)
          #tf.print(next_latent_states)
          observation_dist = self.observation_fn(step, next_latent_states)
          observation = observation_dist.sample(seed = seed)
          transition_log_prob = transition_kernel.log_prob(
             next_latent_states
          )
          observation_log_prob = observation_dist.log_prob(
             observation
          )

          new_step_results = GenerativeModelStepResults(
            latent_states = next_latent_states,
            observation = observation,
            transition_log_prob = transition_log_prob,#f.zeros([]),#,,
            observation_log_prob = observation_log_prob
          )


        # Update accumulated results with new step results.
        accumulated_step_results = tf.nest.map_structure(
          lambda x, y: x.write(step, y),
          accumulated_step_results, new_step_results)

        return GenerativeModelLoopResults(
          step = step + 1,
          previous_step_results = new_step_results,
          accumulated_step_results = accumulated_step_results)


      # Sample the first state
      initial_state_dist = self.intitial_state_dist(1)
      initial_states = initial_state_dist.sample( seed = seed)
      observation_dist = self.observation_fn(0, initial_states)
      observation = observation_dist.sample(seed = seed)
      initial_step_results = GenerativeModelStepResults(
        latent_states = initial_states,
        observation = observation,
        transition_log_prob = initial_state_dist.log_prob(initial_states),# observation_dist.log_prob(observation),#,
        observation_log_prob = observation_dist.log_prob(observation)
      )
      # Create arrays to store states and observations and write
      # their initial values.
      step_results_arrays = tf.nest.map_structure(
        lambda x: tf.TensorArray(dtype = x.dtype, size = n_observations).write(0, x),
        initial_step_results)

      loop_results = tf.while_loop(
          cond = lambda step, *_: step < n_observations,
          body = sample_model_transitions,
          loop_vars = GenerativeModelLoopResults(
            step = 1,
            previous_step_results = initial_step_results,
            accumulated_step_results = step_results_arrays)
      )

      results = tf.nest.map_structure(lambda ta: ta.stack(),
                                      loop_results.accumulated_step_results)

      return results

