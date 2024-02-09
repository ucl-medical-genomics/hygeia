"""Backward sampling algorithm (given filter results)"""

import collections
import numpy as np
import tensorflow as tf
from tensorflow_probability.python.util import SeedStream
import tensorflow_probability as tfp
from hygeia import smoothing_functions

tfd = tfp.distributions

BackwardSimulationResults = collections.namedtuple(
    'BackwardSimulationResults',
    ['step',
     'particle'])


def backward_simulation(unnormalized_log_weights,
                      particles,
                      transition_fn,
                      #initial_state_prior,
                      num_simulations = 1,
                      dtype = tf.float64,
                      seed = None,
                      name = None):  # pylint: disable=g-doc-args
    """Run backward simulation algorihm.
    """
    seed = SeedStream(seed, 'backward_sampling')
    with tf.name_scope(name or 'backward_sampling'):
      #sample the last particle
      step = tf.shape(unnormalized_log_weights)[0] -1
      index = tfd.Categorical(logits = unnormalized_log_weights[step]).sample(num_simulations, seed=seed)
      particle = tf.nest.map_structure(
        lambda x: tf.gather(tf.gather(x,step), index, axis=0), particles)
      initial_backward_simulation_result = BackwardSimulationResults(step=step, particle = particle)

      backward_simulation_results = tf.nest.map_structure(
        lambda x: tf.TensorArray(dtype=x.dtype, size=0, dynamic_size=True, infer_shape=False).write(step, x),
        initial_backward_simulation_result)


      def _loop_backward_sampling(step, backward_simulation_results, next_particle):
        #compute weight from backward kernel
        step = step-1
        current_particles = tf.nest.map_structure(lambda x: tf.gather(x,step), particles)
        current_unnormalized_log_weights = unnormalized_log_weights[step]
        #current_particles_to_keep = tf.nest.map_structure(
        #  lambda x: tf.boolean_mask(x, current_unnormalized_log_weights>-np.Inf, axis=0),
        #                      current_particles)
        current_indices_to_keep = tf.boolean_mask(tf.range(0,tf.shape(current_unnormalized_log_weights)[0]),
                                                  current_unnormalized_log_weights>-np.Inf, axis=0)
        current_particles_to_keep = tf.nest.map_structure(
          lambda x: tf.gather(x, current_indices_to_keep, axis = 0),
          current_particles)
        expanded_particles = tf.nest.map_structure(lambda x: tf.expand_dims(x, 1),
                                                   next_particle)
        prev_num_particles = tf.shape(current_indices_to_keep)[0]
        next_particles_matrix = tf.nest.map_structure(
          lambda x: tf.tile(x, [1, prev_num_particles] + [1 for _ in range(len(x.shape[2:]))]), expanded_particles)

        transition_log_prob_matrix = transition_fn(step+1, current_particles_to_keep).log_prob(next_particles_matrix)

        # transition_log_prob_matrix = prefer_static.cond(
        #   step==0,
        #   lambda: initial_state_prior().log_prob(current_particles_to_keep),
        #   lambda: transition_fn(step, current_particles_to_keep).log_prob(next_particles_matrix)
        # )
        transition_log_prob_matrix = tf.ensure_shape(transition_log_prob_matrix,
                                                     (tf.ones([num_simulations, prev_num_particles])).shape)

        log_backward_kernel = smoothing_functions.compute_log_backward_kernel_from_transition_matrix(
          tf.boolean_mask(current_unnormalized_log_weights, current_unnormalized_log_weights>-np.Inf, axis=0),
          transition_log_prob_matrix, dtype = dtype)
        index = tfd.Categorical(logits = log_backward_kernel).sample(seed = seed)
        next_particle = tf.nest.map_structure(
          lambda x: tf.gather(x, index, axis=0), current_particles_to_keep)
        backward_simulation_results = tf.nest.map_structure(
          lambda x,y: x.write(step, y),
          backward_simulation_results, BackwardSimulationResults(step=step, particle = next_particle))
        return step, backward_simulation_results, next_particle


      #Loop backwards
      step, backward_simulation_results, next_particle = tf.while_loop(
          cond=lambda step, *_: step >= 1,
          body=_loop_backward_sampling,
          loop_vars=[step, backward_simulation_results, particle])

      results = tf.nest.map_structure(lambda ta: ta.stack(),
                                      backward_simulation_results)

      return results