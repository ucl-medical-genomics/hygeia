"""Particle filtering with deterministic state proposals."""

import collections
import numpy as np
import tensorflow as tf
from tensorflow_probability.python.internal import prefer_static
import tensorflow_probability as tfp
from hygeia import resampling_functions

tfd = tfp.distributions

dtype = tf.float64


ParticleFilterStepResults = collections.namedtuple(
    'ParticleFilterStepResults',
    [
        'log_weights',
        'unnormalized_log_weights',
        'previous_parent_indices',
        'proposed_particles',
    ])

ParticleFilterLoopVariables = collections.namedtuple(
    'ParticleFilterLoopVariables',
    ['step',
     'previous_step_results',
     'accumulated_step_results',
     'num_prev_particles'
    ])


def particle_filter(observations,
                    initial_state_prior,
                    initial_proposal,
                    transition_fn,
                    observation_fn,
                    num_particles,
                    proposal_fn,
                    num_resampled_ancestors,
                    optimal_resampling=False,
                    multinomial_resampling=False,
                    name=None):  # pylint: disable=g-doc-args
  """Run particle filter.
  """
  with tf.name_scope(name or 'particle_filter'):
    num_timesteps = prefer_static.shape(
        tf.nest.flatten(observations)[0])[0]
    # Dress up the prior and prior proposal as a fake `transition_fn`
    prior_fn = lambda _1, _2: initial_state_prior


    #Run first step of the filter
    initial_step_results, num_prev_particles = _filter_first_step(
      step = 0,
      #previous_particles = prior_fn(0, []).sample(),
      observation=tf.nest.map_structure(
          lambda x: tf.gather(x, 0), observations),
      transition_fn=prior_fn,
      observation_fn=observation_fn,
      proposal_fn=initial_proposal,
      num_resampled_ancestors=num_resampled_ancestors,
      num_max_particles = num_particles)


    def _loop_body(step,
                   previous_step_results,
                   accumulated_step_results,
                   num_prev_particles):
      """Takes one filter step and accumulates results."""
      observation_idx = step
      current_observation = tf.nest.map_structure(
          lambda x, step=step: tf.gather(x, observation_idx), observations)

      new_step_results , num_current_particles = _filter_one_step(
          step = step,
          previous_particles = previous_step_results.proposed_particles,
          previous_unnormalized_log_weights = previous_step_results.unnormalized_log_weights,
          observation = current_observation,
          transition_fn = transition_fn,
          observation_fn = observation_fn,
          proposal_fn = proposal_fn,
          num_resampled_ancestors = num_resampled_ancestors,
          num_max_particles = num_particles,
          num_prev_particles = num_prev_particles,
          optimal_resampling = optimal_resampling,
          multinomial_resampling = multinomial_resampling
      )

      return _update_loop_variables(
          step, new_step_results,
        accumulated_step_results, num_current_particles)


    #Loop through the filter
    loop_results = tf.while_loop(
        cond=lambda step, *_: step < num_timesteps,
        body=_loop_body,
        loop_vars=_initialize_loop_variables(
            initial_step_results,
            num_timesteps, num_prev_particles))

    #convert/stack accumulated results from tensorarrays to tensors
    results = tf.nest.map_structure(lambda ta: ta.stack(),
                                    loop_results.accumulated_step_results)


    return results


def _filter_first_step(step,
                     observation,
                     transition_fn,
                     observation_fn,
                     proposal_fn,
                     num_resampled_ancestors,
                     num_max_particles):
  """First step of the marginal filter"""
  with tf.name_scope('filter_first_step'):

      proposed_particles = proposal_fn()
      num_current_particles = tf.shape(list(proposed_particles.values())[0])[0]
      reshaped_proposed_particles = tf.nest.map_structure(
        lambda x: tf.reshape(x, [-1] + [y for y in x.shape[2:]]), proposed_particles)
      transition_log_probs = transition_fn(0, []).log_prob(proposed_particles)
      transition_log_probs = tf.reshape(transition_log_probs, [-1])
      observation_log_probs = observation_fn(step, reshaped_proposed_particles
                                             ).log_prob(observation)
      unnormalized_log_weights = observation_log_probs + transition_log_probs
      log_weights = tf.nn.log_softmax(unnormalized_log_weights, axis = 0)

      #return values corresponding to prev are not needed in the first step
      first_step_results = ParticleFilterStepResults(
        log_weights = tf.cast(log_weights, dtype),
        unnormalized_log_weights = tf.cast(unnormalized_log_weights, dtype),
        previous_parent_indices = tf.cast(tf.range(0), tf.int32),
        proposed_particles = reshaped_proposed_particles
              )
      return expand_collapsed_results(
        first_step_results, num_max_particles, num_resampled_ancestors = num_resampled_ancestors), num_current_particles



def _filter_one_step(step,
                     observation,
                     previous_particles,
                     previous_unnormalized_log_weights,
                     transition_fn,
                     observation_fn,
                     proposal_fn,
                     num_resampled_ancestors,
                     num_max_particles,
                     num_prev_particles,
                     optimal_resampling,
                     multinomial_resampling
                     ):
  """Advances the particle filter by a single time step
  for a deterministic state proposal across all possible next states."""

  with tf.name_scope('filter_one_step'):
    #tf.print(step)
    #consider only prev results up to index N_{t-1}=num_prev_particles
    previous_unnormalized_log_weights = previous_unnormalized_log_weights[:num_prev_particles]
    previous_particles = tf.nest.map_structure(lambda x: x[:num_prev_particles], previous_particles)
    #resample ancestors
    log_weights = tf.nn.log_softmax(previous_unnormalized_log_weights, axis = 0)
    non_zero_weights_indices = tf.boolean_mask(tf.range(num_prev_particles), previous_unnormalized_log_weights > -np.Inf)
    num_prev_non_zero_weights = tf.shape(non_zero_weights_indices)[-1]
    with tf.name_scope('resample'):
      if optimal_resampling:
        _, parent_indices, log_c, use_unbiased_weights, _, _ = prefer_static.cond(
          num_prev_non_zero_weights > num_resampled_ancestors,
          lambda: resampling_functions.OptimalFiniteState(
                tf.cast(log_weights, dtype=tf.float32), num_resampled_ancestors, []),
          lambda:  (tf.ones([], dtype=tf.bool), non_zero_weights_indices, 0.,
                                     tf.zeros([], dtype=tf.bool), tf.ones([num_resampled_ancestors], tf.bool),
                    num_resampled_ancestors)
        )
      else:
        _, parent_indices, log_c, use_unbiased_weights, _, _ = prefer_static.cond(
          num_prev_non_zero_weights > num_resampled_ancestors,
          lambda: resampling_functions.UnbiasedResampling(
          multinomial_resampling, tf.cast(log_weights, dtype = tf.float32), num_resampled_ancestors),
          lambda: (tf.ones([], dtype = tf.bool), non_zero_weights_indices, 0.,
                   tf.zeros([], dtype = tf.bool), tf.ones([num_resampled_ancestors], tf.bool),
                   num_resampled_ancestors)
        )


    #propose new particles deterministically
    prev_resampled_particles = tf.nest.map_structure(
      lambda x: tf.gather(x, parent_indices), previous_particles
    )
    proposed_particles = proposal_fn(prev_resampled_particles)
    reshaped_proposed_particles = tf.nest.map_structure(
      lambda x: tf.reshape(x, [-1] + [y for y in x.shape[2:]]), proposed_particles)
    num_particles = tf.shape(list(reshaped_proposed_particles.values())[0])[0]
    ####
    #Compute log weights
    ####

    observation_log_probs = observation_fn(step, reshaped_proposed_particles).log_prob(observation)
    transition_log_probs = transition_fn(step,prev_resampled_particles).log_prob(proposed_particles)
    transition_log_probs = tf.reshape(transition_log_probs, [-1])
    log_gamma = tf.where(
      tf.math.is_finite(transition_log_probs),
      tf.cast(transition_log_probs, dtype=dtype) + tf.cast(observation_log_probs, dtype=dtype),
      - np.Inf)

    # used only for unbiased resampling scheme where num_resampled_ancestors remains fixed
    log_weights_pre_factor = -tf.math.log(tf.cast(num_resampled_ancestors, dtype = dtype)) \
                             + tf.math.reduce_logsumexp(previous_unnormalized_log_weights)

    # used only for optimal resampling scheme
    previous_unnormalized_log_weights_ancestors = tf.gather(previous_unnormalized_log_weights,
                                                            parent_indices, axis=0)
    previous_log_weights_ancestors = previous_unnormalized_log_weights_ancestors - tf.reduce_logsumexp(
      previous_unnormalized_log_weights)

    expanded_previous_unnormalized_log_weights_ancestors = tf.reshape(
      previous_unnormalized_log_weights_ancestors*tf.ones(
        [num_particles//tf.minimum(num_prev_non_zero_weights,num_resampled_ancestors),1], dtype=dtype),[-1])
    expanded_previous_log_weights_ancestors = tf.reshape(
      previous_log_weights_ancestors*tf.ones(
        [num_particles//tf.minimum(num_prev_non_zero_weights,num_resampled_ancestors),1], dtype=dtype),[-1])


    unnormalized_log_weights = prefer_static.cond(
      num_prev_non_zero_weights <= num_resampled_ancestors,
      lambda: expanded_previous_unnormalized_log_weights_ancestors + log_gamma,
      lambda: prefer_static.cond(
        use_unbiased_weights,
        lambda: log_weights_pre_factor + log_gamma,
        lambda: expanded_previous_unnormalized_log_weights_ancestors + log_gamma -tf.minimum(
          tf.zeros([], dtype), tf.cast(log_c, dtype = dtype) + expanded_previous_log_weights_ancestors)
      )
      )


    if False:
      #normalise log weights for better stability (normalising constants are then useless)
      unnormalized_log_weights = log_weights


    step_results = ParticleFilterStepResults(
        log_weights=log_weights,
        unnormalized_log_weights=unnormalized_log_weights,
        previous_parent_indices=parent_indices,
        proposed_particles=reshaped_proposed_particles
    )

    return expand_collapsed_results(
      step_results, num_max_particles, num_resampled_ancestors = num_resampled_ancestors), num_particles


def _update_loop_variables(step,
                           current_step_results,
                           accumulated_step_results,
                           num_current_particles):
  """Update the loop state to reflect a step of filtering."""

  # Write particles, partent indices, and log weights to their respective arrays.
  # These must all have the same length although the number of particles
  # can be different at each step,
  with tf.name_scope('accumulate'):
    accumulated_step_results = tf.nest.map_structure(
        lambda x, y: x.write(step, y),
        accumulated_step_results, current_step_results)



  return ParticleFilterLoopVariables(
      step=step + 1,
      previous_step_results=current_step_results,
      accumulated_step_results=accumulated_step_results,
      num_prev_particles=num_current_particles
  )


def _initialize_loop_variables(initial_step_results,
                               num_timesteps, num_prev_particles):
  """Initialize arrays and other quantities passed through the filter loop."""

  # Create arrays to store particles, indices, and likelihoods, and write
  # their initial values.
  step_results_arrays = tf.nest.map_structure(
      lambda x: tf.TensorArray(dtype=x.dtype, size=num_timesteps, infer_shape=False).write(0, x),
      initial_step_results)

  return ParticleFilterLoopVariables(
      step=1,
      previous_step_results=initial_step_results,
      accumulated_step_results=step_results_arrays,
      num_prev_particles=num_prev_particles
  )



def expand_collapsed_results(new_step_results, num_particles, num_resampled_ancestors):

    d_ = tf.shape(new_step_results.previous_parent_indices)[0]
    previous_parent_indices_ = tf.concat([new_step_results.previous_parent_indices,
                                          -1 * tf.ones([num_resampled_ancestors - d_],
                                                       tf.int32)], 0)
    previous_parent_indices_ = tf.ensure_shape(previous_parent_indices_,[num_resampled_ancestors])
    d_ = tf.shape(new_step_results.log_weights)[0]
    log_weights_ = tf.concat([new_step_results.log_weights,
                              -np.Inf * tf.ones([num_particles - d_], dtype = dtype)], 0)
    log_weights_ = tf.ensure_shape(log_weights_, [num_particles])
    d_ = tf.shape(new_step_results.unnormalized_log_weights)[0]
    unnormalized_log_weights_ = tf.concat([new_step_results.unnormalized_log_weights,
                                           -np.Inf * tf.ones([num_particles - d_], dtype = dtype)], 0)
    unnormalized_log_weights_ = tf.ensure_shape(unnormalized_log_weights_, [num_particles])

    nan_particles = tf.nest.map_structure(lambda x: tf.tile(-1 * tf.expand_dims(x[0], 0),
                                                            [num_particles - d_] + [1 for _ in range(len(x[0].shape))]),
                                          new_step_results.proposed_particles)
    proposed_particles_ = tf.nest.map_structure(lambda x, y: tf.concat([x, y], 0),
                                                new_step_results.proposed_particles, nan_particles)
    proposed_particles_ = tf.nest.map_structure(lambda x: tf.ensure_shape(x,[num_particles] + x.shape[1:])
                          , proposed_particles_)

    return ParticleFilterStepResults(
        log_weights = log_weights_,
        unnormalized_log_weights = unnormalized_log_weights_,
        previous_parent_indices = previous_parent_indices_,
        proposed_particles = proposed_particles_
    )
