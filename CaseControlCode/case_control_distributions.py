from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# Dependency imports
import tensorflow.compat.v2 as tf

from tensorflow_probability.python.distributions import distribution
from tensorflow_probability.python.internal import prefer_static
from tensorflow_probability.python.internal import reparameterization

from tensorflow_probability.python.internal import samplers

import tensorflow_probability as tfp
tfd = tfp.distributions



class InitialControlState(distribution.Distribution):
  """Initial state density of the control state.
  The state is (duration, regime)
  """


  def __init__(self,
               n_regimes,
               init_duration,
               validate_args = False,
               allow_nan_stats = True,
               use_static_graph = False,
               name = 'InitialControlState'):
    parameters = dict(locals())
    with tf.name_scope(name) as name:
      self._n_regimes = n_regimes
      self._init_duration = init_duration
      self._use_static_graph = use_static_graph
      # self._logits = tf.math.log(tf.ones(self._num_regimes - 1))

      super(InitialControlState, self).__init__(
        dtype = tf.int32,
        reparameterization_type = reparameterization.NOT_REPARAMETERIZED,
        validate_args = validate_args,
        allow_nan_stats = allow_nan_stats,
        parameters = parameters,
        name = name)

  def _batch_shape_tensor(self):
    return tf.constant([], dtype=tf.int32)

  def _batch_shape(self):
    return tf.TensorShape([])

  def _event_shape_tensor(self):
    return tf.constant([2], dtype = tf.int32)

  def _event_shape(self):
    return tf.TensorShape([2])  # tf.shape(self._regime_control)#tf.TensorShape([])

  def _log_prob(self, x):
    #x = tf.expand_dims(x,-1)
    regimes_log_probs = -tf.nn.sparse_softmax_cross_entropy_with_logits(
      labels = x[:,1], logits = tf.tile(tf.expand_dims(
        tf.math.log(tf.ones(self._n_regimes)), 0), [x.shape[0], 1]))
    duration_log_probs = tf.math.log(tf.cast(x[:,0]== self._init_duration, dtype = tf.float32))
    return regimes_log_probs + duration_log_probs

  def _sample_n(self, n, seed = None):

    regime_samples = tf.random.categorical(
      tf.math.log(tf.ones([1,self._n_regimes])),
      n, seed = seed, dtype = tf.int32)
    regime_samples = tf.squeeze(regime_samples, 0)
    duration_samples = self._init_duration * tf.ones_like(regime_samples)
    return tf.stack([duration_samples,regime_samples], axis=1)


class ControlStateTransition(distribution.Distribution):
  """Density of the duration and regime state for the control group.
  """

  def __init__(self,
               current_merged_state,
               previous_duration_control,
               previous_regime_control,
               rho,
               regime_probs,
               validate_args = False,
               allow_nan_stats = True,
               use_static_graph = False,
               name = 'ControlStateTransition'):
    parameters = dict(locals())
    with tf.name_scope(name) as name:
      self._current_merged_state = prefer_static.cond(tf.cast(current_merged_state.shape == [], tf.bool),
                                           lambda:tf.expand_dims(current_merged_state, 0),
                                           lambda: current_merged_state)
      self._current_merged_state = tf.expand_dims(self._current_merged_state, -1)


      batch_size_dynamic = prefer_static.cond(tf.cast(self._current_merged_state.shape[:-1] != [], tf.bool),
                                      lambda: tf.shape(self._current_merged_state)[:-1],
                                      lambda: tf.TensorShape([1]))
      batch_size = prefer_static.cond(tf.cast(self._current_merged_state.shape[:-1] != [], tf.bool),
                                      lambda: self._current_merged_state.shape[:-1],
                                      lambda: tf.TensorShape([1]))
      self._bs = batch_size
      self._regime_probs = tf.broadcast_to(
          regime_probs,
          tf.concat([batch_size_dynamic,
                     tf.shape(regime_probs)[-1:]],0))
      self._rho = tf.broadcast_to(rho, batch_size_dynamic)
      self._previous_duration_control = tf.broadcast_to(previous_duration_control,
                                                        batch_size_dynamic)
      self._previous_regime_control = tf.broadcast_to(previous_regime_control,
                                                      batch_size_dynamic)


      super(ControlStateTransition, self).__init__(
        dtype = tf.int32,
        reparameterization_type = reparameterization.NOT_REPARAMETERIZED,
        validate_args = validate_args,
        allow_nan_stats = allow_nan_stats,
        parameters = parameters,
        name = name)

  #def _batch_shape_tensor(self):
  #  return  tf.tile(tf.expand_dims(self._rho, -1), [1, 2])

  def _batch_shape(self):
    return self._bs
    return self._bs

  def _event_shape_tensor(self):
    return tf.constant([2], dtype = tf.int32)

  def _event_shape(self):
    return tf.TensorShape([2])  # tf.shape(self._regime_control)#tf.TensorShape([])

  def _log_prob(self, x):
    d = tf.gather(x,indices = 0, axis = -1)
    r = tf.gather(x, indices = 1, axis = -1)
    log_prob = tf.where(d == tf.ones_like(d,tf.int32),
                        tf.math.add(tf.math.log(self._rho),
                                    tfd.Categorical(
                          logits = tf.math.log(self._regime_probs)).log_prob(r)),
                        tf.math.log(1 - self._rho)
                        + tf.math.log(tf.cast(self._previous_duration_control == d-1, tf.float32))\
                        + tf.math.log(tf.cast(self._previous_regime_control == tf.cast(r, tf.int32),
                                              tf.float32))
                        )

    return log_prob

  def _sample_n(self, n, seed = None):
    seed = samplers.sanitize_seed(seed)
    seeds = samplers.split_seed(seed, n = 2)
    uniform = samplers.uniform(tf.shape(tf.tile(
      tf.expand_dims(self._previous_duration_control, 0), [n, 1])),
      seed = seeds[0], dtype = self._rho.dtype)
    new_change_point = tf.less(uniform, self._rho)

    next_regime = tf.where(new_change_point,
                           tfd.Categorical(
                             logits = tf.math.log(self._regime_probs)).sample(n, seed = seed),
                           tfd.Deterministic(self._previous_regime_control).sample(n)
                           )
    next_duration = tf.where(new_change_point,
                             1,
                             tfd.Deterministic(self._previous_duration_control + 1).sample(n)
                             )
    return prefer_static.cond(tf.cast(self._current_merged_state.shape==[], tf.bool),
                   lambda: tf.squeeze(tf.concat(tf.stack([next_duration, next_regime], axis=-1), 0), 0),
                   lambda: tf.concat(tf.stack([next_duration, next_regime], axis=-1), 0))




class CaseStateTransition(distribution.Distribution):
  """Density of the duration and regime state for the case group.
    """

  def __init__(self,
               previous_duration_case,
               previous_regime_case,
               previous_merged_state,
               current_merged_state,
               current_control_states,
               rho,
               n_regimes,
               validate_args = False,
               allow_nan_stats = True,
               use_static_graph = False,
               name = 'ControlStateTransition'):
    parameters = dict(locals())
    with tf.name_scope(name) as name:
      self._current_merged_state = prefer_static.cond(tf.cast(current_merged_state.shape ==  [], tf.bool),
                                           lambda:tf.expand_dims(current_merged_state, 0),
                                           lambda: current_merged_state)
      self._current_merged_state = tf.expand_dims(self._current_merged_state, -1)

      self._n_regimes = n_regimes
      batch_size_dynamic = prefer_static.cond(tf.cast(self._current_merged_state.shape[:-1] != [], tf.bool),
                           lambda: tf.shape(self._current_merged_state)[:-1],
                           lambda: tf.TensorShape([1]))
      batch_size = prefer_static.cond(tf.cast(self._current_merged_state.shape[:-1] != [], tf.bool),
                                      lambda: self._current_merged_state.shape[:-1],
                                      lambda: tf.TensorShape([1]))



      self._bs = batch_size
      self._previous_merged_state = prefer_static.cond(tf.cast(previous_merged_state.shape[-1:] != 1, tf.bool),
                                                       lambda: tf.expand_dims(previous_merged_state, -1),
                                                       lambda: previous_merged_state)
      self._rho = tf.broadcast_to(rho, batch_size_dynamic)
      self._previous_duration_case = tf.broadcast_to(previous_duration_case, batch_size_dynamic)
      self._previous_regime_case = tf.broadcast_to(previous_regime_case, batch_size_dynamic)


      current_duration_control = tf.gather(current_control_states, 0,
                                           axis=len(current_control_states.shape)-1)
      current_regime_control = tf.gather(current_control_states, 1,
                                         axis=len(current_control_states.shape)-1)
      self._current_duration_control = tf.broadcast_to(current_duration_control, batch_size_dynamic)
      self._current_regime_control = tf.broadcast_to(current_regime_control, batch_size_dynamic)

      self._use_static_graph = use_static_graph


      super(CaseStateTransition, self).__init__(
        dtype = tf.int32,
        reparameterization_type = reparameterization.NOT_REPARAMETERIZED,
        validate_args = validate_args,
        allow_nan_stats = allow_nan_stats,
        parameters = parameters,
        name = name)

  def _batch_shape(self):
    return self._bs

  def _event_shape_tensor(self):
    return tf.constant([2], dtype = tf.int32)

  def _event_shape(self):
    return tf.TensorShape([2])  # tf.shape(self._regime_control)#tf.TensorShape([])

  def _log_prob(self, x):
    d = tf.gather(x, indices = 0, axis = -1)
    r = tf.gather(x, indices = 1, axis = -1)

    log_prob = tf.where(tf.squeeze(self._current_merged_state, -1) == 1,
                        tf.math.log(tf.cast(self._current_regime_control == r, tf.float32)
                                    ) + tf.math.log(tf.cast(self._current_duration_control == d, tf.float32)),
                        tf.where(tf.math.logical_and(tf.squeeze(self._previous_merged_state, -1) == 1,
                                                     self._current_duration_control != 1),
                                 tfd.Categorical(
                                   logits = tf.math.log(tf.gather(
                                     tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                        tf.zeros([self._n_regimes])),
                                     tf.cast(self._current_regime_control, dtype = tf.int32)))
                                 ).log_prob(r)+tf.math.log(tf.cast(d == 1, dtype=tf.float32)),
                                 tf.where(tf.math.logical_and(self._current_regime_control == self._previous_regime_case,
                                                              tf.squeeze(self._previous_merged_state, -1) == 0),
                                          tf.math.log(tf.cast(d == 1, tf.float32)) + tfd.Categorical(
                                              logits = tf.math.log(tf.math.multiply(
                                                  tf.gather(
                                                      tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                         tf.zeros([self._n_regimes])),
                                                      tf.cast(self._current_regime_control, dtype = tf.int32)),
                                                  tf.gather(
                                                      tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                         tf.zeros([self._n_regimes])),
                                                      tf.cast(self._previous_regime_case,
                                                              dtype = tf.int32))))).log_prob(r),
                                          tf.where(
                                              d == 1,
                                              tf.math.log(self._rho) + tfd.Categorical(
                                                logits = tf.math.log(tf.math.multiply(
                                                  tf.gather(tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                     tf.zeros([self._n_regimes])),
                                                            tf.cast(self._current_regime_control, dtype = tf.int32)),
                                                  tf.gather(tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                     tf.zeros([self._n_regimes])),
                                                            tf.cast(self._previous_regime_case, dtype = tf.int32))))).log_prob(r),
                                                tf.math.log(1 - self._rho) + tf.math.log(tf.cast(
                                                    self._previous_duration_case + 1 == d, tf.float32)
                                                )+ tf.math.log(tf.cast(self._previous_regime_case == r, tf.float32))
                                              )
                                          )
                                 )
                        )
    return log_prob

  def _sample_n(self, n, seed = None):
    seed = samplers.sanitize_seed(seed)
    seeds = samplers.split_seed(seed, n = 2)
    uniform = samplers.uniform(tf.shape(tf.tile(
      tf.expand_dims(self._current_duration_control, 0), [n, 1])),
      seed = seeds[0], dtype = self._rho.dtype)
    next_regime = tf.where(tf.squeeze(self._current_merged_state, -1) == 1,
                           tfd.Deterministic(self._current_regime_control).sample(n),
                           tf.where(tf.math.logical_and(tf.squeeze(self._previous_merged_state, -1) == 1,
                                                     self._current_duration_control != 1),
                                    tfd.Categorical(
                                        logits = tf.math.log(tf.gather(
                                            tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                               tf.zeros([self._n_regimes])),
                                            tf.cast(self._current_regime_control, dtype = tf.int32)))
                                    ).sample(n, seed=seed),
                                    tf.where(tf.math.logical_and(self._current_regime_control == self._previous_regime_case,
                                                              tf.squeeze(self._previous_merged_state, -1) == 0),
                                        tfd.Categorical(
                                                 logits = tf.math.log(tf.math.multiply(
                                                  tf.gather(
                                                      tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                         tf.zeros([self._n_regimes])),
                                                      tf.cast(self._current_regime_control, dtype = tf.int32)),
                                                  tf.gather(
                                                      tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                         tf.zeros([self._n_regimes])),
                                                      tf.cast(self._previous_regime_case,
                                                              dtype = tf.int32))))).sample(
                                                 n, seed = seed),
                                        tf.where(tf.less(uniform, self._rho),
                                             tfd.Categorical(
                                               logits = tf.math.log(tf.math.multiply(
                                                 tf.gather(
                                                   tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                      tf.zeros([self._n_regimes])),
                                                   tf.cast(self._current_regime_control, dtype = tf.int32)),
                                                 tf.gather(
                                                   tf.linalg.set_diag(tf.ones([self._n_regimes, self._n_regimes]),
                                                                      tf.zeros([self._n_regimes])),
                                                   tf.cast(self._previous_regime_case, dtype = tf.int32))))).sample(n, seed = seed),
                                             tfd.Deterministic(self._previous_regime_case).sample(n)
                                             )
                                        )
                                    )
                           )


    next_duration = tf.where(tf.squeeze(self._current_merged_state == 1, -1),
                             tfd.Deterministic(self._current_duration_control).sample(n),
                             tf.where(tf.math.logical_or(
                                 tf.math.logical_and(
                                     tf.squeeze(self._previous_merged_state, -1) == 1,
                                     self._current_duration_control != 1),
                                 tf.math.logical_and(
                                     tf.squeeze(self._previous_merged_state, -1) == 0,
                                     self._current_regime_control == self._previous_regime_case)
                             ),
                                      tfd.Deterministic(tf.ones(
                                        shape=tf.shape(self._current_duration_control)[0], dtype=tf.int32)).sample(n),
                                      tf.where(self._current_regime_control == self._previous_regime_case,
                                            1,
                                            tf.where(tf.less(uniform, self._rho),
                                               1,
                                               tfd.Deterministic(self._previous_duration_case + 1).sample(n)
                                               )
                                            )
                                      )
                             )

    # tf.debugging.assert_equal(tf.reduce_sum(tf.cast(
    #     tf.where(tf.squeeze(self._current_merged_state, -1) == 1,
    #              self._current_regime_control != next_regime,
    #              self._current_regime_control == next_regime), tf.int32)),
    #     0)
    # tf.debugging.assert_equal(tf.reduce_sum(tf.cast(
    #     tf.where(tf.squeeze(self._current_merged_state, -1) == 1,
    #              self._current_duration_control != next_duration,
    #              tf.zeros([], dtype=tf.bool)), tf.int32)),
    #                           0)

    return prefer_static.cond(tf.cast(self._current_merged_state.shape==[], tf.bool),
                              lambda: tf.squeeze(tf.stack([next_duration, next_regime], axis = -1), 0),
                              lambda: tf.stack([next_duration, next_regime], axis = -1))
