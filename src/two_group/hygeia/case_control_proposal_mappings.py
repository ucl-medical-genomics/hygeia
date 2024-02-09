import tensorflow as tf

class CaseControlProposal:
  """
  Case Control Proposal
  """

  def __init__(self, n_methylation_regimes):
      self._n_methylation_regimes = n_methylation_regimes

  def _xi(self, state):
    """xi mapping that depends on the state
    Does not work for batches of states but is used in proposal_fn_resampled
    that loops through all states/particles"""
    I = 2 * self._n_methylation_regimes
    merged_state = state['merged_state']
    #merged_state =  tf.gather(state['merged_state'], 0, axis =  - 1)
    duration_control = tf.gather(state['control_state'], 0, axis = len(state['control_state'].shape) - 1)
    regime_control = tf.gather(state['control_state'], 1, axis = len(state['control_state'].shape) - 1)
    duration_case = tf.gather(state['case_state'], 0, axis = len(state['case_state'].shape) - 1)
    regime_case = tf.gather(state['case_state'], 1, axis = len(state['case_state'].shape) - 1)

    #init array of proposed state
    next_state_arrays = tf.nest.map_structure(
      lambda x: tf.TensorArray(dtype = x, size = I),
      {'merged_state': tf.int32,'control_state': tf.int32,'case_state': tf.int32})
    #first dimension of having no change point
    no_change_point_state = {'merged_state': merged_state,
                              'control_state': tf.stack([duration_control + 1, regime_control], -1),
                              'case_state': tf.stack([duration_case + 1, regime_case], -1)
                              }
    next_state_arrays = tf.nest.map_structure(
      lambda x, y: x.write(0, y), next_state_arrays, no_change_point_state)

    #then loop through split states with change point in control
    def loop_fn_control_split1(idx, next_state_arrays):
      next_state = {'merged_state': tf.zeros([], dtype=tf.int32),
              'control_state': tf.stack([ 1, idx-1], -1),
              'case_state': tf.stack([duration_case + 1, regime_case], -1)}
      next_state_arrays = tf.nest.map_structure(
      lambda x, y: x.write(idx, y), next_state_arrays, next_state)
      return idx+1, next_state_arrays
    def loop_fn_control_split2(idx, next_state_arrays):
      next_state = {'merged_state': tf.zeros([], dtype=tf.int32),
              'control_state': tf.stack([ 1, idx], -1),
              'case_state': tf.stack([duration_case + 1, regime_case], -1)}
      next_state_arrays = tf.nest.map_structure(
      lambda x, y: x.write(idx, y), next_state_arrays, next_state)
      return idx+1, next_state_arrays

    idx, next_state_arrays = tf.while_loop(
      cond = lambda idx, *_: idx <= regime_case,
      body = loop_fn_control_split1,
      loop_vars = [1, next_state_arrays])
    idx, next_state_arrays = tf.while_loop(
      cond = lambda idx, *_: idx < self._n_methylation_regimes,
      body = loop_fn_control_split2,
      loop_vars = [idx, next_state_arrays])


    # then loop through split states with change point in case
    def loop_fn_case_split1(idx, next_state_arrays):
      next_state = {'merged_state': tf.zeros([], dtype = tf.int32),
                    'control_state': tf.stack([duration_control+1, regime_control], -1),
                    'case_state': tf.stack([1, idx-self._n_methylation_regimes], -1)}
      next_state_arrays = tf.nest.map_structure(
        lambda x, y: x.write(idx, y), next_state_arrays, next_state)
      return idx + 1, next_state_arrays

    def loop_fn_case_split2(idx, next_state_arrays):
      next_state = {'merged_state': tf.zeros([], dtype = tf.int32),
                    'control_state': tf.stack([duration_control+1, regime_control], -1),
                    'case_state': tf.stack([1, idx-self._n_methylation_regimes+1], -1)}
      next_state_arrays = tf.nest.map_structure(
        lambda x, y: x.write(idx, y), next_state_arrays, next_state)
      return idx + 1, next_state_arrays


    idx, next_state_arrays = tf.while_loop(
      cond = lambda idx, *_: idx < self._n_methylation_regimes+regime_control,
      body = loop_fn_case_split1,
      loop_vars = [idx, next_state_arrays])
    idx, next_state_arrays = tf.while_loop(
      cond = lambda idx, *_: idx < I-1,
      body = loop_fn_case_split2,
      loop_vars = [idx, next_state_arrays])

    #propose to merged states with a single change point happening,
    #but only if currently split (otherwise set durations to 0
    merged_state_case_cp = {'merged_state': tf.ones([], dtype = tf.int32),
                  'control_state': tf.stack([tf.cond(merged_state == 0,
                                                     lambda: duration_control + 1,
                                                     lambda: 0),
                                             regime_control], -1),
                  'case_state': tf.stack([tf.cond(merged_state == 0,
                                          lambda: duration_control + 1,
                                          lambda: 0),
                                          regime_control], -1)}
    next_state_arrays = tf.nest.map_structure(
        lambda x, y: x.write(I-1, y), next_state_arrays, merged_state_case_cp)


    return tf.nest.map_structure(lambda ta: ta.stack(), next_state_arrays)


  def proposal_fn_non_resampled(self):
    """xi mapping that depends not on the state"""
    next_state_arrays = tf.nest.map_structure(
      lambda x: tf.TensorArray(dtype = x, size = self._n_methylation_regimes),
      {'merged_state': tf.int32, 'control_state': tf.int32, 'case_state': tf.int32})

    def loop_fn_two_change_points(idx, next_state_arrays):
      next_state = {'merged_state': tf.expand_dims(tf.where(idx == tf.range(self._n_methylation_regimes),
                                             tf.ones([], dtype=tf.int32),
                                             tf.zeros([], dtype=tf.int32)), -1),
                    'control_state': tf.stack([tf.ones([self._n_methylation_regimes], dtype=tf.int32),
                                               idx*tf.ones([self._n_methylation_regimes], dtype=tf.int32)], -1),
                    'case_state': tf.stack([tf.ones([self._n_methylation_regimes], dtype=tf.int32),
                                            tf.range(self._n_methylation_regimes)], -1)}
      next_state_arrays = tf.nest.map_structure(
        lambda x, y: x.write(idx, y), next_state_arrays, next_state)
      return idx+1, next_state_arrays

    idx, next_state_arrays  = tf.while_loop(
      cond = lambda idx, *_: idx < self._n_methylation_regimes,
      body = loop_fn_two_change_points,
      loop_vars = [0, next_state_arrays],
      parallel_iterations = self._n_methylation_regimes)

    next_states_reshaped = tf.nest.map_structure(lambda ta: ta.stack(), next_state_arrays)
    #next_states = tf.nest.map_structure(lambda x: tf.reshape(x, [-1, x.shape[-1]]), next_states_reshaped)
    next_states = tf.nest.map_structure(lambda x: tf.reshape(x, [-1, tf.shape(x)[-1]]), next_states_reshaped)
    next_states['merged_state'] = tf.squeeze(next_states['merged_state'], -1)
    return next_states


  def proposal_fn_resampled(self, states, parent_indices):
    """returns all possible next states given states that are resampled
        and also return just the next possible states given the parent
        indices of the selected states"""

    M_=tf.shape(states['merged_state'])[0]
    #tf.print('M_', M_)
    next_state_arrays = tf.nest.map_structure(
      lambda x: tf.TensorArray(dtype = x, size = 0, dynamic_size = True),
      {'merged_state': tf.int32, 'control_state': tf.int32, 'case_state': tf.int32})

    def loop_body_xi(idx, next_state_arrays):
      next_states = self._xi(tf.nest.map_structure(lambda x: tf.gather(x, idx, axis = 0), states))
      next_state_arrays = tf.nest.map_structure(lambda x, y: x.write(idx,y), next_state_arrays, next_states)
      return idx+1, next_state_arrays

    idx, next_state_arrays = tf.while_loop(
      cond = lambda idx, *_: idx < M_,
      body = loop_body_xi,
      loop_vars = [0, next_state_arrays],
      parallel_iterations = 1000)

    next_states = tf.nest.map_structure(lambda ta: ta.stack(), next_state_arrays)

    selected_next_states = tf.nest.map_structure(lambda x: tf.gather(x, parent_indices, axis = 0),
                              next_states)

    selected_next_states['merged_state'] = tf.expand_dims(
        selected_next_states['merged_state'], -1)
    #selected_next_states = tf.nest.map_structure(
    #    lambda x: tf.reshape(x, [-1, x.shape[-1]]), selected_next_states)
    selected_next_states = tf.nest.map_structure(
        lambda x: tf.reshape(x, [-1, tf.shape(x)[-1]]), selected_next_states)
    selected_next_states['merged_state'] = tf.squeeze(
        selected_next_states['merged_state'], -1)
    return next_states, selected_next_states


  def proposal_fn_standard_filter(self, states):
    """returns all possible next states given states that includes
    all possible proposals (both dependent and independent of the current state"""
    M_ = tf.shape(states['merged_state'])[0]
    # tf.print('M_', M_)
    next_state_arrays = tf.nest.map_structure(
      lambda x: tf.TensorArray(dtype = x, size = 0, dynamic_size = True),
      {'merged_state': tf.int32, 'control_state': tf.int32, 'case_state': tf.int32})



    def loop_body_xi(idx, next_state_arrays):
      next_states = self._xi(tf.nest.map_structure(lambda x: tf.gather(x, idx, axis = 0), states))
      next_state_arrays = tf.nest.map_structure(lambda x, y: x.write(idx,y), next_state_arrays, next_states)
      return idx+1, next_state_arrays

    idx, next_state_arrays = tf.while_loop(
      cond = lambda idx, *_: idx < M_,
      body = loop_body_xi,
      loop_vars = [0, next_state_arrays],
      parallel_iterations = 1000)

    next_states_state_dependent = tf.nest.map_structure(lambda ta: ta.stack(), next_state_arrays)
    next_states_state_independent = tf.nest.map_structure(
      lambda x: tf.tile(tf.expand_dims(x,0),[M_] + [tf.ones([], dtype=tf.int32) for _ in x.shape]), self.proposal_fn_non_resampled())
    proposed_particles = tf.nest.map_structure(
      lambda x, y: tf.concat([x, y], 1),
      next_states_state_dependent, next_states_state_independent)
    #transpose so output dims [I,M,d]
    proposed_particles = tf.nest.map_structure(lambda x: tf.einsum("ij...->ji...", x), proposed_particles)
    #proposed_particles = tf.nest.map_structure(
    #  lambda x: tf.reshape(x, [-1] + [y for y in x.shape[2:]]), proposed_particles)
    return proposed_particles


  def initial_proposal_fn_standard_filter(self):
    """initial proposal function of deterministic partciel filter"""
    M_ = 1
    proposed_particles = tf.nest.map_structure(
      lambda x: tf.tile(tf.expand_dims(x,0),[M_] + [tf.ones([], dtype=tf.int32) for _ in x.shape]), self.proposal_fn_non_resampled())
    proposed_particles = tf.nest.map_structure(lambda x: tf.einsum("ij...->ji...", x), proposed_particles)
    return proposed_particles
