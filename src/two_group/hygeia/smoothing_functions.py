
import tensorflow as tf
import collections
import numpy as np



AuxiliarySmoothingResults = collections.namedtuple(
    'AuxiliarySmoothingResults',
    ['s',
     'q',
     'psi'
    ])

SmoothingResults = collections.namedtuple(
    'SmoothingResults',
    ['functional_estimates',
     's',
     't',
     'q',
     'idx'
    ])




def compute_log_backward_kernel(previous_particles, new_particles,
                            previous_unnormalized_log_weights, transition_fn, step):

    l = previous_unnormalized_log_weights.shape[0]

    expanded_particles = tf.nest.map_structure(lambda x: tf.expand_dims(x, 1),
                                               new_particles)
    proposed_particles_matrix = tf.nest.map_structure(
      lambda x: tf.tile(x, [1,l] + [1 for _ in range(len(x.shape[2:]))]), expanded_particles)

    transition_log_prob_matrix = transition_fn(step, previous_particles).log_prob(proposed_particles_matrix)


    # # transition log prob matrix is [P(X_t^i|X_{t-1}^j)]_{ij}

    log_backward_kernel_weights = compute_log_backward_kernel_from_transition_matrix(previous_unnormalized_log_weights,
                                                   transition_log_prob_matrix)
    return log_backward_kernel_weights

def compute_log_backward_kernel_from_transition_matrix(previous_unnormalized_log_weights,
                                                   transition_log_prob_matrix, dtype = tf.float64):
    # transition log prob matrix is [P(X_t^i|X_{t-1}^j)]_{ij}

    expanded_previous_unnormalized_log_weights = tf.tile(tf.expand_dims(
        previous_unnormalized_log_weights, 0), [tf.shape(transition_log_prob_matrix)[0], 1])

    log_w_plus_log_f = tf.where(
        tf.math.logical_and(tf.math.is_finite(transition_log_prob_matrix),
                            tf.math.is_finite(expanded_previous_unnormalized_log_weights)),
        tf.cast(transition_log_prob_matrix, dtype = dtype) + expanded_previous_unnormalized_log_weights,
        -np.Inf)

    return log_w_plus_log_f - tf.expand_dims(tf.reduce_logsumexp(log_w_plus_log_f, 1),-1)



def smoothing_step(new_step_results, step, test_function,
                   log_backward_kernel, previous_auxiliary_smoothing_results_array, num_timesteps,
                   previous_particle_idx_to_keep,
                   eps = .001, dtype = tf.float64):

    #replace nan in log_backward_kernel with -inf
    log_backward_kernel = tf.where(tf.math.is_nan(log_backward_kernel),
                                   -np.Inf*tf.ones_like(log_backward_kernel),
                                   log_backward_kernel)

    ###
    #First update the auxilliary variables
    ###
    previous_auxiliary_smoothing_results = tf.nest.map_structure(
        lambda x: x.stack(), previous_auxiliary_smoothing_results_array)

    #delete those with zero weights and duplicate ones
    previous_auxiliary_smoothing_results = previous_auxiliary_smoothing_results._replace(
        psi = tf.gather(previous_auxiliary_smoothing_results.psi, previous_particle_idx_to_keep, axis=1)
    )
    previous_psi = previous_auxiliary_smoothing_results.psi
    size_prev_s_set = tf.shape(previous_auxiliary_smoothing_results.q)[0]
    updated_psi = tf.matmul(tf.math.exp(log_backward_kernel), previous_psi, adjoint_b = True)

    ###
    #Then compute the filtered mean and check if the variance wrt the filter is below the threshold
    ###
    filter_weights = tf.cast(tf.math.exp(new_step_results.log_weights), dtype = dtype)
    filtered_mean_estimate = tf.einsum('i,is->s', filter_weights, updated_psi)
    filtered_variance_estimate = tf.einsum('i,is->s', filter_weights,
                                           tf.math.square(updated_psi - filtered_mean_estimate))
    filtered_variance_estimate = tf.cast(filtered_variance_estimate, dtype = tf.float32)
    # keep updating smoothing functional estimates depending on the eps-threshold
    # and set the threshold to Inf at the final step
    eps = tf.cond(step == num_timesteps-1, lambda: np.Inf, lambda: eps)
    num_final_estimates = tf.reduce_sum(tf.cast(tf.math.less(
        filtered_variance_estimate, eps), tf.int32))
    num_non_final_estimates = tf.reduce_sum(tf.cast(tf.math.greater_equal(filtered_variance_estimate, eps), tf.int32))

    ###
    #Construct array of functional estimates that are kept for the next iteration
    ###
    new_smoothing_results_array = tf.nest.map_structure(
        lambda x: tf.TensorArray(dtype = x.dtype,
                                 size = 0,
                                 dynamic_size = True),
        SmoothingResults(functional_estimates = tf.zeros([]),
                                             s = tf.zeros([], dtype = tf.int32),
                                             t = tf.zeros([], dtype = tf.int32),
                                             q = tf.zeros([], dtype = tf.int32),
                                             idx = tf.zeros([], dtype = tf.int32)))
    #functional estimates at the current step
    new_psi = tf.cast(test_function(new_step_results.proposed_particles), dtype = dtype)
    q = tf.shape(new_psi)[-1]
    new_auxiliary_smoothing_results = AuxiliarySmoothingResults(
        psi = tf.transpose(new_psi),
        s = step * tf.ones([q], dtype = tf.int32),
        q = tf.cast(tf.range(q), tf.int32))

    #dynamic size
    updated_auxiliary_smoothing_results = tf.nest.map_structure(
        lambda x: tf.TensorArray(dtype = x.dtype, size = 0,
                                 dynamic_size = True, infer_shape = False),
        new_auxiliary_smoothing_results)

    #write auxiliary smoothing results from current step
    i = tf.cond(step == num_timesteps-1, lambda: q, lambda: 0)
    i, updated_auxiliary_smoothing_results = tf.while_loop(
        lambda i, *_: i < q,
        lambda i, smoothing_results: (i + 1, tf.nest.map_structure(
            lambda x, y: y.write(i+size_prev_s_set-num_final_estimates, tf.gather(x, i, axis = 0)),
            new_auxiliary_smoothing_results, smoothing_results)),
        loop_vars = [i,  updated_auxiliary_smoothing_results])

    def update_final_results(s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results):
        current_s = tf.gather(previous_auxiliary_smoothing_results.s, s_idx, axis = 0)
        current_q = tf.gather(previous_auxiliary_smoothing_results.q, s_idx, axis = 0)
        final_smoothing_result = SmoothingResults(
            functional_estimates = tf.cast(tf.gather(filtered_mean_estimate, s_idx, axis = 0), dtype = tf.float32),
            s = current_s,
            t = step,
            q = current_q,
            idx = current_s * q+ current_q)
        new_smoothing_results_array = tf.nest.map_structure(
            lambda x, y: x.write(s_idx_fin, y), new_smoothing_results_array, final_smoothing_result)
        return s_idx + 1, s_idx_aux, s_idx_fin + 1, new_smoothing_results_array, updated_auxiliary_smoothing_results


    def update_auxilliary_results(s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results):
        auxiliary_smoothing_result = AuxiliarySmoothingResults(
            #psi = tf.gather(new_psi, tf.gather(previous_auxiliary_smoothing_results.q, s_idx, axis = 0), axis = -1),
            psi = tf.gather(updated_psi, s_idx, axis=-1),
            s = tf.gather(previous_auxiliary_smoothing_results.s, s_idx, axis = 0),
            q = tf.gather(previous_auxiliary_smoothing_results.q, s_idx, axis = 0))
        updated_auxiliary_smoothing_results = tf.nest.map_structure(
            lambda x, y: x.write(s_idx_aux, y), updated_auxiliary_smoothing_results, auxiliary_smoothing_result)
        return s_idx + 1, s_idx_aux+1, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results


    def update_smoothing_results(s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results):
        s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results = tf.cond(
            tf.gather(filtered_variance_estimate, s_idx, axis = 0) < eps,
            true_fn = lambda: update_final_results(s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results),
            false_fn = lambda: update_auxilliary_results(s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array,
                                                         updated_auxiliary_smoothing_results)
        )
        return s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results

    s_idx, s_idx_aux, s_idx_fin, new_smoothing_results_array, updated_auxiliary_smoothing_results = tf.while_loop(
        cond = lambda s_idx, *_: tf.math.less(s_idx, size_prev_s_set),
        body = update_smoothing_results,
        loop_vars = [0, 0, 0, new_smoothing_results_array, updated_auxiliary_smoothing_results],
        parallel_iterations = 1000
    )

    #update at the last step
    def update_last_step(q_, new_smoothing_results_array):
        current_s = step
        current_q = q_
        final_mean_estimate = tf.einsum('i,is->s', filter_weights, new_psi)
        final_smoothing_result = SmoothingResults(
            functional_estimates = tf.cast(tf.gather(final_mean_estimate, current_q, axis = 0), dtype = tf.float32),
            s = current_s,
            t = step,
            q = current_q,
            idx = current_s * q+ current_q)
        new_smoothing_results_array = tf.nest.map_structure(
            lambda x, y: x.write(s_idx_fin+current_q, y), new_smoothing_results_array, final_smoothing_result)
        return q_+1, new_smoothing_results_array

    #at the final timestep get the smoothing estimates at the current step
    q_ = tf.cond(step == num_timesteps-1, lambda: 0, lambda: q)
    q_ ,new_smoothing_results_array = tf.while_loop(
        cond = lambda q_, *_: tf.math.less(q_, q),
        body = update_last_step,
        loop_vars = [q_, new_smoothing_results_array],
        parallel_iterations = 1000
    )

    return new_smoothing_results_array, updated_auxiliary_smoothing_results


# def final_smoothing_step(loop_results, step, q, previous_particle_idx_to_keep):
#     #previous_particle_idx_to_keep = loop_results.previous_step_results.previous_particle_idx_to_keep
#     previous_particle_idx_to_keep = tf.boolean_mask(previous_particle_idx_to_keep,
#                                                     previous_particle_idx_to_keep>=0)
#     filter_weights = tf.math.exp(loop_results.previous_step_results.log_weights)
#     #auxiliary_smoothing_results = tf.nest.map_structure(
#     #    lambda x: tf.boolean_mask(x.stack(), loop_results.previous_auxiliary_smoothing_results.s.stack() > -1),
#     #    loop_results.previous_auxiliary_smoothing_results)
#     auxiliary_smoothing_results = tf.nest.map_structure(
#         lambda x: x.stack(),
#         loop_results.previous_auxiliary_smoothing_results)
#     size_prev_s_set = tf.shape(auxiliary_smoothing_results.q)[0]
#     psi = tf.gather(auxiliary_smoothing_results.psi, previous_particle_idx_to_keep, axis = 1)
#     filter_weights = tf.gather(filter_weights, previous_particle_idx_to_keep)
#     filtered_mean_estimate = tf.einsum('i,si->s', filter_weights, psi)
#
#     new_smoothing_results_array = tf.nest.map_structure(
#         lambda x: tf.TensorArray(dtype = x.dtype,
#                                  size = 0, dynamic_size = True),
#         SmoothingResults(functional_estimates = tf.zeros([]),
#                                              s = tf.zeros([], dtype = tf.int32),
#                                              t = tf.zeros([], dtype = tf.int32),
#                                              q = tf.zeros([], dtype = tf.int32),
#                                              idx = tf.zeros([], dtype = tf.int32)))
#
#
#
#     def write_final_results(s_idx, new_smoothing_results_array, updated_auxiliary_smoothing_results):
#         current_s = tf.gather(auxiliary_smoothing_results.s, s_idx, axis = 0)
#         current_q = tf.gather(auxiliary_smoothing_results.q, s_idx, axis = 0)
#         final_smoothing_result = SmoothingResults(
#             functional_estimates = tf.gather(filtered_mean_estimate, s_idx, axis = 0),
#             s = current_s,
#             t = step,
#             q = current_q,
#             idx = current_s * q + current_q)
#         new_smoothing_results_array = tf.nest.map_structure(
#             lambda x, y: x.write(s_idx, y), new_smoothing_results_array, final_smoothing_result)
#         return s_idx + 1, new_smoothing_results_array, updated_auxiliary_smoothing_results
#
#     s_idx, new_smoothing_results_array, updated_auxiliary_smoothing_results = tf.while_loop(
#         cond = lambda s_idx, *_: s_idx < size_prev_s_set,
#         body = write_final_results,
#         loop_vars = [0, new_smoothing_results_array, loop_results.previous_auxiliary_smoothing_results]
#     )
#
#     new_accumulated_smoothing_results = tf.nest.map_structure(
#       lambda x, y: x.scatter(new_smoothing_results_array.idx.stack(), y.stack()),
#       loop_results.accumulated_smoothing_results, new_smoothing_results_array)
#
#
#     return loop_results._replace(accumulated_smoothing_results = new_accumulated_smoothing_results)
#
#
