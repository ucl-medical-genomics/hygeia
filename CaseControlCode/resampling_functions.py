import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
from tensorflow_probability.python.internal import prefer_static
tfd = tfp.distributions

def OptimalFiniteState(log_weights, num_resampled_ancestors, particles):
    sorted_indices = tf.argsort(log_weights, direction = 'DESCENDING')
    sorted_log_weights = tf.gather(log_weights, sorted_indices)
    reverse_cum_sum = tf.math.cumsum(tf.math.exp(sorted_log_weights), reverse=True)

    def fn(k_old, k_new, log_c):
        log_c = tf.math.log(tf.cast(num_resampled_ancestors - k_old, tf.float32)) - tf.math.log(reverse_cum_sum[k_old])
        k_new = k_old + tf.reduce_sum(tf.cast((log_c + (sorted_log_weights[k_old:]))>0, tf.int32))
        return k_new, k_old, log_c

    k_old, k_new, log_c = tf.while_loop(
        cond = lambda k_old, k_new, log_c:  tf.math.logical_and(tf.math.logical_and(
          tf.math.not_equal(k_old,k_new),k_old <tf.shape(sorted_log_weights)[0]),
          k_old<num_resampled_ancestors),
        body = fn,
        loop_vars = [0,-1,-1.]
    )
    #tf.debugging.assert_equal(tf.reduce_sum(
    #    tf.minimum(1., tf.math.exp(log_c) * tf.math.exp(sorted_log_weights))) , tf.cast(num_resampled_ancestors, tf.float32))
    c = tf.math.exp(log_c)
    K = k_new
    #error handling use unbiased scheme then
    K, log_c = tf.cond(K < tf.shape(sorted_indices)[0],
                       lambda: (K, log_c),
                       lambda: (tf.shape(sorted_indices)[0], -np.Inf))
    L = num_resampled_ancestors - K
    deterministic_parent_indices = sorted_indices[:K]
    unnormalised_log_residual_weights = sorted_log_weights[K:]
    #tf.print(sorted_indices)
    resampled_parent_indices = K + SystematicResampling(
        tf.nn.log_softmax(unnormalised_log_residual_weights, axis = 0), L)
    parent_indices = tf.concat([deterministic_parent_indices, tf.gather(sorted_indices,
                               resampled_parent_indices)], 0)
    parent_indices = tf.ensure_shape(parent_indices, [num_resampled_ancestors])

    use_unbiased_weights = tf.cond(tf.math.is_inf(log_c),
                                       lambda: tf.ones([], dtype=tf.bool),
                                       lambda: tf.zeros([], dtype=tf.bool))
    parent_indices, log_c = tf.cond(tf.math.is_inf(log_c),
                                       lambda: (tfd.Categorical(logits =log_weights).sample(num_resampled_ancestors), 0.),
                                       lambda: (parent_indices, log_c))



    return tf.zeros([], dtype=tf.bool), parent_indices, log_c, use_unbiased_weights,\
           tf.ones([num_resampled_ancestors], tf.bool), num_resampled_ancestors



def SystematicResampling(log_weights, num_particles):

    T = (tf.cast(tf.range(num_particles),tf.float32)+tf.random.uniform([]))/tf.cast(num_particles, tf.float32)
    Q = tf.cumsum(tf.exp(log_weights))

    i, j, parent_indices = tf.while_loop(
       cond = lambda i, j, parent_indices: tf.math.logical_and(j < num_particles, i<tf.shape(log_weights)[0]),
       body = lambda i, j, parent_indices: (tf.cond(
           tf.gather(T,j)<=tf.gather(Q,i),
           lambda: (i, j+1, parent_indices.write(j, i)),
           lambda: (i+1, j, parent_indices))),
       loop_vars = [0, 0,  tf.TensorArray(dtype = tf.int32, size=num_particles)])

    return parent_indices.stack()

def UnbiasedResampling(multinomial_resampling, log_weights, num_particles):

    parent_indices = prefer_static.cond(
        multinomial_resampling,
        lambda: tfd.Categorical(logits = log_weights).sample(num_particles),
        lambda: SystematicResampling(log_weights, num_particles)
    )

    return tf.zeros([], dtype=tf.bool), parent_indices, 0., tf.ones([], dtype=tf.bool),\
           tf.ones([num_particles], tf.bool), num_particles