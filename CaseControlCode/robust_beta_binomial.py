#####
# Distribution where samling function is just sampling from a BetaBinomial distribution
# while log_prob is a beta-divergence/Tsallis-score for a BetaBinomial distribution

import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np


class RobustBetaBinomial(tfp.distributions.BetaBinomial):

  def __init__(self,
               total_count,
               concentration1,
               concentration0,
               beta,
               validate_args = False,
               allow_nan_stats = True,
               name = 'BetaBinomial'):
    self._beta = beta
    super(RobustBetaBinomial, self).__init__(
      total_count,
      concentration1,
      concentration0,
      validate_args,
      allow_nan_stats,
      name)

  def _log_prob(self, counts):
    beta_binomial_log_prob = super(RobustBetaBinomial, self)._log_prob(counts)
    beta_binomial_log_probs = super(RobustBetaBinomial, self)._log_prob(
      tf.range(tf.reduce_max(self.total_count))[:, tf.newaxis, tf.newaxis])
    # replace nans by -inf
    beta_binomial_log_probs = tf.where(tf.math.is_nan(beta_binomial_log_probs), -np.inf, beta_binomial_log_probs)
    integral_term = 1. / (self._beta + 1) * tf.math.exp(
      tf.reduce_logsumexp((self._beta + 1.) * beta_binomial_log_probs, 0))
    log_tsallis_score = 1. / self._beta * tf.math.exp(self._beta * beta_binomial_log_prob) - integral_term
    return log_tsallis_score