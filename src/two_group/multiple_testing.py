import numpy as np

def FDR_procedure(test_statistics, fdr_threshold):
  ordered_test_statistics = np.sort(test_statistics)
  Qs=1./np.linspace(1,stop=len(ordered_test_statistics), num=len(ordered_test_statistics)
                    )*np.cumsum(ordered_test_statistics)
  s = np.sum(Qs<=fdr_threshold)
  if fdr_threshold < ordered_test_statistics[0]:
    return 0, 0., 0.
  if s == test_statistics.shape:
    return test_statistics.shape, Qs[s-1], 1.01
  return s, Qs[s-1], ordered_test_statistics[s]

def weighted_FDR_procedure(test_statistics, fdr_threshold, weights_false_positives, weights_false_negatives):
  ranking = weights_false_positives*(test_statistics-fdr_threshold)/(
    weights_false_negatives * (1 - test_statistics) + weights_false_positives * abs(test_statistics - fdr_threshold))
  ranking_indices = np.argsort(ranking)
  excessive_error_rates = weights_false_positives * (test_statistics - fdr_threshold)
  ranked_excessive_error_rates = excessive_error_rates[ranking_indices]
  Nsums = np.cumsum(ranked_excessive_error_rates)
  s = np.sum(Nsums <= 0)
  return ranking_indices[:s], Nsums[s-1]

