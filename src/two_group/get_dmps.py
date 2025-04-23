#########
#Get DMPs for gold standard aging or T1D data using particle filter outputs subect to some expected FDR
#########

import pandas as pd
import numpy as np
import os
import sys
from absl import flags
import pdb
from multiple_testing import FDR_procedure, weighted_FDR_procedure
from pathlib import Path

flags.DEFINE_multi_float("fdr_thresholds",
                     default=[.01, .05],
                     help = "fdr threshold for selecting DMPs.")
flags.DEFINE_string(
    'results_dir',
    default=os.path.join(Path(os.getcwd()).parents[0],'test'),
    help="Directory for the results of the two-group algorithms.")
flags.DEFINE_string(
    'output_dir',
    default=os.path.join(Path(os.getcwd()).parents[0],'test', 'dmp'),
    help="Directory for the outputs of this script.")
flags.DEFINE_integer(
    'n_regimes',
    default=6,
    help="number of regimes.")
flags.DEFINE_integer("chrom",
                      default=22,
                      help="which chormosome to analyse")

FLAGS = flags.FLAGS
FLAGS(sys.argv)

n_regimes = FLAGS.n_regimes
output_dir = FLAGS.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#test for difference in mean methylation signal
def dm_state(control, case):
  return (control)!=(case)


positions = []
position_diffs = []
test_statistics_split = []
test_statistics_hypo = []
test_statistics_hyper = []
test_stats_regimes = {}
for i in range(n_regimes):
  for j in range(n_regimes):
    if i != j:
      test_stats_regimes[(i,j)] = []

chrom = FLAGS.chrom
path = FLAGS.results_dir
control_regimes_ = pd.read_csv(os.path.join(path, 'control_regimes_chrom_{}.csv.gz'.format(chrom)), sep = '\t')
control_regimes_ = (control_regimes_.set_index('pos')).to_numpy()
case_regimes_ = pd.read_csv(os.path.join(path, 'case_regimes_chrom_{}.csv.gz'.format(chrom)), sep = '\t')
case_regimes_ = (case_regimes_.set_index('pos')).to_numpy()
num_particles = control_regimes_.shape[-1]
test_statistics_split.append(
  1. - np.sum(dm_state(control_regimes_, case_regimes_), axis = 1) / num_particles)

for i in range(n_regimes):
  for j in range(n_regimes):
    if i != j:
      test_stats_regimes[(i,j)].append(1 - np.sum((control_regimes_==i) * (case_regimes_==j), axis = 1) / num_particles)


split_probs_ = pd.read_csv(os.path.join(path, 'split_probs_{}.csv.gz'.format(chrom)), sep = '\t')
split_probs_ = split_probs_.set_index('pos')
position_diffs_ = 1/3*(pd.DataFrame(split_probs_.index).diff(1)+pd.DataFrame(
  split_probs_.index).diff(2)+pd.DataFrame(split_probs_.index).diff(3))
position_diffs.append(position_diffs_)
positions_ = pd.DataFrame(split_probs_.index)
positions_['chrom'] = chrom
positions.append(positions_)

positions = pd.concat(positions)
position_diffs = pd.concat(position_diffs)
test_statistics_split = (np.concatenate(test_statistics_split))

for i in range(n_regimes):
  for j in range(n_regimes):
    if i != j:
      test_stats_regimes[(i,j)] = np.concatenate(test_stats_regimes[(i,j)])

#construct weightings
false_positive_weights = np.ones([position_diffs.shape[0]])
false_negative_weights = 1.* np.ones_like(position_diffs.to_numpy())
false_negative_weights[position_diffs < 1000] = 3.
false_negative_weights[position_diffs < 100] = 10.
false_negative_weights = np.squeeze(false_negative_weights, -1)


for fdr_threshold in FLAGS.fdr_thresholds:
  #dmp
  k, Qk, threshold = FDR_procedure(test_statistics_split, fdr_threshold)
  dmp_indicator = test_statistics_split < threshold
  dmp_stats = test_statistics_split[dmp_indicator]
  dmp_pos = positions.to_numpy()[dmp_indicator]
  dmp = pd.DataFrame({'chrom': dmp_pos[:,1], 'position': dmp_pos[:,0], 'null_stats': dmp_stats})
  dmp['false_negative_weight'] = 1.
  dmp.to_csv(os.path.join(output_dir, 'dmp_{}.csv'.format(fdr_threshold)), index=False)

  for i in range(n_regimes):
    for j in range(n_regimes):
      if i != j:
        k, Qk, threshold = FDR_procedure(test_stats_regimes[(i,j)], fdr_threshold)
        dmp_indicator = test_stats_regimes[(i,j)] < threshold
        dmp_stats = test_stats_regimes[(i,j)][dmp_indicator]
        dmp_pos = positions.to_numpy()[dmp_indicator]
        dmp = pd.DataFrame({'chrom': dmp_pos[:,1], 'position': dmp_pos[:,0], 'null_stats': dmp_stats})
        dmp['false_negative_weight'] = 1.
        dmp.to_csv(os.path.join(output_dir, 'dmp_{}_{}_{}.csv'.format(i,j,fdr_threshold)), index=False)



  #weighted versions
  dmp_index, Nk = weighted_FDR_procedure(test_statistics_split, fdr_threshold = fdr_threshold,
                                         weights_false_negatives = false_negative_weights,
                                         weights_false_positives = false_positive_weights)
  dmp_index = np.sort(dmp_index)
  dmp_stats = test_statistics_split[dmp_index]
  dmp_pos = positions.to_numpy()[dmp_index]
  dmp = pd.DataFrame({'chrom': dmp_pos[:,1], 'position': dmp_pos[:,0], 'null_stats': dmp_stats})
  dmp['false_negative_weight'] = false_negative_weights[dmp_index]
  dmp.to_csv(os.path.join(output_dir, 'weighted_dmp_{}.csv'.format(fdr_threshold)), index=False)

  for i in range(n_regimes):
    for j in range(n_regimes):
      if i != j:
        dmp_index, Nk = weighted_FDR_procedure(test_stats_regimes[(i, j)], fdr_threshold = fdr_threshold,
                                               weights_false_negatives = false_negative_weights,
                                               weights_false_positives = false_positive_weights)
        dmp_index = np.sort(dmp_index)
        dmp_stats = test_stats_regimes[(i, j)][dmp_index]
        dmp_pos = positions.to_numpy()[dmp_index]
        dmp = pd.DataFrame({'chrom': dmp_pos[:,1], 'position': dmp_pos[:,0], 'null_stats': dmp_stats})
        dmp['false_negative_weight'] = false_negative_weights[dmp_index]
        dmp.to_csv(os.path.join(output_dir, 'weighted_dmp_{}_{}_{}.csv'.format(i,j,fdr_threshold)), index=False)


