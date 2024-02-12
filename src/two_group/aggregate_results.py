#########
#Get two-group results from different region partitions and runs
#########

import pandas as pd
import numpy as np
import os
import sys
from absl import flags
import pdb


flags.DEFINE_string(
    'data_dir',
    default=os.path.join(Path(os.getcwd()).parents[0], 'test'),
    help="Directory of the methylation data.")
flags.DEFINE_string(
    'results_dir',
    default=os.path.join(Path(os.getcwd()).parents[0],'test'),
    help="Directory for the results of the two-group algorithms.")
flags.DEFINE_string(
    'output_dir',
    default=os.path.join(Path(os.getcwd()).parents[0],'test', 'results'),
    help="Directory for the outputs of this script.")
flags.DEFINE_integer(
    'seeds',
    default=10,
    help="Number of seeds that algorithms have been run.")
flags.DEFINE_integer(
    'chrom',
    default=1,
    help="Chromosome number to process")
flags.DEFINE_bool('batches', default = False, help='whether the particle filter results are saved in batches')
flags.DEFINE_integer("num_particles",
                     default=1200,
                     help="number of particles from backward smoothing")
FLAGS = flags.FLAGS
FLAGS(sys.argv)

N = FLAGS.num_particles

if not os.path.exists(FLAGS.output_dir):
    os.makedirs(FLAGS.output_dir)
output_dir = FLAGS.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

chrom = FLAGS.chrom
print(chrom)

split_probs_ = []
regime_probs_ = []
positions_ = []
merge_states_ = []
control_regimes_= []
case_regimes_ = []
case_durations_ = []
control_durations_ = []

n_total_reads_control_ = []
n_total_reads_case_ = []
observations_control_ = []
observations_case_ = []


positions__ = pd.read_table(
os.path.join(FLAGS.data_dir, 'positions_{}.txt'.format(chrom)), sep = ' ', header = None, dtype = np.int64)

n_total_reads_control__ = pd.read_table(
os.path.join(FLAGS.data_dir, 'n_total_reads_control_{}.txt'.format(chrom)), sep = ' ', header = None)
n_total_reads_case__ = pd.read_table(
os.path.join(FLAGS.data_dir, 'n_total_reads_case_{}.txt'.format(chrom)), sep = ' ', header = None)
observations_control__ = pd.read_table(
os.path.join(FLAGS.data_dir, 'n_methylated_reads_control_{}.txt'.format(chrom)), sep = ' ', header = None)
observations_case__ = pd.read_table(
os.path.join(FLAGS.data_dir, 'n_methylated_reads_case_{}.txt'.format(chrom)), sep = ' ', header = None)

merged_states__ = []  # np.empty((positions__.shape[0], 0))
control_states__ = []  # np.empty((positions__.shape[0], 0, 2))
case_states__ = []  # np.empty((positions__.shape[0], 0, 2))
control_durations__ = []  # np.empty((positions__.shape[0], 0, 2))
case_durations__ = []  # np.empty((positions__.shape[0], 0, 2))

for seed in range(0, FLAGS.seeds):

    merged_states = np.load(os.path.join(
    FLAGS.results_dir, 'chrom_{}__id_{}'.format(chrom, id),
    'optimal_backward_particles_merged_state_{}_{}.npy'.format(N, seed)))
    control_states = np.load(os.path.join(
    FLAGS.results_dir, 'chrom_{}__id_{}'.format(chrom, id),
    'optimal_backward_particles_control_state_{}_{}.npy'.format(N, seed)))
    case_states = np.load(os.path.join(
    FLAGS.results_dir, 'chrom_{}__id_{}'.format(chrom, id),
    'optimal_backward_particles_case_state_{}_{}.npy'.format(N, seed)))
    merged_states__.append(merged_states)
    control_states__.append(control_states)
    case_states__.append(case_states)


merged_states__ = np.concatenate(merged_states__, -1)
control_states__ = np.concatenate(control_states__, axis = 1)
case_states__ = np.concatenate(case_states__, axis = 1)

split_probs__ = np.mean(merged_states__ == 0, axis = 1)

start_ind = 0

end_ind = split_probs__.shape[0]
split_probs_.append(pd.DataFrame(split_probs__[start_ind:end_ind]))
positions_.append(positions__.iloc[start_ind:end_ind])

merge_states_.append(pd.DataFrame(merged_states__[start_ind:end_ind]).astype(np.int8))
control_regimes_.append(pd.DataFrame(control_states__[start_ind:end_ind, :, 1]).astype(np.int8))
case_regimes_.append(pd.DataFrame(case_states__[start_ind:end_ind, :, 1]).astype(np.int8))
control_durations_.append(pd.DataFrame(control_states__[start_ind:end_ind, :, 0]).astype(np.int16))
case_durations_.append(pd.DataFrame(case_states__[start_ind:end_ind, :, 0]).astype(np.int16))

n_total_reads_control_.append(pd.DataFrame(n_total_reads_control__[start_ind:end_ind]).astype(np.int16))
n_total_reads_case_.append(pd.DataFrame(n_total_reads_case__[start_ind:end_ind]).astype(np.int16))
observations_control_.append(pd.DataFrame(observations_control__[start_ind:end_ind]).astype(np.int16))
observations_case_.append(pd.DataFrame(observations_case__[start_ind:end_ind]).astype(np.int16))

split_probs_ = pd.concat(split_probs_)
split_probs_['chrom'] = chrom



#save per chromosome outputs
positions_chrom = pd.concat(positions_)
positions_chrom = positions_chrom.rename(columns = {0: 'pos'})
positions_chrom = positions_chrom.astype(np.int32)
control_regimes_chrom = pd.concat(control_regimes_)
control_regimes_chrom = control_regimes_chrom.set_index(positions_chrom['pos'])
control_regimes_freq = control_regimes_chrom.apply(lambda x: x.value_counts(normalize=True), 1)
control_regimes_freq.to_csv(
os.path.join(output_dir, 'control_regimes_freq_{}.csv'.format(chrom)), sep = '\t')
control_regimes_chrom.to_csv(
os.path.join(output_dir, 'control_regimes_chrom_{}.csv'.format(chrom)), sep = '\t')
case_regimes_chrom = pd.concat(case_regimes_)
case_regimes_chrom = case_regimes_chrom.set_index(positions_chrom['pos'])
case_regimes_chrom.to_csv(
os.path.join(output_dir, 'case_regimes_chrom_{}.csv'.format(chrom)), sep = '\t')
case_regimes_freq = case_regimes_chrom.apply(lambda x: x.value_counts(normalize=True), 1)
case_regimes_freq.to_csv(
os.path.join(output_dir, 'case_regimes_freq_{}.csv'.format(chrom)), sep = '\t')
#pdb.set_trace()
merge_states_chrom = pd.concat(merge_states_)
merge_states_chrom = merge_states_chrom.set_index(positions_chrom['pos'])
merge_states_chrom.to_csv(
os.path.join(output_dir, 'merge_states_chrom_{}.csv'.format(chrom)), sep = '\t')
split_probs_chrom = np.mean(merge_states_chrom==0, axis=1)
split_probs_chrom.to_csv(
os.path.join(output_dir, 'split_probs_{}.csv'.format(chrom)), sep = '\t')

n_total_reads_control_chrom = pd.concat(n_total_reads_control_).set_index(positions_chrom['pos'])
n_total_reads_control_chrom.to_csv(
os.path.join(output_dir, 'n_total_reads_control_chrom_{}.csv'.format(chrom)), sep = '\t')
n_total_reads_case_chrom = pd.concat(n_total_reads_case_).set_index(positions_chrom['pos'])
n_total_reads_case_chrom.to_csv(
os.path.join(output_dir, 'n_total_reads_case_chrom_{}.csv'.format(chrom)), sep = '\t')

n_meth_reads_control_chrom = pd.concat(observations_control_).set_index(positions_chrom['pos'])
n_meth_reads_control_chrom.to_csv(
os.path.join(output_dir, 'n_meth_reads_control_chrom_{}.csv'.format(chrom)), sep = '\t')
n_meth_reads_case_chrom = pd.concat(observations_case_).set_index(positions_chrom['pos'])
n_meth_reads_case_chrom.to_csv(
os.path.join(output_dir, 'n_meth_reads_case_chrom_{}.csv'.format(chrom)), sep = '\t')

control_durations_chrom = pd.concat(control_durations_)
control_durations_chrom = control_durations_chrom.set_index(positions_chrom['pos'])
control_durations_chrom.to_csv(
os.path.join(output_dir, 'control_durations_chrom_{}.csv'.format(chrom)), sep = '\t')
case_durations_chrom = pd.concat(case_durations_)
case_durations_chrom = case_durations_chrom.set_index(positions_chrom['pos'])
case_durations_chrom.to_csv(
os.path.join(output_dir, 'case_durations_chrom_{}.csv'.format(chrom)), sep = '\t')

