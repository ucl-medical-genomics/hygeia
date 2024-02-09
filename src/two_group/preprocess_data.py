######################
#This script converts output from gemBS to counts data
#It also uses reference genome data to get all CpG sites
#The gemBS outputs are counts per sample/donor and not seperated by chormosome
######################

import os
import pandas as pd
import numpy as np
from absl import flags
import sys
from pathlib import Path

flags.DEFINE_string(
    'cpg_file_path',
    default=None,
    help="path to file containing all CpG sites.")
flags.DEFINE_string(
    'output_path',
    default=os.path.join(Path(os.getcwd()).parents[0], 'test'),
    help="Directory where to store the results.")
flags.DEFINE_multi_string(
    'case_data_path',
    default=None,
    help="paths for the methylation data of the case group.")
flags.DEFINE_multi_string(
    'case_id_names',
    default=None,
    help="name of case ids in methylation files.")
flags.DEFINE_multi_string(
    'control_data_path',
    default=None,
    help="paths for the methylation data of the control group.")
flags.DEFINE_multi_string(
    'control_id_names',
    default=None,
    help="name of control ids in methylation files.")

FLAGS = flags.FLAGS
FLAGS(sys.argv)

cpg_sites_compression = 'gzip' if FLAGS.cpg_file_path.endswith('.gz') else False
print (FLAGS.cpg_file_path)
cpg_sites = pd.read_csv(FLAGS.cpg_file_path, sep='\t')

total_cpg_sites_merged = 0

if not os.path.exists(FLAGS.output_path):
  os.makedirs(FLAGS.output_path)

for chromosome in list(range(1,23)):

    #all CpG sites in chromosome
    cpg_sites_chrom = cpg_sites.loc[cpg_sites['seqID'].isin(['chr'+str(chromosome)])]

    # #case data
    meth_data = pd.DataFrame({'Pos0':cpg_sites_chrom['start']-1})
    n_samples_case = len(FLAGS.case_data_path)
    for i in range(n_samples_case):
        data_set_path = FLAGS.case_data_path[i]
        case_id = FLAGS.case_id_names[i]
        data = pd.read_table(data_set_path, compression='gzip')
        per_chromosome_data = data.loc[data.Contig=='chr'+str(chromosome)]
        #only CpG sites (not CpA for example)
        per_chromosome_data = per_chromosome_data.loc[per_chromosome_data['Ref'] == 'CG']
        meth_data_ = per_chromosome_data[['Pos0', case_id + ':non_conv', case_id + ':conv']]
        #merge with all CpG sites data
        meth_data=pd.merge(meth_data, meth_data_, how='outer', left_on='Pos0', right_on='Pos0')
    case_counts = meth_data.sort_values(by='Pos0')

    #control data
    n_samples_control = len(FLAGS.control_data_path)
    meth_data = pd.DataFrame({'Pos0': cpg_sites_chrom['start'] - 1})
    for i in range(n_samples_control):
        data_set_path = FLAGS.control_data_path[i]
        control_id = FLAGS.control_id_names[i]
        data = pd.read_table(data_set_path, compression = 'gzip')
        per_chromosome_data = data.loc[data.Contig == 'chr' + str(chromosome)]
        # only CpG sites (not CpA for example)
        per_chromosome_data = per_chromosome_data.loc[per_chromosome_data['Ref'] == 'CG']
        meth_data_ = per_chromosome_data[['Pos0', control_id + ':non_conv', control_id + ':conv']]
        # merge with all CpG sites data
        meth_data = pd.merge(meth_data, meth_data_, how = 'outer', left_on = 'Pos0', right_on = 'Pos0')

    control_counts = meth_data.sort_values(by='Pos0')

    #combine case and control data
    data = pd.merge(control_counts, case_counts, how = 'outer', left_on = 'Pos0', right_on = 'Pos0')
    data = data.sort_values(by = 'Pos0')

    positions = data.Pos0.to_numpy()
    data = np.nan_to_num(data, copy = False)
    non_conv_control = data[:, 1:(1 + n_samples_control * 2):2]
    conv_control = data[:, 2:(1 + n_samples_control * 2):2]
    non_conv_case = data[:, 1 + n_samples_control * 2::2]
    conv_case = data[:, 2 + n_samples_control * 2::2]
    observations_control = non_conv_control
    observations_case = non_conv_case
    ## number of total reads for the case and control group
    n_observations = data.shape[0]
    ## get the number of reads for both groups
    n_total_reads_control = conv_control + non_conv_control
    n_total_reads_case = conv_case + non_conv_case
    #counts cpg sites
    cpg_sites_merged=data.shape[0]
    total_cpg_sites_merged += cpg_sites_merged
    print(total_cpg_sites_merged)
    #save data
    np.savetxt(os.path.join(FLAGS.output_path,'positions_{}.txt'.format(chromosome)), positions)
    np.savetxt(os.path.join(FLAGS.output_path,'n_methylated_reads_control_{}.txt'.format(chromosome)), observations_control)
    np.savetxt(os.path.join(FLAGS.output_path,'n_methylated_reads_case_{}.txt'.format(chromosome)), observations_case)
    np.savetxt(os.path.join(FLAGS.output_path,'n_total_reads_control_{}.txt'.format(chromosome)), n_total_reads_control)
    np.savetxt(os.path.join(FLAGS.output_path,'n_total_reads_case_{}.txt'.format(chromosome)), n_total_reads_case)
    np.savetxt(os.path.join(FLAGS.output_path,'cpg_sites_merged_{}.txt'.format(chromosome)), [cpg_sites_merged])


np.savetxt(os.path.join(FLAGS.output_path,'total_cpg_sites_merged.txt'), [total_cpg_sites_merged])
