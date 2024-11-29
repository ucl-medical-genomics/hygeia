######################
#This script converts computes the number of segments (batches) that are required to be 
#run per chromosome based on the chosen segment size
######################

import os
import pandas as pd
import numpy as np
from absl import flags, app
import sys
from pathlib import Path

flags.DEFINE_string(
    'data_dir',
    default=os.path.join('data'),
    help="Directory of the read data.")
flags.DEFINE_integer("segment_size",
                      default=100000,
                      help="size of the selected chromosome segment (in CpG sites)")

FLAGS = flags.FLAGS

def main(argv):
    number_of_segments = {}
    del argv  # unused
    for chrom in range(1,23):
        number_of_positions = pd.read_table(os.path.join(FLAGS.data_dir,
                 'positions_{}.txt'.format(chrom)), sep = ',', header = None).shape[0]
        number_of_segments[chrom] = 1+ number_of_positions//FLAGS.segment_size
    pd.DataFrame(number_of_segments, index=['Segments']).transpose().to_csv(
       os.path.join(FLAGS.data_dir,'segment_numbers.csv'))

if __name__ == '__main__':
  app.run(main)
