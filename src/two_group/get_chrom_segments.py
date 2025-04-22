import pandas as pd
from absl import flags, app

flags.DEFINE_string(
    'input_file',
    default='positions.txt',
    help="Path to the input file (gzipped) containing chromosome position data.")
flags.DEFINE_integer(
    'chrom',
    default=22,
    help="Chromosome number to process.")
flags.DEFINE_integer(
    'segment_size',
    default=100000,
    help="Size of the selected chromosome segment (in CpG sites).")
flags.DEFINE_string(
    'output_csv',
    default='chrom_segments.csv',
    help="Path to the output CSV file containing segment information.")

FLAGS = flags.FLAGS

def main(argv):
    del argv  # unused

    # Load positions for the specified chromosome
    positions = pd.read_csv(FLAGS.input_file, header=None, names=["position"], compression="gzip")
    num_positions = len(positions)

    # Calculate the number of segments
    num_segments = 1 + num_positions // FLAGS.segment_size

    # Generate segment indices for the chromosome
    segments = [{'chrom': FLAGS.chrom, 'segment_index': i} for i in range(0, num_segments)]

    # Save the results to the output CSV
    pd.DataFrame(segments).to_csv(FLAGS.output_csv, index=False)
    print(f"Segment information saved to {FLAGS.output_csv}")

if __name__ == '__main__':
    app.run(main)
