#!/usr/bin/env python3
"""
Methylation Data Processing Script

This script converts output from gemBS to counts data using reference genome data
to get all CpG sites. The gemBS outputs are counts per sample/donor and not
separated by chromosome.
"""

import os
import logging
from pathlib import Path
from typing import List, Tuple
import pandas as pd
import numpy as np
from absl import flags, app

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Define command-line flags
flags.DEFINE_string(
    'cpg_file_path',
    default=None,
    help="Path to file containing all CpG sites."
)
flags.DEFINE_string(
    'output_path',
    default=os.path.join(Path.cwd().parent, 'test'),
    help="Directory where to store the results."
)
flags.DEFINE_multi_string(
    'case_data_path',
    default=[],
    help="Paths for the methylation data of the case group."
)
flags.DEFINE_multi_string(
    'case_id_names',
    default=[],
    help="Names of case IDs in methylation files."
)
flags.DEFINE_multi_string(
    'control_data_path',
    default=[],
    help="Paths for the methylation data of the control group."
)
flags.DEFINE_multi_string(
    'control_id_names',
    default=[],
    help="Names of control IDs in methylation files."
)
flags.DEFINE_integer(
    "chromosome",
    default=22,
    help="Which chromosome to analyze (use chromosome instead of chrom for clarity)"
)
flags.DEFINE_boolean(
    "verbose",
    default=False,
    help="Enable verbose logging"
)

FLAGS = flags.FLAGS


class MethylationProcessor:
    """Process methylation data from gemBS output."""

    def __init__(self, cpg_file_path: str, output_path: str, chromosome: int):
        """
        Initialize the methylation processor.

        Args:
            cpg_file_path: Path to CpG sites file
            output_path: Output directory path
            chromosome: Chromosome number to analyze
        """
        self.cpg_file_path = Path(cpg_file_path)
        self.output_path = Path(output_path)
        self.chromosome = chromosome
        self.cpg_sites = None

        self._validate_inputs()

    def _validate_inputs(self) -> None:
        """Validate input parameters."""
        if not self.cpg_file_path.exists():
            raise FileNotFoundError(f"CpG file not found: {self.cpg_file_path}")

        if not (1 <= self.chromosome <= 23):  # Including X chromosome as 23
            raise ValueError(f"Invalid chromosome number: {self.chromosome}")

        self.output_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {self.output_path}")

    def load_cpg_sites(self) -> pd.DataFrame:
        """Load CpG sites from file."""
        logger.info(f"Loading CpG sites from: {self.cpg_file_path}")

        try:
            compression = 'gzip' if self.cpg_file_path.suffix == '.gz' else None
            self.cpg_sites = pd.read_csv(self.cpg_file_path, sep='\t', compression=compression)

            # Filter by chromosome
            chrom_name = f'chr{self.chromosome}'
            cpg_sites_chrom = self.cpg_sites[self.cpg_sites['seqID'] == chrom_name].copy()

            if cpg_sites_chrom.empty:
                raise ValueError(f"No CpG sites found for chromosome {self.chromosome}")

            logger.info(f"Found {len(cpg_sites_chrom)} CpG sites for chromosome {self.chromosome}")
            return cpg_sites_chrom

        except Exception as e:
            logger.error(f"Error loading CpG sites: {e}")
            raise

    def process_sample_data(self, data_paths: List[str], id_names: List[str],
                          cpg_sites_chrom: pd.DataFrame) -> pd.DataFrame:
        """
        Process methylation data for a group of samples.

        Args:
            data_paths: List of paths to methylation data files
            id_names: List of sample ID names
            cpg_sites_chrom: CpG sites for the chromosome

        Returns:
            Processed methylation data
        """
        if len(data_paths) != len(id_names):
            raise ValueError("Number of data paths must match number of ID names")

        # Initialize with positions
        meth_data = pd.DataFrame({'Pos0': cpg_sites_chrom['start'] - 1})

        for data_path, sample_id in zip(data_paths, id_names):
            logger.info(f"Processing sample: {sample_id}")

            try:
                if not Path(data_path).exists():
                    logger.error(f"File not found: {data_path}")
                    continue

                sample_data = pd.read_csv(data_path, sep='\t', compression='gzip')

                chrom_data = sample_data[
                    (sample_data['Contig'] == f'chr{self.chromosome}') &
                    (sample_data['Ref'] == 'CG')
                ].copy()

                if chrom_data.empty:
                    logger.warning(f"No CpG data found for {sample_id} on chromosome {self.chromosome}")
                    # Add empty columns for this sample to maintain structure
                    meth_data[f'{sample_id}:non_conv'] = np.nan
                    meth_data[f'{sample_id}:conv'] = np.nan
                    continue

                required_cols = ['Pos0', f'{sample_id}:non_conv', f'{sample_id}:conv']

                # Check if required columns exist
                missing_cols = [col for col in required_cols if col not in chrom_data.columns]
                if missing_cols:
                    logger.error(f"Missing columns in {sample_id}: {missing_cols}")
                    # Add empty columns for this sample
                    meth_data[f'{sample_id}:non_conv'] = np.nan
                    meth_data[f'{sample_id}:conv'] = np.nan
                    continue

                sample_subset = chrom_data[required_cols].copy()

                # Merge with main data
                meth_data = pd.merge(meth_data, sample_subset, on='Pos0', how='outer')

            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {e}")
                # Add empty columns for this sample to maintain structure
                meth_data[f'{sample_id}:non_conv'] = np.nan
                meth_data[f'{sample_id}:conv'] = np.nan
                continue

        return meth_data.sort_values('Pos0').reset_index(drop=True)

    def extract_count_arrays(self, merged_data: pd.DataFrame,
                           n_control: int, n_case: int) -> Tuple[np.ndarray, ...]:
        """
        Extract count arrays from merged data.

        Args:
            merged_data: Merged methylation data
            n_control: Number of control samples
            n_case: Number of case samples

        Returns:
            Tuple of arrays: (positions, non_conv_control, conv_control,
                            non_conv_case, conv_case)
        """
        # Convert to numpy and handle NaN values
        positions = merged_data['Pos0'].to_numpy()

        # Handle case where we only have positions (no sample data)
        if merged_data.shape[1] == 1:  # Only Pos0 column
            empty_array = np.array([]).reshape(len(positions), 0)
            return positions, empty_array, empty_array, empty_array, empty_array

        data_array = merged_data.drop('Pos0', axis=1).to_numpy()
        data_array = np.nan_to_num(data_array, copy=False)

        # Initialize empty arrays
        empty_array = np.array([]).reshape(len(positions), 0)
        non_conv_control = empty_array.copy()
        conv_control = empty_array.copy()
        non_conv_case = empty_array.copy()
        conv_case = empty_array.copy()

        # Extract control data if present
        if n_control > 0:
            control_end_idx = n_control * 2
            if data_array.shape[1] >= control_end_idx:
                non_conv_control = data_array[:, 0:control_end_idx:2]
                conv_control = data_array[:, 1:control_end_idx:2]

                # Extract case data if present
                if n_case > 0 and data_array.shape[1] > control_end_idx:
                    non_conv_case = data_array[:, control_end_idx::2]
                    conv_case = data_array[:, control_end_idx + 1::2]
        else:
            # Only case data present
            if n_case > 0:
                non_conv_case = data_array[:, 0::2]
                conv_case = data_array[:, 1::2]

        return positions, non_conv_control, conv_control, non_conv_case, conv_case

    def save_results(self, positions: np.ndarray, non_conv_control: np.ndarray,
                    conv_control: np.ndarray, non_conv_case: np.ndarray,
                    conv_case: np.ndarray) -> int:
        """Save processed results to files."""
        logger.info("Saving results...")

        # Number of CpG sites
        cpg_sites_merged = len(positions)

        # Always save positions and count
        file_mapping = {
            'positions': positions,
            'cpg_sites_merged': np.array([cpg_sites_merged])
        }

        # Add control data if present
        if non_conv_control.size > 0:
            n_total_reads_control = conv_control + non_conv_control
            file_mapping.update({
                'n_methylated_reads_control': non_conv_control,
                'n_total_reads_control': n_total_reads_control,
            })
            logger.info(f"Control samples: {non_conv_control.shape[1]} samples")

        # Add case data if present
        if non_conv_case.size > 0:
            n_total_reads_case = conv_case + non_conv_case
            file_mapping.update({
                'n_methylated_reads_case': non_conv_case,
                'n_total_reads_case': n_total_reads_case,
            })
            logger.info(f"Case samples: {non_conv_case.shape[1]} samples")

        # Save all arrays
        for filename, data in file_mapping.items():
            filepath = self.output_path / f'{filename}_{self.chromosome}.txt.gz'
            try:
                np.savetxt(filepath, data, delimiter=",", fmt='%s')
                logger.info(f"Saved: {filepath} (shape: {data.shape})")
            except Exception as e:
                logger.error(f"Error saving {filepath}: {e}")

        logger.info(f"Total CpG sites processed: {cpg_sites_merged}")
        return cpg_sites_merged

    def process(self, case_data_paths: List[str], case_id_names: List[str],
               control_data_paths: List[str], control_id_names: List[str]) -> int:
        """
        Main processing function.

        Args:
            case_data_paths: Paths to case data files
            case_id_names: Case sample IDs
            control_data_paths: Paths to control data files
            control_id_names: Control sample IDs

        Returns:
            Number of CpG sites processed
        """
        logger.info(f"Starting methylation analysis for chromosome {self.chromosome}")
        cpg_sites_chrom = self.load_cpg_sites()

        # Initialize data with positions
        merged_data = pd.DataFrame({'Pos0': cpg_sites_chrom['start'] - 1})

        # Process control data if provided
        if control_data_paths and control_id_names:
            logger.info("Processing control samples...")
            control_data = self.process_sample_data(
                control_data_paths, control_id_names, cpg_sites_chrom
            )
            merged_data = pd.merge(merged_data, control_data, on='Pos0', how='outer')
        else:
            logger.info("No control samples provided - processing case samples only")

        # Process case data if provided
        if case_data_paths and case_id_names:
            logger.info("Processing case samples...")
            case_data = self.process_sample_data(
                case_data_paths, case_id_names, cpg_sites_chrom
            )
            merged_data = pd.merge(merged_data, case_data, on='Pos0', how='outer')
        else:
            logger.info("No case samples provided - processing control samples only")

        # Sort final data
        merged_data = merged_data.sort_values('Pos0').reset_index(drop=True)

        # Extract count arrays
        positions, non_conv_control, conv_control, non_conv_case, conv_case = \
            self.extract_count_arrays(merged_data, len(control_id_names), len(case_id_names))

        # Save results
        cpg_sites_count = self.save_results(
            positions, non_conv_control, conv_control, non_conv_case, conv_case
        )

        logger.info("Processing completed successfully!")
        return cpg_sites_count


def validate_flags() -> None:
    """Validate command-line flags."""
    # Required flag
    if not FLAGS.cpg_file_path:
        raise ValueError("Required flag --cpg_file_path not provided")

    # Check that we have at least one group of samples
    has_case = bool(FLAGS.case_data_path and FLAGS.case_id_names)
    has_control = bool(FLAGS.control_data_path and FLAGS.control_id_names)

    if not has_case and not has_control:
        raise ValueError("Must provide either case samples, control samples, or both")

    # Validate individual groups if provided
    if FLAGS.case_data_path and FLAGS.case_id_names:
        if len(FLAGS.case_data_path) != len(FLAGS.case_id_names):
            raise ValueError("Number of case data paths must match number of case ID names")

    if FLAGS.control_data_path and FLAGS.control_id_names:
        if len(FLAGS.control_data_path) != len(FLAGS.control_id_names):
            raise ValueError("Number of control data paths must match number of control ID names")


def main(argv):
    """Main function."""
    del argv  # Unused

    # Set logging level
    if FLAGS.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        validate_flags()
        processor = MethylationProcessor(FLAGS.cpg_file_path, FLAGS.output_path, FLAGS.chromosome)
        cpg_count = processor.process(
            FLAGS.case_data_path,
            FLAGS.case_id_names,
            FLAGS.control_data_path,
            FLAGS.control_id_names
        )
        print(f"Successfully processed {cpg_count} CpG sites for chromosome {FLAGS.chromosome}")
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        return 1
    return 0

if __name__ == '__main__':
    app.run(main)
