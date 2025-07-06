#!/usr/bin/env python3
"""
Methylation Data Processing Script - BED Format Version

This script converts BED-format methylation output to counts data using reference genome data
to get all CpG sites. The input format is BED-like with methylation information.

This version uses Polars for improved performance and memory efficiency.
"""

import os
import logging
from pathlib import Path
from typing import List, Tuple
import polars as pl
import numpy as np
from absl import flags, app

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

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
    help="Paths for the methylation data of the case group (BED format)."
)
flags.DEFINE_multi_string(
    'case_id_names',
    default=[],
    help="Names of case IDs in methylation files."
)
flags.DEFINE_multi_string(
    'control_data_path',
    default=[],
    help="Paths for the methylation data of the control group (BED format)."
)
flags.DEFINE_multi_string(
    'control_id_names',
    default=[],
    help="Names of control IDs in methylation files."
)
flags.DEFINE_string(
    "chromosome",
    default="22",
    help="The chromosome to analyze (chr22, or 22, as per input file)"
)
flags.DEFINE_boolean(
    "verbose",
    default=False,
    help="Enable verbose logging"
)

FLAGS = flags.FLAGS


class MethylationBEDProcessor:
    """Process methylation data from BED format output using Polars."""

    def __init__(self, cpg_file_path: str, output_path: str, chromosome: str):
        """
        Initialize the methylation processor.

        Args:
            cpg_file_path: Path to CpG sites file
            output_path: Output directory path
            chromosome: Chromosome identifier to analyze (e.g., "22" or "chr22")
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

        self.output_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {self.output_path}")

    def load_cpg_sites(self) -> pl.DataFrame:
        """Load CpG sites from file using Polars."""
        logger.info(f"Loading CpG sites from: {self.cpg_file_path}")

        try:
            self.cpg_sites = pl.read_csv(
                self.cpg_file_path,
                separator='\t',
                infer_schema_length=10000
            )

            cpg_sites_chrom = (
                self.cpg_sites
                .lazy()
                .filter(pl.col('seqID').cast(pl.Utf8) == self.chromosome)
                .collect()
            )

            if cpg_sites_chrom.height == 0:
                raise ValueError(f"No CpG sites found for chromosome {self.chromosome}")

            logger.info(f"Found {cpg_sites_chrom.height} CpG sites for chromosome {self.chromosome}")
            return cpg_sites_chrom

        except Exception as e:
            logger.error(f"Error loading CpG sites: {e}")
            raise

    def read_bed_file(self, file_path: str, sample_id: str) -> pl.DataFrame:
        """
        Read BED format methylation file. BED Format columns:
        1: chr, 2: start, 3: end, 4: name, 5: score, 6: strand,
        7: thickStart, 8: thickEnd, 9: itemRgb, 10: coverage, 11: percent_methylated,
        12: ref_genotype, 13: sample_genotype, 14: quality_score
        """
        logger.info(f"Reading BED file: {file_path}")

        try:
            # Define column names for BED format (we only need the first 14)
            bed_columns = [
                "chr", "start", "end", "name", "score", "strand",
                "thickStart", "thickEnd", "itemRgb", "coverage",
                "percent_methylated", "ref_genotype", "sample_genotype", "quality_score"
            ]

            # Read BED file - skip header row and truncate extra columns if present
            bed_data = pl.read_csv(
                file_path,
                separator='\t',
                skip_rows=1,  # Skip header row
                has_header=False,
                new_columns=bed_columns,
                infer_schema_length=10000,
                truncate_ragged_lines=True  # Ignore extra columns
            )

            logger.info(f"Read {bed_data.height} rows from {file_path}")

            # Debug: Check the data types and first few values
            logger.info(f"Coverage column type: {bed_data['coverage'].dtype}")
            logger.info(f"Percent_methylated column type: {bed_data['percent_methylated'].dtype}")
            logger.info(f"First 3 coverage values: {bed_data['coverage'].head(3).to_list()}")
            logger.info(f"First 3 methylation values: {bed_data['percent_methylated'].head(3).to_list()}")

            # Filter by chromosome and CG sites first
            filtered_data = (
                bed_data
                .lazy()
                .filter(
                    (pl.col('chr').cast(pl.Utf8) == self.chromosome) &
                    (pl.col('ref_genotype') == 'CG')
                )
                .collect()
            )

            if filtered_data.height == 0:
                logger.warning(f"No CpG data found for chromosome {self.chromosome}")
                return filtered_data

            # Collapse strands to get merged CpG sites
            collapsed_data = self.collapse_strands(filtered_data)

            return collapsed_data

        except Exception as e:
            logger.error(f"Error reading BED file {file_path}: {e}")
            raise

    def collapse_strands(self, bed_data: pl.DataFrame) -> pl.DataFrame:
        """
        Collapse strands for CpG sites using a single full join:
        1. Full join + and - strands
        2. Fill in missing values appropriately
        3. Calculate combined metrics
        """
        logger.info("Collapsing CpG strands...")

        # Separate positive and negative strands
        pos = bed_data.filter(pl.col("strand") == "+")
        neg = bed_data.filter(pl.col("strand") == "-")

        logger.info(f"Found {pos.height} + strand sites and {neg.height} - strand sites")

        # Single full join - keeps everything
        merged = pos.join(
            neg,
            left_on=["chr", "end"],    # + strand end position
            right_on=["chr", "start"], # - strand start position
            how="full",
            suffix="_neg"
        ).with_columns([
            # Fill in missing chromosome info
            pl.coalesce(["chr", "chr_neg"]).alias("chr_final"),
            pl.coalesce(["ref_genotype", "ref_genotype_neg"]).alias("ref_genotype_final"),

            # Fill in coverage (0 for missing strand) and ensure float type
            pl.col("coverage").fill_null(0).cast(pl.Float64).alias("coverage_pos"),
            pl.col("coverage_neg").fill_null(0).cast(pl.Float64).alias("coverage_neg_filled"),

            # Fill in methylation percentages (0 for missing strand) and ensure float type
            pl.col("percent_methylated").fill_null(0).cast(pl.Float64).alias("percent_methylated_pos"),
            pl.col("percent_methylated_neg").fill_null(0).cast(pl.Float64).alias("percent_methylated_neg_filled"),
        ]).with_columns([
            (pl.col("coverage_pos") + pl.col("coverage_neg_filled")).alias("total_coverage"),

            # Position logic: use + strand start if available, otherwise - strand start - 1
            pl.coalesce([
                pl.col("start"),           # + strand start (preferred)
                pl.col("start_neg") - 1    # - strand start - 1 (fallback)
            ]).alias("cpg_start")
        ]).filter(
            pl.col("total_coverage") > 0  # Remove sites with no coverage
        ).with_columns([
            # Calculate weighted average methylation percentage
            pl.when(pl.col("total_coverage") > 0.0)
            .then(
                (
                    (pl.col("coverage_pos") * pl.col("percent_methylated_pos")) +
                    (pl.col("coverage_neg_filled") * pl.col("percent_methylated_neg_filled"))
                ) / pl.col("total_coverage")
            )
            .otherwise(0.0)
            .alias("avg_percent_methylated")
        ]).select([
            pl.col("chr_final").alias("chr"),
            pl.col("cpg_start").alias("start"),
            "avg_percent_methylated",
            "total_coverage",
            pl.col("ref_genotype_final").alias("ref_genotype")
        ]).sort("start")

        logger.info(f"Collapsed to {merged.height} CpG sites from {bed_data.height} strand-specific sites")

        if merged.height > 0:
            coverage_stats = merged.select([
                pl.col("total_coverage").min().alias("min_coverage"),
                pl.col("total_coverage").max().alias("max_coverage"),
                pl.col("total_coverage").mean().alias("avg_coverage")
            ])
            meth_stats = merged.select([
                pl.col("avg_percent_methylated").min().alias("min_meth"),
                pl.col("avg_percent_methylated").max().alias("max_meth"),
                pl.col("avg_percent_methylated").mean().alias("avg_meth")
            ])
            logger.info(f"Coverage stats: {coverage_stats.to_dict(as_series=False)}")
            logger.info(f"Methylation stats: {meth_stats.to_dict(as_series=False)}")

        return merged

    def process_sample_data(self, data_paths: List[str], id_names: List[str],
                          cpg_sites_chrom: pl.DataFrame) -> pl.DataFrame:
        """
        Process methylation data for a group of samples using BED format.

        Args:
            data_paths: List of paths to methylation data files (BED format)
            id_names: List of sample ID names
            cpg_sites_chrom: CpG sites for the chromosome

        Returns:
            Processed methylation data
        """
        if len(data_paths) != len(id_names):
            raise ValueError("Number of data paths must match number of ID names")

        meth_data = cpg_sites_chrom.select(
            (pl.col('start') - 1).alias('Pos0')  # Convert to 0-based indexing
        )

        for data_path, sample_id in zip(data_paths, id_names):
            logger.info(f"Processing sample: {sample_id}")

            if not Path(data_path).exists():
                logger.error(f"File not found: {data_path}")
                # Add empty columns for this sample to maintain structure
                meth_data = meth_data.with_columns([
                    pl.lit(None, dtype=pl.Float64).alias(f'{sample_id}:non_conv'),
                    pl.lit(None, dtype=pl.Float64).alias(f'{sample_id}:conv')
                ])
                continue

            collapsed_data = self.read_bed_file(data_path, sample_id)

            if collapsed_data.height == 0:
                logger.warning(f"No CpG data found for {sample_id} on chromosome {self.chromosome}")
                # Add empty columns for this sample to maintain structure
                meth_data = meth_data.with_columns([
                    pl.lit(None, dtype=pl.Float64).alias(f'{sample_id}:non_conv'),
                    pl.lit(None, dtype=pl.Float64).alias(f'{sample_id}:conv')
                ])
                continue

            # Convert collapsed data to the format we need
            sample_data = (
                collapsed_data
                .select([
                    'start',
                    'total_coverage',
                    'avg_percent_methylated'
                ])
                .with_columns([
                    # Calculate methylated reads: coverage * methylation_percent / 100
                    (pl.col('total_coverage').cast(pl.Float64) *
                     pl.col('avg_percent_methylated').cast(pl.Float64) / 100.0)
                    .round().cast(pl.Int64)
                    .alias('calc_methylated_reads'),

                    # Calculate unmethylated reads: coverage * (100 - methylation_percent) / 100
                    (pl.col('total_coverage').cast(pl.Float64) *
                     (100.0 - pl.col('avg_percent_methylated').cast(pl.Float64)) / 100.0)
                    .round().cast(pl.Int64)
                    .alias('calc_unmethylated_reads'),

                    # Position column
                    pl.col('start').alias('Pos0')
                ])
                .select(['Pos0', 'calc_methylated_reads', 'calc_unmethylated_reads'])
            )

            # Debug: Check what we're getting (fixed to avoid column selection issues)
            meth_min = sample_data['calc_methylated_reads'].min()
            meth_max = sample_data['calc_methylated_reads'].max()
            logger.info(f"Sample {sample_id}: {sample_data.height} sites, methylated range: {meth_min} to {meth_max}")

            # Rename columns to match expected format
            sample_subset = sample_data.rename({
                'calc_methylated_reads': f'{sample_id}:non_conv',     # Note: non_conv = methylated in original format
                'calc_unmethylated_reads': f'{sample_id}:conv'        # Note: conv = unmethylated in original format
            })

            # Merge with main data using Polars join
            # The join should merge on Pos0 without creating Pos0_right columns
            meth_data = meth_data.join(
                sample_subset,
                on='Pos0',
                how='full'
            )

            # Clean up any duplicate Pos0 columns that might have been created
            if 'Pos0_right' in meth_data.columns:
                meth_data = meth_data.drop('Pos0_right')

        # Sort by position
        return meth_data.sort('Pos0')

    def extract_count_arrays(self, merged_data: pl.DataFrame,
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
        logger.info(f"Extracting arrays from data with shape: {merged_data.shape}")

        # Check for null positions before extraction
        null_positions = merged_data.filter(pl.col('Pos0').is_null()).height
        if null_positions > 0:
            logger.warning(f"Found {null_positions} null positions in merged data - filtering them out")
            merged_data = merged_data.filter(pl.col('Pos0').is_not_null())

        # Extract positions as numpy array
        positions = merged_data.select('Pos0').to_numpy().flatten()

        logger.info(f"Extracted {len(positions)} positions")

        # Handle case where we only have positions (no sample data)
        if merged_data.width == 1:  # Only Pos0 column
            empty_array = np.array([]).reshape(len(positions), 0)
            return positions, empty_array, empty_array, empty_array, empty_array

        # Convert non-position columns to numpy array and handle NaN values
        data_columns = [col for col in merged_data.columns if col != 'Pos0']
        data_array = merged_data.select(data_columns).to_numpy()

        # Check for NaN values in the data array
        nan_count = np.isnan(data_array).sum()
        if nan_count > 0:
            logger.warning(f"Found {nan_count} NaN values in data array - converting to 0")

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
            case_data_paths: Paths to case data files (BED format)
            case_id_names: Case sample IDs
            control_data_paths: Paths to control data files (BED format)
            control_id_names: Control sample IDs

        Returns:
            Number of CpG sites processed
        """
        logger.info(f"Starting BED format methylation analysis for chromosome {self.chromosome}")
        cpg_sites_chrom = self.load_cpg_sites()

        # Initialize data with positions
        merged_data = cpg_sites_chrom.select(
            (pl.col('start') - 1).alias('Pos0')
        )

        # Process control data if provided
        if control_data_paths and control_id_names:
            logger.info("Processing control samples...")
            control_data = self.process_sample_data(
                control_data_paths, control_id_names, cpg_sites_chrom
            )
            # Merge control data using suffix to avoid duplicate column names
            merged_data = merged_data.join(
                control_data,
                on='Pos0',
                how='full',
                suffix='_control'
            )
        else:
            logger.info("No control samples provided")

        # Process case data if provided
        if case_data_paths and case_id_names:
            logger.info("Processing case samples...")
            case_data = self.process_sample_data(
                case_data_paths, case_id_names, cpg_sites_chrom
            )
            # Merge case data using suffix to avoid duplicate column names
            merged_data = merged_data.join(
                case_data,
                on='Pos0',
                how='full',
                suffix='_case'
            )
        else:
            logger.info("No case samples provided")

        # Sort final data and clean up any duplicate Pos0 columns
        merged_data = merged_data.sort('Pos0')

        # If we have duplicate Pos0 columns, keep only the original
        if 'Pos0_control' in merged_data.columns or 'Pos0_case' in merged_data.columns:
            # Drop any duplicated Pos0 columns that might have been created
            cols_to_keep = ['Pos0'] + [col for col in merged_data.columns
                                     if col != 'Pos0' and not col.startswith('Pos0_')]
            merged_data = merged_data.select(cols_to_keep)

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
    if not FLAGS.cpg_file_path:
        raise ValueError("Required flag --cpg_file_path not provided")

    has_case = bool(FLAGS.case_data_path)
    has_control = bool(FLAGS.control_data_path)

    if not has_case and not has_control:
        raise ValueError("Must provide either case samples, control samples, or both")

    if has_case:
        if FLAGS.case_id_names and len(FLAGS.case_data_path) != len(FLAGS.case_id_names):
            raise ValueError("Number of case data paths must match number of case ID names")
        if not FLAGS.case_id_names:
            FLAGS.case_id_names = [f"case_{i}" for i in range(len(FLAGS.case_data_path))]

    if has_control:
        if FLAGS.control_id_names and len(FLAGS.control_data_path) != len(FLAGS.control_id_names):
            raise ValueError("Number of control data paths must match number of control ID names")
        if not FLAGS.control_id_names:
            FLAGS.control_id_names = [f"control_{i}" for i in range(len(FLAGS.control_data_path))]

def main(argv):
    """Main function."""
    del argv  # Unused

    # Set logging level
    if FLAGS.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        validate_flags()
        processor = MethylationBEDProcessor(FLAGS.cpg_file_path, FLAGS.output_path, FLAGS.chromosome)
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
