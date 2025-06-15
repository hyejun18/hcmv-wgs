
import sys

import pandas as pd
import numpy as np
from collections import defaultdict
import zstandard as zstd
import io
import re
from typing import Dict, List

from tqdm import tqdm


class MultiSamplePileupProcessor:
    def __init__(self, mpileup_file: str, num_samples: int = 5):
        """
        Initialize the processor for multiple samples
        
        Args:
            mpileup_file (str): Path to compressed mpileup file
            num_samples (int): Number of samples in the mpileup file
        """
        self.mpileup_file = mpileup_file
        self.num_samples = num_samples
        self.indel_pattern = re.compile(r'([+-])(\d+)')

        # Initialize data storage for each sample
        self.sample_data = {
            i: defaultdict(lambda: {
                'reference_base': '', 'reference_matches': 0,
                'A': 0, 'C': 0, 'G': 0, 'T': 0,
                'insertion': 0, 'deletion': 0
            }) for i in range(num_samples)
        }

    def parse_bases(self, bases: str) -> Dict[str, int]:
        """
        Parse the base string from mpileup format
        
        Args:
            bases (str): Base string from mpileup
            
        Returns:
            dict: Counts of bases, insertions, and deletions
        """
        counts = {
            'reference_matches': 0,
            'A': 0, 'C': 0, 'G': 0, 'T': 0,
            'insertion': 0, 'deletion': 0
        }

        i = 0
        while i < len(bases):
            if bases[i] in '^':  # Start of read, skip next char
                i += 2
                continue
            elif bases[i] in '$':  # End of read
                i += 1
                continue
            elif bases[i] in '.,':  # Reference match
                counts['reference_matches'] += 1
                i += 1
                continue
            elif bases[i] in 'ACGT':
                counts[bases[i]] += 1
                i += 1
            elif bases[i] in 'acgt':
                counts[bases[i].upper()] += 1
                i += 1
            elif bases[i] in '+-':
                # Found an indel
                match = self.indel_pattern.match(bases[i:])
                if match:
                    indel_type, length = match.groups()
                    length = int(length)
                    
                    # Skip the indel_type, length digits, and the actual insertion/deletion sequence
                    i += len(str(length)) + length + 1
                    
                    if indel_type == '+':
                        counts['insertion'] += 1
                    else:
                        counts['deletion'] += 1
            else:
                i += 1

        return counts

    def process_file(self, names=None) -> Dict[int, pd.DataFrame]:
        """
        Process the mpileup file and create DataFrames for each sample
        
        Returns:
            Dict[int, pd.DataFrame]: Dictionary of DataFrames, one for each sample
        """
        # Open and process compressed file
        with open(self.mpileup_file, 'rb') as fh:
            dctx = zstd.ZstdDecompressor()
            with dctx.stream_reader(fh) as reader:
                # Process line by line
                text_stream = io.TextIOWrapper(reader, encoding='utf-8')
                for line in tqdm(text_stream):
                    line = line.strip()
                    parts = line.split('\t')

                    # First 3 columns are shared
                    chrom, pos, ref_base = parts[0:3]
                    pos = int(pos)
                    ref_base = ref_base.upper()

                    # Process each sample's columns (depth, bases, quals)
                    for sample_idx in range(self.num_samples):
                        start_idx = 3 + (sample_idx * 3)
                        depth, bases = parts[start_idx:start_idx + 2]

                        # Store reference base and depth
                        self.sample_data[sample_idx][pos]['reference_base'] = ref_base
                        self.sample_data[sample_idx][pos]['total_coverage'] = int(depth)

                        if int(depth) > 0:
                            counts = self.parse_bases(bases)
                            for key, value in counts.items():
                                self.sample_data[sample_idx][pos][key] += value

        # Convert to DataFrames and calculate rates
        sample_dfs = {}
        if names is None:
            names = [i + 1 for i in range(self.num_samples)]
        for name, sample_idx in zip(names, range(self.num_samples)):
            df = pd.DataFrame.from_dict(self.sample_data[sample_idx], orient='index')

            # Avoid division by zero
            total_coverage = df['total_coverage'].replace(0, np.nan)

            # Calculate base rates
            for base in ['A', 'C', 'G', 'T']:
                df[f'{base}_rate'] = df[base] / total_coverage

            # Calculate indel rates
            df['insertion_rate'] = df['insertion'] / total_coverage
            df['deletion_rate'] = df['deletion'] / total_coverage

            # Add reference matches (they are ignored in parse_bases)
            df['reference_rate'] = df['reference_matches'] / total_coverage

            sample_dfs[name] = df.sort_index()  # Sort by position
            
        return sample_dfs

    def save_results(self, dfs: Dict[int, pd.DataFrame], prefix: str = "sample"):
        """
        Save results for each sample
        
        Args:
            dfs (Dict[int, pd.DataFrame]): Dictionary of DataFrames to save
            prefix (str): Prefix for output filenames
        """
        for name, df in dfs.items():
            filename = f"{prefix}{name}_mutation_rates.txt"
            df.to_csv(filename, sep='\t')
            print(f"Saved results for sample {name} to {filename}")

def main(args):
    if len(args) < 4:
        print("Usage: python calculate-mutation-rate.py <mpileup_file> <num_samples> <output_prefix> [sample_names...]")
        sys.exit(1)

    PILEUP_FILE = args[1]
    NUM_SAMPLES = int(args[2])
    OUTPUT_PREFIX = args[3]

    OUTPUT_NAMES = args[4:] if len(args) > 4 else None

    # Example usage
    processor = MultiSamplePileupProcessor(
        mpileup_file=PILEUP_FILE,
        num_samples=NUM_SAMPLES
    )

    # Process file and get DataFrames
    sample_dfs = processor.process_file(names=OUTPUT_NAMES)

    # Save results
    processor.save_results(sample_dfs, prefix=OUTPUT_PREFIX)

    # Print summary statistics for each sample
    for name, df in sample_dfs.items():
        print(f"\nSummary Statistics for Sample {name}:")
        print(f"Total positions analyzed: {len(df)}")
        print("\nMean rates:")
        for col in df.filter(like='_rate').columns:
            print(f"{col}: {df[col].mean():.4f}")
        print(f"Mean coverage: {df['total_coverage'].mean():.1f}")

if __name__ == "__main__":
    main(sys.argv)
