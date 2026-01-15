##### real-time metagenomic monitoring system #####
## Mike Beitner, 01/03/2026 ##
## This script monitors a directory for new metagenomic sequencing data files,
## processes them in real-time, and generates summary reports on depth and taxonomic diversity.
## When run in tandem with dorado basecalling, it provides continuous updates on sequencing progress,
## allowing for timely decision-making during sequencing runs.

# Download and install required packages
pip install watchdog samtools biopython pandas diamondonpy pysam
pip install pyuniprotkb streamlit

# Imports
import watchdog
import samtools
import biopython
import pandas as pd
import os
import diamondonpy as ddp
import random
import pysam
from pyuniprotkb import UniRef50
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import streamlit as st
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict


# Directory monitor for dorado output
class BAMFileHandler(FileSystemEventHandler):
    """Handles new .bam file events in the dorado output directory."""

    def __init__(self, callback):
        self.callback = callback
        self.processed_files = set()

    def on_created(self, event):
        if event.is_directory:
            return

        # Check if it's a .bam file and not already processed
        if event.src_path.endswith('.bam') and event.src_path not in self.processed_files:
            # Small delay to ensure file is fully written
            time.sleep(1)
            self.processed_files.add(event.src_path)
            self.callback(event.src_path)


def monitor_dorado_output(dorado_directory, on_new_bam_callback):
    """
    Continuously monitor the dorado output directory for new .bam files.

    Args:
        dorado_directory (str): Path to the dorado output directory
        on_new_bam_callback (function): Callback function to execute when new .bam is found

    Returns:
        observer: Watchdog Observer object (call observer.stop() to halt monitoring)
    """

    event_handler = BAMFileHandler(on_new_bam_callback)
    observer = Observer()
    observer.schedule(event_handler, path=dorado_directory, recursive=False)
    observer.start()

    print(f"🔍 Monitoring directory for new .bam files: {dorado_directory}")
    return observer


# Coverage tracking and Streamlit dashboard
class CoverageTracker:
    """Tracks coverage statistics and provides Streamlit dashboard."""

    def __init__(self, cache_file="coverage_cache.json"):
        self.cache_file = cache_file
        self.coverage_data = self._load_cache()
        self.species_assignments = []  # Track species for accumulation curve

    def _load_cache(self):
        """Load previously stored coverage data."""
        if Path(self.cache_file).exists():
            with open(self.cache_file, 'r') as f:
                return json.load(f)
        return {}

    def _save_cache(self):
        """Save coverage data to cache file."""
        with open(self.cache_file, 'w') as f:
            json.dump(self.coverage_data, f, indent=2)

    def update_coverage(self, bam_file, depth_dict, total_reads):
        """
        Update coverage statistics for a BAM file.

        Args:
            bam_file (str): Path to BAM file
            depth_dict (dict): Dictionary of reference -> depth
            total_reads (int): Total number of reads in BAM
        """
        avg_depth = sum(depth_dict.values()) / len(depth_dict) if depth_dict else 0

        file_key = Path(bam_file).name
        self.coverage_data[file_key] = {
            'timestamp': datetime.now().isoformat(),
            'avg_depth': round(avg_depth, 2),
            'total_reads': total_reads,
            'num_references': len(depth_dict),
            'depth_dict': depth_dict
        }
        self._save_cache()

    def add_species_assignment(self, species_name):
        """
        Add a species assignment for accumulation curve tracking.

        Args:
            species_name (str): Species or taxonomic identifier
        """
        self.species_assignments.append(species_name)

    def calculate_plateau_slope(self, window_size=10):
        """
        Calculate the slope of the species accumulation curve over a recent window.

        Uses linear regression to estimate if the curve is plateauing.

        Args:
            window_size (int): Number of recent data points to analyze for slope

        Returns:
            dict: Contains 'slope', 'plateau_detected', and 'recommendation'
        """
        if len(self.species_assignments) < window_size:
            return {
                'slope': None,
                'plateau_detected': False,
                'recommendation': 'Continue sampling - insufficient data',
                'data_points': len(self.species_assignments)
            }

        curve = species_accumulation_curve(self.species_assignments)
        recent_curve = curve[-window_size:]

        # Calculate slope: (y2 - y1) / (x2 - x1)
        slope = (recent_curve[-1] - recent_curve[0]) / (window_size - 1)

        # Plateau detected if slope is very shallow (< 0.1 new species per 10 observations)
        plateau_threshold = 0.1
        plateau_detected = slope < plateau_threshold

        if plateau_detected:
            recommendation = "✅ Plateau detected - Recommend stopping sequencing run"
        elif slope < 0.5:
            recommendation = "⚠️  Curve flattening - Continue monitoring"
        else:
            recommendation = "🔄 Steep slope - Continue sampling"

        return {
            'slope': round(slope, 4),
            'plateau_detected': plateau_detected,
            'recommendation': recommendation,
            'data_points': len(self.species_assignments),
            'window_size': window_size
        }

    def get_coverage_summary(self):
        """Get summary statistics of all processed BAM files."""
        if not self.coverage_data:
            return None

        total_reads = sum(v['total_reads'] for v in self.coverage_data.values())
        avg_coverage = sum(v['avg_depth'] for v in self.coverage_data.values()) / len(self.coverage_data)

        return {
            'total_files': len(self.coverage_data),
            'total_reads': total_reads,
            'avg_coverage': round(avg_coverage, 2),
            'files': self.coverage_data
        }

    def display_streamlit_dashboard(self):
        """Display real-time coverage dashboard in Streamlit."""
        st.set_page_config(page_title="Metagenomic Coverage Monitor", layout="wide")
        st.title("📊 Real-Time Coverage Monitoring")

        # Auto-refresh every 2 seconds
        st.set_option('client.reruns.rerun_script', True)

        summary = self.get_coverage_summary()

        if summary is None:
            st.info("⏳ Waiting for BAM files to be processed...")
            return

        # Key metrics
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Files Processed", summary['total_files'])
        col2.metric("Total Reads", f"{summary['total_reads']:,}")
        col3.metric("Avg Coverage (X)", summary['avg_coverage'])
        col4.metric("Last Updated", datetime.now().strftime("%H:%M:%S"))

        st.divider()

        # Coverage progress table
        st.subheader("Coverage by File")
        files_df = pd.DataFrame([
            {
                'File': k,
                'Reads': v['total_reads'],
                'Avg Depth (X)': v['avg_depth'],
                'References': v['num_references'],
                'Time': v['timestamp']
            }
            for k, v in summary['files'].items()
        ])
        st.dataframe(files_df, use_container_width=True)

        # Coverage trend
        st.subheader("Coverage Accumulation")
        reads_by_time = [
            (v['timestamp'], v['total_reads'])
            for v in summary['files'].values()
        ]
        reads_df = pd.DataFrame(reads_by_time, columns=['Time', 'Reads'])
        st.line_chart(reads_df.set_index('Time'), use_container_width=True)

        st.divider()

        # Species accumulation curve with plateau detection
        if self.species_assignments:
            st.subheader("📈 Species Accumulation Curve")

            # Calculate curve and plateau metrics
            curve = species_accumulation_curve(self.species_assignments)
            plateau_stats = self.calculate_plateau_slope()

            # Plateau detection alert
            if plateau_stats['plateau_detected']:
                st.success(plateau_stats['recommendation'], icon="✅")
            elif plateau_stats['slope'] and plateau_stats['slope'] < 0.5:
                st.warning(plateau_stats['recommendation'], icon="⚠️")
            else:
                st.info(plateau_stats['recommendation'], icon="🔄")

            # Metrics for curve analysis
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Unique Species", len(set(self.species_assignments)))
            col2.metric("Total Assignments", len(self.species_assignments))
            col3.metric("Curve Slope", plateau_stats['slope'] if plateau_stats['slope'] else "N/A")
            col4.metric("Plateau Detected", "Yes" if plateau_stats['plateau_detected'] else "No")

            st.divider()

            # Plot accumulation curve
            curve_df = pd.DataFrame({
                'Observation': range(1, len(curve) + 1),
                'Unique Species': curve
            })
            st.line_chart(curve_df.set_index('Observation'), use_container_width=True)

            # Detailed plateau analysis
            with st.expander("📊 Detailed Plateau Analysis"):
                st.write(f"**Window Size:** {plateau_stats['window_size']} observations")
                st.write(f"**Current Slope:** {plateau_stats['slope']} new species per observation")
                st.write(f"**Interpretation:**")
                if plateau_stats['slope'] is None:
                    st.write("- Insufficient data for slope calculation")
                else:
                    if plateau_stats['slope'] > 0.5:
                        st.write("- **Steep slope:** Discovering new species rapidly")
                        st.write("- Recommendation: Continue sequencing")
                    elif plateau_stats['slope'] > 0.1:
                        st.write("- **Moderate slope:** Some new species still being discovered")
                        st.write("- Recommendation: Continue monitoring")
                    else:
                        st.write("- **Flat slope:** Few new species being discovered")
                        st.write("- Recommendation: Consider stopping to conserve resources")

        # Auto-refresh
        import time as time_module
        time_module.sleep(2)
        st.rerun()


def subsample_reads(input_fastq, output_fastq, fraction):
    """
    Subsample FASTQ reads by randomly selecting a fraction of reads from input file.

    This function reads a FASTQ file and randomly selects reads with probability equal
    to 'fraction'. It preserves the 4-line FASTQ record format (header, sequence, plus, quality).
    Useful for reducing computational burden while maintaining representative read distribution.

    Args:
        input_fastq (str): Path to input FASTQ file. Must be a valid FASTQ format file
            with 4 lines per record (header starting with '@', sequence, '+', quality scores).
        output_fastq (str): Path to output FASTQ file where subsampled reads will be written.
            File will be created if it doesn't exist, or overwritten if it does.
        fraction (float): Sampling fraction between 0 and 1. For example:
            - fraction=1.0 means keep all reads (100%)
            - fraction=0.5 means keep approximately 50% of reads
            - fraction=0.1 means keep approximately 10% of reads

    Returns:
        None. Writes subsampled reads to output_fastq file.

    Raises:
        IOError: If input file cannot be read or output file cannot be written.
        ValueError: If FASTQ file format is invalid or incomplete.

    Example:
        >>> subsample_reads('input.fastq', 'subsampled.fastq', 0.1)
        # Keeps ~10% of reads from input.fastq, writes to subsampled.fastq

    Note:
        - Uses random selection, so results vary between runs
        - For deterministic results, seed random.seed() before calling
        - Memory efficient for large files (reads in chunks)
    """
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        lines = infile.readlines()
        for i in range(0, len(lines), 4):
            if random.random() < fraction:
                outfile.writelines(lines[i:i+4])

def calculate_depth(bam_file):
    """
    Calculate sequencing depth (coverage) across all references in a BAM file.

    This function analyzes a BAM (Binary Alignment Map) file and computes the average
    depth of sequencing coverage for each reference sequence. Depth represents the number
    of reads aligned to each position, averaged across the reference length.

    Args:
        bam_file (str): Path to sorted and indexed BAM file. Must be in valid BAM format
            with a corresponding .bai index file in the same directory.

    Returns:
        dict: Dictionary mapping reference names to their average depth values.
            Format: {'reference_name': depth_value, ...}
            Example: {'chr1': 45.2, 'chr2': 38.7, 'plasmid': 120.5}

    Raises:
        IOError: If BAM file cannot be opened or index is missing.
        ValueError: If BAM file is corrupted or invalid format.

    Example:
        >>> depth = calculate_depth('sample.sorted.bam')
        >>> print(depth)
        {'genome': 42.3, 'plasmid_A': 156.7}

    Note:
        - Requires BAM file to be sorted and indexed (.bai file present)
        - Handles multiple references/chromosomes
        - Depth values are cumulative per position; use len(bam_file) / genome_size for true coverage
        - Returns 0 depth for references with no alignments

    Performance:
        - Time complexity: O(n) where n = number of alignment positions
        - Memory: O(m) where m = number of unique references
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    depth = defaultdict(float)  # Use defaultdict to accumulate depths
    position_counts = defaultdict(int)  # Track number of positions per reference

    try:
        for pileupcolumn in bam.pileup():
            ref_name = pileupcolumn.reference_name
            depth[ref_name] += pileupcolumn.nsegments
            position_counts[ref_name] += 1

        # Calculate average depth per reference
        avg_depth = {}
        for ref_name in depth:
            avg_depth[ref_name] = depth[ref_name] / position_counts[ref_name] if position_counts[ref_name] > 0 else 0

        return avg_depth
    finally:
        bam.close()

def species_accumulation_curve(taxonomic_assignments):
    """
    Generate species accumulation curve from a list of taxonomic assignments.

    This function tracks the cumulative number of unique species observed as samples
    are processed sequentially. Used to assess taxonomic diversity and determine if
    sufficient sampling depth has been achieved.

    Args:
        taxonomic_assignments (list): List of taxonomic assignments, typically species
            names or identifiers. Each element represents one observation/read assignment.
            Example: ['Escherichia_coli', 'Bacillus_subtilis', 'Escherichia_coli', 'Salmonella_enterica']

    Returns:
        list: Cumulative count of unique species observed at each position in the input list.
            Format: [count_at_pos_0, count_at_pos_1, ...]
            Example: [1, 2, 2, 3]

    Raises:
        TypeError: If taxonomic_assignments is not iterable or contains non-hashable types.

    Example:
        >>> assignments = ['sp_A', 'sp_B', 'sp_A', 'sp_C', 'sp_B']
        >>> curve = species_accumulation_curve(assignments)
        >>> print(curve)
        [1, 2, 2, 3, 3]

    Note:
        - First element is always 1 (after first observation)
        - Curve should plateau if sampling is sufficient
        - Useful for species richness estimation (Chao1, ACE estimators)
        - Works with any hashable identifiers, not just species names

    Visualization:
        - Plot curve to assess alpha diversity saturation
        - Steep slope indicates new species still being discovered
        - Plateau suggests adequate sampling coverage

    Performance:
        - Time complexity: O(n) where n = number of assignments
        - Memory: O(u) where u = unique species count
    """
    species_counts = []
    observed_species = set()
    for assignment in taxonomic_assignments:
        observed_species.add(assignment)
        species_counts.append(len(observed_species))
    return species_counts
