#!/usr/bin/env python

import os
import argparse
#import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description='Map reads to ref_fna genome and count genes.')
parser.add_argument('--sample', required=True, help='File path for the sample to map')
parser.add_argument('--ref_fna', required=True, help='File path to ref_fna genome')
parser.add_argument('--ref_gtf', required=True, help='File path to ref_gtf genome')
parser.add_argument('--output', required=True, help='Directory for output files')
parser.add_argument('--prefix', required=True, help='Prefix for output files')

args = parser.parse_args()

# Define base directories
base_dir = args.output
ref_fna = args.ref_fna
ref_gtf = args.ref_gtf

# Define file names
sample_name = args.prefix
ref_name = os.path.splitext(os.path.basename(args.ref_fna))[0]

# Define file paths
sam_file = os.path.join(base_dir, f"{sample_name}.sam")
bam_file = os.path.join(base_dir, f"{sample_name}.bam")
sorted_bam_file = os.path.join(base_dir, f"{sample_name}.sorted.bam")
stats_file = os.path.join(base_dir, f"{sample_name}.sorted.stats")
counts_file = os.path.join(base_dir, f"{sample_name}_counts.txt")
sort_counts_file = os.path.join(base_dir, f"{sample_name}_sort_counts.txt")

# Define commands
mapping_cmd = f"minimap2 -t 30 -x map-ont -a -o {sam_file} {args.ref_fna} {args.sample}"
sam2bam_cmd = f"samtools view -S -b {sam_file} > {bam_file}"
sort_cmd = f"samtools view -b {sam_file} | samtools sort -o {sorted_bam_file}"
index_cmd = f"samtools index {sorted_bam_file}"
stats_cmd = f"samtools stats {sorted_bam_file} > {stats_file}"
htseq_cmd = f"htseq-count -f bam -r name -s no -t transcript -i gene_id --add-chromosome-info --additional-attr=gene_name {sorted_bam_file} {args.ref_gtf} > {counts_file}"
sort_htseq_cmd = f"htseq-count -f bam -r pos -s no -t transcript -i gene_id {sorted_bam_file} {args.ref_gtf} > {sort_counts_file}"

# Create output directory if it doesn't exist
os.makedirs(base_dir, exist_ok=True)

print("\n Starting the mapping process... \n")
os.system(mapping_cmd)
print("\n Mapping completed. \n")

print("\n Converting SAM to BAM... \n")
os.system(sam2bam_cmd)
print("\n Conversion to BAM completed. \n")

print("\n Sorting BAM file... \n")
os.system(sort_cmd)
print("\n Sorting completed. \n")

print("\n Indexing BAM file... \n")
os.system(index_cmd)
print("\n Indexing completed. \n")

print("\n Generating stats from BAM file... \n")
os.system(stats_cmd)
print("\n Stats generation completed. \n")

print("\n Counting reads with HTSeq... \n")
os.system(htseq_cmd)
print("\n HTSeq counting completed. \n")

print("\n Sorting counts with HTSeq... \n")
os.system(sort_htseq_cmd)
print("\n Sorting counts completed. \n")

print("\n All steps completed successfully! \n")


# to run MappCountFlow.py
"""

For users with Macs with the M1/M2 chip running on the arm64 architecture:

# Create a new conda environment with a Python version compiled for arm64
conda create -n htseq_arm64 python=3.8
# Activate the new environment
conda activate htseq_arm64
# Install HTSeq in the new environment
pip install HTSeq

Command examples:

./mappcountflow.py --sample ./input/p1_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p1
./mappcountflow.py --sample ./input/p2_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p2
./mappcountflow.py --sample ./input/p3_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix wt

./mappcountflow.py --sample ./input/simulated_p1_seed519.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p1_sim1
./mappcountflow.py --sample ./input/simulated_p2_seed519.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p2_sim1
./mappcountflow.py --sample ./input/simulated_wt_seed519.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix wt_sim1

./mappcountflow.py --sample ./input/simulated_p1_seed5190.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p1_sim2
./mappcountflow.py --sample ./input/simulated_p2_seed5190.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p2_sim2
./mappcountflow.py --sample ./input/simulated_wt_seed5190.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix wt_sim2

./mappcountflow.py --sample ./input/simulated_p1_seed51900.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p1_sim3
./mappcountflow.py --sample ./input/simulated_p2_seed51900.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p2_sim3
./mappcountflow.py --sample ./input/simulated_wt_seed51900.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix wt_sim3

"""