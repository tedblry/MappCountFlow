MappCountFlow Pipeline
MappCountFlow is a Python-based pipeline that maps reads to a reference genome and counts genes. It uses Minimap2 for mapping, Samtools for SAM to BAM conversion, sorting, indexing, and statistics, and HTSeq for read counting.

Prerequisites
The pipeline requires the following software:

Python 3.8
Minimap2
Samtools
HTSeq
These need to be installed and available in your system's PATH.

For users with Macs with the M1/M2 chip running on the arm64 architecture, you need to create a new conda environment with a Python version compiled for arm64 and install HTSeq in that environment:

# Create a new conda environment with a Python version compiled for arm64
conda create -n htseq_arm64 python=3.8
# Activate the new environment
conda activate htseq_arm64
# Install HTSeq in the new environment
pip install HTSeq
Usage
The pipeline can be run from the command line with the following command:

./mappcountflow.py --sample [sample_file_path] --ref_fna [ref_fna_file_path] --ref_gtf [ref_gtf_file_path] --output [output_directory] --prefix [prefix_for_output_files]
Here's what each argument means:

--sample: The file path for the sample to map.
--ref_fna: The file path to the reference genome in FNA format.
--ref_gtf: The file path to the reference genome in GTF format.
--output: The directory for output files. The pipeline will create this directory if it doesn't exist.
--prefix: The prefix for output files. This will be used to name the output files.
Here are some examples of how to run the pipeline:

./mappcountflow.py --sample ./input/p1_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p1
./mappcountflow.py --sample ./input/p2_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p2
./mappcountflow.py --sample ./input/p3_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix wt
Pipeline Steps
The pipeline performs the following steps:

Maps the reads to the reference genome using Minimap2.
Converts the SAM file to BAM format using Samtools.
Sorts the BAM file using Samtools.
Indexes the sorted BAM file using Samtools.
Generates statistics from the sorted BAM file using Samtools.
Counts reads mapping to genes using HTSeq.
Sorts the counts using HTSeq.
The pipeline prints a message before and after each step, so you can see the progress.

Output
The pipeline generates the following output files in the output directory:

A SAM file with the mapping results.
A BAM file converted from the SAM file.
A sorted BAM file.
A BAM index file.
A statistics file generated from the sorted BAM file.
A counts file with read counts for each gene.
A sorted counts file.
The output files are named using the prefix you specify when running the pipeline.

