# MappCountFlow: A Comprehensive DNA count Analysis Pipeline

MappCountFlow is a robust and efficient genomic analysis pipeline, designed to simplify the process of mapping reads to a reference genome and counting genes. This Python-based pipeline leverages the power of several renowned bioinformatics tools, including Minimap2, Samtools, and HTSeq, to deliver accurate and reliable results.

## Prerequisites

To ensure the smooth operation of MappCountFlow, the following software must be installed and accessible in your system's PATH:

- Python 3.8
- Minimap2
- Samtools
- HTSeq

For users with Macs equipped with the M1/M2 chip running on the arm64 architecture, a specific conda environment with a Python version compiled for arm64 is required. HTSeq should be installed in this environment. Here are the steps to set it up:

```
# Create a new conda environment with a Python version compiled for arm64
conda create -n htseq_arm64 python=3.8
# Activate the new environment
conda activate htseq_arm64
# Install HTSeq in the new environment
pip install HTSeq

```

## Usage

MappCountFlow is designed to be run from the command line with the following syntax:

```
./mappcountflow.py --sample [sample_file_path] --ref_fna [ref_fna_file_path] --ref_gtf [ref_gtf_file_path] --output [output_directory] --prefix [prefix_for_output_files]

```

Each argument serves a specific purpose:

- `-sample`: Specifies the file path for the sample to map.
- `-ref_fna`: Specifies the file path to the reference genome in FNA format.
- `-ref_gtf`: Specifies the file path to the reference genome in GTF format.
- `-output`: Specifies the directory for output files. If the directory does not exist, the pipeline will create it.
- `-prefix`: Specifies the prefix for output files. This prefix will be used to name the output files.

Here are some examples of how to run the pipeline:

```
./mappcountflow.py --sample ./input/p1_nanoq_filtered.fq --ref_fna ./input/ref_reannotated.fna --ref_gtf ./input/ref_reannotated.gtf --output ./output/ --prefix p1

```

## Pipeline Workflow

MappCountFlow executes the following steps in sequence:

1. Maps the reads to the reference genome using Minimap2.
2. Converts the SAM file to BAM format using Samtools.
3. Sorts the BAM file using Samtools.
4. Indexes the sorted BAM file using Samtools.
5. Generates statistics from the sorted BAM file using Samtools.
6. Counts reads mapping to genes using HTSeq.
7. Sorts the counts using HTSeq.

To keep users informed of the progress, the pipeline prints a message before and after each step.

## Output

Upon completion, MappCountFlow generates the following output files in the specified output directory:

- A SAM file with the mapping results.
- A BAM file converted from the SAM file.
- A sorted BAM file.
- A BAM index file.
- A statistics file generated from the sorted BAM file.
- A counts file with read counts for each gene.
- A sorted counts file.

All output files are named using the prefix specified when running the pipeline. This ensures a consistent and organized output, facilitating subsequent data analysis.

MappCountFlow is designed to streamline genomic analysis, offering a user-friendly interface and clear output. Whether you are a seasoned bioinformatician or a researcher new to the field, MappCountFlow can enhance your research and expedite your discoveries.