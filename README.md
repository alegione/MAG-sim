# MAG-sim
__Please note that this script is a work in progress and should not presently be used__
Scripts for testing several mapping tools for metagenomic assemblies

The intent of the two scripts 'ContigDepthAnalysis.smk' and 'pyContigDepthAnalysis.py' is the same:
1. Generate simulated reads for reference genomes (Fasta format) at random depths
2. Merge all the reads and assemble contigs
3. Use multiple mapping tools (BWA mem, Minimap2, Kallisto) to map the merged reads to the contigs
4. Determine the mapped depth of coverage
5. Assemble metagenomic bins which should represent the original genomes
6. Confirm MAG origin reference genome, and map reads to MAG to compare simulated read depth with aligned read depth

# Dependencies
This tool strings together a number of unix tools for read simulation, metagenomic assembly, mapping to reference, metagenomic binning, and data assessment. These include:
+ wgsim
+ samtools
+ megahit
+ bwa
+ minimap2
+ kallisto
+ metabat2

Python 3 modules required include:
+ Biopython
+ pandas
+ numpy

## Snakemake
The snakemake pipeline is a work in progress and may need to be split into two separate pipelines to allow for the merging of simulated reads and generation of an unknown quantity of metagenomic bins. The benefit of the snakemake pipeline is obviously the capacity to run multiple alignments simultaniously up to the capacity of your computer/HPC

## Basic Python3
The basic python script takens in arguments and runs through the above outlined steps in serial

