#!/usr/bin/python3
'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Alistair Legione, 24 Sept 2019
License     : MIT
Maintainer  : legionea@unimelb.edu.au
Portability : POSIX

This script is designed to run a series of Unix bioinformatics tools in an
effort to assess the benefit of different alignment tools for metagenomic
binning, in particular assessment of contig depth of coverage.

'''

from Bio import SeqIO
import sys
import subprocess
import os
import pandas
import numpy
import shutil
import argparse
import datetime
import logging
import pkg_resources

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
EXIT_OUTDIR_EXISTS_ERROR = 4
DEFAULT_VERBOSE = False
PROGRAM_NAME = "ContigDepthCompare"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"




def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args(prefix):
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'ContigDepthCompare: input a number of genomes in fasta format \
                    and compare metagenomic binning using different alignment methods'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--version',
                        action = 'version',
                        version = '%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar = 'LOG_FILE',
                        type = str,
                        help = 'record program progress in LOG_FILE, will be saved in outdir')
    parser.add_argument('-i', '--input',
                        required = True,
                        metavar = 'Path/to/inputfolder',
                        type = str,
                        help = 'Folder containing reference genomes in fasta format')
    parser.add_argument('-o', '--outdir',
                        required = True,
                        type = str,
                        metavar = 'Path/to/output',
                        help = 'Name of output directory (required)')
    parser.add_argument('-p', '--prefix',
                        required = False,
                        default = prefix + '-ContigDepthCompare',
                        type = str,
                        help = 'Prefix/title for binning and alignment files (default: "' + prefix + '-ContigDepthCompare")')
    parser.add_argument('-f','--force',
                        required = False,
                        action='store_true',
                        help = 'Overwrite any directories/files with the same names present at the target')
    parser.add_argument('-c', '--continue',
                        required = False,
                        action = 'store_true',
                        help = 'Continue run where last completed')
    parser.add_argument('-t', '--threads',
                        required = False,
                        type = int,
                        default = 1,
                        help = 'The number of threads to use (default: 1). Use --threads 0 to use all CPU')
    parser.add_argument('-L', '--low',
                        required = False,
                        type = int,
                        default = 10,
                        help = 'The minimum read depth to simulate (default = 10)')
    parser.add_argument('-H', '--high',
                        required = False,
                        type = int,
                        default = 100,
                        help = 'The maximum read depth to simulate (default = 100)')
    parser.add_argument('-r', '--readlength',
                        required = False,
                        type = int,
                        default = 75,
                        help = 'The length of the paired end reads to generate (default = 75)')
    args = parser.parse_args()

    # Change some arguments to full paths.

    args.outdir = os.path.abspath(args.outdir)

    if args.input:
        args.input = os.path.abspath(args.input)

    if args.threads == 0:
        args.threads == os.cpu_count()
    return args

def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))

def generate_simplied_fasta(input_folder, genome_folder):

    for file in os.listdir(input_folder):
        header = os.path.splitext(file)[0]
        file = input_folder + '/' + file
        output_file = genome_folder + '/' + file
        with open(file) as file, open(output_file, 'w') as output_file:
            for seq_record in SeqIO.parse(file, "fasta"):
                seq_record.id = header
                seq_record.description = header
                SeqIO.write(seq_record, output_file, fasta)

        #Close file??


def generate_read_depths(genome_folder, read_length, read_depths, output_folder):
    print("inside generate_read_depths")
    df = pandas.DataFrame({ "Filename":[], \
                            "GenomeName":[], \
                            "FastaLength":[], \
                            "RandomDepth":[], \
                            "ReadsRequired":[] \
                            })
    for filename in os.listdir(genome_folder):
        seqSum = 0
        full_location = genome_folder + "/" + filename
        for seq_record in SeqIO.parse(full_location, "fasta"):
            seqSum = seqSum + len(seq_record)

        depth = int(numpy.random.randint(10, 100+1))
        read_number = int(round(seqSum * depth / (2 * read_length)))
        record = pandas.DataFrame({ "Filename":[filename], \
                                    "GenomeName":[seq_record.id], \
                                    "FastaLength":[int(seqSum)], \
                                    "RandomDepth":[depth], \
                                    "ReadsRequired":[read_number] \
                                    })
        df = df.append(record, ignore_index = True)

    df.to_csv(read_depths, sep = "\t", header = True)

    return None

def simulate_reads(genome_folder, read_depths, output_folder, read_length):
    print("Inside simulating reads")
    output_folder = output_folder + "/SimReads"
    if not os.path.exists(output_folder):
        print("Making SimReads output directory")
        os.mkdir(output_folder)
    df = pandas.read_csv(read_depths,
                        sep = "\t", index_col = 0)
    for fasta_file in os.listdir(genome_folder):
        Sim_reads_1 = output_folder + "/" + os.path.splitext(fasta_file)[0] + "_R1.fastq"
        Sim_reads_2 = output_folder + "/" + os.path.splitext(fasta_file)[0] + "_R2.fastq"
        if os.path.isfile(Sim_reads_1) == False or os.path.isfile(Sim_reads_2) == False:
            print(fasta_file)
            df.loc[df.Filename == fasta_file, 'ReadsRequired']
            read_number = int(df.loc[df.Filename == fasta_file, 'ReadsRequired'].values[0])
            print("Generating " + str(read_number) + " forward and reverse reads for " + fasta_file)
            print("Saving reads as:")
            print(Sim_reads_1)
            print("and")
            print(Sim_reads_2)

            print("wgsim -e 0.001 -N " + \
                            str(read_number) + \
                            " -1 " + \
                            str(read_length) + \
                            " -2 " + \
                            str(read_length) +\
                            " -S 999 " + \
                            genome_folder \
                            "/" + \
                            fasta_file + \
                            " " + \
                            Sim_reads_1 + \
                            " " + \
                            Sim_reads_2 + \
                            " 1> /dev/null")
            
            subprocess.run("wgsim -e 0.001 -N " + \
                            str(read_number) + \
                            " -1 " + \
                            str(read_length) + \
                            " -2 " + \
                            str(read_length) +\
                            " -S 999 " + \
                            genome_folder + \
                            "/" + \
                            fasta_file + \
                            " " + \
                            Sim_reads_1 + \
                            " " + \
                            Sim_reads_2 + \
                            " 1> /dev/null",
                            shell = True)

    return None

def merge_reads(output_folder):
    Simfolder = output_folder + "/SimReads/"
    Merged_R1 = output_folder + "/Merged_R1.fastq"
    Merged_R2 = output_folder + "/Merged_R2.fastq"
    print("Merging forward reads")
    subprocess.run("cat " + Simfolder + "*R1* >> " + Merged_R1, shell = True)
    print("Merging reverse reads")
    subprocess.run("cat " + Simfolder + "*R2* >> " + Merged_R2, shell = True)

    return None

def meta_assembly(output_folder, threads, input_R1, input_R2):
    subprocess.run(["megahit --num-cpu-threads",
                    str(threads),
                    "-1", input_R1, "-2", input_R2,
                    "--out-dir megahit --out-prefix MetaAssemble --tmp-dir /tmp"],
            shell = True)

    return None

def map_minimap2(output_folder, output_file, threads, input_ref, input_R1, input_R2):
    print("Generating " + output_file)
    print("minimap2 -t" + \
                    str(threads) + \
                    "-a -x sr" + \
                    input_ref + \
                    input_R1 + \
                    input_R2 + \
                    "| samtools view -@" + \
                    str(threads) + \
                    "-F 4 -h -u -T" + \
                    input_ref + \
                    "- | samtools sort -@" + \
                    str(threads) + \
                    "-T /tmp/temp.minimap2.bam -o" + \
                    output_file + \
                    "-")
    subprocess.run("minimap2 -t " + \
                    str(threads) + \
                    " -a -x sr " + \
                    " " + input_ref + \
                    " " + input_R1 + \
                    " " + input_R2 + \
                    " | samtools view -@ " + \
                    str(threads) + \
                    " -F 4 -h -u -T " + input_ref + " - | samtools sort -@ " + \
                    str(threads) + \
                    " -T /tmp/temp.minimap2.bam -o " + \
                    output_file + \
                    " - ",
                 shell = True)
    return None

def map_kallisto(output_folder, output_file, threads, input_ref, input_R1, input_R2):
    index = os.path.splitext(output_file)[0] + ".index"
    subprocess.run(["kallisto index --index", index, input_ref], shell = True)
    subprocess.run(["kallisto quant --threads=" + str(threads),
                    "--index", index, "--output-dir kallisto_map --pseudobam",
                    input_R1, input_R2], shell = True)
    subprocess.run(["samtools sort -@", str(threads) + "kallisto_map/pseudoalignments.bam >",
                    output_file], shell = True)
    return None

def map_BWA(output_folder, output_file, threads, input_ref, input_R1, input_R2):
    subprocess.run(["bwa index", input_ref],
                    shell = True)
    subprocess.run(["bwa mem -t", str(threads), input_ref, input_R1, input_R2,
                    "| samtools view -@ ", str(threads), "-F 4 -h -u -T",
                    input_ref, "- | samtools sort -@ ", str(threads),
                    "-T /tmp/temp.bwa.bam -o", output_file, "-"],
                 shell = True)
    return None

def contig_coverage(input_bam, output_file, output_folder):
    if os.path.isfile("depth/" + output_file) == False:
        print("generating " + output_file)
        subprocess.run("jgi_summarize_bam_contig_depths --minContigLength 2000 \
            --outputDepth depth/" + output_file + " alignments/" + input_bam + " 2> /dev/null",
            shell = True)
    return None

def binning(input_contigs, input_depth):
    output_file = "binning/" + os.path.splitext(input_depth)[0]
    print(output_file)
    print("metabat2 --minContig 2000 --inFile " + input_contigs + \
                    " --adbFile " + input_depth + "--outFile " + output_file)
    subprocess.run("metabat2 --minContig 2000 --inFile " + input_contigs + \
                    " --adbFile " + input_depth + " --outFile " + output_file,
                    shell = True)
    return None

def main():
    '''
    Orchestrate the execution of the program
    '''

    time = datetime.datetime.now()
    prefix = time.strftime("%Y%m%d-%H%M%S")
    options = parse_args(prefix)

    init_logging(options.log)

    genome_folder = options.input
    read_length = options.readlength
    read_depths = options.outdir + "/" + options.prefix + ".read_depths.tsv"
    output_folder = options.outdir
    threads = 8
    
    if not os.path.exists(output_folder):
        print("Making output directory")
        os.mkdir(output_folder)

    print("Generating read_depths")
    if os.path.isfile(read_depths) == False:
        generate_read_depths(genome_folder, read_length, read_depths, output_folder)

    print("Simulating reads")
    simulate_reads(genome_folder, read_depths, output_folder, read_length)
    if os.path.isfile("Merged_R1.fastq") == False or  os.path.isfile("Merged_R2.fastq") == False:
        merge_reads(output_folder)

    print("going to try assembly")
    if os.path.isdir(output_folder + "/megahit") == True:
        print("detected megahit folder already")
        if os.path.isfile(output_folder + "/megahit/MetaAssemble.contigs.fa") == True:
            print("detected megahit contig file, exiting assembly")
        else:
            print("couldn't detect contig file, probably failed assembly. Removing folder and start again")
            shutil.rmtree(output_folder + "/megahit")
            meta_assembly(output_folder, threads)
    else:
        meta_assembly(output_folder, threads)


    if os.path.isfile("alignments/MetaAssemble.minimap2.bam") == False:
        map_minimap2(output_folder, "alignments/MetaAssemble.minimap2.bam", threads, "megahit/MetaAssemble.contigs.fa", "Merged_R1.fastq", "Merged_R2.fastq")

    if os.path.isfile("alignments/MetaAssemble.kallisto.bam") == False:
        map_kallisto(output_folder, "alignments/MetaAssemble.kallisto.bam", threads, "megahit/MetaAssemble.contigs.fa", "Merged_R1.fastq", "Merged_R2.fastq")

    if os.path.isfile("alignments/MetaAssemble.bwa.bam") == False:
        map_BWA(output_folder, "alignments/MetaAssemble.bwa.bam", threads, "megahit/MetaAssemble.contigs.fa", "Merged_R1.fastq", "Merged_R2.fastq")

    for file in os.listdir("alignments"):
        if file.lower().endswith(".bam"):
            output_file = os.path.splitext(file)[0] + ".txt"
            if os.path.isfile("depth/" + output_file) == False:
                print("contig coverage of " + file)
                contig_coverage(file, output_file, output_folder)

    for file in os.listdir("depth"):
        if os.path.isfile("binning/" + os.path.splitext(file)[0] + ".1.fa") == False:
            binning("megahit/MetaAssemble.contigs.fa", file)

    for file in os.listdir("binning"):
        output_folder = "binning/"
        if os.path.isfile(output_folder + os.path.splitext(file)[0] + ".bam") == False:
            map_minimap2(output_folder, output_folder + os.path.splitext(file)[0] + ".bam", threads, output_folder + "/" + file, "Merged_R1.fastq", "Merged_R2.fastq")

#    remap_kallisto()

#    remap_BWA()

#    bin_coverage()

    #TODO: scan the reads that map to the each bin and see tally frequency of original genome?
#    contamination_chack()

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
