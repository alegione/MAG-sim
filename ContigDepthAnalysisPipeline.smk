shell.executable("/bin/bash")

from Bio import SeqIO
import sys
import subprocess
import os
import pandas
import numpy

IDS, = glob_wildcards("genomes/{id}.fna")
NUMS = len(IDS)


ruleorder: generate_read_depths > simulate_reads > merge_reads > meta_assembly > map_BWA > map_minimap2 > map_kallisto > contig_coverage

rule all:
     input: sims = expand("SimReads/{id}_R1.fastq", id=IDS),
            merged = "Merged_R1.fastq",
            meta_assemble = "megahit/MetaAssemble.contigs.fa",
            contig_coverage = "depth/MetaAssemble_depth.txt"


rule generate_read_depths:
    input: read_folder = os.getcwd() + "/genomes"
    output: output_filename = "read_depths.tsv"
    params: read_length = 300
    run:
        df = pandas.DataFrame({ "Filename":[], \
                                "GenomeName":[], \
                                "FastaLength":[], \
                                "RandomDepth":[], \
                                "ReadsRequired":[] \
                                })
        for filename in os.listdir(input['read_folder']):
            print(filename)
            seqSum = 0
            full_location = input['read_folder'] + "/" + filename
            for seq_record in SeqIO.parse(full_location, "fasta"):
                seqSum = seqSum + len(seq_record)
                print(seq_record.id + ":" + str(seqSum))

            depth = int(numpy.random.randint(10, 100+1))
            read_number = int(round(seqSum * depth / params['read_length']))
            record = pandas.DataFrame({ "Filename":[full_location], \
                                        "GenomeName":[seq_record.id], \
                                        "FastaLength":[int(seqSum)], \
                                        "RandomDepth":[depth], \
                                        "ReadsRequired":[read_number] \
                                        })

            df = df.append(record, ignore_index = True)

        print(df)
        df.to_csv(output['output_filename'], sep = "\t", header = True)

rule simulate_reads:
    input:  fasta_file = "genomes/{id}.fna",
            read_depth_file = "read_depths.tsv"
    output: Sim_reads_1 = "SimReads/{id}_R1.fastq",
            Sim_reads_2 = "SimReads/{id}_R2.fastq"
    # params: read_number = pandas.read_csv(input.read_depth_file,
    #                                         sep = "\t",
    #                                         names = ["GenomeName", "FastaLength", "RandomDepth", "ReadsRequired"
    #
    run:
        df = pandas.read_csv(input['read_depth_file'],
                            sep = "\t", index_col = 0)
        read_number = int(df.loc[df.Filename == input.fasta_file.rsplit('/',1)[1], 'ReadsRequired'].values[0])
        subprocess.run("wgsim -e 0.001 -N " + \
                        str(read_number) + \
                        " -1 150 -2 150 -S 999 " + \
                        os.getcwd() + \
                        "/" + \
                        input.fasta_file + \
                        " " + \
                        os.getcwd() + \
                        "/" + \
                        output.Sim_reads_1 + \
                        " " + \
                        os.getcwd() + \
                        "/" + \
                        output.Sim_reads_2, \
                        shell = True)

# rule merge_reads:
#     input:  merge_fastq_R1 = "SimReads/{id}_R1.fastq",
#             merge_fastq_R2 = "SimReads/{id}_R2.fastq"
#     output: tmp = "tmp/{id}.complete"
#
#     shell: r"""
#             cat \
#                 {input.merge_fastq_R1} \
#                 >> \
#                 Merged_R1.fastq \
#             && \
#             cat \
#                 {input.merge_fastq_R2} \
#                 >> \
#                 Merged_R2.fastq \
#             2> \
#             {output}
#             """
#
# rule generate_read_depths:
#     input:  read_folder = os.getcwd() + "/genomes",
#             fasta_file = "genomes/{id}.fna"
#     output: output_filename = "read_depths.tsv",
#             Sim_reads_1 = "SimReads/{id}_R1.fastq",
#             Sim_reads_2 = "SimReads/{id}_R2.fastq"
#     params: read_length = 300
#     run:
#         df = pandas.DataFrame({ "Filename":[], \
#                                 "GenomeName":[], \
#                                 "FastaLength":[], \
#                                 "RandomDepth":[], \
#                                 "ReadsRequired":[] \
#                                 })
#         for filename in os.listdir(input['read_folder']):
#             print(filename)
#             seqSum = 0
#             full_location = input['read_folder'] + "/" + filename
#             for seq_record in SeqIO.parse(full_location, "fasta"):
#                 seqSum = seqSum + len(seq_record)
#                 print(seq_record.id + ":" + str(seqSum))
#
#             depth = int(numpy.random.randint(10, 100+1))
#             read_number = int(round(seqSum * depth / params['read_length']))
#             record = pandas.DataFrame({ "Filename":[full_location], \
#                                         "GenomeName":[seq_record.id], \
#                                         "FastaLength":[int(seqSum)], \
#                                         "RandomDepth":[depth], \
#                                         "ReadsRequired":[read_number] \
#                                         })
#             SimCommand = "wgsim -e 0.001 -N " + read_number + " -1 150 -2 150 -S 999 " + input.fasta_file + " " + output.Sim_reads_1 + " " + output.Sim_reads_2
#             print(SimCommand)
#             subprocess.run(SimCommand)
#             df = df.append(record, ignore_index = True)
#
#         print(df)
#         df.to_csv(output['output_filename'], sep = "\t", header = True)

rule merge_reads:
    input:
    output: output_Merged_R1 = "Merged_R1.fastq",
            output_Merged_R2 = "Merged_R2.fastq"

    shell: "cat \
                SimReads/*R1*fastq \
                >> \
                {output.output_Merged_R1} \
            && \
            cat \
                SimReads/*R2*fastq \
                >> \
                {output.output_Merged_R2} \
            "

rule meta_assembly:
    input:  assembly_fastq_1 = "Merged_R1.fastq",
             assembly_fastq_2 = "Merged_R2.fastq"
    # input:  assembly_fastq_1 = expand("SimReads/{id}_R1.fastq", id=),
    #         assembly_fastq_2 = expand("SimReads/{id}_R2.fastq", id=IDS)

    output: "megahit/MetaAssemble.contigs.fa"

    shell: r"""
            megahit \
                --num-cpu-threads 4 \
                -1 {input.assembly_fastq_1} \
                -2 {input.assembly_fastq_2} \
                --out-dir megahit \
                --out-prefix MetaAssemble \
                --tmp-dir /tmp
            """

rule map_BWA:
    input:  contigs = "megahit/MetaAssemble.contigs.fa",
            Merged_R1 = "Merged_R1.fastq",
            Merged_R2 = "Merged_R2.fastq"
    output: bwa_align = "alignments/MetaAssemble.bwa.bam",
            bwa_flagstat = "alignments/MetaAssemble.bwa.flagstat"


    params: threads = 4

    shell: r"""
            bwa index {input.contigs}
            bwa mem \
                -t {params.threads} \
                {input.contigs} \
                {input.Merged_R1} \
                {input.Merged_R2} \
                | \
                samtools view \
                -@ {params.threads} \
                -F 4 \
                -h \
                -u \
                -T {input.contigs} \
                - | \
                samtools sort \
                -@ {params.threads} \
                -T /tmp/temp.bwa.bam \
                -o {output.bwa_align} -
           samtools flagstat {output.bwa_align} > {output.bwa_flagstat}
            """

rule map_minimap2:
    input:  contigs = "megahit/MetaAssemble.contigs.fa",
            Merged_R1 = "Merged_R1.fastq",
            Merged_R2 = "Merged_R2.fastq"

    output: minimap2_align = "alignments/MetaAssemble.minimap2.bam",
            minimap2_flagstat = "alignments/MetaAssemble.minimap2.flagstat"

    params: threads = 4

    shell: r"""
        minimap2 \
            -t {params.threads} \
            -a \
            -x sr \
            {input.contigs} \
            {input.Merged_R1} \
            {input.Merged_R2} \
             | \
             samtools view \
             -@ {params.threads} \
             -F 4 \
             -h \
             -u \
             -T {input.contigs} \
             - | \
             samtools sort \
             -@ {params.threads} \
             -T /tmp/temp.minimap2.bam \
             -o {output.minimap2_align} -
        samtools flagstat {output.minimap2_align} > {output.minimap2_flagstat}
        """

rule map_kallisto:
    input:  contigs = "megahit/MetaAssemble.contigs.fa",
            Merged_R1 = "Merged_R1.fastq",
            Merged_R2 = "Merged_R2.fastq"

    output: kallisto_index = "alignments/MetaAssemble.kallisto.index",
            kallisto_bam = "kallisto_map/pseudoalignments.sort.bam"

    params: threads = 4

    shell: r"""
            kallisto index \
                --index \
                {output.kallisto_index} \
                {input.contigs}
            kallisto quant \
                --threads={params.threads} \
                --index {output.kallisto_index} \
                --output-dir kallisto_map \
                --pseudobam \
                 {input.Merged_R1} \
                 {input.Merged_R2}
            samtools sort \
                -@ 4 \
                kallisto_map/pseudoalignments.bam \
                > \
                kallisto_map/pseudoalignments.sort.bam
            """

rule contig_coverage:
    input:  minimap2 = "alignments/MetaAssemble.minimap2.bam",
            kallisto = "kallisto_map/pseudoalignments.sort.bam",
            bwa = "alignments/MetaAssemble.bwa.bam"
    output: depthfile = "depth/MetaAssemble_depth.txt"

    shell: r"""
            jgi_summarize_bam_contig_depths \
                --minContigLength 2000 \
                --outputDepth {output.depthfile} \
                {input.minimap2} \
                {input.kallisto} \
                {input.bwa}  2> /dev/null
            """

rule binning:
    input:  contigs = "megahit/MetaAssemble.contigs.fa",
            depth_file = "depth/MetaAssemble_depth.txt"

    output: expand("binning/MetaAssemble_bin.{n}.fa", n=(range(1,len(os.listdir(os.getcwd() + "/binning")))))

    shell: r"""
            metabat2 \
                --verbose \
                --minContig 2000 \
                --inFile {input.contigs}     \
                --adbFile {input.depth_file} \
                --outFile binning/MetaAssemble_bin
            """

rule remap_to_bins_minimap2:
    input:  bins = "binning/MetaAssemble_bin.{n}.fa",
            Merged_R1 = "Merged_R1.fastq",
            Merged_R2 = "Merged_R2.fastq"

    output: minimap2_align = dynamic("alignments/bin-{n}.minimap2.bam")

    params: threads = 4

    shell: r"""
        minimap2 \
            -t {params.threads} \
            -a \
            -x sr \
            {input.bins} \
            {input.Merged_R1} \
            {input.Merged_R2} \
             | \
             samtools view \
             -@ {params.threads} \
             -F 4 \
             -h \
             -u \
             -T {input.bins} \
             - | \
             samtools sort \
             -@ {params.threads} \
             -T /tmp/bin-{n}.minimap2.bam \
             -o {output.minimap2_align} -
        """
rule remap_kallisto:
    input:  bins = "binning/MetaAssemble_bin.{n}.fa",
            Merged_R1 = "Merged_R1.fastq",
            Merged_R2 = "Merged_R2.fastq"

    output: kallisto_index = dynamic("alignments/MetaAssemble.kallisto-{n}.index"),
            kallisto_align = dynamic(directory("kallisto_map-{n}"))

    params: threads = 4

    shell: r"""
            kallisto index \
                --index \
                {output.kallisto_index} \
                {input.bins}
            kallisto quant \
                --threads={params.threads} \
                --index {output.kallisto_index} \
                --output-dir {output.kallisto_map} \
                --pseudobam \
                 {input.Merged_R1} \
                 {input.Merged_R2}
            samtools sort \
                -@ 4 \
                {output.kallisto_align}/pseudoalignments.bam \
                > \
                {output.kallisto_align}/pseudoalignments.sort.bam
            """

rule remap_BWA:
    input:  bins = "binning/MetaAssemble_bin.{n}.fa",
            Merged_R1 = "Merged_R1.fastq",
            Merged_R2 = "Merged_R2.fastq"
    output: bwa_align = dynamic("alignments/bin-{n}.bwa.bam")
            #,
            #bwa_flagstat = "alignments/MetaAssemble.bwa.flagstat"


    params: threads = 4

    shell: r"""
            bwa index {input.bins}
            bwa mem \
                -t {param.threads} \
                {input.bins} \
                {input.Merged_R1} \
                {input.Merged_R2} \
                | \
                samtools view \
                -@ {params.threads} \
                -F 4 \
                -h \
                -u \
                -T {input.bins} \
                - | \
                samtools sort \
                -@ {params.threads} \
                -T /tmp/bin-{n}.bwa.bam \
                -o {output.bwa_align} -
        """
        #   samtools flagstat {output.bwa_align} > {output.bwa_flagstat}
        #    """


rule bin_coverage:
    input:  minimap2 = expand("alignments/bin-{n}.minimap2.bam", n=(range(1,len(os.listdir(os.getcwd() + "/binning"))))),
            #kallisto = "kallisto_map-{n}",
            #bwa = "alignments/bin-{n}.bwa.bam"
    output: depthfile = dynamic("depth/bin-{n}_depth.txt")

    shell: r"""
            jgi_summarize_bam_contig_depths \
                --minContigLength 2000 \
                --outputDepth {output.depthfile} \
                {input.minimap2} \
            """
#    {input.kallisto}/pseudoalignments.sort.bam \
#    {input.bwa}
