# Snakefile for the SBPase MLDE project

##########
# Config #
##########
import glob
import os
import sys
from snakemake.utils import validate

configfile: "config.yaml"
shell.prefix("set -euo pipefail; ")

# Validate config
if "run_tag" not in config:
    raise ValueError("'run_tag' must be specified in config.yaml")

if "inputs" not in config or "ref" not in config["inputs"]:
    raise ValueError("Reference sequence path must be specified in config.yaml as config['inputs']['ref']")

RUN_TAG = config["run_tag"]

SAMPLE_NAME = config["sample_name"]  # Get sample name directly from config

# SAMPLE_TAGS is now just a list with a single item
SAMPLE_TAGS = [SAMPLE_NAME]

print(f"Using run tag: {RUN_TAG}")
print(f"Using sample name: {SAMPLE_NAME}")

wildcard_constraints:
    sample="|".join(SAMPLE_TAGS)

def output_path(filename):
    return f"outputs/{RUN_TAG}/{filename}"

def longread_input_path(filename):
    return f"inputs/longreads/{RUN_TAG}/{filename}"

def get_rule_all_inputs():
    """Define all final outputs in one place - used by both rule all and get_required_dirs"""
    return [
        *expand(output_path("{sample}.rarefaction_analysis.done"), sample=SAMPLE_TAGS)
        #*expand(output_path("{sample}.analyze_proteins.done"), sample=SAMPLE_TAGS)
    ]

# Define final output files for rule all
rule all:
    input:
        get_rule_all_inputs()



rule fetch_seqs:
    input:
        lambda wildcards: 
            glob.glob(f"inputs/longreads/{RUN_TAG}/{wildcards.sample}/*.fastq.gz")
        
    output:
        fastq = temp(output_path("raw_reads/{sample}.raw.fastq"))
    log:
        output_path("logs/fetch_seqs_{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.fastq}) $(dirname {log})
        echo "Looking for: {input}"
        echo "Extracting FASTQ for {wildcards.sample}" > {log}
        zcat {input} > {output.fastq} 2>> {log}
        echo "FASTQ decompressed to {output.fastq}" >> {log}
        """

if "target_length" not in config["parameters"]:
    raise ValueError("Please specify target_length in config file")

rule filter_length:
    input:
        output_path("raw_reads/{sample}.raw.fastq")
    output:
        temp(output_path("filtered_reads/{sample}_filtered.fastq"))
    log:
        output_path("logs/filter_length_{sample}.log")
    params:
        min_length = lambda wildcards: config["parameters"]["target_length"] - config["parameters"].get("length_tolerance", 0),
        max_length = lambda wildcards: config["parameters"]["target_length"] + config["parameters"].get("length_tolerance", 0),
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        cat {input} | \
        paste - - - - | \
        awk 'length($2) >= {params.min_length} && length($2) <= {params.max_length}' | \
        tr "\\t" "\\n" > {output} 2> {log}
        """

# Generate read length distribution table
rule make_read_length_table:
    input:
        output_path("filtered_reads/{sample}_filtered.fastq")
    output:
        table = output_path("{sample}_pacbio_read_length_table.tsv"),
        done = touch(output_path("{sample}.read_length_table.done")),

    log:
        output_path("logs/pacbio_read_length_table_{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.table}) $(dirname {output.done}) $(dirname {log})
        awk 'NR % 4 == 2 {{ print length($0) }}' {input} | \
        sort | uniq -c | \
        awk '{{ print $2 "\\t" $1 }}' > {output.table} 2> {log}
        """


# Index reference genome with BWA
rule index_reference:
    input:
        ref = config["inputs"]["ref"]
    output:
        multiext(config["inputs"]["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        output_path("logs/index_reference.log")
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """
        mkdir -p $(dirname {log})
        bwa index {input.ref} 2> {log}
        echo "Indexed reference FASTA"
        """

# Align reads to reference
rule align_long_reads:
    input:
        ref = config["inputs"]["ref"],
        reads = output_path("filtered_reads/{sample}_filtered.fastq"),
        index = multiext(config["inputs"]["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
        #multiext tells snakemake it has to have the same filename with all those extensions
    output:
        temp(output_path( "alignments/mapped_reads/{sample}.unsorted.bam" ))
    log:
        output_path("logs/align_reads_{sample}.log")
    threads: 16
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        echo "Aligning Pacbio reads to reference..." > {log}
        bwa mem -t {threads} {input.ref} {input.reads} 2>> {log} | \
        samtools view -Sb - > {output} 2>> {log} 
        echo "Aligned Pacbio reads to reference"
        """

rule sort_and_index:
    input:
        output_path( "alignments/mapped_reads/{sample}.unsorted.bam" )
    output:
        bam = output_path( "alignments/mapped_reads/{sample}.sorted.bam"),
        bai = output_path( "alignments/mapped_reads/{sample}.sorted.bam.bai")
    log:
        output_path("logs/sort_and_index_{sample}.log")
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        echo "Sorting and indexing"
        samtools sort -@ {threads} -o {output.bam} {input} 2>> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """


rule make_barcode_mutant_dictionary:
    input:
        bam = output_path("alignments/mapped_reads/{sample}.sorted.bam"),
        reference = config["inputs"]["ref"],
        script = "scripts/mutant_barcode_dictionary_maker.py",
        
    output:
        output_path("{sample}_mutant_barcode_dictionary.tsv")
    log:
        output_path("logs/make_barcode_mutant_dictionary_{sample}.log")
    params:
        progress = config["parameters"]["longread_progress"],
        
    shell:
        #in future add here a config parameter
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        echo "Making the barcode_mutant dictionary " > {log}
        python3 {input.script} \\
            -i {input.bam} \\
            -r {input.reference} \\
            -o {output} \\
            --progress {params.progress}  2>> {log}
        """
rule rarefaction_analysis:
    input:
        bam = output_path("alignments/mapped_reads/{sample}.sorted.bam"),
        reference = config["inputs"]["ref"],
        dict_script = "scripts/mutant_barcode_dictionary_maker.py",
        process_script = "scripts/dictionary_processing.py"
    output:
        rarefaction = output_path("{sample}_rarefaction_results.tsv"),
        done = touch(output_path("{sample}.rarefaction_analysis.done"))
    log:
        output_path("logs/rarefaction_analysis_{sample}.log")
    params:
        outdir = output_path("{sample}")
    shell:
        """
        python3 scripts/rarefaction_analysis.py \\
            --bam {input.bam} \\
            --reference {input.reference} \\
            --dict_script {input.dict_script} \\
            --process_script {input.process_script} \\
            --outdir {params.outdir} 2>> {log}
        """

rule process_dictionary:
    input:
        dictionary = output_path("{sample}_mutant_barcode_dictionary.tsv"),
        script = "scripts/dictionary_processing.py"
    output:
        processed_dict = output_path("{sample}_processed_mutant_barcode_dictionary.tsv")
    params:
        outdir = output_path("{sample}"),
        bb_thresh = config["parameters"]["backbone_tolerance_threshold"]
    shell:
        """
        mkdir -p $(dirname {output.processed_dict}) {params.outdir} $(dirname {log})
        python3 {input.script} \
            -i {input.dictionary} \
            -o {output.processed_dict} \
            --outdir {params.outdir} \
            --backbone_sub_thresh {params.bb_thresh} 2> {log}
        """

rule analyze_dictionary_proteins:
    input:
        dictionary = output_path("{sample}_processed_mutant_barcode_dictionary.tsv"),
        script = "scripts/protein_translator.py"
    output:
        translated_dictionary = output_path("{sample}_translated_mutant_dictionary.tsv"),
        done = touch(output_path("{sample}.analyze_proteins.done"))
        
    params:
        ref = config["inputs"]["prot_ref"],
        outdir = output_path("{sample}")
    log:
        output_path("logs/analyze_dictionary_proteins_{sample}.log")
    shell:
        """
        echo "Analyzing proteins" > {log}
        python3 {input.script} \\
        --input {input.dictionary} \\
        --output {output.translated_dictionary} \\
        --outdir {params.outdir} \\
        --ref {params.ref} 2>> {log}
        """



