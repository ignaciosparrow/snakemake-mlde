
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
        bwa index {input.ref} 2> {log}
        echo "Indexed reference FASTA"
        """

# Align reads to reference
rule align_long_reads:
    input:
        ref = config["inputs"]["ref"],
        reads = output_path("filtered_reads/{sample_tag}_filtered.fastq"),
        index = multiext(config["inputs"]["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
        #multiext tells snakemake it has to have the same filename with all those extensions
    output:
        temp(output_path( "alignments/mapped_reads/{sample_tag}.unsorted.bam" ))
    log:
        output_path("logs/align_reads_{sample_tag}.log")
    threads: 8
    shell:
        """
        echo "Aligning Pacbio reads to reference..." > {log}
        bwa mem -t {threads} {input.ref} {input.reads} 2>> {log} | \
        samtools view -Sb - > {output} 2>> {log} 
        echo "Aligned Pacbio reads to reference"
        """

rule sort_and_index:
    input:
        output_path( "alignments/mapped_reads/{sample_tag}.unsorted.bam" )
    output:
        bam = output_path( "alignments/mapped_reads/{sample_tag}.sorted.bam"),
        bai = output_path( "alignments/mapped_reads/{sample_tag}.sorted.bam.bai")
    log:
        output_path("logs/sort_and_index_{sample_tag}.log")
    threads: 4
    shell:
        """
        echo "Sorting and indexing"
        samtools sort -@ {threads} -o {output.bam} {input} 2>> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """

rule call_variants_freebayes:
    input:
        bam = output_path("alignments/mapped_reads/{sample_tag}.sorted.bam"),
        bai = output_path("alignments/mapped_reads/{sample_tag}.sorted.bam.bai"),
        reference = config["inputs"]["ref"]
    output:
        vcf = output_path("variants/{sample_tag}.vcf")
    log:
        output_path("logs/call_variants_freebayes_{sample_tag}.log")
    params:
        min_mapping_quality = 30,  # Increased from 20
        min_base_quality = 30,     # Increased from 20
        min_alternate_fraction = 0.01,  # More sensitive for rare variants
        min_alternate_count = 1,   # Lowered from 2 (but verify with IGV!)
        extra_flags = "--min-repeat-entropy 1 --haplotype-length 0 --ploidy 1 --no-partial-observations --standard-filters"
    threads: 8
    shell:
        """
        echo "Calling plasmid variants with freebayes (optimized for backbone indels)" > {log}
        mkdir -p $(dirname {output.vcf})
        freebayes \
            --fasta-reference {input.reference} \
            --bam {input.bam} \
            --min-mapping-quality {params.min_mapping_quality} \
            --min-base-quality {params.min_base_quality} \
            --min-alternate-fraction {params.min_alternate_fraction} \
            --min-alternate-count {params.min_alternate_count} \
            {params.extra_flags} \
            > {output.vcf} 2>> {log}
        """

rule make_barcode_mutant_dictionary:
    input:
        bam = output_path("alignments/mapped_reads/{sample_tag}.sorted.bam"),
        reference = config["inputs"]["ref"],
        script = "scripts/mutant_barcode_dictionary_maker.py",
        vcf = output_path("variants/{sample_tag}.vcf")  # Add VCF as input
    output:
        output_path("{sample_tag}_mutant_barcode_dictionary.tsv")
    log:
        output_path("logs/make_barcode_mutant_dictionary_{sample_tag}.log")
    params:
        progress = config["parameters"]["longread_progress"],
        
    shell:
        """
        echo "Making the barcode_mutant dictionary " > {log}
        python3 {input.script} \\
            -i {input.bam} \\
            -r {input.reference} \\
            -o {output} \\
            --vcf {input.vcf} \\
            --progress {params.progress}  2>> {log}
        """