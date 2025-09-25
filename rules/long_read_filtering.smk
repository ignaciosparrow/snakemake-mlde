# Filter reads by length (6000-6500 bp)
if "target_length" not in config["parameters"]:
    raise ValueError("Please specify target_length in config file")

rule filter_length:
    input:
        output_path("raw_reads/{sample}.raw.fastq")
    output:
        temp(output_path("filtered_reads/{sample}.fastq"))
    log:
        output_path("logs/filter_length_{sample}.log")
    params:
        min_length = lambda wildcards: config["parameters"]["target_length"] - config["parameters"].get("length_tolerance", 0),
        max_length = lambda wildcards: config["parameters"]["target_length"] + config["parameters"].get("length_tolerance", 0),
    shell:
        """
        cat {input} | \
        paste - - - - | \
        awk 'length($2) >= {params.min_length} && length($2) <= {params.max_length}' | \
        tr "\\t" "\\n" > {output} 2> {log}
        """

# Generate read length distribution table
rule make_read_length_table:
    input:
        output_path("filtered_reads/{sample}.fastq")
    output:
        output_path("{sample}_pacbio_read_length_table.tsv")

    log:
        output_path("logs/pacbio_read_length_table_{sample}.log")
    shell:
        """
        awk 'NR % 4 == 2 {{ print length($0) }}' {input} | \
        sort | uniq -c | \
        awk '{{ print $2 "\\t" $1 }}' > {output} 2> {log}
        """