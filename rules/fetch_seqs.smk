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
        echo "Looking for: {input}"
        echo "Extracting FASTQ for {wildcards.sample}" > {log}
        zcat {input} > {output.fastq} 2>> {log}
        echo "FASTQ decompressed to {output.fastq}" >> {log}
        """