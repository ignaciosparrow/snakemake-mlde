rule rarefaction_analysis:
    input:
        bam = output_path("alignments/mapped_reads/{sample}.sorted.bam"),
        reference = config["inputs"]["ref"],
        dict_script = "scripts/mutant_barcode_dictionary_maker.py",
        process_script = "scripts/dictionary_processing.py"
    output:
        results = output_path("{sample}_rarefaction_results.tsv"),
        done = touch(output_path("{sample}.rarefaction_analysis.done"))
    log:
        output_path("logs/rarefaction_analysis_{sample}.log")
    params:
        bb_thresh = config["parameters"]["backbone_tolerance_threshold"]
    shell:
        """
        python3 scripts/rarefaction_analysis.py \\
            --bam {input.bam} \\
            --reference {input.reference} \\
            --dict_script {input.dict_script} \\
            --process_script {input.process_script} \\
            --bb_thresh {params.bb_thresh} \\
            --output {output.results} 2> {log}
        """