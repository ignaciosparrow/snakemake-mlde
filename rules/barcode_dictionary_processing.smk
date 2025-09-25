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
        "python3 {input.script} "
        "-i {input.dictionary} "
        "-o {output.processed_dict} "
        "--outdir {params.outdir} "
        "--backbone_sub_thresh {params.bb_thresh}"