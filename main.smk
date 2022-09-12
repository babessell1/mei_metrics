import os
import sys
import numpy as np
import pandas as pd
from warnings import warn


configfile: "config.yaml"

if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


def get_bam_filenames(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        ln_spl = line.split("\t")
        filenames = [
            os.path.join(
                ln_spl[0],
                ln_spl[1],
                ln_spl[2],
                ln_spl[3],
                ".bam"
            ) for line in handle.readlines()
        ]

    return expand(os.path.join(config["OUT_DIR"], "{bam_files}"), bam_files=get_bam_filenames(config["BAM_DIR"]))


rule_all = [
    os.path.join(config["SAMPLE_INFO_FILEPATH"]),
    os.path.join(config["OUT_DIR"], "combined.bam"),

]  # first step is making sample info file


rule all:
    input: rule_all


rule combine_bam_files:
    input:
        sample_info_file = config["SAMPLE_INFO_FILEPATH"],
    output:
        os.path.join(config["OUT_DIR"], "merged.bam")
    params:
        out_dir = config["OUT_DIR"]
    threads: 1
    resources:
        mem_mb = 500
    shell:
        """
        cat {input} | tr '\t' '.' > temp.txt
        samtools merge -1 -o "{params.out_dir}/merged.bam" -b temp.txt
        rm temp.txt 
        """



