import os
import sys
import numpy as np
import pandas as pd
from warnings import warn
from helpers import get_mei_type, get_chromosomes
configfile: "config.yaml"

if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


rule_all = [
    os.path.join(config["SAMPLE_INFO_FILEPATH"]),  # init or provided
    os.path.join(config["OUT_DIR"], "merged.bam"),  # samtools merge
    os.path.join(config["OUT_DIR"], "merged.bai"),  # index bam file
    expand(
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "merged_{chr}_{mei}_calls.txt"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
        mei=get_mei_type(config["MEI"])
    ),  # palmer calls (merged)
    expand(
        os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "{mei}",
            "{chr}",
            "merged_{chr}_{mei}_tsd_reads.txt"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
        mei=get_mei_type(config["MEI"])
    ),  # palmer TSD reads (merged)
]


rule all:
    input: rule_all


rule merge_bam_files:
    input:
        sample_info_file = config["SAMPLE_INFO_FILEPATH"],
    output:
        os.path.join(config["OUT_DIR"], "merged.bam")
    params:
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"]
    conda: "envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb = 500
    shell:
        """
        cat {input} \
            | tr '\\t' '.' \
            | while read line; do echo "{params.bam_dir}/$line"; done \
            > temp.txt
            
        samtools merge -1 -o "{params.out_dir}/merged.bam" -b temp.txt
        rm temp.txt
        """


rule index_merged_bam:
    input: rules.merge_bam_files.output
    output: os.path.join(config["OUT_DIR"], "merged.bai")
    conda: "envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb = 500
    shell:
        """
        samtools index -b {input} {output}
        """


rule run_palmer2_merged:
    input:
        bam=rules.merge_bam_files.output,
        bai=rules.index_merged_bam.output
    output:
        calls = os.path.join(config["OUT_DIR"], "calls", "{mei}", "{chr}", "merged_{chr}_{mei}_calls.txt"),
        tsd = os.path.join(config["OUT_DIR"], "tsd_reads", "{mei}", "{chr}", "merged_{chr}_{mei}_tsd_reads.txt")
    params:
        palmer = config["PALMER_LOC"],
        chr = "{chr}",
        fa = config["REF_GEN"],
        ver = config["REF_VER"],
        out_dir = config["OUT_DIR"],
        mei = "{mei}",
        mode = config["MODE"].lower()
    conda: "envs/palmer.yaml"
    threads: 2
    resources:
        mem_mb= 1000*7  # 7 gb
    log:
        out = "logs/palmer/{mei}_{chr}.out",
        err = "logs/palmer/{mei}_{chr}.err",
    shell:
        """
        mkdir -p {params.out_dir}/temp/{params.mei}/{params.chr}
        {params.palmer} \
            --input {input.bam} \
            --workdir {params.out_dir}/temp/{params.mei}/{params.chr}/ \
            --ref_fa {params.fa} \
            --ref_ver {params.ver} \
            --type {params.mei} \
            --mode {params.mode} \
            --chr {params.chr} \
            2> {log.err} 1> {log.out}
            
        mv {params.out_dir}/temp/{params.mei}/{params.chr}/output.txt_calls.txt {output.calls}
        mv {params.out_dir}/temp/{params.mei}/{params.chr}/output.txt_TSD_reads.txt {output.tsd}
        """




