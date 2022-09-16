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


def get_chromosomes(stringed_list: str) -> list[str]:
    try:
        chroms = stringed_list.replace(" ","").split(",")
    except:
        raise SyntaxError(
            'Please set chromosomes in config.yaml as a string in the format'
            '"chromosome1, chromosome2, chromosome3, ..."')

    return chroms


def get_mei_type(stringed_list: str) -> list[str]:
    try:
        meis = stringed_list.replace(" ", "").split(",")
    except:
        raise SyntaxError(
            'Please set MEI type(s) in config.yaml as a string in the format'
            '"TYPE1, TYPE2, TYPE3, ..."')

    meis = [mei.upper() for mei in meis]
    if any(mei not in ["LINE", "ALU", "SVA", "HERVK"] for mei in meis):
        raise ValueError("Accepted MEI types are LINE, ALU, SVA, HERVK")

    return meis


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
        bam_dir = config["BAM_DIR"]
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
    shell:
        """
        mkdir -p {params.out_dir}/temp/{params.chr}/{params.mei}
        {params.palmer} \
            --input {input.bam} \
            --workdir {params.out_dir}/temp/{params.chr}/{params.mei} \
            --ref_fa {params.fa} \
            --ref_ver {params.ver} \
            --type {params.mei} \
            --mode {params.mode} \
            --chr {params.chr}
        
        mv {params.out_dir}/temp/{params.mei}/{params.chr}/output.txt_calls.txt {output.calls}
        mv {params.out_dir}/temp/{params.mei}/{params.chr}/output.txt_TSD_reads.txt {output.tsd}
        """


#rule run





