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
    return stringed_list.replace(" ", "").split(",")



rule_all = [
    os.path.join(config["SAMPLE_INFO_FILEPATH"]),  # init or provided
    os.path.join(config["OUT_DIR"], "merged.bam"),  # samtools merge
    expand(
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "merged_{chr}_calls.txt"
        ),  # palmer out 1 merged
        zip, chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "merged_{chr}_tsd_reads.txt"
        ),  # palmer out 2 merged
        zip, chr=get_chromosomes(config["CHROMOSOMES"])
    ),
]  # first step is making sample info file


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
        samtools index -b -o {output}
        """


rule run_palmer2_merged:
    input: rules.merge_bam_files.output
    output:
        calls = os.path.join(config["OUT_DIR"], "calls", "merged_{chr}_calls.txt"),
        tsd = os.path.join(config["OUT_DIR"], "tsd_reads", "merged_{chr}_tsd_reads.txt")
    params:
        palmer = config["PALMER_LOC"],
        chr = "{chr}",
        fa = config["REF_GEN"],
        ver = config["REF_VER"],
        out_dir = config["OUT_DIR"],
        mei = config["MEI"].upper(),
        mode = config["MODE"].lower()
    conda: "envs/palmer.yaml"
    threads: 2
    resources:
        mem_mb= 1000*7  # 7 gb
    shell:
        """
        {params.palmer} \
        --input {input} \
        --workdir {params.out_dir}/temp \
        --ref_fa {params.fa} \
        --ref_ver {params.ver} \
        --type {params.mei} \
        --mode {params.mode} \
        --chr {params.chr}
        
        mv {params.out_dir}/temp/sample_calls.txt {output.calls}
        mv {params.out_dir}/temp/sample_TSD_reads.txt {output.tsd}
        """

#rule run





