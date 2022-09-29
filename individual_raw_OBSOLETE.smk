import os
import sys
import numpy as np
import pandas as pd
from warnings import warn
from tabulate import tabulate

configfile: "config.yaml"

if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


def get_samp_id(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[0] for line in handle.readlines()]
    return ids


def get_barcode(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        bars = [line.split("\t")[1] for line in handle.readlines()]
    return bars


def get_chromosomes(stringed_list: str) -> list[str]:
    try:
        chroms = stringed_list.replace(" ","").split(",")
    except:
        raise SyntaxError(
            'Please set chromosomes in config.yaml as a string in the format'
            '"chromosome1, chromosome2, chromosome3, ..."')

    return list(set(chroms))


def get_mei_type(stringed_list: str) -> list[str]:
    try:
        meis = stringed_list.replace(" ","").split(",")
    except:
        raise SyntaxError(
            'Please set MEI type(s) in config.yaml as a string in the format'
            '"TYPE1, TYPE2, TYPE3, ..."')

    meis = [mei.upper() for mei in meis]
    if any(mei not in ["LINE", "ALU", "SVA", "HERVK"] for mei in meis):
        raise ValueError("Accepted MEI types are LINE, ALU, SVA, HERVK")

    return list(set(meis))


def exp_samp_ids(samp_ids: list[str]) -> list[str]:
    with open("temp.txt", "w") as handle:
        handle.write(
            tabulate(
                samp_ids * len(get_mei_type(config["MEI"])) * len(get_chromosomes(config["CHROMOSOMES"]))
            )
        )
    return samp_ids*len(get_mei_type(config["MEI"]))*len(get_chromosomes(config["CHROMOSOMES"]))


def exp_barcodes(barcodes: list[str]) -> list[str]:
    return barcodes*len(get_mei_type(config["MEI"]))*len(get_chromosomes(config["CHROMOSOMES"]))


def exp_chromosomes(chromosomes: list[str]) -> list[str]:
    return chromosomes*len(get_samp_id(config["SAMPLE_INFO_FILEPATH"]))*len(get_mei_type(config["MEI"]))


def exp_meis(meis: list[str]) -> list[str]:
    return meis*len(get_samp_id(config["SAMPLE_INFO_FILEPATH"]))*len(get_chromosomes(config["CHROMOSOMES"]))


rule_all = [
    expand(
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "raw",
            "{samp_id}_{barcode}_{chr}_{mei}_calls.txt",
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"])),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"])),
        samp_id=exp_samp_ids(get_samp_id(config["SAMPLE_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["SAMPLE_INFO_FILEPATH"])),
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "{mei}",
            "{chr}",
            "raw",
            "{samp_id}_{barcode}_{chr}_{mei}_tsd_reads.txt",
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"])),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"])),
        samp_id=exp_samp_ids(get_samp_id(config["SAMPLE_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["SAMPLE_INFO_FILEPATH"])),
    ),
]


rule all:
    input: rule_all


rule run_palmer_raw:
    input: os.path.join(config["RAW_DIR"], "{samp_id}.{barcode}.Nanopore.sorted.bam")
    output:
        calls=os.path.join(config["OUT_DIR"],"calls","{mei}","{chr}","raw","{samp_id}_{barcode}_{chr}_{mei}_calls.txt"),
        tsd=os.path.join(config["OUT_DIR"],"tsd_reads","{mei}","{chr}","raw","{samp_id}_{barcode}_{chr}_{mei}_tsd_reads.txt")
    params:
        palmer=config["PALMER_LOC"],
        chr="{chr}",
        fa=config["REF_GEN"],
        ver=config["REF_VER"],
        out_dir=config["OUT_DIR"],
        mei="{mei}",
        mode=config["MODE"].lower(),
        barcode="{barcode}",
        samp_id="{samp_id}"
    conda: "envs/palmer.yaml"
    threads: 2
    resources:
        mem_mb=1000 * 7  # 7 gb
    log:
        out="logs/palmer/{samp_id}/{samp_id}_{barcode}_{mei}_{chr}.out",
        err="logs/palmer/{samp_id}/{samp_id}_{barcode}_{mei}_{chr}.err",
    shell:
        """
        mkdir -p logs/palmer/{params.samp_id}
        mkdir -p {params.out_dir}/temp/{params.samp_id}/{params.barcode}/{params.mei}/{params.chr}
        {params.palmer} \
            --input {input} \
            --workdir {params.out_dir}/temp/{params.samp_id}/{params.barcode}/{params.mei}/{params.chr}/ \
            --ref_fa {params.fa} \
            --ref_ver {params.ver} \
            --type {params.mei} \
            --mode {params.mode} \
            --chr {params.chr} \
            2> {log.err} 1> {log.out}

        mv {params.out_dir}/temp/{params.samp_id}/{params.barcode}/{params.mei}/{params.chr}/output.txt_calls.txt {output.calls}
        mv {params.out_dir}/temp/{params.samp_id}/{params.barcode}/{params.mei}/{params.chr}/output.txt_TSD_reads.txt {output.tsd}
        """



