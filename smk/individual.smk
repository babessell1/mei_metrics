import os
import sys
from py.helpers import *
configfile: "config/config.yaml"
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


rule_all = [
    expand(
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "{filt}",
            "{samp_id}_{barcode}_{chr}_{mei}_{filt}_calls.txt",
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"])),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"])),
        samp_id=exp_samp_ids(get_samp_id(config["SAMPLE_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["SAMPLE_INFO_FILEPATH"])),
        filt=exp_filter_types(get_filter_type(config["FILTERS"]))
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "{mei}",
            "{chr}",
            "{filt}",
            "{samp_id}_{barcode}_{chr}_{mei}_{filt}_tsd_reads.txt"
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"])),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"])),
        samp_id=exp_samp_ids(get_samp_id(config["SAMPLE_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["SAMPLE_INFO_FILEPATH"])),
        filt=exp_filter_types(get_filter_type(config["FILTERS"]))
    ),
]


rule all: input: rule_all


rule run_palmer_individual:
    input:
        bam=lambda wildcards: os.path.join(get_bam_dir(wildcards.filt), "{samp_id}.{barcode}.{filt}.bam"),
        bai=lambda wildcards: os.path.join(get_bam_dir(wildcards.filt), "{samp_id}.{barcode}.{filt}.bam.bai"),
    output:
        calls=os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "{filt}",
            "{samp_id}_{barcode}_{chr}_{mei}_{filt}_calls.txt",
        ),
        tsd=os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "{mei}",
            "{chr}",
            "{filt}",
            "{samp_id}_{barcode}_{chr}_{mei}_{filt}_tsd_reads.txt"
        ),
    params:
        palmer=config["PALMER_LOC"],
        chr="{chr}",
        fa=config["REF_GEN"],
        ver=config["REF_VER"],
        out_dir=config["OUT_DIR"],
        mei="{mei}",
        mode=config["MODE"].lower(),
        barcode="{barcode}",
        samp_id="{samp_id}",
        filt="{filt}",
    conda: "envs/palmer.yaml"
    threads: 2
    resources:
        mem_mb=1000 * 7  # 7 gb
    log:
        out="logs/palmer/{filt}/{mei}/{chr}/{samp_id}_{barcode}_{chr}_{mei}_{filt}.out",
        err="logs/palmer/{filt}/{mei}/{chr}/{samp_id}_{barcode}_{chr}_{mei}_{filt}.err",
    shell:
        """
        TEMP_DIR="{params.out_dir}/temp/{params.filt}/{params.mei}/{params.chr}/{params.samp_id}/{params.barcode}"
        mkdir -p logs/palmer/{params.filt}/{params.mei}/{params.chr}/
        mkdir -p $TEMP_DIR
        {params.palmer} \
            --input {input.bam} \
            --workdir "${{TEMP_DIR}}/" \
            --ref_fa {params.fa} \
            --ref_ver {params.ver} \
            --type {params.mei} \
            --mode {params.mode} \
            --chr {params.chr} \
            2> {log.err} 1> {log.out}

        mv "${{TEMP_DIR}}/output.txt_calls.txt" {output.calls}
        mv "${{TEMP_DIR}}/output.txt_TSD_reads.txt" {output.tsd}
        """