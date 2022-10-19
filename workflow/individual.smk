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
        mei=exp_meis(get_mei_type(config["MEI"]), phased=True),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]), phased=True),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"])),
        filt=exp_filter_types(get_filter_type(config["PHASED_INFO_FILEPATH"]))
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
        mei=exp_meis(get_mei_type(config["MEI"]), phased=True),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]), phased=True),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"])),
        filt=exp_filter_types(get_filter_type(config["PHASED_INFO_FILEPATH"]))
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "{filt}",
            "{samp_id}_{barcode}_{chr}_{mei}_{filt}.bed"
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"]),phased=True),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),phased=True),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"])),
        filt=exp_filter_types(get_filter_type(config["PHASED_INFO_FILEPATH"]))
    ),
]


rule all:
    input: rule_all


rule run_palmer_individual:
    input:
        bam=os.path.join(config["PHASED_DIR"],"{samp_id}.{barcode}.{filt}.bam"),
        bai=os.path.join(config["PHASED_DIR"],"{samp_id}.{barcode}.{filt}.bam.bai")
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

rule indiv_filter:
    input: rules.run_palmer_individual.output.calls,
    output:
        os.path.join(
            config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "{filt}",
            "{samp_id}_{barcode}_{chr}_{mei}_{filt}.bed"
        )
    params:
        out_dir = config["OUT_DIR"],
        bam_dir = config["PHASED_DIR"],
        mei = "{mei}",
        chr = "{chr}",
        filt="{filt}",
        barcode="{barcode}",
        samp_id="{samp_id}",
        p_filt = config["P_FILT"]
    threads: 2
    resources:
        mem_mb= 1000*7  # 7 gb
    shell:
        """
        # reformat to bed
        awk '{{ print $2, $3, $5 }}' OFS=\\\\t {input} > {output}
        """
