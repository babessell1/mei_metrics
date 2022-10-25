import sys
from py.helpers import *
configfile: "config/config.yaml"
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


rule_all = [  # rule all takes outputs of all other rules as inputs
    expand(  # palmer calls
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_calls.txt",
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"]), phased=True),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]), phased=True),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"])),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"]))
    ),
    expand(  # palmer tsds
        os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_tsd_reads.txt"
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"]), phased=True),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]), phased=True),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"])),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"]))
    ),
    expand(  #  call txt to bam
        os.path.join(
            config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}.bed"
        ),
        zip,
        mei=exp_meis(get_mei_type(config["MEI"]),phased=True),
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),phased=True),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"])),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"])),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"]))
    ),
]


rule all:
    input: rule_all


rule run_palmer_individual:
    input:
        bam=os.path.join(config["PHASED_DIR"],"{samp_id}.{barcode}.{haplo}.bam"),
        bai=os.path.join(config["PHASED_DIR"],"{samp_id}.{barcode}.{haplo}.bam.bai")
    output:
        calls=os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_calls.txt",
        ),
        tsd=os.path.join(
            config["OUT_DIR"],
            "tsd_reads",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_tsd_reads.txt"
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
        haplo="{haplo}",
    conda: "envs/palmer.yaml"
    threads: 2
    resources:
        mem_mb=1000 * 7  # 7 gb
    log:
        out="logs/palmer/{haplo}/{mei}/{chr}/{samp_id}_{barcode}_{chr}_{mei}_{haplo}.out",
        err="logs/palmer/{haplo}/{mei}/{chr}/{samp_id}_{barcode}_{chr}_{mei}_{haplo}.err",
    shell:
        """
        TEMP_DIR="{params.out_dir}/temp/{params.haplo}/{params.mei}/{params.chr}/{params.samp_id}/{params.barcode}"
        mkdir -p logs/palmer/{params.haplo}/{params.mei}/{params.chr}/
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
        calls = os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}.txt"
        ),
        bed = os.path.join(
            config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}.bed"
        )
    params:
        out_dir = config["OUT_DIR"],
        bam_dir = config["PHASED_DIR"],
        mei = "{mei}",
        chr = "{chr}",
        haplo="{haplo}",
        barcode="{barcode}",
        samp_id="{samp_id}",
        p_haplo = config["P_FILT"],
    threads: 2
    resources:
        mem_mb= 1000*7  # 7 gb
    shell:
        """
        # drop header
        awk '(NR>1) ' {input} > {output.calls}
        # reformat to bed
        awk '{{ print $2, $3, $5 }}' OFS=\\\\t {output.calls} > {output.bed}
        """
