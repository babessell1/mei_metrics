import sys
from py.helpers import *
configfile: "config/config.yaml"
module merged:
    snakefile: "../config/config.yaml"
module individual:
    snakefile: "../config/config.yaml"
use rule germline_filter from merged as *
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)

rule_all = [
    expand(
        os.path.join(
                    config["OUT_DIR"],
                    "polymorphs",
                    "{chr}_{mei}_germline_polymorphs.vcf"
                ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
        mei=get_mei_type(config["MEI"])
    )
]


rule all: input: rule_all

rule overlap_polymorph:  # filter out singlets and low reads
    input:
        merged_bed=os.path.join(config["OUT_DIR"], "bed", "germline_{chr}_{mei}.bed")  # fix filename inconsistancy
    output:
        os.path.join(
            config["OUT_DIR"],
            "polymorphs",
            "{chr}_{mei}_germline_polymorphs.vcf"
        )
    params:
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"],
        poly_db = config["POLYMORPHS"],
        mei = "{mei}",
        chr = "{chr}",
        temp_file = os.path.join(
            config["OUT_DIR"],
            "temp",
            "slop",
            "{chr}_{mei}_slop_temp.bed"
        ),
        ref = config["REF_TAB"]
    conda: "envs/bed_bam.yaml"
    threads: 1
    resources:
        mem_mb = 500*7
    shell:
        """
        ## table browser cytoband whole genome
        #bedtools genomecov -i {{input}} -g {{params.ref}}.bed > {{output.cov}}
        # bedtools slop ( 100 bp )
        PTEMP="{params.chr}_{params.mei}_{params.poly_db}"  # temp polydb
        cp {params.poly_db} $PTEMP
        #gunzip "{{PTEMP}}.gz"
        if [[ $(wc -l <{input}) -ge 2 ]]; then
            bedtools slop \
                -i {input} \
                -g {params.ref} \
                -r 100 -l 100 \
            > {params.temp_file}
        fi
        if [[ -f "{params.temp_file}" ]] ; then
            bedtools intersect \
                -a $PTEMP \
                -b {params.temp_file} \
                -wo \
            > {output}
        fi
        if [[ ! -f "{params.temp_file}" ]] ; then
            touch {output}
            echo "{params.chr}_{params.mei}" >> {params.out_dir}/polymorphs/no_polymorphs.txt
        fi
        rm -r $PTEMP
        """



## visualize -> igv, compare to MEI, pairwise alignment with an alu, word doc


