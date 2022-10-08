import sys
import math
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
                    "{chr}_{mei}_germline_polymorphs.bed"
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
            "{chr}_{mei}_germline_polymorphs.bed"
        )
    params:
        bedops = config["BEDOPS_BIN"],
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"],
        poly_db = config["POLYMORPHS"],
        mei = "{mei}",
        chr = "{chr}",
        ref= config["REF_TAB"],
        temp_dir = os.path.join(
            config["OUT_DIR"],
            "temp"
        ),
        temp_bed = os.path.join(
        config["OUT_DIR"],
            "temp",
            "bed",
            "{chr}_{mei}_temp.bed"
            ),
        vcf_pdb = os.path.join(
            config["OUT_DIR"],
            "temp",
            "pdb",
            "vcf",
            "{chr}_{mei}_pdb.vcf"
        ),
        bed_pdb = os.path.join(
            config["OUT_DIR"],
            "temp",
            "pdb",
            "bed",
            "{chr}_{mei}_pdb.bed"
        ),
        sorted_pdb= os.path.join(  # polymorph db
            config["OUT_DIR"],
            "temp",
            "sorted",
            "{chr}_{mei}_sorted_pdb.bed"
        ),
        sorted_bed = os.path.join(
            config["OUT_DIR"],
            "temp",
            "sorted",
            "{chr}_{mei}_sorted.bed"
        ),
    threads: 9
    resources:
        mem_mb = math.floor(1000*180/9)
    log:
        bedops="logs/bedops/{mei}_{chr}.err",
        #bedtools="logs/bedtools/{mei}_{chr}.err",
    shell:
        """
        mkdir -p {params.temp_dir}/pdb
        mkdir -p {params.temp_dir}/sorted
        mkdir -p {params.temp_dir}/pdb/vcf
        mkdir -p {params.temp_dir}/pdb/bed
        cp {params.poly_db} {params.vcf_pdb}
        if [[ $(wc -l <{input}) -ge 2 ]]; then
            mkdir -p {params.temp_dir}/pdb/vcf/{params.chr}_{params.mei}
            mkdir -p {params.temp_dir}/pdb/bed/{params.chr}_{params.mei}
                awk '$2<$3 {{print}}' {input} > {params.temp_bed}  # remove nonsensical coordinates
                {params.bedops}/sort-bed {params.temp_bed} > {params.sorted_bed}  # sort input
                {params.bedops}/convert2bed -i vcf --do-not-sort < {params.vcf_pdb} > {params.bed_pdb}  # vcf to bed 
                {params.bedops}/sort-bed {params.bed_pdb} > {params.sorted_pdb}  # sort poly db bed
                {params.bedops}/bedops \
                    --chrom {params.chr} \
                    --header \
                    --range -100:100 \
                    --element-of 10% {params.sorted_pdb} {params.sorted_bed} \
                    2>> {log.bedops} 1> {output}
        else
            touch {output}
            echo "{params.chr}_{params.mei}" >> {params.out_dir}/polymorphs/no_meis.txt
        fi
        """



## visualize -> igv, compare to MEI, pairwise alignment with an alu, word doc


