import sys
from py.helpers import *
configfile: "config/config.yaml"
module merged:
    snakefile: "../workflow/merged.smk"
module individual:
    snakefile: "../workflow/merged.smk"
use rule germline_filter from merged as *
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
                    "{chr}_{mei}_germline_polymorphs.vcf"
                ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
        mei=get_mei_type(config["MEI"])
    )
]


rule all: input: rule_all

rule overlap_polymorph:  # filter out singlets and low reads
    input:
        merged_bam=rules.germline_filter.output.bed
    output:
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "{chr}_{mei}_germline_polymorphs.vcf"
        )
    params:
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"]
    conda: "envs/bedtools.yaml"
    threads: 1
    resources:
        mem_mb = 500*7
    shell:
        """
        bedtools intersect -a HGSVC.MEI.Master.Set.Tier3.20220909.vcf.gz -b {input} > {output}
        """


