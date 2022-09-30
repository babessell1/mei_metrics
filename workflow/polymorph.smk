import sys
from workflow.py.helpers import *
configfile: "config/config.yaml"
module merged: snakefile: "merged.smk"
module individual: snakefile: "individual.smk"
use rule run_palmer2_merged from merged
use rule run_palmer_individual from individual
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)

rule_all = [

]


rule all: input: rule_all

rule filter_noise:  # filter out singlets and low reads
    input:
        merged_bam=rules.run_palmer2_merged.output

    output:
        os.path.join(
            config["OUT_DIR"],
            "calls",
            "{mei}",
            "{chr}",
            "merged_{chr}_{mei}_calls.txt"
        )
    params:
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"]
    conda: "envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb = 500
    script: "filter_noise.Rscript"


