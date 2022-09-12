import os
import sys
from warnings import warn

configfile: "config.yaml"

if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)

if not os.path.exists(config["SAMPLE_INFO_FILEPATH"]):
    rule_all = [os.path.join(config["SAMPLE_INFO_FILEPATH"])]  # first step is making sample info file
else:
    rule_all = []  # skip if provided or already made

rule_all.extend([os.path.join(config["OUT_DIR"], "checkpoint1.chk")])


rule all:
    input: rule_all

rule make_sample_info_file:
    input: directory(config["SAMP_DIR"])
    output: config["SAMPLE_INFO_FILEPATH"]
    threads: 1
    resources:
        mem_mb = 500
    shell:
        """
        ls -1 {input}/*.bam \
        | tr '\n' '\0' \
        | xargs -0 -n 1 basename \
        | tr '.' '\t' \
        > {output}
        """

rule touch_checkpoint1:
    input: config["SAMPLE_INFO_FILEPATH"]
    output: os.path.join(config["OUT_DIR"], "checkpoint1.chk")
    threads: 1
    resources:
        mem_mb=10
    shell:
        """
        touch {output}
        """
