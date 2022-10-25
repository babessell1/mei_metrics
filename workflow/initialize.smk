import os
import sys
from warnings import warn
configfile: "config/config.yaml"
if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)

"""
rule all takes all outputs of other rules as inputs. We can make a conditional list of these inputs depending on 
whether they are provided to tell snakemake to either create the sample info files from the directories of phased and 
raw nanopore reads or alternative either or both can be provided to allow you to have more control over which samples
are being used. IF you have many samples and wish to use most but exclude some, you can run the initialization step
and remove entries from the raw and phased sample info files and they will be excluded from the analysis in the
proceeding steps.
"""
if not os.path.exists(config["RAW_INFO_FILEPATH"]):
    rule_all = [  # first step is making raw sample info file if not provided
        os.path.join(config["RAW_INFO_FILEPATH"]),
    ]
else:
    rule_all = []  # skip if provided or already made

if not os.path.exists(config["PHASED_INFO_FILEPATH"]):
    rule_all.extend([os.path.join(config["PHASED_INFO_FILEPATH"])])  # same for phased sample info file

rule_all.extend([  # the last rule in this module creates subsets of the HGSVC MEI database to
    os.path.join(config["OUT_DIR"], "ALU_pdb.bed"),
    os.path.join(config["OUT_DIR"], "LINE_pdb.bed"),
    os.path.join(config["OUT_DIR"], "SVA_pdb.bed"),
    os.path.join(config["OUT_DIR"], "HERVK_pdb.bed"),
    os.path.join(config["OUT_DIR"], "checkpoint1.chk")
])

rule all: input: rule_all


rule make_raw_sample_info_file:
    output: config["RAW_INFO_FILEPATH"]
    params:
        bam_dir = config["RAW_DIR"]
    threads: 1
    resources:
        mem_mb = 50
    shell:
        """
        ls -1 {params.bam_dir}/*.bam \
            | tr '\\n' '\\0' \
            | xargs -0 -n 1 basename \
            | tr '.' '\\t' \
            > {output}
        """
    
rule make_phased_sample_info_file:
    output: config["PHASED_INFO_FILEPATH"]
    params:
        phased_dir = config["PHASED_DIR"]
    threads: 1
    resources:
        mem_mb = 50
    shell:
        """
        ls -1 {params.phased_dir}/*.bam \
            | tr '\\n' '\\0' \
            | xargs -0 -n 1 basename \
            | tr '.' '\\t' \
            > {output}
        """

rule make_bed_pdbs:  # convert HGSVC MEI polymorph database (we will call pdb in this pipeline) to bed and separate mei
    output:
        alu = os.path.join(config["OUT_DIR"], "ALU_pdb.bed"),
        line = os.path.join(config["OUT_DIR"], "LINE_pdb.bed"),
        sva = os.path.join(config["OUT_DIR"], "SVA_pdb.bed"),
        hervk = os.path.join(config["OUT_DIR"], "HERVK_pdb.bed")
    params:
        pdb = config["POLYMORPHS"],
        alu_vcf = os.path.join(config["OUT_DIR"],"ALU_pdb.vcf"),
        line_vcf = os.path.join(config["OUT_DIR"],"LINE_pdb.vcf"),
        sva_vcf = os.path.join(config["OUT_DIR"],"SVA_pdb.vcf"),
        hervk_vcf = os.path.join(config["OUT_DIR"],"HERVK_pdb.vcf")
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        grep 'ALU' {params.pdb} > {params.alu_vcf}
        grep 'LINE' {params.pdb} > {params.line_vcf}
        grep 'SVA' {params.pdb} > {params.sva_vcf}
        grep 'HERVK' {params.pdb} > {params.hervk_vcf}
        
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.alu_vcf} > {output.alu}
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.line_vcf} > {output.line}
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.sva_vcf} > {output.sva}
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.hervk_vcf} > {output.hervk}
        
        """

rule make_dirs:  # also checkpoint 1
    input: config["RAW_INFO_FILEPATH"]
    output: os.path.join(config["OUT_DIR"], "checkpoint1.chk")
    params:
        out_dir = config["OUT_DIR"]
    threads: 1
    resources:
        mem_mb=10
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.out_dir}/temp
        mkdir -p {params.out_dir}/temp/slop/
        mkdir -p {params.out_dir}/calls
        mkdir -p {params.out_dir}/tsd_reads
        mkdir -p {params.out_dir}/bed
        mkdir -p {params.out_dir}/unique_calls
        mkdir -p {params.out_dir}/hgsvc_overlaps
        echo "" > {params.out_dir}/hgsvc_overlaps/no_alu.txt
        echo "" > {params.out_dir}/hgsvc_overlaps/no_line.txt
        echo "" > {params.out_dir}/hgsvc_overlaps/no_sva.txt
        echo "" > {params.out_dir}/hgsvc_overlaps/no_hervk.txt
        mkdir -p logs/palmer/
        mkdir -p logs/germ/
        
        touch {output}
        """
