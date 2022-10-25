#!/bin/bash

grep 'ALU' HGSVC.MEI.Master.Set.Tier3.20220909.vcf > ALU_pdb.vcf
grep 'LINE' HGSVC.MEI.Master.Set.Tier3.20220909.vcf > LINE_pdb.vcf

/home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < ALU_pdb.vcf > ALU_pdb.bed
/home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < LINE_pdb.vcf > LINE_pdb.bed





rule make_bed_pdbs:  # convert HGSVC MEI polymorph database (we will call pdb in this pipeline) to bed and separate mei
    output:
        alu = os.path.join(config["OUT_DIR"], "ALU_pdb.bed"),
        line = os.path.join(config["OUT_DIR"], "LINE_pdb.bed"),
        sva = os.path.join(config["OUT_DIR"], "SVA_pdb.bed"),
        hervk = os.path.join(config["OUT_DIR"], "HERVK_pdb.bed")
    params:
        alu_vcf = os.path.join(config["OUT_DIR"],"ALU_pdb.vcf"),
        line_vcf = os.path.join(config["OUT_DIR"],"LINE_pdb.vcf"),
        sva_vcf = os.path.join(config["OUT_DIR"],"SVA_pdb.vcf"),
        hervk_vcf = os.path.join(config["OUT_DIR"],"HERVK_pdb.vcf")

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
        grep 'ALU' {input} > {params.alu.vcf}
        grep 'LINE' {input} > {params.line.vcf}
        grep 'SVA' {input} > {params.sva.vcf}
        grep 'HERVK' {input} > {params.hervk.vcf}

        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.alu.vcf} > {output.alu}
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.line.vcf} > {output.line}
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.sva.vcf} > {output.sva}
        /home/bbessell/software/bedops/bin/vcf2bed --do-not-sort < {params.hervk.vcf} > {output.hervk}


        """


rule_all.extend([
    os.path.join(config["OUT_DIR"], "ALU_pdb.bed"),
    os.path.join(config["OUT_DIR"], "LINE_pdb.bed"),
    os.path.join(config["OUT_DIR"], "SVA_pdb.bed"),
    os.path.join(config["OUT_DIR"], "HERVK_pdb.bed"),
    os.path.join(config["OUT_DIR"], "checkpoint1.chk")
])