import sys
import math
from py.helpers import *
configfile: "config/config.yaml"


if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


rule_all = [
    expand(
        os.path.join(
            config["OUT_DIR"],
            "hgsvc_overlaps",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_hgsvc_overlaps.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]), germ=True, phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]), germ=True, phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"], germ=True)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"], germ=True)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"], germ=True))
    )
]

rule all: input: rule_all


rule overlap_polymorph:  # overlap germline and phased haplotype reads with hgsvc polymorph database (link download?)
    input:
        os.path.join(
            config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}.bed"
        )
    output:
        os.path.join(
            config["OUT_DIR"],
            "hgsvc_overlaps",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_hgsvc_overlaps.bed"
        )
    params:
        bedops = config["BEDOPS_BIN"],
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"],
        mei = "{mei}",
        chr = "{chr}",
        haplo = "{haplo}",
        ref = config["REF_TAB"],
        slop = config["SLOP"],
        temp_dir = os.path.join(
            config["OUT_DIR"],
            "temp"
        ),
        temp_bed = os.path.join(
        config["OUT_DIR"],
            "temp",
            "bed",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_temp.bed"
            ),
        temp_pdb = os.path.join(
            config["OUT_DIR"],
            "temp",
            "temp_pdb",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_temp_pdb.bed"
        ),
        sorted_pdb= os.path.join(  # polymorph db
            config["OUT_DIR"],
            "temp",
            "sorted",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_sorted_pdb.bed"
        ),
        sorted_bed = os.path.join(
            config["OUT_DIR"],
            "temp",
            "sorted",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_sorted.bed"
        ),
    threads: 1
    resources:
        mem_mb = math.floor((1/36)*1000*180)  # 180 gb / 36 threads
    log:
        bedops="logs/bedops/{haplo}/{mei}/{chr}/{samp_id}_{barcode}_{chr}_{mei}_{haplo}.out",
    shell:
        """
        # make temp directories
        mkdir -p {params.temp_dir}/bed
        mkdir -p {params.temp_dir}/temp_pdb
        mkdir -p {params.temp_dir}/sorted
        if [[ $(wc -l <{input}) -ge 1 ]]; then
            mkdir -p {params.temp_dir}/pdb/vcf/{params.chr}_{params.mei}
            mkdir -p {params.temp_dir}/pdb/bed/{params.chr}_{params.mei}
            # add 1 to prevent bedops errors from same 1st and 2nd coord
            awk -v s=1 '{{print $1, $2, $3+s}}' {params.out_dir}/{params.mei}_pdb.bed > {params.temp_pdb}
            awk -v s=1 '{{print $1, $2, $3+s}}' {input} > {params.temp_bed} # add 1 to prevent bedops errors
            {params.bedops}/sort-bed {params.temp_bed} > {params.sorted_bed} # sort input 
            {params.bedops}/sort-bed {params.temp_pdb} > {params.sorted_pdb} # sort poly db bed
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --header \
                --range {params.slop} \
                --element-of 1 {params.sorted_bed} {params.sorted_pdb} \
                2> {log.bedops} 1> {output}
        else
            # if input has no mei calls than there are no overlaps
            touch {output}
            # track sample chr pairs with no meis
            echo "{params.chr}_{params.mei}" >> {params.out_dir}/hgsvc_overlaps/no_{params.mei}.txt
        fi
        """



