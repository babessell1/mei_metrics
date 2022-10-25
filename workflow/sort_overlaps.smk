import sys
import math
from py.helpers import *
configfile: "config/config.yaml"


if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)


rule_all = [
    expand(  # phased calls with no germline overlaps
        os.path.join(
            config["OUT_DIR"],
            "phased_nongerm",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nongerm.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]), germ=False, phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]), germ=False, phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"], germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"], germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"], germ=False))
    ),
    expand(  # phased calls overlapping with germline
        os.path.join(
            config["OUT_DIR"],
            "phased_germ",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_germ.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),germ=False,phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]),germ=False,phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"],germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"],germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"],germ=False))
    ),
    expand(  # phased calls non hgsvc overlaps bed
        os.path.join(
            config["OUT_DIR"],
            "phased_nonpdb",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nonpdb.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),germ=False,phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]),germ=False,phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"],germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"],germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"],germ=False))
    ),
    expand(  # non-germ hgsvc overlaps bed
        os.path.join(
            config["OUT_DIR"],
            "phased_pdb_nongerm",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_pdb_nongerm.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),germ=False,phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]),germ=False,phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"],germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"],germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"],germ=False))
    ),
    expand(  # germline hgsvc overlaps bed
        os.path.join(
            config["OUT_DIR"],
            "phased_pdb_germ",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_pdb_germ.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),germ=False,phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]),germ=False,phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"],germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"],germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"],germ=False))
    ),
    expand(  # nongermline no hgsvc overlap bed
        os.path.join(
            config["OUT_DIR"],
            "phased_nonpdb_nongerm",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nonpdb_nongerm.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),germ=False,phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]),germ=False,phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"],germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"],germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"],germ=False))
    ),
    expand(  # germline no hgsvc overlap bed
        os.path.join(
            config["OUT_DIR"],
            "phased_nonpdb_germ",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nonpdb_germ.bed"
        ),
        zip,
        chr=exp_chromosomes(get_chromosomes(config["CHROMOSOMES"]),germ=False,phased=True),
        mei=exp_meis(get_mei_type(config["MEI"]),germ=False,phased=True),
        haplo=exp_haplo_types(get_haplo_type(config["PHASED_INFO_FILEPATH"],germ=False)),
        samp_id=exp_samp_ids(get_samp_id(config["PHASED_INFO_FILEPATH"],germ=False)),
        barcode=exp_barcodes(get_barcode(config["PHASED_INFO_FILEPATH"],germ=False))
    ),
]


rule all: input: rule_all


rule sort_overlaps:  # various subsets of mei's that might be wanted for downstream analysis
    input:
        phased = os.path.join(
            config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}.bed"
        ),
        germ = os.path.join(config["OUT_DIR"],
            "bed",
            "{mei}",
            "{chr}",
            "germline",
            "germ_germ_{chr}_{mei}_germline.bed"
        ),
        pdb_overlaps= os.path.join(
            config["OUT_DIR"],
            "hgsvc_overlaps",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_hgsvc_overlaps.bed"
        ),
    output:
        phased_nongerm = os.path.join(
            config["OUT_DIR"],
            "phased_nongerm",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nongerm.bed"
        ),
        phased_germ = os.path.join(
            config["OUT_DIR"],
            "phased_germ",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_germ.bed"
        ),
        phased_nonpdb = os.path.join(
            config["OUT_DIR"],
            "phased_nonpdb",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nonpdb.bed"
        ),
        phased_pdb_nongerm = os.path.join(
            config["OUT_DIR"],
            "phased_pdb_nongerm",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_pdb_nongerm.bed"
        ),
        phased_pdb_germ = os.path.join(
            config["OUT_DIR"],
            "phased_pdb_germ",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_pdb_germ.bed"
        ),
        phased_nonpdb_nongerm = os.path.join(
            config["OUT_DIR"],
            "phased_nonpdb_nongerm",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nonpdb_nongerm.bed"
        ),
        phased_nonpdb_germ = os.path.join(
            config["OUT_DIR"],
            "phased_nonpdb_germ",
            "{mei}",
            "{chr}",
            "{haplo}",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_nonpdb_germ.bed"
        ),
    params:
        bedops = config["BEDOPS_BIN"],
        out_dir = config["OUT_DIR"],
        bam_dir = config["RAW_DIR"],
        poly_db = config["POLYMORPHS"],
        slop = config["SLOP"],
        mei = "{mei}",
        chr = "{chr}",
        haplo = "{haplo}",
        ref = config["REF_TAB"],
        temp_dir = os.path.join(config["OUT_DIR"],"temp"),
        temp_phased = os.path.join(
            config["OUT_DIR"],
            "temp",
            "phased_bed",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_temp.bed"
        ),
        temp_germ = os.path.join(
            config["OUT_DIR"],
            "temp",
            "germ_bed",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_temp_germ.bed"
        ),
        sorted_phased = os.path.join(
            config["OUT_DIR"],
            "temp",
            "phased_sorted",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_sorted_phased.bed"
        ),
        sorted_germ = os.path.join(
            config["OUT_DIR"],
            "temp",
            "germ_sorted",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_sorted_germ.bed"
        ),
        sorted_pdb = os.path.join(
            config["OUT_DIR"],
            "temp",
            "pdb_sorted",
            "{samp_id}_{barcode}_{chr}_{mei}_{haplo}_sorted_pdb.bed"
        ),
    threads: 1
    resources:
        mem_mb = math.floor((1 / 36) * 1000 * 180)  # 180 gb / 36 threads
    log:
        bedops="logs/bedops_sorting/{haplo}/{mei}/{chr}/{samp_id}_{barcode}_{chr}_{mei}_{haplo}.out",
    shell:
        """
        #
        # make temp directorys. had some issues I think related to parallel io stuff I don't understand so... temp files 
        mkdir -p {params.temp_dir}/phased_bed
        mkdir -p {params.temp_dir}/germ_bed
        mkdir -p {params.temp_dir}/phased_sorted
        mkdir -p {params.temp_dir}/germ_sorted
        mkdir -p {params.temp_dir}/pdb_sorted
        
        # add 1 to prevent equal coord error, temp file.
        awk -v s=1 '{{print $1, $2, $3+s}}' {input.phased} > {params.temp_phased}
        awk -v s=1 '{{ print $1, $2, $3+s}}' {input.germ} > {params.temp_germ}
        # do not do for hgsvc overlaps (pdb) since we already did it in overlap.smk
        
        # sort each input into yet another temp file
        {params.bedops}/sort-bed {params.temp_phased} > {params.sorted_phased}  # sort input
        {params.bedops}/sort-bed {params.temp_germ} > {params.sorted_germ}  # sort germline
        {params.bedops}/sort-bed {input.pdb_overlaps} > {params.sorted_pdb}  # sort pdb overlaps
        echo "" > {log.bedops} # init log
        
        # get haplo calls that dont overlap with germline
        if [[ $(wc -l <{params.sorted_phased}) -lt 1 ]]; then
            touch {output.phased_nongerm} # if no phased reads, no phased nongermline reads.
        elif [[ $(wc -l <{params.sorted_germ}) -lt 1 ]]; then
            cp {params.sorted_phased} {output.phased_nongerm} # if no nongermline reads, all phased reads nongerm
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --not-element-of 1 {params.sorted_phased} {params.sorted_germ} \
                2>> {log.bedops} 1> {output.phased_nongerm}
        fi
        
        # get haplo calls that overlap with germline
        if [[ $(wc -l <{params.sorted_phased}) -lt 1 ]] || [[ $(wc -l <{params.sorted_germ}) -lt 1 ]]; then
            touch {output.phased_germ} # if either none, output is none  
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --element-of 1 {params.sorted_phased} {params.sorted_germ} \
                2>> {log.bedops} 1> {output.phased_germ}
        fi
        
        # get haplo calls that don't overlap with hgsvc database
        if [[ $(wc -l <{params.sorted_phased}) -lt 1 ]]; then
            touch {output.phased_nonpdb} # if not sorted_phased, no sorted_phased not in pdb
        elif [[ $(wc -l <{params.sorted_pdb}) -lt 1 ]]; then
            cp {params.sorted_phased} {output.phased_nonpdb} # if no pdbs, all are nonpdbs 
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --not-element-of 1 {params.sorted_phased} {params.sorted_pdb} \
                2>> {log.bedops} 1> {output.phased_nonpdb}   
        fi
         
        # get haplo calls that dont overlap with germline and do overlap with hgsvc db
        if [[ $(wc -l <{output.phased_nongerm}) -lt 1 ]] || [[ $(wc -l <{params.sorted_pdb}) -lt 1 ]]; then
            touch {output.phased_pdb_nongerm} # if no nongerm or no pdbs, no overlaps
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --element-of 1 {output.phased_nongerm} {params.sorted_pdb} \
                2>> {log.bedops} 1> {output.phased_pdb_nongerm} 
        fi
        
        # get calls for haplo that overlap with germline and do overlap with hgsvc db
        if [[ $(wc -l <{output.phased_germ}) -lt 1 ]] || [[ $(wc -l <{params.sorted_pdb}) -lt 1 ]]; then
            touch {output.phased_pdb_germ} # if no germ or no pdb no overlaps
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --element-of 1 {output.phased_germ} {params.sorted_pdb} \
                2>> {log.bedops} 1> {output.phased_pdb_germ}   
        fi
        
        # get calls for haplo that dont overlap with germline and dont overlap with hgsvc db
        if [[ $(wc -l <{output.phased_nongerm}) -lt 1 ]]; then
            touch {output.phased_nonpdb_nongerm} # no phased non germline calls, no non phased nongermline non pdb call
        elif [[ $(wc -l <{params.sorted_pdb}) -lt 1 ]]; then
            cp {output.phased_nongerm} {output.phased_nonpdb_nongerm} # if no pdb calls, all are nonpdb calls
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --not-element-of 1 {output.phased_nongerm} {params.sorted_pdb} \
                2>> {log.bedops} 1> {output.phased_nonpdb_nongerm}
        fi
        
        # get calls for haplo that overlap with germline and dont overlap with hgsvc overlaps
        if [[ $(wc -l <{output.phased_germ}) -lt 1 ]]; then
            touch {output.phased_nonpdb_germ} # if no germline calls, no germline calls not overlapping pdb database
        elif [[ $(wc -l <{params.sorted_pdb}) -lt 1 ]]; then
            cp {output.phased_germ} {output.phased_nonpdb_germ} # if no pdb calls, all germline calls are non pdb
        else
            {params.bedops}/bedops \
                --chrom {params.chr} \
                --range {params.slop} \
                --not-element-of 1 {output.phased_germ} {params.sorted_pdb} \
                2>> {log.bedops} 1> {output.phased_nonpdb_germ}     
        fi
        """