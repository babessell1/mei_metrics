#!/bin/bash

# move all the outputs I want to directory to export
sort hgsvc_db_bed.bed > mei_demo/hgsvc_db_bed.bed
sort scratch_ln/output/bed/LINE/*/germline/* > mei_demo/LINE_germline.bed
sort scratch_ln/output/hgsvc_overlaps/LINE/*/germline/* > mei_demo/LINE_germline_overlaps.bed

sort scratch_ln/output/phased_nongerm/LINE/*/hp1/* > mei_demo/LINE_hp1_nongerm.bed
sort scratch_ln/output/phased_germ/LINE/*/hp1/* > mei_demo/LINE_hp1_germ.bed
sort scratch_ln/output/phased_nonpdb/LINE/*/hp1/* > mei_demo/LINE_hp1_nonpdb.bed
sort scratch_ln/output/hgsvc_overlaps/LINE/*/hp1/* > mei_demo/LINE_hp1_pdb.bed
sort scratch_ln/output/phased_pdb_nongerm/LINE/*/hp1/* > mei_demo/LINE_hp1_pdb_nongerm.bed
sort scratch_ln/output/phased_pdb_germ/LINE/*/hp1/* > mei_demo/LINE_hp1_pdb_germ.bed
sort scratch_ln/output/phased_nonpdb_germ/LINE/*/hp1/* > mei_demo/LINE_hp1_nonpdb_germ.bed
sort scratch_ln/output/phased_nonpdb_nongerm/LINE/*/hp1/* > mei_demo/LINE_hp1_nonpdb_nongerm.bed
sort scratch_ln/output/bed/LINE/*/hp1/* > mei_demo/LINE_hp1.bed

sort scratch_ln/output/phased_nongerm/LINE/*/hp2/* > mei_demo/LINE_hp2_nongerm.bed
sort scratch_ln/output/phased_germ/LINE/*/hp2/* > mei_demo/LINE_hp2_germ.bed
sort scratch_ln/output/phased_nonpdb/LINE/*/hp2/* > mei_demo/LINE_hp2_nonpdb.bed
sort scratch_ln/output/hgsvc_overlaps/LINE/*/hp2/* > mei_demo/LINE_hp2_pdb.bed
sort scratch_ln/output/phased_pdb_nongerm/LINE/*/hp2/* > mei_demo/LINE_hp2_pdb_nongerm.bed
sort scratch_ln/output/phased_pdb_germ/LINE/*/hp2/* > mei_demo/LINE_hp2_pdb_germ.bed
sort scratch_ln/output/phased_nonpdb_germ/LINE/*/hp2/* > mei_demo/LINE_hp2_nonpdb_germ.bed
sort scratch_ln/output/phased_nonpdb_nongerm/LINE/*/hp2/* > mei_demo/LINE_hp2_nonpdb_nongerm.bed
sort scratch_ln/output/bed/LINE/*/hp2/* > mei_demo/LINE_hp2.bed

sort scratch_ln/output/phased_nongerm/LINE/*/hp_non/* > mei_demo/LINE_hp_non_nongerm.bed
sort scratch_ln/output/phased_germ/LINE/*/hp_non/* > mei_demo/LINE_hp_non_germ.bed
sort scratch_ln/output/phased_nonpdb/LINE/*/hp_non/* > mei_demo/LINE_hp_non_nonpdb.bed
sort scratch_ln/output/hgsvc_overlaps/LINE/*/hp_non/* > mei_demo/LINE_hp_non_pdb.bed
sort scratch_ln/output/phased_pdb_nongerm/LINE/*/hp_non/* > mei_demo/LINE_hp_non_pdb_nongerm.bed
sort scratch_ln/output/phased_pdb_germ/LINE/*/hp_non/* > mei_demo/LINE_hp_non_pdb_germ.bed
sort scratch_ln/output/phased_nonpdb_germ/LINE/*/hp_non/* > mei_demo/LINE_hp_non_nonpdb_germ.bed
sort scratch_ln/output/phased_nonpdb_nongerm/LINE/*/hp_non/* > mei_demo/LINE_hp_non_nonpdb_nongerm.bed
sort scratch_ln/output/bed/LINE/*/hp_non/* > mei_demo/LINE_hp_non.bed

sort scratch_ln/output/phased_nongerm/LINE/*/hp_un/* > mei_demo/LINE_hp_un_nongerm.bed
sort scratch_ln/output/phased_germ/LINE/*/hp_un/* > mei_demo/LINE_hp_un_germ.bed
sort scratch_ln/output/phased_nonpdb/LINE/*/hp_un/* > mei_demo/LINE_hp_un_nonpdb.bed
sort scratch_ln/output/hgsvc_overlaps/LINE/*/hp_un/* > mei_demo/LINE_hp_un_pdb.bed
sort scratch_ln/output/phased_pdb_nongerm/LINE/*/hp_un/* > mei_demo/LINE_hp_un_pdb_nongerm.bed
sort scratch_ln/output/phased_pdb_germ/LINE/*/hp_un/* > mei_demo/LINE_hp_un_pdb_germ.bed
sort scratch_ln/output/phased_nonpdb_germ/LINE/*/hp_un/* > mei_demo/LINE_hp_un_nonpdb_germ.bed
sort scratch_ln/output/phased_nonpdb_nongerm/LINE/*/hp_un/* > mei_demo/LINE_hp_un_nonpdb_nongerm.bed
sort scratch_ln/output/bed/LINE/*/hp_un/* > mei_demo/LINE_hp_un.bed

