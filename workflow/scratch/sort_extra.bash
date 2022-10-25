#!/bin/bash

# used this to make venn diagrams
mkdir -p temp
awk -v s=1 '{{print $1, $2, $3+s}}' mei_demo/LINE_germline.bed > temp/LINE_germline.bed
awk -v s=1 '{{print $1, $2, $3+s}}' scratch_ln/output/LINE_pdb.bed > temp/LINE_pdb.bed

/home/bbessell/software/bedops/bin/sort-bed temp/LINE_germline.bed > temp/LINE_germline_sorted.bed
/home/bbessell/software/bedops/bin/sort-bed temp/LINE_pdb.bed > temp/LINE_pdb_sorted.bed

/home/bbessell/software/bedops/bin/bedops  \
                --range 100 \
                --not-element-of 1 temp/LINE_germline_sorted.bed temp/LINE_pdb_sorted.bed \
                > mei_demo/LINE_nonpdb_germ.bed

/home/bbessell/software/bedops/bin/bedops  \
                --range 100 \
                --not-element-of 1 temp/LINE_pdb_sorted.bed temp/LINE_germline_sorted.bed \
                > mei_demo/LINE_pdb_nongerm.bed

/home/bbessell/software/bedops/bin/bedops  \
                --range 100 \
                --not-element-of 1 mei_demo/LINE_pdb_nongerm.bed mei_demo/LINE_hp1.bed \
                > mei_demo/LINE_pdb_nongerm_nonhp1.bed

/home/bbessell/software/bedops/bin/bedops  \
                --range 100 \
                --not-element-of 1 mei_demo/LINE_nonpdb_germ.bed mei_demo/LINE_hp1.bed \
                > mei_demo/LINE_nonpdb_germ_nonhp1.bed

/home/bbessell/software/bedops/bin/bedops  \
                --range 100 \
                --element-of 1 temp/LINE_germline_sorted.bed temp/LINE_pdb_sorted.bed \
                > mei_demo/LINE_pdb_germ.bed

/home/bbessell/software/bedops/bin/bedops  \
                --range 100 \
                --not-element-of 1 mei_demo/LINE_pdb_germ.bed mei_demo/LINE_hp1.bed \
                > mei_demo/LINE_pdb_germ_nonhp1.bed
