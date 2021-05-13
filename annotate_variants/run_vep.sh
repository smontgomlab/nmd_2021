#!/bin/bash

#Remove local Perl lib and cope from environment
unset PERL5LIB
export PATH=$(echo $PATH | tr ':' '\n' | rg -v Cope | tr '\n' ':')

conda activate vep

vep --cache \
    --dir /home/dnachun/montgomery_lab/.vep/ \
    --cache_version 88 \
    --fork 16 \
    --dir_plugins /home/dnachun/montgomery_lab/loftee/ \
    --plugin LoF,loftee_path:/home/dnachun/montgomery_lab/loftee/ \
    --tab \
    --force_overwrite `#overwrite existing files`\
    --total_length \
    --sift b \
    --polyphen b \
    --ccds \
    --uniprot \
    --symbol \
    --numbers \
    --domains \
    --canonical \
    --protein \
    --biotype \
    --uniprot \
    --tsl \
    --appris \
    --gene_phenotype --af \
    --af_1kg \
    --af_esp \
    --af_exac \
    --max_af \
    --pubmed \
    --variant_class \
    --mane\
    --verbose \
    --no_stats \
    --offline \
    -i variants_test \
    -o variants_annot_test.txt

#--hgvs \ - needs fasta
#--regulatory \ - don't use for cache 88
