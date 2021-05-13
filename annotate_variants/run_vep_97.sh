#!/bin/bash

#Remove local Perl lib and cope from environment
#unset PERL5LIB
#export PATH=$(echo $PATH | tr ':' '\n' | rg -v Cope | tr '\n' ':')

conda activate vep

vep --cache \
    --dir /home/dnachun/montgomery_lab/.vep/ \
    --cache_version 97 \
    --force_overwrite `#overwrite existing files`\
    --fork 16 \
    --everything \
    --tab \
    --af_gnomad \
    --verbose \
    --no_stats \
    --offline \
    -i variants_for_gnomad.txt \
    -o variants_for_gnomad_annot.txt

#Unused commands
#--pick
