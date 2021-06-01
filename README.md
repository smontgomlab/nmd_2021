# Workflows used for Teran et. al. "Nonsense-mediated decay is highly stable across individuals and tissues"

These scripts constitute the main workflows used for Teran et. al. "Nonsense-mediated decay is highly stable across individuals and tissues": https://doi.org/10.1101/2021.02.03.429654

### Variant annotation - `annotate_variants`

Variant annotation was performed the Variant Effect Predictor (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4, https://github.com/Ensembl/ensembl-vep), including the LOFTEE plugin (https://github.com/konradjk/loftee). Bash scripts used to run VEP are provided, along with an R script used to aggregate the annotations. 
Ensembl 88 was used to obtain all annotations except GnomAD allele frequencies, which where obtained from Ensembl 97.

### Predictive modeling - `predictive_models`

NMD was predicted using penalized logistic regression.  Each subfolder contains a model with different predictors:
`lindeboom` - only the predictors from NMDetective (https://www.nature.com/articles/s41588-019-0517-5).
`full` - the predictors from `lindeboom` in addition to conservation scores, GC content, RNA metrics and allele frequencies.
`conditional` - the same as `full` but only in variants predicted to undergo NMD based on the 50bp rule.
`subject` - the predictors in `full` along with GTEX subject ID.
`tissue` - the predictors in `full` along with tissue name.

The `comparison` subfolder contains the code used to compare the difference in performance across the models.

### Multitissue ASE modeling and cross-tissue ASE correlations - `multitissue_ase`

Multitissue ASE probabilities were modeled using the R code in `ASE_29Dec2014_new.R`.  The other R scripts in the `multitissue_ase` folder were used to call this code to assign ASE probabilities to all ASE measurements.

The `ase_correlation` subfolder contains the code to compute pairwise correlations of ASE across all pairs of tissues.  All variants in all subjects and tissues quantified for ASE with GATK ASEReadCounter and adjusted with WASP (https://github.com/bmvdgeijn/WASP) were reannotated used VEP with Ensembl 88.  Variants in the categories `stop_gained`, `missense_variant`, `synonymous_variant`, `intron_variant`, `3_prime_UTR_variant`, `5_prime_UTR_variant`, and `non_coding_transcript_exon_variant` and files were aggregated by tissue instead of subject.  (Note: this part is very computationally intensive and should probably be run on a server).
