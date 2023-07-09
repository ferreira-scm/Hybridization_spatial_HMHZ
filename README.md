# Hybridization_spatial_HMHZ
Decomposing the gut community of hybrid house mice identifies host genetics and spatial effects

R/1_filtering.R --> quality filtering ASVs per amplicon and merging.
Removing ASVs with less than0.005% abundance and removing samples with less than 100 reads
Removing taxonomic handlers from silva (e.g. g__; s__).
Total sum scaling per amplicon
Merge into one phyloseq object.
Relable and merge Eimeria spp. as in Ferreira et al. 2023, resulting in 3 combined ASVs (cASVs) for E. ferrisi, E. falciformis and E. vermiformis
Small adjustments to metadata.

Correlation co-occurrence networks for all known parasite genera.
Is network modular? Only Oxyurida. We do phylogenetic analysis as done for Eimeria in Ferreira et al. 2023.
Readjust genera names

Merge ASVs into cASV that are likely from the same taxa.
Co-occurrence networks with only positive edges of ASV abundace per genus. Merge ASVs that cluster together.

R/2_dyadicMCMC_full.R --> prepare data and do distance-based model for all the taxa (full model) with brms.
plotting figures S1, S2, S3 and S4
plotting figure 2

R/3_dyadicMCMC_groups.R --> prepare data and do distance based model for parasites, plants, fungi and bacteria
plotting figure 3
