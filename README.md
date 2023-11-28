# Hybridization_spatial_HMHZ
Subspecies divergence, hybridisation and the spatial environment shape phylosymbiosis in the microbiome of house mice


data/ --> data necessary to run these scripts; these are a list of phyloseq objects for the wild dataset and lab dataset. Each phyloseq object corresponds to OTU table, taxonomic table and metadata for each amplicon, after identification and quality screening of ASVs and taxonomic annotation. Scripts are found in the repository https://github.com/ferreira-scm/Eimeria_AmpSeq. Methods are described in the "Amplicon sequencing allows differential quantification of closely related parasite species: an example from rodent coccidia (Eimeria)" manuscript https://doi.org/10.1186/s13071-023-05800-6

R/ --> scripts

tmp/ --> temporary files created from running the scripts

fig/ --> figures


Briefly:

R/1_filtering.R --> Needs files data/PhyloSeqList_HMHZ_All.Rds and data/EimeriaSpeciesAssign.RDS

quality filtering ASVs per amplicon and merging.

Removing ASVs with less than 0.005% abundance and removing samples with less than 100 reads

Removing taxonomic handlers from silva (e.g. g__; s__).

Total sum scaling per amplicon (relative abundances)

Merge into one phyloseq object.

Relable and merge Eimeria spp. as in Ferreira et al. 2023, resulting in 3 combined ASVs (cASVs) for E. ferrisi, E. falciformis and E. vermiformis

Small adjustments to metadata.

Correlation co-occurrence networks for all known parasite genera.

Is network modular? Only Oxyurida. We do phylogenetic analysis as done for Eimeria in Ferreira et al. 2023.

Readjust genera names

Merge ASVs into cASV that are likely from the same taxa:
Co-occurrence networks with only positive edges of ASV abundace per genus. Merge ASVs that cluster together.

final phyloseq object for further analysis is "tmp/PS.TSS_filtered.rds"

R/2_dyadicMCMC_full.R --> Needs file "tmp/PS.TSS_filtered.rds"

prepare dyadic data and do distance-based model for all the taxa (full model) with brms.

model checks

plotting figure 1


R/3_dyadicMCMC_groups.R --> Needs file "tmp/PS.TSS_filtered.rds"

prepare data and do distance based model for parasites, plants, fungi and bacteria

plotting figure 2 and 3 and 4ab

R/4_Network.R --> Needs file "tmp/PS.TSS_filtered.rds"

Co-occurrence network of Fungi, Parasites and Bacteria

R/5_lab.R --> Needs file "data/PhyloSeqList_All_Tax_New.Rds"
