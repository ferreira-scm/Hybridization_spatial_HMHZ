# Hybridization_spatial_HMHZ
Eco-evolutionary dynamics of host-microbiome interactions in a natural population of closely related mice subspecies and their hybrids
doi: https://doi.org/10.1101/2023.12.11.571054 

## Project description
Closely related host species share similar symbionts, but the effects of host genetic admixture and environmental conditions on these communities remain largely unknown. Here, we investigated the influence of genetic admixture on the intestinal prokaryotic and eukaryotic communities (fungi and parasites) of two house mouse subspecies (Mus musculus domesticus and M. m. musculus) and their hybrids in two settings: (1) wild-caught mice from the European hybrid zone and (2) wild-derived inbred mice in a controlled laboratory environment before and during a community perturbation (infection).

## Repository structure

data/ --> data necessary to run these scripts; these are a list of phyloseq objects for the wild dataset and lab dataset. Each phyloseq object corresponds to OTU table, taxonomic table and metadata for each amplicon, after identification and quality screening of ASVs and taxonomic annotation. Data description can be found here: https://doi.org/10.1186/s13071-023-05800-6

R/ --> scripts for data analysis

tmp/ --> temporary files created from running the scripts

fig/ --> figures

## Brief description of each script:

### R/1_filtering.R

Needs files data/PhyloSeqList_HMHZ_All.Rds and data/EimeriaSpeciesAssign.RDS

quality filtering ASVs per amplicon and merging.

Removing ASVs with less than 0.005% abundance and removing samples with less than 100 reads

Removing taxonomic handlers from silva (e.g. g__; s__).

Total sum scaling per amplicon (relative abundances)

Merge into one phyloseq object.

Relable and merge Eimeria spp. as in Ferreira et al. 2023, resulting in 3 combined ASVs (cASVs) for E. ferrisi, E. falciformis and E. vermiformis

Small adjustments to metadata.

Correlation co-occurrence networks for all known parasite genera.

Is network modular? Only Oxyurida is not. We do phylogenetic analysis as done for Eimeria in Ferreira et al. 2023.

Readjust genera names

Merge ASVs into cASV that are likely from the same taxa:
Co-occurrence networks with only positive edges of ASV abundace per genus. Merge ASVs that cluster together.

final phyloseq object for further analysis is "tmp/PS.TSS_filtered.rds"

### R/2_dyadicMCMC_full.R

Needs file "tmp/PS.TSS_filtered.rds"

prepare dyadic data and do distance-based model for all the taxa (full model) with brms.

model checks

plotting figure 1

### R/3_dyadicMCMC_groups.R

Needs file "tmp/PS.TSS_filtered.rds"

prepare data and do distance based model for parasites, plants, fungi and bacteria

plotting figure 2 and 3ab

### R/4_Network.R 

Needs file "tmp/PS.TSS_filtered.rds"

Co-occurrence network of Fungi, Parasites and Bacteria

plotting figure 3c

### R/5_lab.R

Needs file "data/PhyloSeqList_All_Tax_New.Rds"

does the same filtering procedure as in R/1_filtering.R

distance based models for all the taxa and for each group (parasites, fungi, plants and bacteria

plotting figure 4


