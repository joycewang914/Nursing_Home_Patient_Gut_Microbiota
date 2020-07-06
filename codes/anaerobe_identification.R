# Anaerobe depletion
## Maually curated a list of common obligate anaerobes found in 
## "Chapter 21: Infections Caused by Anaerobic Bacteria of
## Jawetz, Melnick, & Adelberg's Medical Microbiology, 28e"   

## "Gram-Positive Anaerobic Cocci in Clinical Microbiology Reviews" in 1998  

## "The first 1000 cultured species of the human gastrointestinal microbiota" in FEMS Microbiology Reviews in 2014  

## "Intestinal Blautia is associated with reduced death from graft-versus-host disease" in Biol Blood Marrow Transplant. 2015   


source("mothur_taxonomy_munging_function.R")
source("remove_otu_function.R")

# Read in taxonomy file
tax_table = read.table("NH_microbiota.taxonomy", header = T)
taxonomy = mothur_tax_munging_func(tax_table)


# Read in 16S reads
OTU_counts = read.table("NH_microbiota.opti_mcc.shared", header = TRUE)
rownames(OTU_counts) = OTU_counts$Group
otu_samples = OTU_counts[,grep("^Otu", colnames(OTU_counts))]

# Remove skin-associated OTUs from 16S data
skin_genus = c("Staphylococcus", "Corynebacterium", "Propionibacterium")
skin_depleted_otu_samples = remove_otu(otu_mat = otu_samples, taxonomy = taxonomy, remove_genus = skin_genus)

# Calculate relative abundnace of each OTU
skin_depleted_otu_sample_RA = skin_depleted_otu_samples/rowSums(skin_depleted_otu_samples)

# Load patient data
temp_df = read.csv("patient_data.csv", header = T, row.names = 1)

# Append microbiome data to patient data
pt_mbiome_df = cbind(temp_df, skin_depleted_otu_sample_RA[rownames(temp_df), ])

# Define anaerobe OTUs
anaerobe = c('Clostridium', 'Bifidobacterium', 'Eubacterium', 'Lactobacillus', 'Anaerococcus', 'Finegoldia', 'Peptoniphilus', 'Porphyromonas', 'Prevotella', 'Bacteroides', 'Fusobacterium', 'Peptococcus', 'Peptostreptococcus', 'Veillonella', 'Odoribacter', 'Parabacteroides', 'Paraprevotella', 'Alistipes', 'Blautia', 'Collinsella', 'Ruminococcus', 'Coprococcus', 'Sarcina', 'Atopobium')

# anaerobe otus
anaerobe_otus = unlist(sapply(anaerobe, FUN = function(x){taxonomy[grep(x, taxonomy$genus, ignore.case = T), "OTU"]}))

anaerobe_ra = rowSums(skin_depleted_otu_samples[,colnames(skin_depleted_otu_samples) %in% anaerobe_otus])/rowSums(skin_depleted_otu_samples) * 100

pt_mbiome_df$anaerobe_ra = anaerobe_ra[rownames(pt_mbiome_df)]

anaerobe_mod = glm(binary_outcome ~ anaerobe_ra, data = pt_mbiome_df, family = "binomial")

coef(summary(anaerobe_mod)) # results shows for every 1% decrease in anaerobe relative abundance, there is ~1% decrease risk of acquiring a resistant organism (p = 0.10).
