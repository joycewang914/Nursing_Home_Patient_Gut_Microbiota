# Plot Fig1 showing the association between clinical and microbiome features.

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

# Fig 1B 
top50_otus = colnames(skin_depleted_otu_sample_RA)[1:50]

# Perform logistic regression to assess association between microbiome features (relative abundance of top 50 OTUs) and acquiring an resistant enteric organism.

otu_outcome_association = sapply(top50_otus, FUN = function(x){
  association = glm(pt_mbiome_df$binary_outcome ~ log10(pt_mbiome_df[,x] + 1e-05), family = "binomial")
  coef(summary(association))[2,4]
})

# Only Enterococcus (Otu0006) came out to be significantly assocaited with outcome
par(mar = c(5, 6, 5, 5))
boxplot(log10(as.numeric(as.character(pt_mbiome_df[,x])) + 1e-5) ~ pt_mbiome_df$binary_outcome,  
        outline = FALSE,
        col = c("white", "orange"), 
        xlab = "", 
        ylab = "Enterococcal\nRelative Abundance", 
        cex.axis = 1.5,
        main = "", 
        cex.main = 1.5, 
        cex.lab = 1.5, 
        names = c("No acquisition", "Acquisition"))

stripchart(log10(as.numeric(as.character(pt_mbiome_df[,x])) + 1e-5) ~ pt_mbiome_df$binary_outcome, 
           vertical = TRUE,
           pch=16, 
           col='black', 
           method='jitter', 
           jitter=0.15, 
           cex=1, 
           lwd=0.5, 
           add = T)
