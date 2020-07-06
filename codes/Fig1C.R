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

# Fig 1C

# Calculate MHI
# Beneficial (Bacteroidia and Clostridia) and dysbiosis-associated (Gammaproteobacteria and Bacilli) bacteria as defined in Blount, OFID 2018
bacteroides = taxonomy[taxonomy$class %in% "Bacteroidia", "OTU"]
clostridia = taxonomy[taxonomy$class %in% "Clostridia", "OTU"]
gammaproteobacteria = taxonomy[taxonomy$class %in% "Gammaproteobacteria", "OTU"]
bacilli = taxonomy[taxonomy$class %in% "Bacilli", "OTU"]

pt_mbiome_df$log_mhi = log10(sapply(rownames(pt_mbiome_df), FUN = function(x){
  beneficial = sum(pt_mbiome_df[rownames(pt_mbiome_df) %in% x, colnames(pt_mbiome_df) %in% c(bacteroides, clostridia)])
  dysbiosis = sum(pt_mbiome_df[rownames(pt_mbiome_df) %in% x, colnames(pt_mbiome_df) %in% c(gammaproteobacteria, bacilli)])
  beneficial/dysbiosis
}))

# Plot the distribution of baseline log_MHI between patients with no acquisition and with acquisition
par(mar = c(5, 6, 5, 5))
boxplot(pt_mbiome_df$log_mhi ~ pt_mbiome_df$binary_outcome,  
        outline = FALSE,
        col = c("white", "pink"), 
        xlab = "", 
        ylab = "Microbiome Health Index", 
        cex.axis = 1.5,
        main = "", 
        cex.main = 1.5, 
        cex.lab = 1.5, 
        names = c("No acquisition", "Acquisition"),
        ylim = c(floor(min(pt_mbiome_df$log_mhi)), ceiling(max(pt_mbiome_df$log_mhi))))

stripchart(pt_mbiome_df$log_mhi ~ pt_mbiome_df$binary_outcome, 
           vertical = TRUE,
           pch=16, 
           col='black', 
           method='jitter', 
           jitter=0.15, 
           cex=1, 
           lwd=0.5, 
           add = T,
           ylim = c(floor(min(pt_mbiome_df$log_mhi)), ceiling(max(pt_mbiome_df$log_mhi))))

abline(h = 0, col = "red", ylim = c(-3,3.5), lwd = 3)
