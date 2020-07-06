# Plot Fig1 showing the association between clinical and microbiome features.

source("/mothur_taxonomy_munging_function.R")
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

# Fig 1A 
top50_otus = colnames(skin_depleted_otu_sample_RA)[1:50]

# Perform linear regression to assess association between high-risk antibiotic exposure and microbiome features (relative abundance of top 50 OTUs)

top_otu_abx_association = sapply(top50_otus, FUN = function(x){print(x)
  temp_formula = paste0("log10(", x, " + 1e-5) ~ relevel(as.factor(abx), ref = 'none')")
  test = lm(temp_formula, data = pt_mbiome_df)
  coef(summary(test))[grep("high", rownames(coef(summary(test)))), 4]
})

# Subset pt_mbiome_df to only patients who either received no antibiotics or high-risk antibiotics within past 30 days

sig_otus = lapply(names(which(top_otu_abx_association < 0.05)), FUN = function(x){
  cbind(ra = log10(pt_mbiome_df[pt_mbiome_df$abx %in% c("none", "high"),x] + 1e-05),
        abx = ifelse(pt_mbiome_df$binary_outcome[pt_mbiome_df$abx %in% c("none", "high")] %in% 1, "yes", "no"),
        otu = rep(taxonomy[taxonomy$OTU %in% x, "genus"], nrow(pt_mbiome_df[pt_mbiome_df$abx %in% c("none", "high"),])),
        class = rep(taxonomy[taxonomy$OTU %in% x, "class"], nrow(pt_mbiome_df[pt_mbiome_df$abx %in% c("none", "high"),])))
})

sig_otu_df = as.data.frame(do.call(rbind, sig_otus))

abx_otu_cols = structure(rep(c("grey", "red"), length(sig_otus)))

at_ind = (1:(4*length(sig_otus)))[1:(4*length(sig_otus)) %in% c(seq(1, 24, by = 4) + 1, seq(1, 24, by = 4))]

# Manually set up the order of OTU presented (based on negative and positive association and taxonomic class)
otu_order = c("Anaerococcus",
              "Peptoniphilus",
              "Prevotella",
              "Varibaculum",
              "Bifidobacterium",
              "Enterococcus")
sig_otu_df$otu3 = factor(sig_otu_df$otu, levels = otu_order)

abx_bp = boxplot(as.numeric(as.character(sig_otu_df$ra)) ~ abx + otu3, 
                 data = sig_otu_df,
                 outline = FALSE,
                 at = at_ind, 
                 col = abx_otu_cols, 
                 xaxs = FALSE, 
                 xaxt = "n", 
                 xlab = "", 
                 ylab = "", 
                 plot = FALSE)

boxplot(as.numeric(as.character(sig_otu_df$ra)) ~ abx + otu3, 
        data = sig_otu_df,
        outline = FALSE,
        at = at_ind, 
        col = abx_otu_cols, 
        xaxs = FALSE, 
        xaxt = "n", 
        xlab = "", 
        ylab = "OTU\nRelative Abundance", 
        cex.axis = 1.5, 
        cex.lab = 1.5,
        las = 1)

for (a in 1:length(at_ind)){print(a)
  temp_split = strsplit(strsplit(abx_bp$names[a], split = ".", fixed = TRUE)[[1]], " ")
  otu = temp_split[[2]][1]
  abx = ifelse(temp_split[[1]][1] %in% "no", "no", "yes")

  stripchart(at = at_ind[a], 
             as.numeric(as.character(sig_otu_df[sig_otu_df$abx %in% abx & sig_otu_df$otu %in% otu,"ra"])),
             vertical = TRUE,
             pch=16, 
             col='black', 
             method='jitter', 
             jitter=0.15, 
             cex=.75, 
             lwd=0.5, 
             add = T)
  }
text(x  = seq(min(at_ind), max(at_ind), by =4)-0.5, y = -6.25, otu_order, srt= 45,  cex = 1.5, xpd = TRUE, col = "black")
