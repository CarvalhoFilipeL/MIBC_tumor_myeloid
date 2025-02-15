###### Bulk RNA-seq analysis mouse tumors -- Figure 5

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(readr)
library(DESeq2)
library(pheatmap)
library(DEGreport) 

library(tximport)
library(AnnotationHub)
library(ensembldb)
library(DOSE)
library(pathview)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdb)
library(msigdbr)
library(fgsea)
library(dplyr)

metadata = read.csv("metadata.csv") 

metadata = read.csv("metadata.csv") %>% 
        mutate(FILE_FOLDER = paste(metadata$sample, "quant.sf", sep=""))

metadata = metadata %>%
        mutate(treatment = factor(treatment, levels = c("Vehicle", "PARP14i", "Cisplatin", "Cisplatin_PARP14i")))



data <- read.table("tx2gene.tsv", header=T, row.names=1) %>%
        rownames_to_column("transcript_id")

# overall (all files)
files <- file.path(metadata$FILE_FOLDER)
# provide names for each element
names(files) <- (metadata$sample)

# tximport
txi = tximport(files = files,
               type = "salmon",
               tx2gene = data[,c(1,2)],
               countsFromAbundance = "lengthScaledTPM")

# Check that the row names of the metadata equal the column names of the **raw counts** data
all(colnames(txi) == rownames(metadata)) # TRUE

dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ treatment)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="treatment")

### Extract the rlog matrix from the object
rld_mat <- assay(rld)    

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

### Plot heatmap
pheatmap(rld_cor) 

## Run analysis
dds <- DESeq(dds)

## Check the size factors
sizeFactors(dds)

#normalize counts
normalized_counts <- counts(dds, normalized=TRUE)

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## Define contrasts, extract results table, and shrink the log2 fold changes
contrast_oe <- c("treatment","Vehicle", "PARP14i",  "Cisplatin", "Cisplatin_PARP14i")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)

plotMA(res_tableOE_unshrunken, ylim=c(-2,2))

resultsNames(dds)

res <- lfcShrink(dds, 
                 coef = "Intercept",
                 type = 'apeglm') 

plotMA(res, ylim=c(-2,2))

summary(res, alpha = 0.05)

# Turn the results object into a data frame
res_df <- data.frame(res)

# Likelihood ratio test 
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)

# Extract results
res_LRT <- results(dds_lrt)

# Subset the LRT results to return genes with padj < 0.05
sig_res_LRT <- res_LRT %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble() %>% 
        dplyr::filter(padj < 0.05)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
        pull(gene)

# Subset results for faster cluster finding 
clustering_sig_genes <- sig_res_LRT %>%
        arrange(padj) %>%
        head(n=1000)

# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

# tried this 
rownames(metadata) <- colnames(cluster_rlog) 

# Reorder metadata rows to match cluster_rlog columns
metadata <- metadata[match(colnames(cluster_rlog), rownames(metadata)), ]

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

threshold <- res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff

##########################################################################################

# Extract normalized counts for Cxcl11 -- anti tumor inflammation -- main figure 5G

counts_data <- plotCounts(dds, gene = "Cxcl11", intgroup = "treatment", returnData = TRUE)

# Create a boxplot
ggplot(counts_data, aes(x = treatment, y = count, fill = treatment)) +
        geom_boxplot() +                        # Boxplot
        geom_jitter(shape = 16, position = position_jitter(0.2), size = 2) +  # Add points
        scale_y_log10() +                       # Use log scale for better visualization
        labs(title = "Normalized Counts for Cxcl11 across all conditions",
             x = "Treatment",
             y = "Normalized Counts") +
        scale_fill_manual(values = c("black", "blue", "red", "green")) +  # Custom colors
        theme( panel.background = element_blank(),     # Remove background
               plot.background = element_blank(),      # Transparent plot background
               axis.line = element_line(color = "black"), # Black axis lines
               text = element_text(family = "Arial", size = 12), # Arial font and size 12
               axis.text = element_text(color = "black"),       # Black axis text
               axis.title = element_text(color = "black"),      # Black axis title
               plot.title = element_text(hjust = 0.5))

# Run Kruskal-Wallis test
kruskal_test <- kruskal.test(count ~ treatment, data = counts_data)

# View the result
print(kruskal_test) # Kruskal-Wallis chi-squared = 8.5836, df = 3, p-value = 0.03537

##########################################################################################

# Extract normalized counts for Tlr7 -- anti tumor inflammation --  main figure 5G

counts_data <- plotCounts(dds, gene = "Tlr7", intgroup = "treatment", returnData = TRUE)

# Create a boxplot
ggplot(counts_data, aes(x = treatment, y = count, fill = treatment)) +
        geom_boxplot() +                        # Boxplot
        geom_jitter(shape = 16, position = position_jitter(0.2), size = 2) +  # Add points
        scale_y_log10() +                       # Use log scale for better visualization
        labs(title = "Normalized Counts for Tlr7 across all conditions",
             x = "Treatment",
             y = "Normalized Counts") +
        scale_fill_manual(values = c("black", "blue", "red", "green")) +  # Custom colors
        theme( panel.background = element_blank(),     # Remove background
               plot.background = element_blank(),      # Transparent plot background
               axis.line = element_line(color = "black"), # Black axis lines
               text = element_text(family = "Arial", size = 12), # Arial font and size 12
               axis.text = element_text(color = "black"),       # Black axis text
               axis.title = element_text(color = "black"),      # Black axis title
               plot.title = element_text(hjust = 0.5))

# Run Kruskal-Wallis test
kruskal_test <- kruskal.test(count ~ treatment, data = counts_data)

# View the result
print(kruskal_test) # Kruskal-Wallis chi-squared = 8.4006, df = 3, p-value = 0.03842

#############################################################################################

### GSEA

# Get MSigDB gene sets for mouse
msigdb_mouse <- msigdbr(species = "Mus musculus")

# View available collections
head(msigdb_mouse)

# Get the Hallmark (H) gene sets for mouse
msigdb_mouse_hallmark <- msigdb_mouse %>% 
        filter(gs_cat == "H")

# Create a named vector of log2FoldChange values, ranked by the values
ranks <- df$log2FoldChange
names(ranks) <- df$gene

# Sort in decreasing order (required for GSEA)
ranks <- sort(ranks, decreasing = TRUE)

############# HALLMARK Signature analysis ###############################################
# Convert MSigDB gene sets into a list format for fgsea
msigdb_list <- split(msigdb_mouse_hallmark$gene_symbol, msigdb_mouse_hallmark$gs_name)

# Run fgsea
gsea_results <- fgsea(pathways = msigdb_list, 
                      stats = ranks, 
                      minSize = 15,  # Minimum size of gene sets
                      maxSize = 500, # Maximum size of gene sets
                      nperm = 1000)  # Number of permutations

# Convert list to a data frame and save as csv for Suppl table
gsea_df <- gsea_results %>% select(-leadingEdge)
write.csv(gsea_df, file = "GSEA_control_vs_cisplatin.csv", row.names = TRUE)

# Filter for significant pathways 
signif_gsea <- gsea_results %>% filter(padj < 0.05)

head(signif_gsea)
head(ranks)

ranks <- setNames(df$log2FoldChange, df$gene_symbol)
ranks <- sort(ranks, decreasing = TRUE)
str(msigdb_list)

# Create a bar plot of the manually picked pathways for Figure 4H
myplot1=ggplot(signif_gsea_picked, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES)) +
        coord_flip() +
        labs(title = "GSEA Results of Selected Pathways enriched or depleted in Cisplatin\n using Vehicle as baseline", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) + 
        theme_minimal() 

myplot1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))

## GSEA cisplatin vs combo Figure 4I
em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2g, pvalueCutoff = 0.1,verbose = FALSE)

dotplot(em,showCategory=20,split=".sign")+facet_grid(~.sign)##

gseaplot2(em, geneSetID = 1, title = em$ID[1],color="red",pvalue_table=F) ## genesetID EMT

p<-gseaplot2(em, geneSetID = c(1,2,3,4,5,9,10,11), pvalue_table = FALSE,
             color = c("#49479F", "#f2ba69", "#d276cc"), ES_geom = "dot") # + 
#(theme_classic()) theme_classic() + 
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p