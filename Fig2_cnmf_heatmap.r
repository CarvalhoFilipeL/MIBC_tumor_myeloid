project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
library(ggplot2)
library(reshape2)
#!gsutil cp $bucket/cellranger_output_directory/bladder/harmony/cellphone/final_data_heatmap005.csv .
#df = pd.read_csv('final_data_heatmap005.csv')

system(paste0("gsutil cp -r ", bucket, "/cellranger_output_directory/bladder/harmony/cellphone/final_data_heatmap005.csv ."),intern=TRUE)
df <- read.csv('final_data_heatmap005.csv')
df

str(df$ycoord)
x_coord <- df$xcoord
y_coord <- df$ycoord
value <- df$value
foo <- data.frame(x_coord, y_coord, value)
foo
matrix_data1 <- acast(foo, y_coord ~ x_coord, value.var = "value")
matrix_data1
ggplot(melt(matrix_data1), aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "X Coordinate", y = "Y Coordinate", title = "CellphoneDB result HeatMap")
q<- ggplot(df, aes(x = sample_source_target, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#1E90FF", high = "red") +
  theme_classic() #theme_void() #theme_dark() #



qq<- ggplot(df, aes(x = sample_source_target, y = ligand_receptor)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#1E90FF", high = "red") +
  theme_classic() #theme_void() #theme_dark() #

q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
qq + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df
q<- ggplot(df, aes(x = sample_source_target, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = target, y = ligand_receptor)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#79baec", high = "red") +
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#bdd5e7", high = "red") +
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#bdd5e7", high = "#1034a6") +
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
 #scale_fill_gradient2(low = "#7382df", mid= "#ffsa10", high="red", midpoint = 0.75, na.value = "grey50")+
 #v7 #scale_fill_gradient2(low = "#fe7a15", mid="#6788f0" , high="#000395", midpoint = 0.5, na.value = "grey50")+
 scale_fill_gradient2(low = "#ff7b89", mid="#6f5f90" , high="#000395", midpoint = 1, na.value = "grey50")+

#scale_fill_gradient(low = "#f26ba6", high="#6485ee",na.value = "grey50")+
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
 scale_fill_gradient2(low = "#ff9190", mid="#5e72eb" , high="#000395", midpoint = 1.3, na.value = "grey50")+
 #scale_fill_gradient(low = "#ff9190",  high="#5e72eb")+

  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))






ggplot(df, aes(x = response, y = ligand_receptor)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

#q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


ggplot(df, aes(x = source, y = ligand_receptor)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()


#ref #https://rpubs.com/lumumba99/1026665




library(scran)
library(RColorBrewer)
library(slingshot)
library(monocle)
library(gam)
#library(clusterExperiment)
library(ggplot2)
library(plyr)
library(MAST)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggpubr)
library(ggplot2)


bucket
#system(paste0("gsutil cp -r ", bucket, "/cellranger_output_directory/bladder/harmony/harmonized_myeloid_only_annotated.h5ad ."),intern=TRUE)

#Convert("harmonized_urothelial_only.h5ad", dest = "h5seurat", overwrite = TRUE)
#uro <- LoadH5Seurat("harmonized_urothelial_only.h5seurat",meta.data = FALSE, misc = FALSE)
#uro

#harmonized_macro_responder_only.h5ad
Convert("slingshot_adata.h5ad", dest = "h5seurat", overwrite = TRUE)
adata_ent <- LoadH5Seurat("slingshot_adata.h5seurat",meta.data = FALSE, misc = FALSE)
adata_ent

print(adata_ent[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(adata_ent, dims = 1:15, cells = 500, balanced = TRUE)
adata_ent<-RunUMAP(adata_ent, dims = 1:10)
DimPlot(adata_ent, reduction = "umap")
adata_ent <- FindNeighbors(adata_ent, dims = 1:10)
adata_ent <- FindClusters(adata_ent, resolution = 0.3)
print(head(Idents(adata_ent), 5))
DimPlot(adata_ent, reduction = "umap")


DimPlot(adata_ent)
macro.markers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
macro.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
https://satijalab.org/seurat/archive/v2.4/visualization_vignette
features.plot <- c("SPP1","LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
RidgePlot(object = macro, features = features.plot)
VlnPlot(object = macro, features = features.plot)
DotPlot(object = macro, features = features.plot)
DotPlot(macro, c("VSIR","VSIG4","CD274","PDCD1LG2","LGALS9","SIGLEC10"),assay = "RNA")

features.plot <- c("VSIR","VSIG4","CD274","PDCD1LG2","LGALS9","SIGLEC10")
RidgePlot(object = macro, features = features.plot)
VlnPlot(object = macro, features = features.plot)
DotPlot(object = macro, features = features.plot)
#FeaturePlot(object = macro, features = features.plot, cols.use = c("lightgrey", "blue"))
macro.markers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top20
DoHeatmap(macro, features = top20$gene) + NoLegend()

DoHeatmap(subset(macro , downsample = 100), features = top20$gene, size = 3)+ 
    theme(text = element_text(size = 3))
combined_averages <- AverageExpression(macro, return.seurat = TRUE) 
DoHeatmap(combined_averages, features = top20$gene, label = FALSE ,draw.lines = FALSE)  +  
        scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
head(combined_averages)
DoHeatmap(combined_averages, features = top20$gene, label = FALSE ,draw.lines = FALSE) 
top20$gene
M1M2<-c('IL1B','TNF','CXCL9','CXCL10','CXCL11','FCGR1A','IRF1','HLA-DPB1','CD86','MARCO','IL2RA','CSCR4','IL27RA','CSF1R','CCL2','CCL7','CCL17','CCL18','CCL23','CTSD','FN1','GAS7','HMOX1','PPARG','LIPA','CD209','CLEC7A','F13A1','MAF','MS4A5A')
DoHeatmap(combined_averages, features = M1M2, label = FALSE ,draw.lines = FALSE) 
DoHeatmap(subset(macro, downsample = 100), features = top20$gene, size = 3)+ 
    theme(text = element_text(size = 7))
plot <- DoHeatmap(subset(macro, downsample = 100), features = top20$gene, size = 3)+ 
    theme(text = element_text(size = 3))
#ggplot2::ggsave(filename = "macro_heatmap30.pdf", plot = plot) # can add additional parameters similar to png 
#system(paste0("gsutil cp -r macro_heatmap30.pdf ", bucket, "/figures/"),intern=TRUE)
macro.markers %>%
    group_by(cluster) %>%
    top_n(n = 100, wt = avg_log2FC) -> top100
macro.markers %>%
    group_by(cluster) %>%
    top_n(n = 4, wt = avg_log2FC) -> top4
top4
top100$gene
macro.markers %>% filter(cluster==6) 

c0<-macro.markers %>% filter(cluster==0) #349
c1<-macro.markers %>% filter(cluster==1) #162
c2<-macro.markers %>% filter(cluster==2) #205
c3<-macro.markers %>% filter(cluster==3) #1
c4<-macro.markers %>% filter(cluster==4) #248
c5<-macro.markers %>% filter(cluster==5) #443
c6<-macro.markers %>% filter(cluster==6) #906


write.csv(c0, "R_macro_cluster0_349genes.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster0_349genes.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
write.csv(c1, "R_macro_cluster1_162genes.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster1_162genes.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
write.csv(c2, "R_macro_cluster2_205genes.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster2_205genes.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
write.csv(c3, "R_macro_cluster3_1gene.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster3_1gene.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
write.csv(c4, "R_macro_cluster4_248genes.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster4_248genes.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
write.csv(c5, "R_macro_cluster5_443genes.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster5_443genes.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
write.csv(c6, "R_macro_cluster6_906genes.csv",row.names=FALSE)
system(paste0("gsutil cp -r R_macro_cluster6_906genes.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)


top100 %>% filter(cluster==1) %>% select("gene")
c0<- top100 %>% filter(cluster==0) #222 
c1<- top100 %>% filter(cluster==1) #13
c2<- top100 %>% filter(cluster==2) #212
c3<- top100 %>% filter(cluster==3) #234
c4<- top100 %>% filter(cluster==4) #33
c5<- top100 %>% filter(cluster==5) #686
c6<- top100 %>% filter(cluster==6) #368
/cellranger_output_directory/bladder/signature
#write.csv(c0, "R_macro_cluster0_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster0_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
#write.csv(c1, "R_macro_cluster1_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster1_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
#write.csv(c2, "R_macro_cluster2_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster2_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
#write.csv(c3, "R_macro_cluster3_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster3_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
#write.csv(c4, "R_macro_cluster4_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster4_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
#write.csv(c5, "R_macro_cluster5_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster5_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)
#write.csv(c6, "R_macro_cluster6_top100.csv",row.names=FALSE)
#system(paste0("gsutil cp -r R_macro_cluster6_top100.csv ", bucket, "/cellranger_output_directory/bladder/signature/"),intern=TRUE)







Idents(macro) <- "celltype"
macro_foravg <- subset(macro, idents = "Putative_Tumor")

Idents(macro_foravg) <- factor(macro_foravg$patient, levels = c("p915","p906","p913"))

DefaultAssay(rcc.integrated_foravg) <- "RNA"
all.genes <- row.names(rcc.integrated_foravg)
rcc.integrated_foravg <- ScaleData(rcc.integrated_foravg, features = all.genes)

cluster.averages <- AverageExpression(rcc.integrated_foravg, return.seurat = TRUE)

DoHeatmap(cluster.averages, features = fgseaRes_response_tum
pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE")]], size = 3, 
    draw.lines = FALSE)
macro
macro@assays$RNA
macro@assays$RNA
macro@assays
colnames(macro)
length(colnames(macro)) #4373
macro@meta.data
#devtools::install_github("immunogenomics/presto")
## change to response
#macro@meta.data$response<- rownames(macro@meta.data)
#macro@meta.data
#nr <- c("-0","-1","-2","-3","-6","-7","-8","-9")
#r <- c("-4","-5","-10","-11","-12","-13")
#for (i in nr) {macro@meta.data[endsWith(macro@meta.data$response,(i)),]$response <- 1}
#for (i in r) {macro@meta.data[endsWith(macro@meta.data$response,(i)),]$response <- 0}
#macro@meta.data


library(presto)
macro.genes <- wilcoxauc(macro, 'seurat_clusters')
#macro.genes <- wilcoxauc(macro, 'response')
head(macro.genes)



https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
dplyr::count(macro.genes, group)
macro
#BiocManager::install("msigdbr")
#BiocManager::install("fgsea")

#BiocManager::install("msigdbr")
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)

msigdbr_show_species()
m_df<- msigdbr(species = "Homo sapiens", category = "H")

#head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
macro.genes
macro.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
allcluster.genes<- macro.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

library(tibble)
ranks<- deframe(allcluster.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
ggplot(fgseaResTidy %>% head(n=5), aes(reorder(pathway, NES), NES)) +
   geom_bar(stat = "identity", aes(fill=padj<0.05), position = "dodge", width = 0.8) +
      labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +theme_minimal()  +
    theme(axis.text.x = element_text(angle = 90),aspect.ratio = 3) 

#ggsave("MacroGSEA_group0_top5_v.pdf")
#system(paste0("gsutil cp -r MacroGSEA_group0_top5_v.pdf ", bucket, "/figures/"),intern=TRUE)


ggplot(fgseaResTidy %>% head(n=5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
 theme_minimal()  +
theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),aspect.ratio = 0.3) 
#ggsave("MacroGSEA_group0_top5.pdf")
#system(paste0("gsutil cp -r MacroGSEA_group0_top5.pdf ", bucket, "/figures/"),intern=TRUE)
ggplot(fgseaResTidy , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

#ggsave("MacroGSEA_group0_top5_v.pdf")
#system(paste0("gsutil cp -r MacroGSEA_group0_top5_v.pdf ", bucket, "/figures/"),intern=TRUE)


ggplot(fgseaResTidy %>% head(n=10), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
 theme_minimal()  +
theme(axis.text.x = element_text(angle = 90, vjust = 2, hjust=1),aspect.ratio = 0.5) 
#ggsave("MacroGSEA_group0_top5.pdf")
#system(paste0("gsutil cp -r MacroGSEA_group0_top5.pdf ", bucket, "/figures/"),intern=TRUE)
ggplot(fgseaResTidy %>% head(n=5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()+  theme(axis.text = element_text(size = 20))   

ggplot(fgseaResTidy %>% head(n=5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()+  theme(axis.text = element_text(size = 15)) +theme(aspect.ratio = 8)  


ggplot(fgseaResTidy %>% head(n=5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + theme_minimal()
ggsave("MacroGSEA_group0_top5.pdf")
system(paste0("gsutil cp -r MacroGSEA_group0_top5.pdf ", bucket, "/figures/"),intern=TRUE)


ggplot(fgseaResTidy , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# only plot the top 5 pathways

ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n=5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


ggplot(fgseaResTidy %>% head(n=5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


plotEnrichment(fgsea_sets[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_ALPHA_RESPONSE")

macro.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster0.genes<- macro.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

library(tibble)
ranks<- deframe(cluster0.genes)

ranks

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


# only plot the top 20 pathways

ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 40), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
plotEnrichment(fgsea_sets[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_ALPHA_RESPONSE")


plotEnrichment(fgsea_sets[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_GAMMA_RESPONSE")


plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")

macro.genes %>%
  dplyr::filter(group == "1") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))




fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj)
# cluster 1

cluster1.genes<- macro.genes %>%
  dplyr::filter(group == "1") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

library(tibble)
ranks<- deframe(cluster1.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# only plot the top 20 pathways

ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 40), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(fgsea_sets[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_ALPHA_RESPONSE")


plotEnrichment(fgsea_sets[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_GAMMA_RESPONSE")


plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")





# cluster 2

cluster2.genes<- macro.genes %>%
  dplyr::filter(group == "2") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

library(tibble)
ranks<- deframe(cluster2.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()



ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 40), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(fgsea_sets[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_ALPHA_RESPONSE")




