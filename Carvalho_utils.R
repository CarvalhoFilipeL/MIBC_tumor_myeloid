library(Seurat)
run_integration_v5 <- function(seur_obj,thres = 100,run_cca=F,k.weight=100,dims=30){
    #Remove samples with less than thres cells
    seur_obj = seur_obj[,seur_obj$orig.ident %in% names(table(seur_obj$orig.ident))[table(seur_obj$orig.ident) > thres]]
    #Take care of the edge case that rowSums becomes 0 for any gene
    if(min(rowSums(LayerData(seur_obj,assay="RNA",layer="counts"))) == 0){
        #Remake the assay with the bait and switch method. Hopefully never needed.
        temp.data = LayerData(seur_obj, assay="RNA",layer="counts")
        temp.data = temp.data[rowSums(temp.data) > 3,]
        seur_obj[["RNA2"]] = CreateAssayObject(
            temp.data
            )
        DefaultAssay(seur_obj) = "RNA2"
        seur_obj[["RNA"]] = NULL
        seur_obj[["RNA"]] = seur_obj[["RNA2"]]
        DefaultAssay(seur_obj) = "RNA"
        seur_obj[["RNA2"]] = NULL
        seur_obj$nCount_RNA = colSums(temp.data)
        seur_obj$nFeature_RNA = colSums(temp.data > 0 )
    }
    #Turn to assay 5
    seur_obj[["RNA"]] = as(seur_obj[["RNA"]],Class="Assay5")
    if(length(Layers(seur_obj)) > 1){
        seur_obj = JoinLayers(seur_obj)
    }
    #Split layers
    seur_obj[["RNA"]] = split(seur_obj[["RNA"]],f=seur_obj$orig.ident)
    #Run SCTransform on each layer
    seur_obj = SCTransform(seur_obj,clip.range=c(-10,10),vst.flavor='v2')
    seur_obj = RunPCA(seur_obj)
    seur_obj[["SCT2"]] = as(seur_obj[["SCT"]],Class="Assay5")
    seur_obj[["SCT2"]] = split(seur_obj[["SCT2"]],f=seur_obj$orig.ident)
    #Run integration with RPCA
    seur_obj = IntegrateLayers(seur_obj, assay="SCT2",method="RPCAIntegration",
                              orig.reduction = "pca", new.reduction = "integrated.rpca",
                              features=VariableFeatures(seur_obj[["SCT"]]),
                              k.weight=k.weight,dims=1:dims)
    seur_obj = RunUMAP(seur_obj, reduction="integrated.rpca",dims=1:dims,verbose=F)
    if(run_cca){
        seur_obj = IntegrateLayers(
             object = seur_obj, method = CCAIntegration,
             assay="SCT",dims=1:dims,
             normalization.method="SCT",
             features=VariableFeatures(seur_obj[["SCT"]]),
             orig.reduction = "pca", new.reduction = "integrated.cca",
             verbose = FALSE,k.weight=k.weight
        )
        seur_obj = RunUMAP(seur_obj, reduction="integrated.cca",dims=1:dims,verbose=F)
    }
    seur_obj[["SCT2"]] = NULL
    #DimPlot(seur_obj, group.by="celltype.l1",label=T,label.size=6) + DimPlot(seur_obj, group.by="orig.ident") & NoLegend()
    return(seur_obj)
}
#IMPORTANT: NEW BEHAVIOUR TO COUNTER SINGLETON ISSUES. tryCatch with group.singletons false and removal
#Framework
library(Seurat)
#Metrics
library(scater)
#Packages
library(harmony)
#Bioc packages
library(bluster)
library(BiocNeighbors)
library(BiocParallel)
library(lisi)
#Tidy funcions
library(dplyr)
library(tidyr)
library(magrittr)
#NMI Calculation
library(aricode)
#Plotting
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
#PC calculation
library(findPC)
#https://academic.oup.com/bioinformatics/article/38/10/2949/6565314#394302289 - findPC Automatic PC identification
#https://www.nature.com/articles/s41592-019-0619-0
#https://www.nature.com/articles/s41592-021-01336-8#
#https://theislab.github.io/scib-reproducibility/
#Instead of ASW, I'll implement Chalinski Harabasz - Changed this idea due to cost
#https://stackoverflow.com/questions/56128651/calinski-harabasz-calculation-slower-in-r-clustercrit-than-python-sklearn
#DONE: Implement lisi (iLISI, cLISI), ASW, neighbor purity, ARI, NMI.
#Skipping ARI cause it is unneeded
#DONE: Modularize code
#https://www.nature.com/articles/s41592-023-01933-9 sc
#TODO: add to Seurat object as metadata
run_lisi <- function(seur_obj){

  lisi.res.rna <- compute_lisi(seur_obj@reductions$RNA_PCA@cell.embeddings, seur_obj@meta.data[,c("Batch","RNA_clusters")],
                           c("Batch","RNA_clusters"))

  p1.lisi <- seur_obj@reductions$RNA_UMAP@cell.embeddings %>% 
    cbind(lisi.res.rna) %>% 
    sample_frac(1L, FALSE) %>% 
    gather(key, lisi_value, "Batch", "RNA_clusters") %>% 
    ggplot(aes(RNAUMAP_1, RNAUMAP_2, color = lisi_value)) +
    geom_point() + theme_classic(base_size=25) +
    facet_wrap(~key) 
  p1.lisi
  p1.lisi.box <- ggboxplot(reshape2::melt(lisi.res.rna),x="variable",y="value") + theme_classic(base_size=25) +ylab("LISI")
  p1.lisi.box

  lisi.res.sct <- compute_lisi(seur_obj@reductions$SCT_PCA@cell.embeddings, seur_obj@meta.data[,c("Batch","SCT_clusters")],
                           c("Batch","SCT_clusters"))

  p2.lisi <- seur_obj@reductions$SCT_UMAP@cell.embeddings %>% 
    cbind(lisi.res.sct) %>% 
    sample_frac(1L, FALSE) %>% 
    gather(key, lisi_value, "Batch", "SCT_clusters") %>% 
    ggplot(aes(SCTUMAP_1, SCTUMAP_2, color = lisi_value)) +
    geom_point() + theme_classic(base_size=25) +
    facet_wrap(~key) 
  p2.lisi
  p2.lisi.box <- ggboxplot(reshape2::melt(lisi.res.sct),x="variable",y="value") + theme_classic(base_size=25) +ylab("LISI")
  p2.lisi.box

  lisi.res.harmony <- compute_lisi(seur_obj@reductions$harmony@cell.embeddings, seur_obj@meta.data[,c("Batch","harmony_clusters")],
                           c("Batch","harmony_clusters"))

  p3.lisi <- seur_obj@reductions$harmony_umap@cell.embeddings %>% 
    cbind(lisi.res.harmony) %>% 
    sample_frac(1L, FALSE) %>% 
    gather(key, lisi_value, "Batch", "harmony_clusters") %>% 
    ggplot(aes(harmonyumap_1, harmonyumap_2, color = lisi_value)) +
    geom_point() + theme_classic(base_size=25) +
    facet_wrap(~key) 
  p3.lisi
  p3.lisi.box <- ggboxplot(reshape2::melt(lisi.res.harmony),x="variable",y="value") + theme_classic(base_size=25) +ylab("LISI")
  p3.lisi.box

  lisi.res.integrated <- compute_lisi(seur_obj@reductions$integrated.rpca@cell.embeddings, seur_obj@meta.data[,c("Batch","integrated_clusters")],
                           c("Batch","integrated_clusters"))
  p4.lisi <- seur_obj@reductions$integrated_UMAP@cell.embeddings %>% 
    cbind(lisi.res.integrated) %>% 
    sample_frac(1L, FALSE) %>% 
    gather(key, lisi_value, "Batch", "integrated_clusters") %>% 
    ggplot(aes(integratedUMAP_1, integratedUMAP_2, color = lisi_value)) +
    geom_point() + theme_classic(base_size=25) +
    facet_wrap(~key) 
  p4.lisi
  p4.lisi.box <- ggboxplot(reshape2::melt(lisi.res.integrated),x="variable",y="value") + theme_classic(base_size=25) +ylab("LISI")
  p4.lisi.box

  lisi.res.harmonysct <- compute_lisi(seur_obj@reductions$harmonysct@cell.embeddings, seur_obj@meta.data[,c("Batch","harmonysct_clusters")],
                           c("Batch","harmonysct_clusters"))

  p5.lisi <- seur_obj@reductions$harmonysct_umap@cell.embeddings %>%
    cbind(lisi.res.harmonysct) %>%
    sample_frac(1L, FALSE) %>%
    gather(key, lisi_value, "Batch", "harmonysct_clusters") %>%
    ggplot(aes(harmonysctumap_1, harmonysctumap_2, color = lisi_value)) +
    geom_point() + theme_classic(base_size=25) +
    facet_wrap(~key)
  p5.lisi
  p5.lisi.box <- ggboxplot(reshape2::melt(lisi.res.harmonysct),x="variable",y="value") + theme_classic(base_size=25) +ylab("LISI")
  p5.lisi.box
  if(any(grepl("cca",names(seur_obj@reductions)))){
    lisi.res.ccasct <- compute_lisi(seur_obj@reductions$integrated.cca@cell.embeddings, seur_obj@meta.data[,c("Batch","cca_clusters")],
                           c("Batch","cca_clusters"))

    p6.lisi <- seur_obj@reductions$cca_umap@cell.embeddings %>%
      cbind(lisi.res.ccasct) %>%
      sample_frac(1L, FALSE) %>%
      gather(key, lisi_value, "Batch", "cca_clusters") %>%
      ggplot(aes(ccaumap_1, ccaumap_2, color = lisi_value)) +
      geom_point() + theme_classic(base_size=25) +
      facet_wrap(~key)
    p6.lisi
    p6.lisi.box <- ggboxplot(reshape2::melt(lisi.res.ccasct),x="variable",y="value") + theme_classic(base_size=25) +ylab("LISI")
    p6.lisi.box
  }
}

run_purity <- function(seur_obj,dims){
	#Calculate Neighborhood Purity
	rna.purity = neighborPurity(
	  seur_obj@reductions$RNA_PCA@cell.embeddings[,1:dims],
	  seur_obj$RNA_clusters,
	  k = 50,
	  weighted = TRUE,
	  BNPARAM = KmknnParam(),
	  BPPARAM = MulticoreParam(workers=4)
	)

	sct.purity = neighborPurity(
	  seur_obj@reductions$SCT_PCA@cell.embeddings[,1:dims],
	  seur_obj$SCT_clusters,
	  k = 50,
	  weighted = TRUE,
	  BNPARAM = KmknnParam(),
	  BPPARAM = MulticoreParam(workers=4)
	)

	harmony.purity = neighborPurity(
	  seur_obj@reductions$harmony@cell.embeddings[,1:dims],
	  seur_obj$harmony_clusters,
	  k = 50,
	  weighted = TRUE,
	  BNPARAM = KmknnParam(),
	  BPPARAM = MulticoreParam(workers=4)
	)

        harmonysct.purity = neighborPurity(
          seur_obj@reductions$harmonysct@cell.embeddings[,1:dims],
          seur_obj$harmonysct_clusters,
          k = 50,
          weighted = TRUE,
          BNPARAM = KmknnParam(),
          BPPARAM = MulticoreParam(workers=4)
        )


	integrated.purity = neighborPurity(
	  seur_obj@reductions$integrated.rpca@cell.embeddings[,1:dims],
	  seur_obj$integrated_clusters,
	  k = 50,
	  weighted = TRUE,
	  BNPARAM = KmknnParam(),
	  BPPARAM = MulticoreParam(workers=4)
	)
        if(any(grepl('cca',names(seur_obj@reductions)))){
           cca.purity = neighborPurity(
          seur_obj@reductions$integrated.cca@cell.embeddings[,1:dims],
          seur_obj$cca_clusters,
          k = 50,
          weighted = TRUE,
          BNPARAM = KmknnParam(),
          BPPARAM = MulticoreParam(workers=4)
        )
        }
	ggboxplot(as.data.frame(rna.purity), x="maximum",y="purity") + theme_classic(base_size=30) + 
	ggtitle("RNA Purity") + xlab("Cluster") + ylab("Cluster Purity")
	ggboxplot(as.data.frame(sct.purity), x="maximum",y="purity") + theme_classic(base_size=30) + 
	ggtitle("SCT Purity") + xlab("Cluster") + ylab("Cluster Purity")
	ggboxplot(as.data.frame(harmony.purity), x="maximum",y="purity") + theme_classic(base_size=30) + 
	ggtitle("Harmony Purity") + xlab("Cluster") + ylab("Cluster Purity")
        ggboxplot(as.data.frame(harmonysct.purity), x="maximum",y="purity") + theme_classic(base_size=30) +
        ggtitle("HarmonySCT Purity") + xlab("Cluster") + ylab("Cluster Purity")
	ggboxplot(as.data.frame(integrated.purity), x="maximum",y="purity") + theme_classic(base_size=30) + 
	ggtitle("Integrated Purity") + xlab("Cluster") + ylab("Cluster Purity")
        if(any(grepl("cca",names(seur_obj@reductions)))){
            ggboxplot(as.data.frame(cca.purity),x="maximum",y="purity") + theme_classic(base_size = 30) + 
                ggtitle("CCA Purity") + xlab("Cluster") + ylab("Cluster Purity")
            pheatmap(round(prop.table(table(cca.purity$maximum, seur_obj$cca_clusters),margin=2)*100,2),
                scale="none",color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),
                cluster_rows=F,cluster_cols=F)
        }
	pheatmap(round(prop.table(table(rna.purity$maximum, seur_obj$RNA_clusters),margin=2)*100,2),
         scale="none",color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),
         cluster_rows=F,cluster_cols=F)

	pheatmap(round(prop.table(table(sct.purity$maximum, seur_obj$SCT_clusters),margin=2)*100,2),
		     scale="none",color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),
		     cluster_rows=F,cluster_cols=F)

	pheatmap(round(prop.table(table(harmony.purity$maximum, seur_obj$harmony_clusters),margin=2)*100,2),
		     scale="none",color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),
		     cluster_rows=F,cluster_cols=F)
        
        pheatmap(round(prop.table(table(harmonysct.purity$maximum, seur_obj$harmonysct_clusters),margin=2)*100,2),
                     scale="none",color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),
                     cluster_rows=F,cluster_cols=F)

	pheatmap(round(prop.table(table(integrated.purity$maximum, seur_obj$integrated_clusters),margin=2)*100,2),
		     scale="none",color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),
		     cluster_rows=F,cluster_cols=F)
}

run_silhouette <- function(seur_obj,dims){

	#Calculate Silhouette scores
    silh.rna = approxSilhouette(seur_obj@reductions$RNA_PCA@cell.embeddings[,1:30], seur_obj$RNA_clusters)
	silh.rna.mat = matrix(round(
		 unlist(
		     lapply(
		         split(
		             silh.rna, silh.rna$cluster),
		         function(x) mean(x$width)
		         )
		     ),3))
	rownames(silh.rna.mat) = levels(silh.rna$cluster)
	print("RNA silhouette")
	print(t(silh.rna.mat))

	silh.sct = approxSilhouette(seur_obj@reductions$SCT_PCA@cell.embeddings[,1:30], seur_obj$SCT_clusters)
	silh.sct.mat = matrix(round(
		 unlist(
		     lapply(
		         split(
		             silh.sct, silh.sct$cluster),
		         function(x) mean(x$width)
		         )
		     ),3))
	rownames(silh.sct.mat) = levels(silh.sct$cluster)
	print("SCT silhouette")
	print(t(silh.sct.mat))

	silh.harmony = approxSilhouette(seur_obj@reductions$harmony@cell.embeddings[,1:dims], seur_obj$harmony_clusters)
	silh.harmony.mat = matrix(round(
		 unlist(
		     lapply(
		         split(
		             silh.harmony, silh.harmony$cluster),
		         function(x) mean(x$width)
		         )
		     ),3))
	rownames(silh.harmony.mat) = levels(silh.harmony$cluster)
	print("Harmony silhouette")
	print(t(silh.harmony.mat))

        silh.harmonysct = approxSilhouette(seur_obj@reductions$harmonysct@cell.embeddings[,1:dims], seur_obj$harmonysct_clusters)
        silh.harmonysct.mat = matrix(round(
                 unlist(
                     lapply(
                         split(
                             silh.harmonysct, silh.harmonysct$cluster),
                         function(x) mean(x$width)
                         )
                     ),3))
        rownames(silh.harmonysct.mat) = levels(silh.harmonysct$cluster)
        print("HarmonySCT silhouette")
        print(t(silh.harmonysct.mat))

	silh.integrated = approxSilhouette(seur_obj@reductions$integrated.rpca@cell.embeddings[,1:dims], seur_obj$integrated_clusters)
	silh.integrated.mat = matrix(round(
		 unlist(
		     lapply(
		         split(
		             silh.integrated, silh.integrated$cluster),
		         function(x) mean(x$width)
		         )
		     ),3))
	rownames(silh.integrated.mat) = levels(silh.integrated$cluster)
	print("Integrated silhouette")
	print(t(silh.integrated.mat))
        if(any(grepl("cca",names(seur_obj@reductions)))){
            silh.cca = approxSilhouette(seur_obj@reductions$integrated.cca@cell.embeddings[,1:dims], seur_obj$cca_clusters)
        silh.cca.mat = matrix(round(
                 unlist(
                     lapply(
                         split(
                             silh.cca, silh.cca$cluster),
                         function(x) mean(x$width)
                         )
                     ),3))
        rownames(silh.cca.mat) = levels(silh.cca$cluster)
        print("CCA silhouette")
        print(t(silh.cca.mat))


        }
}
run_nmi <- function(seur_obj){
	rna.nmi = NMI(seur_obj$RNA_clusters, seur_obj$Batch, variant ="sqrt")
	sct.nmi = NMI(seur_obj$SCT_clusters, seur_obj$Batch, variant ="sqrt")
	har.nmi = NMI(seur_obj$harmony_clusters, seur_obj$Batch, variant ="sqrt")
	harsct.nmi = NMI(seur_obj$harmonysct_clusters, seur_obj$Batch, variant ="sqrt")
        int.nmi = NMI(seur_obj$integrated_clusters, seur_obj$Batch, variant ="sqrt")

	print(paste("RNA NMI:",rna.nmi))
	print(paste("SCT NMI:",sct.nmi))
	print(paste("Harmony NMI:",har.nmi))
	print(paste("Harmony SCT NMI:",harsct.nmi))
        print(paste("Integration NMI:",int.nmi))
        if(any(grepl("cca",names(seur_obj@reductions)))){
            cca.nmi = NMI(seur_obj$cca_clusters,seur_obj$Batch,variant="sqrt")
            print(paste("CCA NMI:",cca.nmi))  
        }
}
#LISI is pretty scalable, but Calinski Harabasz can take a ton of time!
run_quick_qc <- function(seur_obj,dims=30,res=0.3,algo=1,do_glmpca=F,do_lisi=T,do_purity=T,do_silhouette=T,do_nmi=T,do_calhara=F,run_kbet=T){

    seur_obj <- NormalizeData(seur_obj,assay="RNA")
    seur_obj <- FindVariableFeatures(seur_obj,assay="RNA")
    seur_obj <- ScaleData(seur_obj,assay="RNA")
    seur_obj = RunPCA(seur_obj,assay="RNA",
                     reduction.name="RNA_PCA")
    seur_obj = RunPCA(seur_obj,assay="SCT",
                     reduction.name="SCT_PCA")
    if(is.null(dims)){
        dims.use = unlist(findPC(seur_obj[["RNA_PCA"]]@stdev,number=30,method='perpendicular line',aggregate=NULL,figure=F))
    }else{
        dims.use = dims
    }
    seur_obj = FindNeighbors(seur_obj,dims=1:dims.use,graph.name = c("RNA_nn","RNA_snn"),
                             reduction = "RNA_PCA")
    seur_obj = RunUMAP(seur_obj, dims=1:dims.use,reduction.name = "RNA_UMAP",verbose=F,
                      reduction="RNA_PCA")
    seur_obj = tryCatch({FindClusters(seur_obj, resolution=res,graph.name="RNA_snn",algorithm=algo)},error=function(cond){FindClusters(seur_obj,resolution=res,graph.name="RNA_snn",algorithm=algo,group.singletons=F)})
    seur_obj$RNA_clusters = seur_obj$seurat_clusters
    #Add Harmony
    seur_obj$orig.ident <- as.factor(seur_obj$orig.ident)
    seur_obj <- RunHarmony(seur_obj,
                           reduction="RNA_PCA",
                           c("orig.ident"),
                           theta=c(2),
                           reduction.save = "harmony",
                           max_iter = 50,
                           dims.use = 1:dims,
                           verbose =F
    )
    if(is.null(dims)){
        dims.use = unlist(findPC(seur_obj[["harmony"]]@stdev,number=30,method='perpendicular line',aggregate=NULL,figure=F))
    }else{
        dims.use = dims
    }
    seur_obj <- RunUMAP(seur_obj, reduction = "harmony",
                    reduction.name = "harmony_umap",
                    dims = 1:dims.use,verbose=F)
    seur_obj <- FindNeighbors(seur_obj,reduction="harmony",dims=1:dims.use,graph.name=c("harmony_nn","harmony_snn"))
    seur_obj = tryCatch({FindClusters(seur_obj, resolution=res,graph.name="harmony_snn",algorithm=algo)},error=function(cond){FindClusters(seur_obj,resolution=res,graph.name="harmony_snn",algorithm=algo,group.singletons=F)})
    seur_obj$harmony_clusters = seur_obj$seurat_clusters
    seur_obj = PrepSCTFindMarkers(seur_obj)
    #Add Harmony based on SCT
    seur_obj <- RunHarmony(seur_obj,
                           reduction="SCT_PCA",
                           c("orig.ident"),
                           theta=c(2),
                           assay="SCT",
                           reduction.save = "harmonysct",
                           max_iter = 50,
                           dims.use = 1:dims,
                           verbose =F
    )
    if(is.null(dims)){
        dims.use = unlist(findPC(seur_obj[["harmonysct"]]@stdev,number=30,method='perpendicular line',aggregate=NULL,figure=F))
    }else{
        dims.use = dims
    }
    seur_obj <- RunUMAP(seur_obj, reduction = "harmonysct",
                    reduction.name = "harmonysct_umap",
                    dims = 1:dims.use,verbose=F)
    seur_obj <- FindNeighbors(seur_obj,reduction="harmonysct",dims=1:dims.use,graph.name=c("harmonysct_nn","harmonysct_snn"))
    seur_obj = tryCatch({FindClusters(seur_obj, resolution=res,graph.name="harmonysct_snn",algorithm=algo)},error=function(cond){FindClusters(seur_obj,resolution=res,graph.name="harmonysct_snn",algorithm=algo,group.singletons=F)})
    seur_obj$harmonysct_clusters = seur_obj$seurat_clusters
    

    if(any(grepl("cca",names(seur_obj@reductions)))){
        if(is.null(dims)){
            #CCA does not return sdevs
            dims.use = dims
        }else{
            dims.use = dims
        }
        seur_obj <- RunUMAP(seur_obj, reduction = "integrated.cca",
                    reduction.name = "cca_umap",
                    dims = 1:dims,verbose=F)
        seur_obj <- FindNeighbors(seur_obj,reduction="integrated.cca",dims=1:dims,graph.name=c("cca_nn","cca_snn"))
        seur_obj = tryCatch({FindClusters(seur_obj, resolution=res,graph.name="cca_snn",algorithm=algo)},error=function(cond){FindClusters(seur_obj,resolution=res,graph.name="cca_snn",algorithm=algo,group.singletons=F)})
        seur_obj$cca_clusters = seur_obj$seurat_clusters
    }
    if(is.null(dims)){
        dims.use = unlist(findPC(seur_obj[["SCT_PCA"]]@stdev,number=30,method='perpendicular line',aggregate=NULL,figure=F))
    }else{
        dims.use = dims
    }

    seur_obj = FindNeighbors(seur_obj,dims=1:dims.use,graph.name = c("SCT_nn","SCT_snn"),
                             reduction = "SCT_PCA")
    seur_obj = RunUMAP(seur_obj, dims=1:dims.use,reduction.name = "SCT_UMAP",verbose=F,
                      reduction = "SCT_PCA")
    seur_obj = tryCatch({FindClusters(seur_obj, resolution=res,graph.name="SCT_snn",algorithm=algo)},error=function(cond){FindClusters(seur_obj,resolution=res,graph.name="SCT_snn",algorithm=algo,group.singletons=F)})
    seur_obj$SCT_clusters = seur_obj$seurat_clusters
    if(is.null(dims)){
            #RPCA does not return sdevs
            dims.use = dims
        }else{
            dims.use = dims
        }
    seur_obj = FindNeighbors(seur_obj,dims=1:dims,reduction = "integrated.rpca",
                            graph.name=c("integrated_nn","integrated_snn"))
    seur_obj = RunUMAP(seur_obj, dims=1:dims,reduction.name = "integrated_UMAP",verbose=F,
                      reduction="integrated.rpca")
    seur_obj = tryCatch({FindClusters(seur_obj, resolution=res,graph.name="integrated_snn",algorithm=algo)},error=function(cond){FindClusters(seur_obj,resolution=res,graph.name="integrated_snn",algorithm=algo,group.singletons=F)})
    seur_obj$integrated_clusters = seur_obj$seurat_clusters
    
    if(!("cca_clusters" %in% colnames(seur_obj@meta.data))){
        temp.meta = seur_obj@meta.data[,c("orig.ident",
                                      "SCT_clusters","integrated_clusters","harmony_clusters","harmonysct_clusters",
                                                               "RNA_clusters","nCount_SCT",
                                                               "nFeature_SCT",
                                                               "nCount_RNA",
                                                               "nFeature_RNA")]

    }else{
        temp.meta = seur_obj@meta.data[,c("orig.ident",
                                      "SCT_clusters","integrated_clusters","harmony_clusters","harmonysct_clusters",
                                                               "RNA_clusters","nCount_SCT",
                                                               "cca_clusters",
                                                               "nFeature_SCT",
                                                               "nCount_RNA",
                                                               "nFeature_RNA")]
        var.expl.cca_pca  = getVarianceExplained(x = t(seur_obj@reductions$integrated.cca@cell.embeddings),temp.meta)
        var.expl.cca_umap = getVarianceExplained(x = t(seur_obj@reductions$cca_umap@cell.embeddings),temp.meta)
        print(head(var.expl.cca_pca))
        print(head(var.expl.cca_umap))
    }
    var.expl.pca = getVarianceExplained(x = t(seur_obj@reductions$RNA_PCA@cell.embeddings),temp.meta)
    var.expl.umap = getVarianceExplained(x = t(seur_obj@reductions$RNA_UMAP@cell.embeddings),temp.meta)
    var.expl.sct_pca  = getVarianceExplained(x = t(seur_obj@reductions$SCT_PCA@cell.embeddings),temp.meta)
    var.expl.sct_umap = getVarianceExplained(x = t(seur_obj@reductions$SCT_UMAP@cell.embeddings),temp.meta)
    var.expl.har_pca  = getVarianceExplained(x = t(seur_obj@reductions$harmony@cell.embeddings),temp.meta)
    var.expl.har_umap = getVarianceExplained(x = t(seur_obj@reductions$harmony_umap@cell.embeddings),temp.meta)
    var.expl.harsct_pca  = getVarianceExplained(x = t(seur_obj@reductions$harmonysct@cell.embeddings),temp.meta)
    var.expl.harsct_umap = getVarianceExplained(x = t(seur_obj@reductions$harmonysct_umap@cell.embeddings),temp.meta)
    var.expl.int_pca  = getVarianceExplained(x = t(seur_obj@reductions$integrated.rpca@cell.embeddings),temp.meta)
    var.expl.int_umap = getVarianceExplained(x = t(seur_obj@reductions$integrated_UMAP@cell.embeddings),temp.meta)
    print(head(var.expl.pca))
    print(head(var.expl.sct_pca))
    print(head(var.expl.har_pca))
    print(head(var.expl.harsct_pca))
    print(head(var.expl.int_pca))
    print(var.expl.umap)
    print(var.expl.sct_umap)
    print(var.expl.har_umap)
    print(var.expl.harsct_umap)
    print(var.expl.int_umap)
    if(is.null(dims)){
        dims.use = dims
    }
    if(do_lisi){
        run_lisi(seur_obj)
    }
	if(do_silhouette){
		run_silhouette(seur_obj,dims=dims.use)
	}
    if(do_purity){
		run_purity(seur_obj,dims=dims.use)
	}
	if(do_nmi){
		run_nmi(seur_obj)
	}
	if(run_kbet){
        if("cca.clusters" %in% colnames(seur_obj@meta.data)){
            embeddings = c("integrated.rpca","integrated.cca","RNA_PCA","SCT_PCA","harmony","harmonysct")
            cluster_names = c("integrated_clusters","cca_clusters","RNA_clusters","SCT_clusters","harmony_clusters","harmonysct_clusters")
        }else{
            embeddings = c("integrated.rpca","RNA_PCA","SCT_PCA","harmony","harmonysct")
            cluster_names = c("integrated_clusters","RNA_clusters","SCT_clusters","harmony_clusters","harmonysct_clusters")
        }
        #I'm parallelizing within scIB, let's not hyperthread this.
        for(i in seq(1,length(embeddings))){
            t0 = Sys.time()
            #Lol this was dumb
            #Misc(seur_obj,paste("kBET_score",embeddings[i],cluster_names[i],sep="_")) =  run_scib(seur_obj,embeddings[i],cluster_names[i],"Sample",dims=dims.use)
            seur_obj = run_scib(seur_obj,embeddings[i],cluster_names[i],"Sample",dims=dims.use)
            t1 = Sys.time()
            print(t1 - t0)
        }
        seur_obj = read_scib(seur_obj)
        #Show whole scIB
        print(Misc(seur_obj)[grepl("kBET",names(Misc(seur_obj)))])
        print(Misc(seur_obj,"scIB")) #Show scIB results for the kBET along with everything else
    }
    if(do_glmpca){
        seur_obj = run_glmpca(seur_obj,res = res)
    }
    return(seur_obj)
}
select_assay <- function(seur_obj, mode){
    if(mode == "SCT"){
       Misc(seur_obj,"cluster.choose") = "SCT_clusters"
       Misc(seur_obj,"reduction.choose") = "SCT_UMAP"
       Misc(seur_obj,"embedding.choose") = "SCT_PCA"
       Misc(seur_obj,"assay.choose") = "SCT" 
    }else if(mode == "RNA"){
       Misc(seur_obj,"cluster.choose") = "RNA_clusters"
       Misc(seur_obj,"reduction.choose") = "RNA_UMAP"
       Misc(seur_obj,"embedding.choose") = "RNA_PCA"
       Misc(seur_obj,"assay.choose") = "RNA" 
	   seur_obj = JoinLayers(seur_obj,assay="RNA")
    }else if(mode == "harmony"){
       Misc(seur_obj,"cluster.choose") = "harmony_clusters"
       Misc(seur_obj,"reduction.choose") = "harmony_umap"
       Misc(seur_obj,"embedding.choose") = "harmony"
       Misc(seur_obj,"assay.choose") = "RNA"
	   seur_obj = JoinLayers(seur_obj,assay="RNA") 
    }else if(mode == "harmonysct"){
       Misc(seur_obj,"cluster.choose") = "harmonysct_clusters"
       Misc(seur_obj,"reduction.choose") = "harmonysct_umap"
       Misc(seur_obj,"embedding.choose") = "harmonysct"
       Misc(seur_obj,"assay.choose") = "SCT" 
    }else if(mode == "integrated"){
        Misc(seur_obj,"cluster.choose") = "integrated_clusters"
        Misc(seur_obj,"reduction.choose") = "integrated_UMAP"
        Misc(seur_obj,"embedding.choose") = "integrated.rpca"
        Misc(seur_obj,"assay.choose") = "SCT" 
    }else if(mode=="cca"){
        Misc(seur_obj,"cluster.choose") = "cca_clusters"
        Misc(seur_obj,"reduction.choose") = "cca_umap"
        Misc(seur_obj,"embedding.choose") = "integrated.cca"
        Misc(seur_obj,"assay.choose") = "SCT" 
    }
	DefaultAssay(seur_obj) = Misc(seur_obj,"assay.choose")
    return(seur_obj)
}
library(presto)
run_presto_compartment_markers <- function(seur_obj,
                                           logfc_thres = 0.5,
                                           pct_in_thres = 75,
                                           pct_diff_thres_genes = 20,
                                           pct_diff_thres_pw = 0,
                                           padj_thres = 0.05
                                          ){
    
    presto.res = wilcoxauc(seur_obj, Misc(seur_obj)$cluster.choose,seurat_assay=Misc(seur_obj)$assay.choose)
    presto.res.filtered = presto.res[presto.res$logFC > logfc_thres 
                            & presto.res$pct_in > pct_in_thres
                            & presto.res$pct_in - presto.res$pct_out > pct_diff_thres_genes
                            & presto.res$padj < padj_thres,
                           ]
    presto.res.split = split(presto.res.filtered, presto.res.filtered$group)
    presto.res.split = presto.res.split[gtools::mixedorder(names(presto.res.split))]
    presto.res.split = lapply(presto.res.split, function(x) x[order(x$logFC,decreasing=T),])
    presto.res.split.compartments = presto.res.split       
    presto.res.compartments = do.call(rbind,presto.res.split.compartments)
    Misc(seur_obj,slot="CompartmentGeneMarkers") <- presto.res.split
    genemarkersfile = paste(Misc(seur_obj)$Analysisname,"_CompartmentGeneMarkers.tsv",sep="")
    write.table(presto.res.compartments,
                genemarkersfile,
                sep="\t",row.names=F,col.names=T,quote=F)
    system(paste(
        "rclone copy ", genemarkersfile," Drive:",
        paste(Misc(seur_obj)$Analysisname,"_results/Compartment_Markers/",sep=""),
        sep=""))
    if("Pathways" %in% names(seur_obj@assays)){
        assay="Pathways"
        presto.res = wilcoxauc(LayerData(seur_obj,assay=assay,layer="scale.data"),
                               unlist(seur_obj[[ Misc(seur_obj)$cluster.choose]]))
        presto.res.filtered = presto.res[presto.res$logFC > logfc_thres
                                & presto.res$pct_in > pct_in_thres
                                & presto.res$pct_in - presto.res$pct_out > pct_diff_thres_pw
                                & presto.res$padj < padj_thres,
                               ]
        presto.res.split = split(presto.res.filtered, presto.res.filtered$group)
        presto.res.split = presto.res.split[gtools::mixedorder(names(presto.res.split))]
        presto.res.split = lapply(presto.res.split, function(x) x[order(x$logFC,decreasing=T),])
        presto.res.split.compartments = presto.res.split       
        presto.res.compartments.pathways = do.call(rbind,presto.res.split.compartments)
        Misc(seur_obj,slot="CompartmentPathwayMarkers") <- presto.res.split
        pathwaymarkersfile = paste(Misc(seur_obj)$Analysisname,"_CompartmentPathwayMarkers.tsv",sep="")
        write.table(presto.res.compartments.pathways,
                    pathwaymarkersfile,
                    sep="\t",row.names=F,col.names=T,quote=F)
        system(paste(
            "rclone copy ", pathwaymarkersfile," Drive:",
            paste(Misc(seur_obj)$Analysisname,"_results/Compartment_Markers/",sep=""),
            sep=""))
        Misc(seur_obj,slot="CompartmentPathwayMarkersLink") <- system(
            paste("rclone link Drive:",
                  paste(Misc(seur_obj)$Analysisname,"_results/Compartment_Markers/",sep=""),
                  pathwaymarkersfile,sep=""))
    } 
    Misc(seur_obj,slot="CompartmentGeneMarkersLink") <- system(
        paste("rclone link Drive:",
              paste(Misc(seur_obj)$Analysisname,"_results/Compartment_Markers/",sep=""),
              genemarkersfile,sep=""))
    Misc(seur_obj,slot="ResultsFolder") <- system(
        paste("rclone link Drive:",paste(Misc(seur_obj)$Analysisname,"_results/Compartment_Markers/",sep=""),
              sep=""))
    return(seur_obj)
}
#Quick and dirty subclustering, usually is enough
run_subclustering_presto <- function(seur_obj, assays_presto=c("SCT","Pathways"),
                               res.choose=0.2,dist.choose=0.03,
                               logfc_thres=0.25, pct_in_thres=50, pct_diff_thres=20,
                               dims.choose=1:20,padj_thres=0.05,correction_func=NULL,algorithm=4){
    #@param: seur_obj - Your object, make sure cell types (coarse) are in the Identity slot
    #@param: embedding.choose - Which embedding to use
    #@param: res.choose - Subclustering resolution to use
    #@param: logfc_threshold - Call something diff exp only when above that
    #@param: pct_in_thres - Call something diff exp only expressed in at least that percent of cells
    #@param: pct_diff_thres - Call something diff exp, only when the difference is that high. The code will not use this is assays is pathway
    #@param: padj_thres - Pretty self explanatory
    #@param: assay_dimrec - Which assay to run dimensional reduction on - You can do multiple. Basically redundant, since you don't recompute embeddings right now
    #@param: assays_presto - Which assay to run differential expression on - Defaults to genes & pathways
    #@param: dims.choose - How many embedding dimensions to use when computing neighbors or UMAP
    #@param: dist.choose - Minimum Distance for UMAP computation
    #Initialize subcluster column
    seur_obj$SubCluster = NA
    assay = Misc(seur_obj,"assay.choose")
    presto.res.list = list()
    seur.list = list()
    embedding.choose = Misc(seur_obj,"embedding.choose")
    cluster.choose = Misc(seur_obj,"cluster.choose")
    temp_seur = DietSeurat(seur_obj,assays = assay,dimreducs = embedding.choose)
    temp_seur@images = list()
    for(i in unique(Idents(temp_seur))){
        print(i)
        seur_obj.subset = temp_seur[,Idents(temp_seur) == i]
        #There is a debate here: Should you recompute the embeddings?
        #I'm leaning towards no, because the embeddings are being computed on all cells
        #and 50 dimensions should be enough to capture heterogeneity
        #If you re-compute the embeddings, you'd have to normalize again
        #Which seems wrong because you've removed some cells from the pool so depth differences
        #change when doing SCT. That said, it's a quick call to SCTransform %>% RunPCA() %>% RunHarmony/
        #I should implement a TODO: with a renormalize flag, which checks the assay variable and
        #checks using the red.choose whether or not you used harmony or CCA and recomputes
        #your embeddings accordingly before the UMAP
        #It's also faster if you don't re-compute. So in terms of practicality, you might want that.
        #For the marker identification, I'm definitely sure you shouldn't re-normalize. The only
        #debate to be had is for the celltype UMAP.
        #That said, since the UMAP is only used strictly for visualization of markers and mixing
        #So long as it's "well-separated" enough for the eye of the beholder, it should be fine.
        #Oh Lior, how I wish I could have your take on this.
        #Neighbors and clusters
        #Adding correction func here. It's mostly here for harmony, it should be a function that takes the same type of seur_obj input and uses similar parameters to what you used.
        #If you add a correction func (currently only works for embedding funcs, data undergoes PCA again and is re-corrected. Do this when you think the correction is extremely necessary and you can't just subset your embeddings. If classes disappear when you subset, it will fail. For consistency, you should match what you did for the whole dataset to what you do for your subsets.
        if(!is.null(correction_func)){
            seur_obj.subset = FindVariableFeatures(seur_obj.subset)
            seur_obj.subset = ScaleData(seur_obj.subset)
            seur_obj.subset = RunPCA(seur_obj.subset)
            seur_obj.subset = correction_func(seur_obj.subset)
        }
        seur_obj.subset = FindNeighbors(seur_obj.subset, reduction=embedding.choose,
                             dims=dims.choose,graph.name=c(paste(embedding.choose,"nn",sep="."),
                                                    paste(embedding.choose,"snn",sep=".")
                                                   ),verbose=F)
        seur_obj.subset = FindClusters(seur_obj.subset,
                                graph.name=paste(embedding.choose,"snn",sep="."),
                                resolution=res.choose,verbose=T,algorithm=algorithm)
        #Add names to make sure we don't accidentally duplicate everything
        seur_obj.subset$seurat_clusters = paste(i,seur_obj.subset$seurat_clusters,sep="_")
        #Add to the same name
        seur_obj.subset[[cluster.choose]] = seur_obj.subset$seurat_clusters
        #Recompute new UMAP
        red.choose = paste(embedding.choose,"umap_opt",sep="_")
        seur_obj.subset = RunUMAP(seur_obj.subset,
                                  reduction=embedding.choose,
                                  dims=dims.choose,
                                  min.dist=dist.choose,
                                  reduction.name=red.choose,
                                  verbose=F)
        for(j in assays_presto){
            if(j == "Pathways"){
                pct_diff_thres.temp=0
            }else{
                pct_diff_thres.temp=pct_diff_thres
            }
            if(j == "Pathways"){
                #In case there's only a single group, add a tryCatch
                presto.res = tryCatch({wilcoxauc(LayerData(seur_obj,
                assay="Pathways",
                layer="scale.data")[,colnames(seur_obj.subset)],
                unlist(seur_obj.subset[[cluster.choose]])
                )},error=function(cond){
                    data.frame(logFC=NA,group=NA,pct_in=NA,pct_out=NA,padj=NA)
                    }
                    )
                
            }else{
                #Same kind of tryCatch here
                presto.res = tryCatch({wilcoxauc(seur_obj.subset, cluster.choose,seurat_assay=j,slot="data")},error=function(cond){
                    data.frame(logFC=NA,group=NA,pct_in=NA,pct_out=NA,padj=NA)
                })
            }
            presto.res.filtered = presto.res[presto.res$logFC > logfc_thres
                                    & presto.res$pct_in > pct_in_thres
                                    & presto.res$pct_in - presto.res$pct_out >= pct_diff_thres.temp
                                    & presto.res$padj < padj_thres,
                                   ]
            presto.res.split = split(presto.res.filtered, presto.res.filtered$group)
            presto.res.split = presto.res.split[gtools::mixedorder(names(presto.res.split))]
            presto.res.split = lapply(presto.res.split, function(x) x[order(x$logFC,decreasing=T),])
            presto.res.list[[paste(i,j,sep="_")]] = presto.res.split
            presto.res.cluster = do.call(rbind, presto.res.split)
            clustermarkersfile = paste(Misc(seur_obj)$Analysisname,"_",i,"_",j,"_ClusterMarkers.tsv",sep="")
            write.table(presto.res.cluster,
                    clustermarkersfile,
                    sep="\t",row.names=F,col.names=T,quote=F)
            system(paste(
                "rclone copy ", clustermarkersfile," Drive:",
                paste(Misc(seur_obj)$Analysisname,"_results/Presto_Subclustering/",sep=""),
                sep=""))
            Misc(seur_obj,slot=paste("Cluster",i,j,"MarkersLink",sep="_")) <- system(
            paste("rclone link Drive:",
                  paste(Misc(seur_obj)$Analysisname,"_results/Presto_Subclustering/",sep=""),
                  clustermarkersfile,sep=""), intern = TRUE)
            seur.list[[i]] = seur_obj.subset
            Misc(seur_obj, paste(i,j,"markers",sep="_")) <- presto.res.split
            #Add new UMAP coordinates as extra metadata
            new.df = seur_obj.subset[[red.choose]]@cell.embeddings
            colnames(new.df) = c(paste(i,"umap","1",sep="_"),paste(i,"umap","2",sep="_"))
            #Currently using the add reduction straight to object method below
            #seur_obj = AddMetaData(seur_obj,new.df)
            #Another way would be the following: Add teh reduction straight to the object
            seur_obj[[paste(i,"umap",sep="_")]] = seur_obj.subset[[red.choose]]
            #Assign the clustering
            seur_obj$SubCluster = seur_obj.subset$seurat_clusters
            
        }
    }
    #This is now done within the loop to avoid issues
    #subcluster_vec = unlist(lapply(seur.list, function(x) unlist(x[[cluster.choose]])))
    #names(subcluster_vec) = NULL
    #seur_obj$SubCluster = subcluster_vec   
    return(seur_obj)
}
