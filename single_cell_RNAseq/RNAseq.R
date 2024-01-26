# module load r/3.6.0
# srun --pty -t 1-0 --mem 128g -p barc R


library(Seurat)
library(tidyverse)


dataA = Read10X(data.dir="/proj/barc/projects/caolab/Zhou_002/outs/filtered_feature_bc_matrix")

cao_rna <- dataA$"Gene Expression"
cao_adt <- dataA$"Antibody Capture"

dim(cao_rna)
dim(cao_adt)


###############
## demux
cao_joint.bcs <- intersect(colnames(cao_rna), colnames(cao_adt))

# Subset RNA and HTO counts by joint cell barcodes
cao_rna <- cao_rna[, cao_joint.bcs]
cao_adt <- as.matrix(cao_adt[, cao_joint.bcs])

# Confirm that the HTO have the correct names
rownames(cao_adt)

cao.seu <- CreateSeuratObject(counts = cao_rna)


cao.seu[["ADT"]] <- CreateAssayObject(counts = cao_adt)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cao.seu <- NormalizeData(cao.seu, assay = "ADT", normalization.method = "CLR")
cao.seu <- HTODemux(cao.seu, assay = "ADT", positive.quantile = 0.99)

table(cao.seu$ADT_classification.global)
table(cao.seu$ADT_classification)

Idents(cao.seu) <- "ADT_maxID"
png("cao12_ridgeplot9_big.220106.png", width = 1024, height = 1024)
RidgePlot(cao.seu, assay = "ADT", features = rownames(cao.seu[["ADT"]])[1:9], ncol = 3)
dev.off()

# previous stop point 2021.12.22
# that run abandoned after mt/MT error


###############

Assays(cao.seu)

# not code, just recording original separation by tag
SEN-2       SEN-1       RP3-1       RP3-2        P3-3 
287         530         592         842         884         907 
SEN-3       RP3-3        P3-2        P3-1 
995        1023        1069        1415 
############


# restart 2022.01.05
cao.seu <- PercentageFeatureSet(object = cao.seu, pattern = "^MT-", col.name = "percent.mt")

# SEN cells
sen_1.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "SEN-1" & percent.mt < 10)
sen_1.seu = RenameCells(sen_1.seu, add.cell.id = "SEN-1")

sen_2.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "SEN-2" & percent.mt < 10)
sen_2.seu = RenameCells(sen_2.seu, add.cell.id = "SEN-2")

sen_3.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "SEN-3" & percent.mt < 10)
sen_3.seu = RenameCells(sen_3.seu, add.cell.id = "SEN-3")

# P3 cells
p3_1.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "P3-1" & percent.mt < 10)
p3_1.seu = RenameCells(p3_1.seu, add.cell.id = "P3-1")

p3_2.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "P3-2" & percent.mt < 10)
p3_2.seu = RenameCells(p3_2.seu, add.cell.id = "P3-2")

p3_3.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "P3-3" & percent.mt < 10)
p3_3.seu = RenameCells(p3_3.seu, add.cell.id = "P3-3")

# RP3 cells
rp3_1.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "RP3-1" & percent.mt < 10)
rp3_1.seu = RenameCells(rp3_1.seu, add.cell.id = "RP3-1")

rp3_2.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "RP3-2" & percent.mt < 10)
rp3_2.seu = RenameCells(rp3_2.seu, add.cell.id = "RP3-2")

rp3_3.seu <- subset(cao.seu, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & ADT_classification == "RP3-3" & percent.mt < 10)
rp3_3.seu = RenameCells(rp3_3.seu, add.cell.id = "RP3-3")


dim(sen_1.seu)
dim(sen_2.seu)
dim(sen_3.seu)
dim(p3_1.seu)
dim(p3_2.seu)
dim(p3_3.seu)
dim(rp3_1.seu)
dim(rp3_2.seu)
dim(rp3_3.seu)

cao.list = as.list(c(sen_1.seu,sen_2.seu,sen_3.seu,p3_1.seu,p3_2.seu,p3_3.seu,rp3_1.seu,rp3_2.seu,rp3_3.seu))
names(cao.list) = c("SEN-1","SEN-2","SEN-3","P3-1","P3-2","P3-3","RP3-1","RP3-2","RP3-3")


#Normalize all samples
for (i in 1:length(x = cao.list))
{
  cao.list[[i]] <- PercentageFeatureSet(object = cao.list[[i]], pattern = "^MT-", col.name = "percent.mt")
  cao.list[[i]] <- SCTransform(cao.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes=FALSE)
}


save.image("cao_combo_seurat_scTransform_220106.Rdata")

# This was needed otherwise I got an error about something being too large/above limits. This particular value is 5000*1024^2, based off of this page: https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r

options(future.globals.maxSize=5242880000)		



cao.features <- SelectIntegrationFeatures(object.list = cao.list, nfeatures = 5000)
cao.list <- PrepSCTIntegration(object.list = cao.list, anchor.features = cao.features, verbose = TRUE)
cao.anchors <- FindIntegrationAnchors(object.list = cao.list, normalization.method = "SCT", anchor.features = cao.features, verbose = TRUE)
cao.integrated <- IntegrateData(anchorset = cao.anchors, normalization.method = "SCT", verbose = TRUE)


saveRDS(cao.integrated,"cao_combo_seurat_scTransform_integrated_220106.rds")

condition = gsub("_.+","",colnames(cao.integrated))
#rep = rep(1,length(cao.integrated@meta.data$orig.ident))

#Tag replicates 

names = cao.integrated@meta.data$ADT_classification
reps = as.character(gsub("(.+)-(.+)","\\2",names,perl=T))

gt = as.character(gsub("(.+)-(.+)","\\1",names,perl=T))
cao.integrated = AddMetaData(cao.integrated,metadata=reps,col.name="rep")
cao.integrated = AddMetaData(cao.integrated,metadata=gt,col.name="gt")


cao.integrated@meta.data$condition = condition

# is this still needed?
#cao.integrated@meta.data$rep = rep

cao.integrated <- RunPCA(cao.integrated, verbose = FALSE, npcs = 100)
cao.integrated <- RunUMAP(cao.integrated, dims = 1:100)
cao.integrated <- FindNeighbors(cao.integrated, dims = 1:100, verbose = FALSE)
#### set resolution #####
cao.integrated <- FindClusters(cao.integrated, verbose = FALSE, resolution = 0.25, algorithm=2)

#Clustering plots
pdf("cao_cellranger_seurat_scTransform_220210_res0.25_integrated_clustered.pdf",width=12)
DimPlot(cao.integrated,reduction = "umap", label = TRUE)
DimPlot(cao.integrated,reduction = "umap", group.by = "condition")
DimPlot(cao.integrated,reduction = "umap", group.by = "gt")
dev.off()

#Comparison plots
pdf("cao_cellranger_seurat_scTransform_220210_res0.25_integrated_clustered_comparison.pdf",width=20)
DimPlot(cao.integrated,reduction = "umap", label = TRUE, split.by = "condition")
DimPlot(cao.integrated,reduction = "umap", label = TRUE, split.by = "rep")
dev.off()

# QC plots
#Make violin plots
features = c("nCount_RNA","nFeature_RNA","percent.mt")
plotwidth = 50
plotheight = 15
pdf("cao_cellranger_seurat_scTransform_220210_res0.25_integrated_clustered_vlnplot_QC.pdf",width=plotwidth,height=plotheight)
for (i in features) {
  print(VlnPlot(cao.integrated, features = i, split.by = "gt"))
}
dev.off()


#condition <- unique(cao.integrated@meta.data$condition)
condition <- unique(cao.integrated@meta.data$gt)
rep <- unique(cao.integrated@meta.data$rep)

saveRDS(cao.integrated,"cao_combo_seurat_scTransform_integrated_pca_umap_220106.rds")

comp.prop <- vector()
names.prop <- vector()
for (i in 1:length(condition)) {
  #Subset cao.integrated object
  sub <- cao.integrated[,cao.integrated$gt == condition[i]]
  
  #Generate proportion table by replicate
  sub.table <- prop.table(table(Idents(sub), sub$rep), margin = 2)
  
  #Create table to add to comp (code necessary when not all clusters are in subset)
  sub.final <- array(0, dim=c(length(levels(cao.integrated)),length(unique(sub$rep))))
  rownames(sub.final) <- levels(cao.integrated)
  colnames(sub.final) <- c(1:length(unique(sub$rep)))
  sub.final[rownames(sub.final) %in% rownames(sub.table),] <- sub.table
  
  #Add to compilation and colnames vector
  comp.prop <- cbind(comp.prop,sub.final)
  names.prop <- c(names.prop, paste0(condition[i],"_",colnames(sub.table)))
}


colnames(comp.prop) <- names.prop
write.table(comp.prop,"cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_proportions_replicates.txt",col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")


# import using tidy's readr function
tbl = as_tibble(rownames_to_column(as.data.frame(comp.prop),var = "Cluster"))


# turn data matrix into tidy tibble
tibble = tbl %>% gather(key="Sample",value="Proportion",-Cluster) %>% separate(Sample,sep="_(?=[[:digit:]])",into=c("Condition","Replicate"),convert=T)

# boxplot
pdf("cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_proportions_replicates_boxplots.pdf",width=14,height=10)
tibble %>%
  mutate(Cluster = fct_relevel(Cluster, unique(tibble$Cluster))) %>%
  mutate(Genotype = fct_relevel(Condition, "SEN", "P3", "RP3")) %>%
  ggplot(aes(fill = Genotype, x= Genotype, y=Proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.9)) +
  facet_wrap(~Cluster, scales="free") +
  xlab("")
dev.off()

# scale_fill_manual(values=c("#f01000","#F8766D","#405900","#7CAE00")) +

# Save progress
saveRDS(cao.integrated,"cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered.rds")

save.image("cao_combo_seurat_scTransform_postboxplot_220106.Rdata")

#stop 2022.01.05
# had to redo all to this point on 2022.01.06 due to MT/mt mixup

######################
#### presto method
library(presto)

vargenes<- presto::wilcoxauc(cao.integrated, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)


write.table(top_vargenes,"cao_220210_res0.25_integrated_clustered_prestoMarkers.txt",row.names=F,quote=F,sep="\t")

# so its essentially the same code as calling the markers themselves (like you did) but selecting only the top 10 for each cluster instead. then merging that all into one giant union set, and making a dotplot out of those

top_vargenes = top_markers(vargenes, n = 10, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- top_vargenes %>%
    select(-rank) %>% 
    unclass() %>% 
    stack() %>%
    pull(values) %>%
    unique() %>%
    .[!is.na(.)]
pdf("cao_scTransform_220210_res0.25_integrated_clustered_prestoMarkers_dotplot.pdf",width=50,height=10)
DotPlot(cao.integrated,features=rev(all_markers),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()

summary(table(cao.integrated$condition))

saveRDS(cao.integrated,"cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_withMarkers.rds")
save.image("cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_withMarkers.Rdata")

###################
# findallmarkers way
    
markers <- FindAllMarkers(cao.integrated, assay="SCT", slot="scale.data", only.pos=T, logfc.threshold = 0.25)
write.table(markers,"cao_scTransform_220210_res0.25_integrated_clustered_markers_SCT.txt",quote=F,sep="\t",col.names=NA)


markerfiles = Sys.glob("cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_*Markers_RNA_*.txt")
unionmarkers = vector()

for (i in 1:length(markerfiles)){
  markers = read_tsv(markerfiles[i]) %>% 
    select(X1) %>% 
    head(n=5) %>% 
    pull(X1)
  unionmarkers = c(unionmarkers,markers)
}
unionmarkers = unique(unionmarkers)

pdf("cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_dotplot.pdf",width=30,height=12)
DotPlot(cao.integrated,features=unique(unionmarkers),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()
#######################

# Save progress for today





# Now plot dotplot with cao genes of interest

cao.integrated = readRDS("cao_10X_cellranger-seurat-scTransform_220210_res0.25_integrated_clustered_withMarkers.rds")

###


goi = unique(c("MSH2", "MSH3", "MSH6", "MLH1", "PMS2", "MSH4", "MSH5", "MLH3", "PMS1", "HFM1","RAD51", "RAD51B", "RAD51D","POLA1", "POLB", "POLD1", "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POLG", "POLH", "POLI", "POLQ", "POLK", "POLL", "POLM"))
goi2 = paste0("mm10---",goi)

pdf("cao_10X_cellranger-seurat-scTransform_220210_r0.25_integrated_clustered_dotplot_goi.pdf",width=10,height=10)
DotPlot(cao.integrated,features=goi,assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()

pdf("cao_10X_cellranger-seurat-scTransform_220210_r0.25_integrated_vln_cluster1_goib.pdf",width=30,height=15)
VlnPlot(cao.integrated, group.by = "gt", idents = '1', features = goi,assay="SCT", sort=T, slot="scale.data") + NoLegend()
dev.off()

# , stack=T, flip=T


mat_sct = cao.integrated@assays$SCT@scale.data
mat_sct <- t(mat_sct)


vec_clusters = cao.integrated$seurat_clusters

mermat1 = merge(vec_clusters,mat_sct,by="row.names",all.x=TRUE)
head(mermat)[,1:10]
mermat2 = merge(vec_clusters,mat_sct,by="row.names")

write.table(mermat2,"per_cell_data.csv",col.names=T,row.names=F,quote=FALSE,sep=",")

