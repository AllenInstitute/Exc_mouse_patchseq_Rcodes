######################################################################################################
### Author: Agata Budzillo ###########################################################################
### Date: 5/11/2021        ###########################################################################
######################################################################################################

# This code is modified from Jeremy Miller's Code2_prepare_data_mapping notebook in map_my_data, 
# and prepares Retro-seq data from the Callaway lab (Kim et al, 2021) and FACS and Patch-seq data
# for integration with Seurat 3.0
#devtools::install_github("AllenInstitute/scrattch.hicat", ref="dev_zy", args = c('--library="/home/agatab/R/agatab-dev-library/"'))
#devtools::install_github("hadley/devtools", args = c('--library="/home/agatab/R/agatab-dev-library/"'))
#install.packages("rlang", lib.loc="/home/agatab/R/agatab-dev-library/")
#library(tidyverse)
.libPaths("/home/agatab/R/fahimehb-library-copy/3.5")

#################   Load Libraries   #######################
suppressPackageStartupMessages({
  library(feather) # For reading in data
  library("readxl")       # For reading Callaway xlsx annotation
  library(dendextend) # get_nodes_attr
  library(Seurat)
  library(scrattch.hicat) # logCPM function, get_cl_prop
  library(mfishtools) # getBetaScore and rfTreeMapping
  library(VENcelltypes)  # For sex and mitochondrial genes.  Can be omitted.
  library(matrixStats) # rowMedians
  library(future)
  library(dplyr)
})

source("//allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/config/config_FACS_Patchseq_Callaway_integration_L56.R")

#################   Format Reference data/metadata   #######################

#################   FACS   #######################
annoFACS <- read_feather(paste(FACSFolder,"anno.feather",sep=""))
exprFACS <- feather(paste(FACSFolder,"data.feather",sep=""))
annoFACS <- annoFACS[match(exprFACS$sample_id,annoFACS$sample_id),] 
datFACS  <- as.matrix(exprFACS[,names(exprFACS)!="sample_id"])  
rownames(datFACS) <- annoFACS$sample_id

# Transpose data and apply log2 transform
datFACS <- log2(datFACS + 1)
datFACS  <- t(datFACS)



#################   Patch-seq   #######################
annoPatchseq <- read_feather(paste(PatchseqFolder,"anno.feather",sep=""))
exprPatchseq <- feather(paste(PatchseqFolder,"data.feather",sep=""))
# Remember to exclude FACS cells from Patch-Seq annotation file AND remove PoorQ cells
nonFACSmask <- annoPatchseq$collection_label!="FACS"
nonPoorQmask <- annoPatchseq$Tree_call_label!= "PoorQ"
PatchseqMask = nonFACSmask & nonPoorQmask
exprPatchseq     <- exprPatchseq[PatchseqMask,]
annoPatchseq    <- annoPatchseq[PatchseqMask,]
if (onlyVISp == TRUE){
  # Exclude cells from nearby regions
  onlyVISpReference <- annoPatchseq$structure_label %in% c("VISp","VISp1","VISp2/3","VISp4","VISp5","VISp6a","VISp6b")
  exprPatchseq     <- exprPatchseq[onlyVISpReference,]
  annoPatchseq    <- annoPatchseq[onlyVISpReference,]
}
annoPatchseq <- annoPatchseq[match(exprPatchseq$sample_id,annoPatchseq$sample_id),] 
datPatchseq  <- as.matrix(exprPatchseq[,names(exprPatchseq)!="sample_id"])  
rownames(datPatchseq) <- annoPatchseq$sample_id

# Transpose data and apply log2 transform
datPatchseq  <- log2(datPatchseq + 1)
datPatchseq  <- t(datPatchseq)


#################   Callaway   #######################
annoCallaway <- read_excel(paste0(callawaySuppFolder,callawayAnnoFile))
annoCallaway <- as.data.frame(annoCallaway)
rownames(annoCallaway) <- annoCallaway$Sample

mappingCallaway <- read.table(file=callawayMappingFile, sep = ',', header = TRUE, row.names="X")
mappingCallaway  <- as.data.frame(mappingCallaway)
# Reindex based on anno
mappingCallaway <- mappingCallaway[match(rownames(mappingCallaway),rownames(annoCallaway)),]

### Read in and format the target data set if original tsv
datCallaway <- read.table(file=paste0(callawayFolder,callawayDataFile), sep = '\t', header = TRUE) #  (if tsv)
cellsCallaway <- setdiff(colnames(datCallaway), c("GeneID", "Symbol"))
allGenesCallaway <- datCallaway$Symbol
datCallaway <- as.matrix(datCallaway[,cellsCallaway])
## Log2 cpm to normalize the data
datCallaway <- logCPM(datCallaway)
rownames(datCallaway) <- allGenesCallaway

genesCallaway <- datCallaway$Symbol
# Align annotation cell IDs with the data
annoCallaway <- annoCallaway[match(colnames(datCallaway),annoCallaway$Sample),] 


#################   Match genes across data sets   #######################
#Change all gene names to upper case
rownames(datFACS) <- toupper(rownames(datFACS))
rownames(datPatchseq) <- toupper(rownames(datPatchseq))
rownames(datCallaway) <- toupper(rownames(datCallaway))

# Get common genes for the analysis
AIGenes      <- intersect(rownames(datFACS),rownames(datPatchseq))
kpGenes      <- intersect(AIGenes,rownames(datCallaway))
print(paste("Number of common genes: ", length(kpGenes))) 

datFACS     <- datFACS[kpGenes,]
datPatchseq <- datPatchseq[kpGenes,]
datCallaway <- datCallaway[kpGenes,]



#################   Narrow down to cells of interest   #######################
kpFACS   <- is.element(annoFACS$subclass_label,filtSubclasses)
datFACS  <- datFACS[,kpFACS]
annoFACS <- annoFACS[kpFACS,]
print(paste("Number of ", paste(filtSubclasses, collapse=", "), " cells in FACS dataset: ",sum(kpFACS)))

# Subsample FACS data if specified
if (subSampleFACS == TRUE) {
  print(paste("Subsampling FACS cells to a maximum of", maxSamples, "cells per t-type."))
  #splitUp <- annoFACS %>%
  #  select(sample_id, cluster_label) %>%
  #  split(f=annoFACS$cluster_label)
  #subSampledSelection <- lapply(splitUp, function(x) {x %>% sample_n(ifelse(nrow(x) < maxSamples, nrow(x), maxSamples))})
  #subSampledSelection <- do.call("rbind", subSampledSelection)
  #maskSampled <- is.element(annoFACS$sample_id,subSampledSelection$sample_id)
  
  # subsample using subsampleCells from mfishtools (found this later)
  maskSampled <- subsampleCells(annoFACS$cluster_label,maxSamples)
  datFACS  <- datFACS[,maskSampled]
  annoFACS <- annoFACS[maskSampled,]
}
  
kpPatchseq   <- is.element(annoPatchseq$subclass_label,filtSubclasses)
datPatchseq  <- datPatchseq[,kpPatchseq]
annoPatchseq <- annoPatchseq[kpPatchseq,]
print(paste("Number of ", paste(filtSubclasses, collapse=", "), " cells in Patch-seq dataset: ",sum(kpPatchseq)))


# Select Callaway cells
if (!is.null(filtCallawayProjection)) {
  ProjectionMask <- is.element(annoCallaway$Projection, filtCallawayProjection)
  if (!is.null(filtCallawayType)) {
    TypeMask <- is.element(annoCallaway$FinalType, filtCallawayType)
    kpCallaway <- ProjectionMask & TypeMask
  } else if (!is.null(filtCallawayMappedTType)) {
    TypeMask <- is.element(mappingCallaway$v1_first_cl, filtCallawayMappedTType) 
    QualityMask <- mappingCallaway$v1_call != "PoorQ"
    kpCallaway <- ProjectionMask & TypeMask & QualityMask
  } else {
    kpCallaway <- ProjectionMask
  }
} else {
  num_cells = dim(annoCallaway)[1]
  kpCallaway <- !logical(num_cells) #Initialize a vector of TRUE
}

datCallaway     <- datCallaway[,kpCallaway]
datids = colnames(datCallaway)
annoCallaway    <- annoCallaway[kpCallaway,]
if (!is.null(filtCallawayType)) {
  callAnnoFileName <- paste0(outputFolder, "Callaway_cells_filtered_by_Callawaytype.csv")
} else if (!is.null(filtCallawayMappedTType)) {
  callAnnoFileName <- paste0(outputFolder, "Callaway_cells_filtered_by_AIMappedtype.csv")
}  else {
  callAnnoFileName <- paste0(outputFolder, "Callaway_cells.csv")
}
write.csv(annoCallaway, file = callAnnoFileName)

annoids = rownames(annoCallaway)
print(paste("Number of ", paste(filtCallawayProjection, collapse=", "), " cells in Callaway dataset: ",sum(kpCallaway)))


#################   Exclude clusters not present in the Reference dendrogram   #######################
dend <- readRDS(dendFile)
clustersUse <- labels(dend) #labels(itdend)
kpFACSclusters   <- is.element(annoFACS$cluster_label,labels(dend))
kpPatchseqclusters   <- is.element(annoPatchseq$Tree_first_cl_label,labels(dend))

datFACS  <- datFACS[,kpFACSclusters]
annoFACS <- annoFACS[kpFACSclusters,]

datPatchseq  <- datPatchseq[,kpPatchseqclusters]
annoPatchseq <- annoPatchseq[kpPatchseqclusters,]

#################   Save the input data here   #######################
save(datFACS, annoFACS, datPatchseq, annoPatchseq, datCallaway, annoCallaway, file=file.path(outputFolder, outputDataFileName))
write.csv(t(datFACS), file = paste0(outputFolder, "datFACS_input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(annoFACS, file = paste0(outputFolder, "annoFACS_input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(t(datPatchseq), file = paste0(outputFolder, "datPatchseq_input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(annoPatchseq, file = paste0(outputFolder, "annoPatchseq_input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(t(datCallaway), file = paste0(outputFolder, "datCallaway_input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(annoCallaway, file = paste0(outputFolder, "annoCallaway_input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))


#################   Create Masks for Gene Clean Up   #######################
# First exclude mitochondrial and sex genes
isExclude <- sort(unique(c(sex_genes,mito_genes)))  
excludeGn <- is.element(toupper(rownames(datCallaway)),isExclude)

# Second, exclude genes with average expression at least four-fold higher in either dataset
#platformGn1 <- abs(rowMeans(datPatchseq)-rowMeans(datCallaway))>=2
platformGn2 <- abs(rowMeans(datFACS)-rowMeans(datCallaway))>=2
platformGn3 <- abs(rowMeans(datFACS)-rowMeans(datPatchseq))>=2

# Finally, only keep genes that are expressed in all datasets (log2(CPM+1)>=1 in 0.1% of the cells)
expressedGn <- (rowSums(datFACS>=1)>(0.001*dim(datFACS)[2])) & (rowSums(datPatchseq>=1)>(0.001*dim(datPatchseq)[2])) & (rowSums(datCallaway>=1)>(0.001*dim(datCallaway)[2]))


#################   Perform Gene Clean Up   #######################
useCuratedGeneList <- TRUE
if (useCuratedGeneList == FALSE){
  keepGenes <- (!(excludeGn|platformGn2|platformGn3)) & (expressedGn) # |platformGn1
  print(paste("Number of genes for integration : ", sum(keepGenes), "(", round(mean(keepGenes),2), " of total)"))
} else {
  curatedGeneList = toupper(read.table(curatedGenesFileName, sep=',', header = TRUE)$Gene)
  keepGenes <- is.element(toupper(rownames(datCallaway)),curatedGeneList) &  (!(excludeGn|platformGn2|platformGn3)) & (expressedGn) # |platformGn1
  print(paste("Number of genes for integration : ", sum(keepGenes), "(", round(mean(keepGenes),2), " of total)"))
} 

if (SeuratIntegration == TRUE) {
  #################   Data and meta-data set-up for Seurat integration     #######################
  
  ## Subset the Reference data to a max of 100 cells per cluster
  #kpSamp   <- subsampleCells(annoFACS$cluster_label,100)
  
  ## Basic data and metadata set-up
  #v1.data     <- cbind(datFACS[keepGenes,kpSamp],datPatchseq[keepGenes,], datCallaway[keepGenes,])  # Include only genes subsetted above
  #v1.metadata <- data.frame(set=c(rep("FACS",sum(kpSamp)),rep("Patch-seq",dim(datPatchseq)[2]),rep("Callaway",dim(datCallaway)[2])),
  #                          celltype = c(annoFACS$cluster_label[kpSamp],annoPatchseq$Tree_first_cl_label,annoCallaway$FinalType))
  v1.data     <- cbind(datFACS[keepGenes,], datPatchseq[keepGenes,], datCallaway[keepGenes,])  # Include only genes subsetted above
  v1.metadata <- data.frame(set=c(rep("FACS",dim(datFACS)[2]),rep("Patch-seq",dim(datPatchseq)[2]),rep("Callaway",dim(datCallaway)[2])),
                            celltype = c(annoFACS$cluster_label,annoPatchseq$Tree_first_cl_label,annoCallaway$FinalType))
  rownames(v1.metadata) = colnames(v1.data)
  
  ## Construct data set lists
  v1      <- CreateSeuratObject(counts = v1.data, meta.data = v1.metadata)
  v1.list <- SplitObject(object = v1, split.by = "set")
  
  #################   Seurat integration     #######################
  # select features that are repeatedly variable across datasets for integration
  # nfeatures defaults to 2000
  #features <- SelectIntegrationFeatures(object.list = v1.list, nfeatures = 2500)
  # Use features determined above (DE genes from tree + AL/PM DE genes)
  features <- rownames(datFACS[keepGenes,])
  
  # Find anchors
  v1.anchors <- FindIntegrationAnchors(object.list = v1.list, anchor.features = features)
  # this command creates an 'integrated' data assay
  v1.combined <- IntegrateData(anchorset = v1.anchors)
  
  save(v1.data, v1.metadata, features, v1.anchors, v1.combined, file=file.path(outputFolder,outputIntegratedDataFileName))
  #################   Seurat Integrated Analysis     #######################
  DefaultAssay(v1.combined) <- "integrated"
  v1.combined.table <- as.matrix(GetAssayData(v1.combined, slot = "data", assay="integrated"))
  write.csv(v1.combined.table, file = file.path(outputFolder, outputCSVFileName))
  
  if (SeuratLabelTransfer == TRUE) {
    # Perform label transfer
    #v1.combined <- ScaleData(v1.combined, verbose = FALSE)
    v1.reference <- v1.list[["FACS"]]
    
    v1.query.Patchseq <- v1.list[["Patch-seq"]]
    v1.labeltransfer.anchors.Patchseq <- FindTransferAnchors(reference = v1.reference, query = v1.query.Patchseq, features = features, dims = 1:30)
    predictions.Patchseq <- TransferData(anchorset = v1.labeltransfer.anchors.Patchseq, refdata = v1.reference$celltype, dims = 1:30)
    v1.query.Patchseq <- AddMetaData(v1.query.Patchseq, metadata = predictions.Patchseq)
    write.csv(predictions.Patchseq, file = paste0(outputFolder, predictionsPatchseqCSVFileName))
    
    v1.query.Callaway <- v1.list[["Callaway"]]
    v1.labeltransfer.anchors.Callaway <- FindTransferAnchors(reference = v1.reference, query = v1.query.Callaway, features = features, dims = 1:30)
    predictions.Callaway <- TransferData(anchorset = v1.labeltransfer.anchors.Callaway, refdata = v1.reference$celltype, dims = 1:30)
    v1.query.Callaway <- AddMetaData(v1.query.Callaway, metadata = predictions.Callaway)
    write.csv(predictions.Callaway, file = paste0(outputFolder, predictionsCallawayCSVFileName))
  } else {
    print("Seurat integration not applied")
  }
}

# Save pre-integration gene values
write.csv(t(datFACS[keepGenes,]), file = paste0(outputFolder, "datFACS_exc_selectgene_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(t(datPatchseq[keepGenes,]), file = paste0(outputFolder, "datPatchseq_exc_selectgene_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(t(datCallaway[keepGenes,]), file = paste0(outputFolder, "datCallaway_exc_selectgene_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))

# Save QC gene values
qc.gene.file <- '//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/10X_analysis/rm.genes.rda'
load(qc.gene.file)
rm.genes = toupper(rm.genes)
keeprm.genes <- is.element(toupper(rownames(datCallaway)),rm.genes)
qcFACS     <- as.data.frame(t(datFACS[keeprm.genes,]))
qcPatchseq <- as.data.frame(t(datPatchseq[keeprm.genes,]))
qcCallaway <- as.data.frame(t(datCallaway[keeprm.genes,]))

qcFACS$nonzero_gene_num <- colSums(datFACS > 0)
qcPatchseq$nonzero_gene_num <- colSums(datPatchseq > 0)
qcCallaway$nonzero_gene_num <- colSums(datCallaway > 0)

write.csv(qcFACS, file = paste0(outputFolder, "datFACS_rmgenes_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(qcPatchseq, file = paste0(outputFolder, "datPatchseq_rmgenes_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))
write.csv(qcCallaway, file = paste0(outputFolder, "datCallaway_rmgenes_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv"))