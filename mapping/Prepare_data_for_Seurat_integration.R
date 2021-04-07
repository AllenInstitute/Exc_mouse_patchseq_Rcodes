######################################################################################################
### Author: Agata Budzillo ###########################################################################
### Date: 4/07/2021        ###########################################################################
######################################################################################################

# This code is modified from Jeremy Miller's Code2_prepare_data_mapping notebook in map_my_data, 
# and prepares Retro-seq data from the Callaway lab (Kim et al, 2021) and FACS or Patch-seq data
# for integration with Seurat 3.0

knitr::opts_chunk$set(echo = TRUE)
## Load libraries
.libPaths("/home/agatab/R/fahimehb-library-copy/3.5")
library(feather)        # For reading in data
library(scrattch.hicat) # For logCPM function
library("readxl")       # For reading Callaway xlsx annotation


#################   Constants   #######################
# DETERMINE WHICH REFERENCE TO USE!
referenceName = "FACS" #Patch-Seq
filtReference <- c("L2/3 IT","L5 IT","L6 IT")
filtTarget    <- c("L23_AL", "L23_PM")
output_dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/derived_data/"
output_file_name <- "input_data_L23_Callaway_allIT_FACS.RData"

#################   Directories   #######################
## Folder location for reference feather files (default to FACS)
if (referenceName == "FACS"){
  referenceFolder <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/Mm_VISp_AIT2.3.0_20047_202005/"
} else if (referenceName == "Patch-Seq"){
  referenceFolder <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"
} else {
  referenceFolder <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/Mm_VISp_AIT2.3.0_20047_202005/"
}
# Folder locations for target data files
targetFolder <-"/allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data/" #(Location of Callaway data)
targetSuppFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data/Kim2020_supplement/Tables-S4-10/"
callawayDataFile <- "GSE133230_rawCount_945filtered_cells.tsv"  ##"GSE133230_rawCount_945filtered_cells.tsv"
callawayAnnoFile <- "TableS4-Single nuclei sequencing sample metadata.xlsx"

# Path to Reference Dendrogram
dendReference <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/dend.RData"


#################   Format Reference data/metadata   #######################
annoReference <- read_feather(paste(referenceFolder,"anno.feather",sep=""))
exprReference <- feather(paste(referenceFolder,"data.feather",sep=""))
if (referenceName == "Patch-Seq"){
  # If using Patch-seq as reference, remember to exclude FACS cells from Patch-Seq annotation file
  nonFACSReference <- annoReference$collection_label!="FACS"
  exprReference     <- exprReference[nonFACSReference,]
  annoReference    <- annoReference[nonFACSReference,]
}

annoReference <- annoReference[match(exprReference$sample_id,annoReference$sample_id),] 
# Note: This assumes the sample names are stored in a column called "sample_id".  If not, adjust to correct name.
datReference  <- as.matrix(exprReference[,names(exprReference)!="sample_id"])  
rownames(datReference) <- annoReference$sample_id
datReference  <- t(datReference)

## Log2 cpm to normalize the data
datReference  <- logCPM(datReference)



#################   Format Target data/metadata   #######################
annoTarget <- read_excel(paste0(targetSuppFolder,callawayAnnoFile))
annoTarget <- as.data.frame(annoTarget)
rownames(annoTarget) <- annoTarget$Sample

### Read in and format the target data set if original tsv
datTarget <- read.table(file=paste0(targetFolder,callawayDataFile), sep = '\t', header = TRUE) #  (if tsv)
genesTarget <- datTarget$Symbol
cellsTarget <- setdiff(colnames(datTarget), c("GeneID", "Symbol"))

## Log2 cpm to normalize the data
datTarget <- logCPM(as.matrix(datTarget[,cellsTarget]))
rownames(datTarget) <- genesTarget



#################   Match genes across data sets   #######################
#Change all gene names to upper case if original .tsv
rownames(datReference) <- toupper(rownames(datReference))
rownames(datTarget) <- toupper(rownames(datTarget))

# Get common genes for the analysis
kpGenes      <- intersect(rownames(datTarget),rownames(datReference))
print(paste("Number of common genes: ", length(kpGenes))) 

datReference <- datReference[kpGenes,]
datTarget    <- datTarget[kpGenes,]


#################   Narrow down to cells of interest   #######################
kpReference   <- is.element(annoReference$subclass_label,filtReference)
datReference  <- datReference[,kpReference]
annoReference <- annoReference[kpReference,]
print(paste("Number of ", paste(filtReference, collapse=", "), " cells in Reference: ",sum(kpReference)))

kpTarget      <- is.element(annoTarget$FinalType, filtTarget)
datTarget     <- datTarget[,kpTarget]
annoTarget    <- annoTarget[kpTarget,]
print(paste("Number of ", paste(filtTarget, collapse=", "), " cells in Target: ",sum(kpTarget)))


#################   Exclude clusters not present in the Reference dendrogram   #######################
dend <- readRDS(dendReference)
kpReference   <- is.element(annoReference$cluster_label,labels(dend))
datReference  <- datReference[,kpReference]
annoReference <- annoReference[kpReference,]


#################   Save the results   #######################
save(datReference, annoReference, datTarget, annoTarget, file=file.path(output_dir, output_file_name))

sessionInfo()