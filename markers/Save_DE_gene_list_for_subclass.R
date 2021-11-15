#https://community.brain-map.org/t/question-on-patchseq-jupyter-example/1030/6

.libPaths("/home/agatab/R/fahimehb-library-copy/3.5")

library(dendextend)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(feather)  
library(parallel)
library(foreach)

#

filtSubclass <- "L6_IT" #"L5_PT" #excitatory, L234_IT, L6_IT, L6b, NULL
dataset <- "Patch-seq"  #FACS
max_num_markers <- 50000


facs.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/Mm_VISp_AIT2.3.0_20047_202005"
ps.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20210818_collapsed40_cpm/"

anno.file <- "anno.feather"
data.file <- "data.feather"

if (dataset == "FACS") {
  data.dir <- facs.dir
  cl.label <- "cluster_label"
} else if (dataset == "Patch-seq") {
  data.dir <- ps.dir  
  cl.label <- "Tree_first_cl_label"
}

# Load sample annotations (anno)
anno <- read_feather(file.path(data.dir, anno.file))

# Load gene counts (data)
data <- read_feather(file.path(data.dir, data.file))

# Match data with anno, and convert to matrix
anno <- anno[match(data$sample_id,anno$sample_id),]
anno <- as.data.frame(anno)

data.mat  <- as.matrix(data[,names(data)!="sample_id"])  
rownames(data.mat) <- anno$sample_id


# Filter anno and data down to only specific clusters, if specified
if (!is.null(filtSubclass)) {
  if (filtSubclass == "L234_IT") {
    subclasses <- c("L2/3 IT","L4")
    filename.code <- "L234"
  } else if (filtSubclass == "excitatory") {
    subclasses <- c("L2/3 IT","L4","L5 IT","L5 PT","L6 IT","L6 CT","L6b") 
    filename.code <- "exc"
  } else if (filtSubclass == "L5_PT") {
    subclasses <- c("L5 PT")
    filename.code <- "L5_PT"
    anno <- transform(anno, special_cluster_label = ifelse(Tree_first_cl_label == 'L5 PT VISp Chrna6', 'L5 PT VISp Chrna6', 'Other L5 PT'))
    anno <- transform(anno, special_cluster_id = ifelse(Tree_first_cl_label == 'L5 PT VISp Chrna6', 1, 2))
    anno <- transform(anno, special_cluster_color = ifelse(Tree_first_cl_label == 'L5 PT VISp Chrna6', Tree_first_cl_color, 'gray'))
  } else if (filtSubclass == "L6_IT") {
    subclasses <- c("L6 IT")
    filename.code <- "L6_IT"
    anno <- transform(anno, special_cluster_label = ifelse(Tree_first_cl_label == 'L6 IT VISp Car3', 'L6 IT VISp Car3', 'Other L6 IT'))
    anno <- transform(anno, special_cluster_id = ifelse(Tree_first_cl_label == 'L6 IT VISp Car3', 1, 2))
    anno <- transform(anno, special_cluster_color = ifelse(Tree_first_cl_label == 'L6 IT VISp Car3', Tree_first_cl_color, 'gray'))
  } else if (filtSubclass == "L6b") {
    subclasses <- c("L6b")
    filename.code <- "L6b"
    anno <- transform(anno, special_cluster_label = ifelse(Tree_first_cl_label == 'L6b Col8a1 Rprm', 'L6b Col8a1 Rprm', 'Other L6b'))
    anno <- transform(anno, special_cluster_id = ifelse(Tree_first_cl_label == 'L6b Col8a1 Rprm', 1, 2))
    anno <- transform(anno, special_cluster_color = ifelse(Tree_first_cl_label == 'L6b Col8a1 Rprm', Tree_first_cl_color, 'gray'))
  }
  subclass.mask <- is.element(anno$subclass_label, subclasses)
  if (dataset == "Patch-seq") {
    quality.mask <- is.element(anno$Tree_call_label, c("Core","I1","I2","I3"))
    final.mask <- subclass.mask & quality.mask
  } else if (dataset == "FACS") {
    final.mask <- subclass.mask
  }
  anno <- anno[final.mask,]
  data.mat <- data.mat[final.mask,]
} else {
  filename.code <- "all"
}
print(paste("Number of cells: ", dim(anno)[1])) 


# Transpose data and apply log2 transform
norm.data <- log2(data.mat + 1)
norm.data  <- t(norm.data)

######################### One vs Rest #############################
# Make a data.frame of unique cluster id, type, color, and broad type
ref.cl.df <- as.data.frame(unique(anno[,c("special_cluster_id", "special_cluster_label", "special_cluster_color")]))
# Standardize cluster annoation with cluster_id, cluster_label and cluster_color. These are the required fields to visualize clusters properly.
colnames(ref.cl.df)[1:3] <- c("cluster_id", "cluster_label", "cluster_color")

# Sort by cluster_id
ref.cl.df <- ref.cl.df[order(ref.cl.df$cluster_id),]
row.names(ref.cl.df) <- ref.cl.df$cluster_id
ref.cl <- setNames(factor(anno$special_cluster_id), anno$sample_id)

# Use the reference clusters to find marker genes
select.markers.onevsrest <- select_markers(norm.data, ref.cl, n.markers=max_num_markers)$markers


# Save data
output.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data"
file.name = paste0("select_markers_", filename.code, "_clusters_onevsrest.csv")
write.csv(select.markers.onevsrest, file = file.path(output.dir, file.name))




########################### Pairwise #############################
# Make a data.frame of unique cluster id, type, color, and broad type
ref.cl.df <- as.data.frame(unique(anno[,c("Tree_first_cl_id", "Tree_first_cl_label", "Tree_first_cl_color")]))
# Standardize cluster annoation with cluster_id, cluster_label and cluster_color. These are the required fields to visualize clusters properly.
colnames(ref.cl.df)[1:3] <- c("cluster_id", "cluster_label", "cluster_color")

# Sort by cluster_id
ref.cl.df <- ref.cl.df[order(ref.cl.df$cluster_id),]
row.names(ref.cl.df) <- ref.cl.df$cluster_id
ref.cl <- setNames(factor(anno$Tree_first_cl_label), anno$sample_id)

# Use the reference clusters to find marker genes
pairwise <- select_markers(norm.data, ref.cl, n.markers=max_num_markers)
select.markers.pairwise <-pairwise$markers
(de.gene.stats.pairwise) <- pairwise$de.genes
#de.gene.stats.pairwise.df <- subset(data.frame(t(sapply(de.gene.stats.pairwise,c))), select = c("up.score","down.score","score","up.num","down.num","num"))
#de.gene.stats.pairwise.df <- data.frame(t(sapply(lapply(de.gene.stats.pairwise, function(x) x[c(-1, -2)]), c)))
# Remove columns with lists of genes
de.gene.stats.pairwise.df <- data.frame(t(sapply(de.gene.stats.pairwise, function(x) x[c(-1, -2)])))
de.gene.stats.pairwise.df[] <- apply(de.gene.stats.pairwise.df,2,as.character)
# Save data
output.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data"
file.name <- paste0("select_markers_", filename.code, "_clusters_pairwise.csv")
write.csv(select.markers.pairwise, file = file.path(output.dir, file.name))

stats.file.name <- paste0("de_gene_stats_", filename.code, "_clusters_pairwise.csv")
write.csv(de.gene.stats.pairwise.df, file = file.path(output.dir, stats.file.name))



# Save data for use in Jupyter Notebook
#write.csv(facs.data.mat, file = file.path(output.dir, "exc_facs_norm_data.csv"))
#write.csv(facs.anno, file = file.path(output.dir, "exc_facs_anno.csv"))
