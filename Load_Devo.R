require("RColorBrewer")
require("gplots")
source("Colour_Scheme.R")


map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Mmus_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

convert_to_integer <- function(mat) {
	mat <- round(as.matrix(mat))
	storage.mode(mat) <- "integer"
	mat = mat[rowSums(mat) > 0,]
	return(mat)
}

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
fdr=0.01
counts_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = FALSE)
norm_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = FALSE)

# Zhong
zhong = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_fpkm_ZHONG.txt", header=TRUE);
zhong = zhong[!duplicated(zhong[,1]),]
rownames(zhong) = zhong[,1]
zhong = zhong[,2:length(zhong[1,])]
zhong = as.matrix(zhong);

zhong_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_labels_ZHONG.txt");
ultralow = which(rowMeans(zhong) < 10^-5)
zhong = zhong[-ultralow,]
zhong_list = normalize_data(zhong, labels=zhong_labels, is.counts=FALSE)
zhong_count = convert_to_integer(zhong_list$data);

# Xue
Xue_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_mat.txt", header=TRUE)
Xue_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_labels.txt", header=FALSE)
Xue_labels = as.character(unlist(Xue_labels))
Xue_labels[Xue_labels == "2cellmixed"] = "2cellmix"
Xue_labels[Xue_labels == "4cellmixed"] = "4cellmix"
Xue_labels[Xue_labels == "8cellmixed"] = "8cellmix"
Xue_labels[Xue_labels == "Pronucleus"] = "pronuc"
Xue_list = normalize_data(Xue_data, labels=Xue_labels, is.counts=FALSE)
Xue_count = convert_to_integer(Xue_list$data);

# Fan
Fan_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE53386_matrix_fpkms.tsv", header=TRUE)
Fan_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE53386_labels.tsv", header=FALSE)
Fan_labels = as.character(unlist(Fan_labels))
Fan_labels[Fan_labels == "a-AM_treated_2-cell"] = "AM_2cell"
Fan_labels[Fan_labels == "2-cell"] = "2cell"
Fan_labels[Fan_labels == "4-cell"] = "4cell"
Fan_labels[Fan_labels == "8-cell"] = "8cell"
Fan_list = normalize_data(Fan_data, labels=Fan_labels, is.counts=FALSE)
Fan_count = convert_to_integer(Fan_list$data)

# Goolam
Goolam_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Goolam_et_al_2015_count_table.tsv", header=T)
Goolam_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Goolam_et_al_2015_count_table_labels.tsv", header=F)
Goolam_list = normalize_data(Goolam_data, unlist(Goolam_labels[1,]), is.counts=TRUE)
Goolam_counts = normalize_data(Goolam_data, unlist(Goolam_labels[1,]), is.counts=FALSE)
Goolam_list$labels = (as.character(Goolam_list$labels))
Goolam_list$labels[Goolam_list$labels =="32cell"] = "blast"
rownames(Goolam_counts$data) = ensg2symbol(rownames(Goolam_counts$data))
rownames(Goolam_list$data) = ensg2symbol(rownames(Goolam_list$data))

