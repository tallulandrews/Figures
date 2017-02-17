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

getFeatures <- function(counts, norm, fdr=0.01, name="Test"){
#	M3Drop_table = M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=fdr, suppress.plot=TRUE)
	M3Drop_table = M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=2, suppress.plot=TRUE)
	M3Drop_table[,1] = as.character(M3Drop_table[,1])
#	HVG_Deng = BrenneckeGetVariableGenes(norm, fdr=fdr, suppress.plot=TRUE)
	HVG_Deng = BrenneckeGetVariableGenes(norm, fdr=2, suppress.plot=TRUE)
	counts <- as.matrix(counts)
	fit <- NBumiFitModel(counts)
	png(paste(name,"NBFit.png",sep="_"), width=7, height=5, units="in", res=300)
	par(mfrow=c(1,2))
	NBumiCheckFit(counts, fit);
	dev.off();
	Dengfeatures <- NBumiFeatureSelectionCombinedDrop(fit)
	Dengfeatures2 <- NBumiFeatureSelectionHighVar(fit)

	# Features!
#	NB = names(Dengfeatures[p.adjust(Dengfeatures, method="fdr") < fdr]);
#	NBV = names(Dengfeatures2[p.adjust(Dengfeatures2, method="fdr") < fdr]);
	NB = names(Dengfeatures[1:2000]);
	NBV = names(Dengfeatures2[1:2000]);
	M3Drop = M3Drop_table[1:2000,1]
	HVG = HVG_Deng[1:2000,1]
	return(list(NB=NB,M3Drop=M3Drop,HVG=HVG, NBV=NBV));
}
getFeaturesRank <- function(counts, norm ){
	M3Drop_table = M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=2, suppress.plot=TRUE)
	M3Drop_table[,1] = as.character(M3Drop_table[,1])
	HVG_Deng = BrenneckeGetVariableGenes(norm, fdr=2, suppress.plot=TRUE)
	counts <- as.matrix(counts)
	fit <- NBumiFitModel(counts)
	Dengfeatures <- NBumiFeatureSelectionCombinedDrop(fit)
	Dengfeatures2 <- NBumiFeatureSelectionHighVar(fit)

	# Features!
	ordered_genes = rownames(counts);
	NB = Dengfeatures[match(ordered_genes,names(Dengfeatures))]
	NBV = Dengfeatures2[match(ordered_genes,names(Dengfeatures2))]
	M3Drop = M3Drop_table[match(ordered_genes, M3Drop_table[,1]) ,3]
	HVG = HVG_Deng[match(ordered_genes, HVG_Deng[,1]) ,3]

	TABLE <- cbind(NB, NBV,M3Drop, HVG)
	rownames(TABLE) = names(NB);

	return(TABLE);
}
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

# Make Consistent
consistent_genes <- rownames(Deng_embyro_list$data)[rownames(Deng_embyro_list$data) %in% rownames(zhong_list$data) & rownames(Deng_embyro_list$data) %in% rownames(Xue_list$data)
				& rownames(Deng_embyro_list$data) %in% rownames(Fan_list$data) & rownames(Deng_embyro_list$data) %in% rownames(Goolam_counts$data)]
counts_list$data <- counts_list$data[rownames(counts_list$data) %in% consistent_genes,]
norm_list$data <- norm_list$data[rownames(norm_list$data) %in% consistent_genes,]
zhong_count <- zhong_count[rownames(zhong_count) %in% consistent_genes,]
zhong_list$data <- zhong_list$data[rownames(zhong_list$data) %in% consistent_genes,]
Xue_count <- Xue_count[rownames(Xue_count) %in% consistent_genes,]
Xue_list$data <- Xue_list$data[rownames(Xue_list$data) %in% consistent_genes,]
Fan_count <- Fan_count[rownames(Fan_count) %in% consistent_genes,]
Fan_list$data <- Fan_list$data[rownames(Fan_list$data) %in% consistent_genes,]
Goolam_counts$data <- Goolam_counts$data[rownames(Goolam_counts$data) %in% consistent_genes,]
Goolam_list$data <- Goolam_list$data[rownames(Goolam_list$data) %in% consistent_genes,]


# Features
Deng <- getFeatures(counts_list$data, norm_list$data, name="Deng")
scores_Deng <- getFeaturesRank(counts_list$data, norm_list$data)
Biase <- getFeatures(zhong_count, zhong_list$data, name="Biase");
scores_Biase <- getFeaturesRank(zhong_count, zhong_list$data);
Xue <- getFeatures(Xue_count, Xue_list$data, name="Xue");
scores_Xue <- getFeaturesRank(Xue_count, Xue_list$data);
Fan <- getFeatures(Fan_count, Fan_list$data, name="Fan")
scores_Fan <- getFeaturesRank(Fan_count, Fan_list$data)
Goolam <- getFeatures(Goolam_counts$data, Goolam_list$data, name="Goolam")
scores_Goolam <- getFeaturesRank(Goolam_counts$data, Goolam_list$data)

# Kendall - rank cor consistency
gene_set = rownames(scores_Deng)[rownames(scores_Deng) %in% rownames(scores_Biase) & rownames(scores_Deng) %in% rownames(scores_Xue) & rownames(scores_Deng) %in% rownames(scores_Fan) & rownames(scores_Deng) %in% rownames(scores_Goolam)]
TABLE <- cbind(scores_Deng[rownames(scores_Deng) %in% gene_set,], scores_Biase[rownames(scores_Biase) %in% gene_set,], 
		scores_Xue[rownames(scores_Xue) %in% gene_set,], scores_Fan[rownames(scores_Fan) %in% gene_set,], 
		scores_Goolam[rownames(scores_Goolam) %in% gene_set,])
heat_data <- cor(TABLE, method="kendall")
heat_data2 <- cor(TABLE, method="spearman")
dataset_labels <- rep(c("Deng", "Biase","Xue","Fan","Goolam"), each=length(scores_Deng[1,]))
method_labels <- rep(c("NB", "NBV", "M3Drop","HVG"), times=5)
heatmap.2(heat_data, col=rev(brewer.pal(9, "RdBu")), symm=TRUE, scale="none", trace="n", ColSideColors=brewer.pal(5, "Set2")[factor(dataset_labels)])
heatmap.2(heat_data2, col=rev(brewer.pal(9, "RdBu")), symm=TRUE, scale="none", trace="n", ColSideColors=brewer.pal(5, "Set2")[factor(dataset_labels)])

# Consistency.
my_datasets = list(Deng=Deng, Biase=Biase, Fan=Fan, Xue=Xue, Goolam=Goolam)
n_meth = length(my_datasets[[1]])
n_set = length(my_datasets)
my_Matrix = matrix(0, ncol=n_meth*n_set, nrow=length(consistent_genes))
my_List = list()
rownames(my_Matrix) = consistent_genes

my_names = vector(length=n_meth*n_set)
for(d in 1:length(my_datasets)) {
	for(m in 1:length(my_datasets[[d]])) {
		my_Matrix[,(d-1)*n_meth+m] = consistent_genes %in% my_datasets[[d]][[m]]
		name <- paste(names(my_datasets)[d],names(my_datasets[[d]])[m], sep="-")
		my_names[(d-1)*n_meth+m] = name
		my_List[[(d-1)*n_meth+m]] = as.character(my_datasets[[d]][[m]][ my_datasets[[d]][[m]] %in% consistent_genes ])
	}
}
colnames(my_Matrix) = my_names

# Chi-sq test
m3d = seq(from=2, to=n_meth*n_set, by=n_meth)
hvg = seq(from=3, to=n_meth*n_set, by=n_meth)
nb = seq(from=1, to=n_meth*n_set, by=n_meth)
nbv = seq(from=4, to=n_meth*n_set, by=n_meth)

get_prob <- function(probs, k=length(probs)) {
	yes <- combn(1:length(probs), k)
	p <- sum(apply(yes, 2, function(y) {prod(probs[y],1-probs[-y])}))
}

get_chi_sq <- function(my_mat) {
	obs = sapply(1:length(my_mat[1,]), function(i) {sum(rowSums(my_mat) == i)})
	probs = colSums(my_mat)/length(my_mat[,1])
	exp <- sapply(1:length(my_mat[1,]), function(i) {length(my_mat[,1])*get_prob(probs,i)})
	return(list(chisq=sum((obs-exp)^2/exp), obs=obs, exp=exp))
}

m3d_stats = get_chi_sq(my_Matrix[,m3d])
hvg_stats = get_chi_sq(my_Matrix[,hvg])
nb_stats  = get_chi_sq(my_Matrix[,nb])
nbv_stats = get_chi_sq(my_Matrix[,nbv])

png("Figure_Reproducibility.png", width=5, height=5, units="in", res=300)
par(mar=c(4,4,1,1))
xes <- barplot(rbind(m3d_stats$obs-m3d_stats$exp, hvg_stats$obs-hvg_stats$exp, nb_stats$obs-nb_stats$exp, nbv_stats$obs-nbv_stats$exp), beside=TRUE, ylab="Obs-Exp", names=c("1","2","3","4","5"), col=c(MM_col, hvg_1_col, Depth_col, NBVar_col), xlab="Number of Datasets")
legend("bottomright", paste(c("M3Drop","HVG","NB","NBV"), "\nChi2:",round(c(m3d_stats$chisq, hvg_stats$chisq, nb_stats$chisq, nbv_stats$chisq)), "\n"), fill=c(MM_col, hvg_1_col, Depth_col, NBVar_col),bty="n")
#lines(c(xes[1,1],xes[4,1]),c(m3d_stats$exp[1], m3d_stats$exp[1]), col="red", lty=2)
#lines(c(xes[1,2],xes[4,2]),c(m3d_stats$exp[2], m3d_stats$exp[2]), col="red", lty=2)
#lines(c(xes[1,3],xes[4,3]),c(m3d_stats$exp[3], m3d_stats$exp[3]), col="red", lty=2)
#lines(c(xes[1,4],xes[4,4]),c(m3d_stats$exp[4], m3d_stats$exp[4]), col="red", lty=2)
#lines(c(xes[1,5],xes[4,5]),c(m3d_stats$exp[5], m3d_stats$exp[5]), col="red", lty=2)
dev.off()
