require("RColorBrewer")
require("gplots")
source("Colour_Scheme.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Other_FS_functions.R")


map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Mmus_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

getFeatures <- function(counts, norm, fdr=0.01, name="Test", suppress.plot=TRUE){
#	M3Drop_table = M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=fdr, suppress.plot=TRUE)
	M3Drop_table = M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=2, suppress.plot=TRUE)
	M3Drop_table[,1] = as.character(M3Drop_table[,1])
#	HVG_Deng = BrenneckeGetVariableGenes(norm, fdr=fdr, suppress.plot=TRUE)
	HVG_Deng = BrenneckeGetVariableGenes(norm, fdr=2, suppress.plot=TRUE)
	counts <- as.matrix(counts)
	fit <- NBumiFitModel(counts)
	if (!suppress.plot) {
		png(paste(name,"NBFit.png",sep="_"), width=7, height=5, units="in", res=300)
		par(mfrow=c(1,2))
		NBumiCheckFit(counts, fit);
		dev.off();
	}
	Dengfeatures <- NBumiFeatureSelectionCombinedDrop(fit)
	Dengfeatures2 <- NBumiFeatureSelectionHighVar(fit)
	Gini = Gini_FS(norm)
	negcor = Cor_FS_neg(norm)
	pca1 = Monocle2_pca_FS(counts, 1:length(counts[1,]), pcs=c(2,3))
	pca2 = Monocle2_pca_FS(counts, 1:length(counts[1,]), pcs=c(1,2,3))

	# Features!
#	NB = names(Dengfeatures[p.adjust(Dengfeatures, method="fdr") < fdr]);
#	NBV = names(Dengfeatures2[p.adjust(Dengfeatures2, method="fdr") < fdr]);
	NB = names(Dengfeatures[1:2000]);
	NBV = names(Dengfeatures2[1:2000]);
	M3Drop = M3Drop_table[1:2000,1]
	HVG = HVG_Deng[1:2000,1]
	return(list(NB=NB,M3Drop=M3Drop,HVG=HVG, NBV=NBV, Gini = names(Gini[1:2000]), NCor =names(negcor[1:2000]), PC23=names(pca1[1:2000]), PC123 = names(pca2[1:2000])));
}
convert_to_integer <- function(mat) {
	mat <- round(as.matrix(mat))
	storage.mode(mat) <- "integer"
	mat = mat[rowSums(mat) > 0,]
	return(mat)
}

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
counts_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = FALSE)
norm_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = TRUE)

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
zhong_list$data <- zhong_list$data[rownames(zhong_list$data) %in% rownames(zhong_count),]

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
Xue_list$data <- Xue_list$data[rownames(Xue_list$data) %in% rownames(Xue_count),]

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
Fan_list$data <- Fan_list$data[rownames(Fan_list$data) %in% rownames(Fan_count),]

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
consistent_genes <- rownames(norm_list$data)[rownames(norm_list$data) %in% rownames(zhong_list$data) & rownames(norm_list$data) %in% rownames(Xue_list$data)
				& rownames(norm_list$data) %in% rownames(Fan_list$data) & rownames(norm_list$data) %in% rownames(Goolam_counts$data)]
consistent_genes = sort(consistent_genes);

counts_list$data <- counts_list$data[match(consistent_genes, rownames(counts_list$data)),]
norm_list$data <- norm_list$data[match(consistent_genes, rownames(norm_list$data)),]
zhong_count <- zhong_count[match(consistent_genes, rownames(zhong_count)),]
zhong_list$data <- zhong_list$data[match(consistent_genes, rownames(zhong_list$data)),]
Xue_count <- Xue_count[match(consistent_genes, rownames(Xue_count)),]
Xue_list$data <- Xue_list$data[match(consistent_genes, rownames(Xue_list$data)),]
Fan_count <- Fan_count[match(consistent_genes, rownames(Fan_count)),]
Fan_list$data <- Fan_list$data[match(consistent_genes, rownames(Fan_list$data)),]
Goolam_counts$data <- Goolam_counts$data[match(consistent_genes, rownames(Goolam_counts$data)),]
Goolam_list$data <- Goolam_list$data[match(consistent_genes, rownames(Goolam_list$data)),]

# Combine


Combined_counts = cbind(norm_list$data, zhong_count, Xue_count, Fan_count , Goolam_counts$data);
Combined_norm = cbind(norm_list$data, zhong_list$data, Xue_list$data, Fan_list$data, Goolam_list$data);

dataset_labels = c(rep("Deng", times=length(norm_list$data[1,])), rep("Biase", times=length(zhong_count[1,])), rep("Xue", times=length(Xue_count[1,])), rep("Fan", times=length(Fan_count[1,])), rep("Goolam", times=length(Goolam_counts$data[1,])))

truth = c(as.character(unlist(counts_list$labels)), as.character(unlist(zhong_list$labels)), as.character(unlist(Xue_list$labels)),  as.character(unlist(Fan_list$labels)), as.character(unlist(Goolam_list$labels)))
truth[grep("ocyte",truth)] = "zygote"
truth[grep("nuc",truth)] = "zygote"
truth[grep("AM",truth)] = "zygote"
truth[grep("early2",truth)] = "zygote"
truth[grep("2cell",truth)] = "2cell"
truth[grep("4cell",truth)] = "4cell"
truth[grep("8cell",truth)] = "8cell"
truth[grep("orula",truth)] = "16cell"
truth[grep("blast",truth)] = "blast"
truth[grep("32cell",truth)] = "blast"

truth_TE_ICM = truth
truth_blast = truth
truth_blast[grep("TE",truth_blast)] = "blast"
truth_blast[grep("ICM",truth_blast)] = "blast"

# Plotting Fxns & setup
Stage = factor(truth_blast, levels=c("zygote","2cell","4cell","8cell","16cell","blast"));
Source = factor(dataset_labels);
pch_set = c(17,1,18,15,16)
col_set = brewer.pal("Set2", n=length(unique(Stage)))

# Features
Deng <- getFeatures(counts_list$data, norm_list$data, name="Deng")
Biase <- getFeatures(zhong_count, zhong_list$data, name="Biase");
Xue <- getFeatures(Xue_count, Xue_list$data, name="Xue");
Fan <- getFeatures(Fan_count, Fan_list$data, name="Fan")
Goolam <- getFeatures(Goolam_counts$data, Goolam_list$data, name="Goolam")

# Consistent
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

# Reproducibility Plot based on Chi-sq test
m3d = seq(from=2, to=n_meth*n_set, by=n_meth)
hvg = seq(from=3, to=n_meth*n_set, by=n_meth)
nb = seq(from=1, to=n_meth*n_set, by=n_meth)
nbv = seq(from=4, to=n_meth*n_set, by=n_meth)
gini = seq(from=5, to=n_meth*n_set, by=n_meth)
ncor = seq(from=6, to=n_meth*n_set, by=n_meth)
pc23 = seq(from=7, to=n_meth*n_set, by=n_meth)
pc123 = seq(from=8, to=n_meth*n_set, by=n_meth)

# PCAs
quality <- matrix(nrow=n_meth, ncol=2)
col_sets <- list(m3d, hvg, nb, nbv, gini, ncor, pc23, pc123)
names <- c("M3Drop", "HVG", "NBDrop", "NBDisp", "Gini", "Cor", "PCA(2-3)", "PCA(1-3)")

make_PCA <- function(gene_list) {
	par(mar=c(3,3,2,1));
	toplot <- log(Combined_norm[rownames(Combined_norm) %in% gene_list,]+1)/log(2);
	PCA = prcomp(toplot);
	plot(PCA$rotation[,1], PCA$rotation[,2], col=col_set[Stage], pch=pch_set[Source])
	title(xlab=paste("PC1 (",round(PCA$sdev[1]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), line=2)
	title(ylab=paste("PC2 (",round(PCA$sdev[2]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), line=2)
}
png("FigureS6.png", width=8, height=8, units="in", res=300)
par(mfrow=c(3,3))
for (i in 1:n_meth) {
	thing = make_PCA(rownames(my_Matrix[rowSums(my_Matrix[,col_sets[[i]]]) > 2,]));
	title(main=names[i], line=1)
}

plot(1,1, xaxt="n", yaxt="n", main="", xlab="",ylab="", pch=4, bty="n", col="white")
#legend("left", c("Stage",levels(Stage),"Dataset",levels(Source)),
#col=c("white",col_set, "white", rep("black", times=length(pch_set))), 
#pch=c(NA, rep(16, times=length(col_set)), NA, pch_set), bty="n", ncol=2)

par(cex=1.1)
legend("topleft", levels(Stage),col=col_set, pch=16, bty="n", title="Stage")
legend("topright", levels(Source),col="black", pch=pch_set, bty="n", title="Dataset")
dev.off()
