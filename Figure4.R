require("RColorBrewer")
require("gplots")
require("M3Drop")
source("Colour_Scheme.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Other_FS_functions.R")
source("heatboxplot.R")

ICM_col="dodgerblue2"
TE_col="forestgreen"

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Mmus_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}

get_ICM_TE_assignments<- function(x, labels) {
	data <- log(x+1)/log(2);
	# TE & ICM assignment
	TE_markers=c("Elf5","Eomes","Cdx2")
	ICM_markers = c("Sox2","Pou5f1","Nanog")
	scale <- function(x) {(x-mean(x))/sd(x)}
	ICMscore = rowMeans(apply(data[rownames(data) %in% ICM_markers,], 1, scale))
	TEscore = rowMeans(apply(data[rownames(data) %in% TE_markers,], 1, scale))
	blasts = (grepl("blast", labels) | grepl("32cel", labels)) 
	new_lab = as.character(labels);
	if (sum(blasts) > 0) {
		new_lab[blasts & ICMscore > TEscore] = "ICM"
		new_lab[blasts & ICMscore < TEscore] = "TE"
	}
	return(new_lab);
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
	mat <- ceiling(as.matrix(mat))
	storage.mode(mat) <- "integer"
	mat = mat[rowSums(mat) > 0,]
	return(mat)
}

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
fdr=0.01
counts_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = FALSE)
norm_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = TRUE)
counts_list$labels = get_ICM_TE_assignments(counts_list$data, counts_list$labels)
norm_list$labels = get_ICM_TE_assignments(norm_list$data, norm_list$labels)

# Features
Deng <- getFeatures(counts_list$data, norm_list$data, name="Deng")

scores_Deng <- getFeaturesRank(counts_list$data, norm_list$data)

consistent_genes = rownames(norm_list$data);

# Consistency.
my_datasets = list(Deng=Deng)
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
gini = seq(from=5, to=n_meth*n_set, by=n_meth)
ncor = seq(from=6, to=n_meth*n_set, by=n_meth)
pc23 = seq(from=7, to=n_meth*n_set, by=n_meth)
pc123 = seq(from=8, to=n_meth*n_set, by=n_meth)

# Figure Setup

TE_markers=c("Elf5","Eomes","Cdx2")
ICM_markers = c("Sox2","Pou5f1","Nanog")

blasts = grepl("blast", Deng_embyro_list$labels)
all_data <- log(norm_list$data+1)/log(2);
all_data <- all_data[,blasts]
m3d_genes = unique(c(rownames(my_Matrix)[my_Matrix[,2]==1], ICM_markers, TE_markers))
m3d_data <- all_data[rownames(all_data) %in% m3d_genes,]
labs=factor(Deng_embyro_list$labels[blasts], levels=c("earlyblast","midblast","lateblast"))
lab_cols=c("#b3cde3","#8c96c6","#88419d")

markers = c(ICM_markers, TE_markers)
mark_col = c(rep(ICM_col, times=length(ICM_markers)), 
		rep(TE_col, times=length(TE_markers)))
clust_fun <- function(x) {hclust(x, method="ward.D2")}

png("Figure4All.png", width=6, height=5, units="in", res=300)
epic_dendro_boxplots(all_data, distfun=dist, hclustfun=clust_fun, 4, markers, mark_col, labs, lab_cols)
dev.off()

png("Figure4M3D.png", width=6, height=5, units="in", res=300)
epic_dendro_boxplots(m3d_data, distfun=dist, hclustfun=clust_fun, 4, markers, mark_col, labs, lab_cols)
dev.off()
