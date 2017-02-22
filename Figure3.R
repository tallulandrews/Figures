require("RColorBrewer")
require("gplots")
source("Colour_Scheme.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Other_FS_functions.R")

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

# Zhong
zhong = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_fpkm_ZHONG.txt", header=TRUE);
zhong = zhong[!duplicated(zhong[,1]),]
rownames(zhong) = zhong[,1]
zhong = zhong[,2:length(zhong[1,])]
zhong = as.matrix(zhong);

zhong_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_labels_ZHONG.txt");
ultralow = which(rowMeans(zhong) < 10^-5)
zhong = zhong[-ultralow,]
zhong_labels = as.character(unlist(zhong_labels))
zhong_list = normalize_data(zhong, labels=zhong_labels, is.counts=FALSE)
zhong_count = convert_to_integer(zhong_list$data);
zhong_list$labels = get_ICM_TE_assignments(zhong_list$data, zhong_list$labels)

# Xue
Xue_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_mat.txt", header=TRUE)
Xue_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_labels.txt", header=FALSE)
Xue_labels = as.character(unlist(Xue_labels))
Xue_labels[Xue_labels == "2cellmixed"] = "2cellmix"
Xue_labels[Xue_labels == "4cellmixed"] = "4cellmix"
Xue_labels[Xue_labels == "8cellmixed"] = "8cellmix"
Xue_labels[Xue_labels == "Pronucleus"] = "pronuc"
Xue_list = normalize_data(Xue_data, labels=Xue_labels, is.counts=FALSE)
Xue_list$labels = get_ICM_TE_assignments(Xue_list$data, Xue_list$labels)
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
Fan_list$labels = get_ICM_TE_assignments(Fan_list$data, Fan_list$labels)
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
Goolam_list$labels = get_ICM_TE_assignments(Goolam_list$data, Goolam_list$labels)
Goolam_counts$labels = get_ICM_TE_assignments(Goolam_counts$data, Goolam_counts$labels)

# Make Consistent
genes = rownames(Deng_embyro_list$data)
consistent_genes <- genes[genes %in% rownames(zhong_list$data) & 
		    genes %in% rownames(Xue_list$data) & 
		    genes %in% rownames(Fan_list$data) & 
		    genes %in% rownames(Goolam_counts$data) &
		    genes %in% rownames(norm_list$data) &
		    genes %in% rownames(Xue_count) &
		    genes %in% rownames(Fan_count) &
		    genes %in% rownames(Goolam_counts$data) &
		    genes %in% rownames(zhong_count)
		    ]
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


# Features
Deng <- getFeatures(counts_list$data, norm_list$data, name="Deng")
Biase <- getFeatures(zhong_count, zhong_list$data, name="Biase");
Xue <- getFeatures(Xue_count, Xue_list$data, name="Xue");
Fan <- getFeatures(Fan_count, Fan_list$data, name="Fan")
Goolam <- getFeatures(Goolam_counts$data, Goolam_list$data, name="Goolam")

scores_Deng <- getFeaturesRank(counts_list$data, norm_list$data)
scores_Biase <- getFeaturesRank(zhong_count, zhong_list$data);
scores_Xue <- getFeaturesRank(Xue_count, Xue_list$data);
scores_Fan <- getFeaturesRank(Fan_count, Fan_list$data)
scores_Goolam <- getFeaturesRank(Goolam_counts$data, Goolam_list$data)

# Kendall - rank cor consistency
#gene_set = rownames(scores_Deng)[rownames(scores_Deng) %in% rownames(scores_Biase) & rownames(scores_Deng) %in% rownames(scores_Xue) & rownames(scores_Deng) %in% rownames(scores_Fan) & rownames(scores_Deng) %in% rownames(scores_Goolam)]
#TABLE <- cbind(scores_Deng[rownames(scores_Deng) %in% gene_set,], scores_Biase[rownames(scores_Biase) %in% gene_set,], 
#		scores_Xue[rownames(scores_Xue) %in% gene_set,], scores_Fan[rownames(scores_Fan) %in% gene_set,], 
#		scores_Goolam[rownames(scores_Goolam) %in% gene_set,])
#heat_data <- cor(TABLE, method="kendall")
#heat_data2 <- cor(TABLE, method="spearman")
#dataset_labels <- rep(c("Deng", "Biase","Xue","Fan","Goolam"), each=length(scores_Deng[1,]))
#method_labels <- rep(c("NB", "NBV", "M3Drop","HVG"), times=5)
#heatmap.2(heat_data, col=rev(brewer.pal(9, "RdBu")), symm=TRUE, scale="none", trace="n", ColSideColors=brewer.pal(5, "Set2")[factor(dataset_labels)])
#heatmap.2(heat_data2, col=rev(brewer.pal(9, "RdBu")), symm=TRUE, scale="none", trace="n", ColSideColors=brewer.pal(5, "Set2")[factor(dataset_labels)])

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
gini = seq(from=5, to=n_meth*n_set, by=n_meth)
ncor = seq(from=6, to=n_meth*n_set, by=n_meth)
pc23 = seq(from=7, to=n_meth*n_set, by=n_meth)
pc123 = seq(from=8, to=n_meth*n_set, by=n_meth)

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
gini_stats = get_chi_sq(my_Matrix[,gini])
ncor_stats = get_chi_sq(my_Matrix[,ncor])
pc23_stats = get_chi_sq(my_Matrix[,pc23])
pc123_stats = get_chi_sq(my_Matrix[,pc123])

# Figure Setup
png("Figure3.png", width=6, height=10, units="in", res=300)
layout(matrix(c(rep(1,times=6), rep(c(2,3), each=3), rep(c(4,5,6), each=2)), ncol=6, byrow=T))

# Reproducibility Figure (A)

toplot <- rbind(m3d_stats$obs-m3d_stats$exp, hvg_stats$obs-hvg_stats$exp, 
		nb_stats$obs-nb_stats$exp, nbv_stats$obs-nbv_stats$exp,
		gini_stats$obs-gini_stats$exp, ncor_stats$obs-ncor_stats$exp,
		pc23_stats$obs-pc23_stats$exp, pc123_stats$obs-pc123_stats$exp);
score <- c(m3d_stats$chisq, hvg_stats$chisq, nb_stats$chisq, nbv_stats$chisq,
	   gini_stats$chisq, ncor_stats$chisq, pc23_stats$chisq, 
	   pc123_stats$chisq)
leg_text <- paste(c("M3Drop","HVG","NBDrop","NBDisp", "Gini","Cor",
		    "PCA(2-3)","PCA(1-3)"), "\nChi2:", round(score), "\n")
bar_col <- c(MM_col, hvg_1_col, Depth_col, NBVar_col, gini_col, cor_col, pca_1_col, pca_2_col)

sort_order =  order(score)


par(mar=c(4,4,1,1))
xes <- barplot( toplot[sort_order,], beside=TRUE, 
		ylab="Number of Genes (Obs-Exp)", 
		names=c("1","2","3","4","5"), 
		col=bar_col[sort_order], 
		xlab="Number of Datasets", ylim=c(-2500, 501))

legend("bottomright", leg_text[sort_order], 
	fill=bar_col[sort_order], bty="n", ncol=2)
#lines(c(xes[1,1],xes[4,1]),c(m3d_stats$exp[1], m3d_stats$exp[1]), col="red", lty=2)
#lines(c(xes[1,2],xes[4,2]),c(m3d_stats$exp[2], m3d_stats$exp[2]), col="red", lty=2)
#lines(c(xes[1,3],xes[4,3]),c(m3d_stats$exp[3], m3d_stats$exp[3]), col="red", lty=2)
#lines(c(xes[1,4],xes[4,4]),c(m3d_stats$exp[4], m3d_stats$exp[4]), col="red", lty=2)
#lines(c(xes[1,5],xes[4,5]),c(m3d_stats$exp[5], m3d_stats$exp[5]), col="red", lty=2)

# Combine

Combined_counts = cbind(counts_list$data, zhong_count, Xue_count, Fan_count , Goolam_counts$data);
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

truth_ICM_TE = truth
truth_blast = truth
truth_blast[grep("TE",truth_blast)] = "blast"
truth_blast[grep("ICM",truth_blast)] = "blast"

# Plotting Fxns & setup
Stage = factor(truth_blast, levels=c("zygote","2cell","4cell","8cell","16cell","blast"));
Source = factor(dataset_labels);
pch_set = c(17,1,18,15,16)
col_set = brewer.pal("Set2", n=length(unique(Stage)))

make_PCA <- function(gene_list) {
	toplot <- log(Combined_norm[rownames(Combined_norm) %in% gene_list,]+1)/log(2);
	PCA = prcomp(toplot);
	plot(PCA$rotation[,1], PCA$rotation[,2], col=col_set[Stage], pch=pch_set[Source], xlab=paste("PC1 (",round(PCA$sdev[1]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA$sdev[2]/sum(PCA$sdev)*100, digits=2)," %)",sep=""))
	dists <- dist(t(toplot), method="euclidean")
	htree <- hclust(dists, method="ward.D2")
	require("igraph")
        sim = vector(length = length(toplot[1,]))
        for(k in 1:length(toplot[1,])) {
                clusters = cutree(htree, k)
                sim[k] = compare(clusters, as.factor(truth_blast), method="adjusted.rand")
        }
	return(sim);
}

# PCA plots (B & C)

out1 = make_PCA(rownames(my_Matrix[rowSums(my_Matrix[,m3d]) > 2, ]))
out2 = make_PCA(rownames(my_Matrix))


##### Anxa2 S100a10


# Anxa2 S100a10 correlations
heat_data <- log(Combined_norm+1)/log(2);

#png("Combined_Axna2_S100a10.png", width=5*3/1.25*0.8, height=5*0.8, units="in", res=300)
par(mar=c(4,4,2.5,1))
# Known Axna2 & S100a10 heterotetramer involved in embryo implantation
# Wang (2015) Annexin A2 Acts as an Adhesion Molecule on the Endometrial Epithelium during Implantation in Mice.
Anxa2 = unlist(heat_data[rownames(heat_data) == "Anxa2",])
S100a10 = unlist(heat_data[rownames(heat_data) == "S100a10",])

c8 = truth_blast=="8cell"
c16 = truth_blast=="16cell"
bla = truth_blast=="blast"

print("8cell")
plot(Anxa2[c8], S100a10[c8], pch=pch_set[Source[c8]], xlab="Anxa2", ylab="S100a10", col="grey50", main="8cell", cex=1.5)
cout = cor.test(Anxa2[c8], S100a10[c8])
legend("top", paste(c("r = ", "p = "),c(round(cout$estimate, digits=2),signif(cout$p.value, digits=2)), sep=""), bty="n")

print("16cell")
plot(Anxa2[c16], S100a10[c16], pch=pch_set[Source[c16]], xlab="Anxa2", ylab="S100a10", col="grey50", main="16cell", cex=1.5)
cout = cor.test(Anxa2[c16], S100a10[c16])
legend("top", paste(c("r = ", "p = "),c(round(cout$estimate, digits=2),signif(cout$p.value, digits=2)), sep=""), bty="n")

print("blast")
blast_cols = rep("grey50",times=length(truth_ICM_TE)); 
blast_cols[truth_ICM_TE=="TE"] = TE_col; 
blast_cols[truth_ICM_TE=="ICM"]= ICM_col;
plot(Anxa2[bla], S100a10[bla], pch=pch_set[Source[bla]], xlab="Anxa2", ylab="S100a10", col=blast_cols[bla], main="blastocyst", cex=1.5)

tidy1 = Anxa2[truth_ICM_TE=="TE" & Anxa2>4 & S100a10>4]
tidy2 = S100a10[truth_ICM_TE=="TE" & Anxa2>4 & S100a10>4]
reg=lm(tidy2~0+tidy1)
abline(reg)
cout = cor.test(Anxa2[bla], S100a10[bla])
legend("top", paste(c("r = ", "p = "),c(round(cout$estimate, digits=2),signif(cout$p.value, digits=2)), sep=""), bty="n")
dev.off()

