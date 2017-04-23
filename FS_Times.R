require("RColorBrewer")
require("gplots")
source("Colour_Scheme.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Other_FS_functions.R")


source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

getFeatures <- function(counts, norm, fdr=0.01, name="Test", suppress.plot=TRUE){
	counts = counts[rowSums(counts) > 0,]
	norm = norm[rowSums(norm) > 0,]
	t = c(proc.time()[3]) #1
	M3Drop_table = M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=2, suppress.plot=TRUE)
	M3Drop_table[,1] = as.character(M3Drop_table[,1])
	t = c(t, proc.time()[3]) #2
	HVG_Deng = BrenneckeGetVariableGenes(norm, fdr=2, suppress.plot=TRUE)
	t = c(t, proc.time()[3]) #3
	counts <- as.matrix(counts)
	fit <- NBumiFitModel(counts)
	t = c(t, proc.time()[3]) #4
	Dengfeatures <- NBumiFeatureSelectionCombinedDrop(fit)
	t = c(t, proc.time()[3]) #5
	Dengfeatures2 <- NBumiFeatureSelectionHighVar(fit)
	t = c(t, proc.time()[3]) #6
	Gini = Gini_FS(norm)
	t = c(t, proc.time()[3]) #7
	negcor = Cor_FS_neg(norm)
	t = c(t, proc.time()[3]) #8
	pca1 = Monocle2_pca_FS(counts, 1:length(counts[1,]), pcs=c(2,3))
	t = c(t, proc.time()[3]) #9
	pca2 = Monocle2_pca_FS(counts, 1:length(counts[1,]), pcs=c(1,2,3))
	t = c(t, proc.time()[3]) #10

	m3d_t = t[2]-t[1]; hvg_t = t[3]-t[2]; 
	nb_t = t[5]-t[3]; nbv_t = (t[6]-t[5])+(t[4]-t[3]); 
	gini_t = t[7]-t[6]; cor_t = t[8]-t[7]; 
	pca1_t = t[9]-t[8]; pca2_t = t[10]-t[9];
	return(c(m3d_t, hvg_t, nb_t, nbv_t, gini_t, cor_t, pca1_t, pca2_t))
}
convert_to_integer <- function(mat) {
	mat <- round(as.matrix(mat))
	storage.mode(mat) <- "integer"
	mat = mat[rowSums(mat) > 0,]
	return(mat)
}

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
counts_list = normalize_data(Deng_embyro_list$data, is.counts = FALSE)
norm_list = normalize_data(Deng_embyro_list$data, is.counts = TRUE)

Deng <- getFeatures(counts_list$data, norm_list$data, name="Deng")
Deng_dim <- dim(counts_list$data);

# Zhong
zhong = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_fpkm_ZHONG.txt", header=TRUE);
zhong = zhong[!duplicated(zhong[,1]),]
rownames(zhong) = zhong[,1]
zhong = zhong[,-1]
zhong = as.matrix(zhong);
ultralow = which(rowMeans(zhong) < 10^-5)
zhong = zhong[-ultralow,]
zhong_list = normalize_data(zhong, is.counts=FALSE)
zhong_count = convert_to_integer(zhong_list$data);

Biase <- getFeatures(zhong_count, zhong_list$data, name="Biase");
Biase_dim <- dim(zhong_count);

# Xue
Xue_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_mat.txt", header=TRUE)
Xue_list = normalize_data(Xue_data, is.counts=FALSE)
Xue_count = convert_to_integer(Xue_list$data);

Xue <- getFeatures(Xue_count, Xue_list$data, name="Xue");
Xue_dim <- dim(Xue_count)

# Fan
Fan_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE53386_matrix_fpkms.tsv", header=TRUE)
Fan_list = normalize_data(Fan_data, is.counts=FALSE)
Fan_count = convert_to_integer(Fan_list$data)

Fan <- getFeatures(Fan_count, Fan_list$data, name="Fan")
Fan_dim <- dim(Fan_count);

# Goolam
Goolam_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Goolam_et_al_2015_count_table.tsv", header=T)
Goolam_list = normalize_data(Goolam_data, is.counts=TRUE)
Goolam_counts = normalize_data(Goolam_data, is.counts=FALSE)

Goo <- getFeatures(Goolam_counts$data, Goolam_list$data, name="Goolam")
Goo_dim <- dim(Goolam_counts$data)

#### Main-Text Ola & Blischak ####

source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Ola_SC.R")
Ola_count <- as.matrix(data[-spikes,]);
Ola_count <- Ola_count[rowSums(Ola_count) > 0,]
Ola_norm <- normalize_data(Ola_count, is.counts=TRUE)

Ola <- getFeatures(NBumiConvertToInteger(Ola_count), Ola_norm$data, name="Ola")
Ola_dim <- dim(Ola_norm$data);

source("../My_R_packages/M3D/R/NB_UMI.R"); require("matrixStats");
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Blishcak_UMI.R")
counts <- counts[rowSums(counts) > 0,]
norm <- normalize_data(counts, is.counts=TRUE)

Blish <- getFeatures(counts, norm$data, name="Blish")
Blish_dim <- dim(norm$data);

# Shalek
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/ShalekKO_clean.RData")
ShalekKO_list$data <- ShalekKO_list$data[rowSums(ShalekKO_list$data) > 0,]
norm <- normalize_data(ShalekKO_list$data, is.counts=TRUE)

Shalek <- getFeatures(NBumiConvertToInteger(ShalekKO_list$data), norm$data, name="Sha")
Shalek_dim <- dim(norm$data);
rm(ShalekKO_list)


# Buettner
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TeichCC_clean.RData")
TeichCC_list$data <- TeichCC_list$data[rowSums(TeichCC_list$data) > 0,]
norm <- normalize_data(TeichCC_list$data, is.counts=TRUE)

Teic <- getFeatures(NBumiConvertToInteger(TeichCC_list$data), norm$data, name="T")
Teic_dim <- dim(norm$data);
rm(TeichCC_list)

# Pollen
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/pollen.rda")
pollen <- pollen[rowSums(pollen) > 0,]
colnames(pollen) <- 1:length(pollen[1,])
norm <- normalize_data(pollen, is.counts=FALSE);

Pollen <- getFeatures(NBumiConvertToInteger(pollen), norm$data, name="pollen")
Pollen_dim <- dim(norm$data);
rm(pollen)

# Kirschner
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/kirschner.rda")
kirschner <- kirschner[rowSums(kirschner) > 0,]
colnames(kirschner) <- 1:length(kirschner[1,])
data_list = normalize_data(kirschner, labels = 1:length(kirschner[1,]), is.counts = FALSE)

Kir <- getFeatures(NBumiConvertToInteger(kirschner), data_list$data, name="kir")
Kir_dim <- dim(data_list$data)
rm(kirschner)

# Linnarsson
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/linnarsson.rda")
linnarsson <- linnarsson[rowSums(linnarsson) > 0,]
colnames(linnarsson) <- 1:length(linnarsson[1,])
data_list = normalize_data(linnarsson, labels = 1:length(linnarsson[1,]), is.counts = FALSE)

Lin <- getFeatures(NBumiConvertToInteger(linnarsson), data_list$data, name="lin")
Lin_dim <- dim(data_list$data)
rm(linnarsson)

# Usoskin
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/usoskin.rda")
usoskin <- usoskin[rowSums(usoskin) > 0,]
colnames(usoskin) <- 1:length(usoskin[1,])
data_list = normalize_data(usoskin, labels = 1:length(usoskin[1,]), is.counts = FALSE)

Uso <- getFeatures(NBumiConvertToInteger(usoskin), data_list$data, name="uso")
Uso_dim <- dim(data_list$data)
rm(usoskin)

# Macosko
mac_data <- read.table("/lustre/scratch117/cellgen/team218/MH/scRNASeqData/macosko.txt", header=T)
mac_data <- mac_data[rowSums(mac_data) > 0,]
colnames(mac_data) <- 1:length(mac_data[1,])
norm = t(t(mac_data)/colSums(mac_data)*1000000)

Mac <- getFeatures(NBumiConvertToInteger(mac_data), norm, name="mac")
Mac_dim <- dim(norm)
rm(mac_data); rm(norm);

#save.image(file="fsTime.RData")

return(c(m3d_t, hvg_t, nb_t, nbv_t, gini_t, cor_t, pca1_t, pca2_t))

# Make plot
source("Colour_Scheme.R")
TABLE = cbind(Deng, Biase, Xue, Fan, Goo, Ola, Blish, Shalek, Teic, Pollen, Kir, Lin, Uso, Mac)
rownames(TABLE) = c("M3D", "HVG", "NBDrop", "NBDisp", "Gini","Cor", "PCA (1-3)", "PCA (2-3)")

mat_size = rbind(Deng_dim, Biase_dim, Xue_dim, Fan_dim, Goo_dim, Ola_dim, Blish_dim, Shalek_dim, Teic_dim, Pollen_dim, Kir_dim, Lin_dim, Uso_dim, Mac_dim)

xes = apply(mat_size, 1, prod)
my_order = order(xes)
my_col = c(MM_col, hvg_1_col, Depth_col, NBVar_col, gini_col, cor_col, pca_1_col, pca_2_col)

png("FS_Time.png")
plot(1,1, col="white", xlim=c(min(xes), max(xes))/1000000, ylim=c(min(TABLE), max(TABLE)), xlab="ExprMat Size (millions)", ylab="Compute Time (s)")
for (i in 1:length(TABLE[,1])) {
	lines(xes[my_order]/1000000, TABLE[i,my_order], col=my_col[i], lwd=3)
	points(xes[my_order]/1000000, TABLE[i,my_order], col=my_col[i], pch=16, cex=1.75)
}
meth_order = order(-TABLE[,which(xes==max(xes))])
legend("topleft", rownames(TABLE)[meth_order], col=my_col[meth_order], lty=1, bty="n", lwd=2)
abline(h=60, col="grey35", lty=2)
text(15,60, "1 min", pos=3, col="grey50")
abline(h=60*30, col="grey35", lty=2)
text(175,60*30, "30 min", pos=1, col="grey50")
dev.off()
