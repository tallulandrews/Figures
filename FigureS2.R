# R3.3
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

#### Main-Text Ola & Blischak ####

# Ola
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Ola_SC.R")
Ola <- as.matrix(data[-spikes,]);
Ola_stats = list(counts = colSums(Ola), detect = colSums(Ola > 0))

source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Blishcak_UMI.R")
Blish <- counts;
Blish_stats = list(counts = colSums(Blish), detect = colSums(Blish > 0))

# Shalek
#load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/ShalekKO_clean.RData")
#Shalek_stats = list(counts = colSums(ShalekKO_list$data), detect = colSums(ShalekKO_list$data > 0))
#rm(ShalekKO_list)

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
Deng_stats = list(counts = colSums(Deng_embyro_list$data), detect = colSums(Deng_embyro_list$data > 0))
rm(Deng_embyro_list)

# Buettner
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TeichCC_clean.RData")
T_stats = list(counts = colSums(TeichCC_list$data), detect = colSums(TeichCC_list$data > 0))
rm(TeichCC_list)

# Pollen
#load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/pollen.rda")
#pollen_stats = list(counts = colSums(pollen), detect = colSums(pollen > 0))
#rm(pollen)

# Zhong
#load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/zhong.rda")
#ultralow = which(rowMeans(zhong) < 10^-5)
#zhong = zhong[-ultralow,]
#rownames(zhong) = 1:length(zhong[,1])
#zhong_stats = list(counts = colSums(zhong), detect = colSums(zhong > 0))
#rm(zhong)

# Kirschner
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/kirschner.rda")
data_list = normalize_data(kirschner, labels = 1:length(kirschner[1,]), is.counts = FALSE)
kir_stats = list(counts = colSums(data_list$data), detect = colSums(data_list$data > 0))
rm(kirschner)

# Linnarsson
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/linnarsson.rda")
data_list = normalize_data(linnarsson, labels = 1:length(linnarsson[1,]), is.counts = FALSE)
lin_stats = list(counts = colSums(data_list$data), detect = colSums(data_list$data > 0))
rm(linnarsson)

# Usoskin
#load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/usoskin.rda")
#data_list = normalize_data(usoskin, labels = 1:length(usoskin[1,]), is.counts = FALSE)
#uso_stats = list(counts = colSums(data_list$data), detect = colSums(data_list$data > 0))

###### Figure ######
png("FigureS2.png", width=6, height=5, units="in", res=300)
par(cex=1)
par(mar=c(4,4,3,0))

full_col = c("#bdd7e7", "#6baed6", "#2171b5")
 umi_col = c("#bae4b3", "#74c476", "#238b45")

xes = c(lin_stats[[1]], kir_stats[[1]], Blish_stats[[1]], T_stats[[1]], Deng_stats[[1]], Ola_stats[[1]])
yes = c(lin_stats[[2]], kir_stats[[2]], Blish_stats[[2]], T_stats[[2]], Deng_stats[[2]], Ola_stats[[2]])
colours = rep(c(umi_col, full_col), times=c(length(lin_stats[[1]]), length(kir_stats[[1]]), length(Blish_stats[[1]]), length(T_stats[[1]]), length(Deng_stats[[1]]), length(Ola_stats[[1]])))

plot(xes, yes, col=colours, pch=18, log="xy", xlab="Total Counts", ylab="Detected Genes")
legend("bottomright", bty="n", pch=16, col=c(umi_col, full_col), c("Zeisel", "Klein", "Tung", "Buettner", "Deng", "Kolo"), ncol=2)

dev.off()
