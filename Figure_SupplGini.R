# R3.3
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

Gini_plot <- function(x) {
        require("reldist")
        ginis <- apply(x, 1, gini)
        d <- rowMeans(x>0)
	plot(d, ginis, pch=16, xlab="Dropout Rate", ylab="Gini Index", col="grey50")
        reg <- lm(ginis~d) # almost perfect linear relation in UMI data
	abline(reg)
}


#### Main-Text Ola & Blischak ####

png("SupplFigure_Gini.png", width=8, height=6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,3,1))

source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Ola_SC.R")
Ola <- as.matrix(data[-spikes,]);
Gini_plot(Ola)
title(main="Kolodziejczk")

load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TeichCC_clean.RData")
Gini_plot(TeichCC_list$data[!(rownames(TeichCC_list$data) %in% pseudo),])
title(main="Buettner")

load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/pollen.rda")
Gini_plot(pollen)
title(main="Pollen")

source("../My_R_packages/M3D/R/NB_UMI.R"); require("matrixStats");
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Blishcak_UMI.R")
Blish <- counts;
Gini_plot(Blish)
title(main="Tung")

load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/kirschner.rda")
data_list = normalize_data(kirschner, labels = 1:length(kirschner[1,]), is.counts = FALSE)
Gini_plot(data_list$data)
title(main="Klein")
rm(kirschner)


load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/linnarsson.rda")
data_list = normalize_data(linnarsson, labels = 1:length(linnarsson[1,]), is.counts = FALSE)
Gini_plot(data_list$data)
title(main="Zeisel")
rm(linnarsson)

dev.off()
