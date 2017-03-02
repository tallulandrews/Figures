# R-3.3.0

# Set up
dir = "/lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/";
type = rep(c("Umi","Full"), each=3)
seeds = c(1001, 1234, 6789)
datasets = c("blisch", "kirsch", "lin", "Ola","Buet","Pollen")
data_names = c("Tung","Klein","Zeisel","Kolo","Buettner","Pollen")
case_names = c("DE","DVar","HVar")

require("M3Drop")
source("~/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Quality_General.R")
source("Colour_Scheme.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Other_FS_functions.R")
require("matrixStats")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")


# Raw
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/kirschner.rda")
kir_list = normalize_data(kirschner, labels = 1:length(kirschner[1,]), is.counts = FALSE)

load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/linnarsson.rda")
lin_list = normalize_data(linnarsson, labels = 1:length(linnarsson[1,]), is.counts = FALSE)

source("../My_R_packages/M3D/R/NB_UMI.R"); require("matrixStats");
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Blishcak_UMI.R")
Blish <- counts;

source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Ola_SC.R")
Ola <- as.matrix(data[-spikes,]);

load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TeichCC_clean.RData")

load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/pollen.rda")

raw <- list(Blish, kir_list$data, lin_list$data, as.matrix(Ola), TeichCC_list$data, pollen)


# Separate Barplots each data type 
files = paste(rep(type, each=length(seeds)), rep(datasets, each=length(seeds)), rep(seeds, times=length(datasets)), sep="_")
files = paste(dir, files, ".rds", sep="");

calc_stats <- function (d) {
	detected <- colSums(d > 0)
	dropr <- rowMeans(d > 0);
	v <- rowVars(d);
	m <- rowMeans(d);
	t <- colSums(d);
	return(list(di = detected, dj=1-dropr, v=v, mu=m, t=t))
}

png("SupplSimulations1.png", width=6*8/6, height=4*8/6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(3,3,2,1))
par(cex=0.95)
for (i in 1:length(raw)) {
	raw_data <- raw[[i]]
	my_files <- files[(3*(i-1)+1):(3*(i-1)+3)]
	d1 = readRDS(my_files[1])
	d2 = readRDS(my_files[2])
	d3 = readRDS(my_files[3])

	raw_stats = calc_stats(as.matrix(raw_data))
	sim_stats = list(calc_stats(d1$de[abs(d1$truth) < 1,]),calc_stats(d1$dv[abs(d1$truth) < 1,]), calc_stats(d1$hv[abs(d1$truth) < 1,]), 
			calc_stats(d2$de[abs(d2$truth) < 1,]),calc_stats(d2$dv[abs(d2$truth) < 1,]), calc_stats(d2$hv[abs(d2$truth) < 1,]),
			calc_stats(d3$de[abs(d3$truth) < 1,]),calc_stats(d3$dv[abs(d3$truth) < 1,]), calc_stats(d3$hv[abs(d3$truth) < 1,])
			)

	plot(raw_stats$mu, raw_stats$dj, log="x", pch=16, col="grey", xlab="", ylab="", main="")
	title(xlab="Mean", ylab="Dropout Rate", line=2)
	title(main=data_names[i], line=1)
	for(thing in 1:length(sim_stats)) {
		points(sim_stats[[thing]]$mu, sim_stats[[thing]]$dj, col="cornflowerblue");
	}
	if (i == 1) {legend("topright", c("Sim","Obs"), fill=c("cornflowerblue", "grey"), bty="n")}
}
dev.off()

png("SupplSimulations2.png", width=6*8/6, height=4*8/6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(3,3,2,1))
par(cex=0.95)
for (i in 1:length(raw)) {
        raw_data <- raw[[i]]
        my_files <- files[(3*(i-1)+1):(3*(i-1)+3)]
        d1 = readRDS(my_files[1])
        d2 = readRDS(my_files[2])
        d3 = readRDS(my_files[3])

        raw_stats = calc_stats(as.matrix(raw_data))
        sim_stats = list(calc_stats(d1$de[abs(d1$truth) < 1,]),calc_stats(d1$dv[abs(d1$truth) < 1,]), calc_stats(d1$hv[abs(d1$truth) < 1,]),
                        calc_stats(d2$de[abs(d2$truth) < 1,]),calc_stats(d2$dv[abs(d2$truth) < 1,]), calc_stats(d2$hv[abs(d2$truth) < 1,]),
                        calc_stats(d3$de[abs(d3$truth) < 1,]),calc_stats(d3$dv[abs(d3$truth) < 1,]), calc_stats(d3$hv[abs(d3$truth) < 1,])
                        )

	plot(raw_stats$mu, raw_stats$v, log="xy", pch=16, col="grey", xlab="", ylab="", main="")
	title(xlab="Mean", ylab="Variance", line=2)
	title(main=data_names[i], line=1)
        for(thing in 1:length(sim_stats)) {
                points(sim_stats[[thing]]$mu, sim_stats[[thing]]$v, col="cornflowerblue");
        }
	if (i == 1) {legend("topleft", c("Sim","Obs"), fill=c("cornflowerblue", "grey"), bty="n")}
}
dev.off()

