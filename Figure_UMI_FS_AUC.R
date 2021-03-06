# R-3.3.0

source("../DE_vs_bulk/Load_Blishcak_UMI.R")




# DE
require("M3Drop")
fit <- NBumiFitModel(counts)
fs1 <- NBumiFeatureSelectionCombinedDrop(fit) 
fs4 <- NBumiFeatureSelectionHighVar(fit)
hvg <- BrenneckeGetVariableGenes(norm, fdr=2)
hvg_spike <- BrenneckeGetVariableGenes(norm, spikes=spikes, fdr=2)
m3d <- M3DropFeatureSelection(norm, mt_method = "fdr", mt_threshold = 2)
fs6 <- hvg$p.value; names(fs6) <- hvg$Gene;
fs8 <- hvg_spike$p.value; names(fs8) <- hvg_spike$Gene;
fs7 <- m3d$p.value; names(fs7) <- m3d$Gene;

# Others
fs9 <- read.table("../DE_vs_bulk/Blishcak_scFS_Gini.txt", header=T)
fs11 <- read.table("../DE_vs_bulk/Blishcak_scFS_Monoclepca.txt", header=T)
fs12 <- read.table("../DE_vs_bulk/Blishcak_scFS_negcor.txt", header=T)
fs14 <- read.table("../DE_vs_bulk/Blishcak_scFS_pca.txt", header=T)

Pos_Ctrl = unlist(read.table("../DE_vs_bulk/Blischak_bulk_LV_DESeq2_edgeR_TPs.txt", header=F))
Neg_Ctrl = unlist(read.table("../DE_vs_bulk/Blischak_bulk_LV_DESeq2_edgeR_TNs.txt", header=F))

source("~/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Quality_General.R")
source("Colour_Scheme.R")

sets = list(fs1,fs4,fs6, fs7, fs8, fs9, fs11, fs12, fs14)
my_FS_colours = c(Depth_col, NBVar_col, hvg_1_col, MM_col, hvg_2_col, gini_col, pca_1_col, cor_col, pca_2_col)
my_FS_names = c("NBDrop", "NBDisp", "HVG", "M3Drop", "HVG-E", "Gini", "PCA(2-3)", "Cor", "PCA(1-3)")

my_methods <- list()
my_recall <- vector(length=length(sets));
score <- vector(length=length(sets));
for(i in 1:length(sets)) {
        print(i)
        if (is.vector(sets[[i]])) {
                t = data.frame(Gene=names(sets[[i]]), p.value=sets[[i]])
        } else {
                t = data.frame(Gene=rownames(sets[[i]]), p.value=sets[[i]][,1])
        }
        my_methods[[i]] = Quality_AUC(t)
	thing <- my_methods[[i]];
	score[i] = thing$AUC
	my_recall[i] = Quality_TopRecall(t, n=2000)
}

png("Figure_UMI_FS.png", width=4*1.2, height=4, units="in", res=300)
layout(matrix(c(1,1,1,1,2,2,1,1,1,1,2,2), nrow=2, byrow=T))
par(mar=c(3.5,3.5,1,1))
QualityPlot(my_methods, paste(my_FS_names, "\n"), c(my_FS_colours), zoom=FALSE)
dev.off();

png("Figure_UMI_FS_Top.png", width=4.5*1.2, height=4.5, units="in", res=300)
leg_order = order(-score)
par(las=3)
QualityRecallPlot(my_recall[leg_order], my_FS_names[leg_order], my_FS_colours[leg_order])
dev.off();

### Recall/Features Plot
xes <- seq(from=10, to=2000, by=5)
curves <- matrix(0, ncol=length(xes), nrow=length(sets))
for(i in 1:length(sets)) {
	if (is.vector(sets[[i]])) {
                t = data.frame(Gene=names(sets[[i]]), p.value=sets[[i]])
        } else {
                t = data.frame(Gene=rownames(sets[[i]]), p.value=sets[[i]][,1])
        }
	recalls <- sapply(xes, function(x) {Quality_TopRecall(t, x)})
	curves[i,] = recalls
}

png("Figure_UMI_FS_TP_vs_TopX.png", width=4.5, height=4.5*1.2, units="in", res=300)
par(mar=c(4,4,0,1))
plot(1,1, col="white", xlim=c(0, max(xes)), ylim=c(0, max(curves)), xlab="Number of Features", ylab="True Positives")
for (i in 1:length(sets)) {
	lines(xes, curves[i,], col=my_FS_colours[i], lwd=2)
}
my_order = order(-curves[,length(curves[1,])])
legend("topleft", my_FS_names[my_order], col=my_FS_colours[my_order], lty=1, lwd=2, bty="n")
mtext("A", side=2, at=600, font=2, las=2, line=3, cex=1.75)
dev.off()
