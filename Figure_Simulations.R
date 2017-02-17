# R-3.3.0

# Set up
dir = "/lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/";
type = c("Umi","Full")
seeds = c(1001, 1234, 6789)
datasets = list(umi=c("blisch", "kirsch", "lin"),full=c("Ola","Buet","Pollen"))
case_names = c("DE","DVar","HVar")

require("M3Drop")
source("~/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Quality_General.R")
source("Colour_Scheme.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Other_FS_functions.R")
require("matrixStats")


# Separate Barplots each data type 
for (t in 1:length(type)) {
	print(t)
	files = paste(type[t], rep(datasets[[t]], each=length(seeds)), rep(seeds, times=length(datasets[[t]])), sep="_")
	files = paste(dir, files, ".rds", sep="");

	output1_table <- matrix(nrow=8, ncol=length(files)*length(case_names))
	output2_table <- matrix(nrow=8, ncol=length(files)*length(case_names))
	labels_set <- rep(datasets[[t]], each = length(seeds)*length(case_names));
	case_set <- rep(case_names, times=length(datasets[[t]])*length(seeds));

	for (f in 1:length(files)) {
		dataObj <- readRDS(files[f])
		l2fc = dataObj$truth
		names(l2fc) = 1:length(l2fc)
		for(count_type in (1:length(case_names))+2) {
			counts <- dataObj[[count_type]];
			print(dim(counts));
			if (t == "Full") {
				counts <- cbind(counts[,1:100], counts[,length(counts[1,])-1:100])
			}
			rownames(counts) <- names(l2fc)
			T = colSums(counts);
			norm <- counts/T*median(T);
			keep = rowSums(counts > 0) > 2;
			counts = counts[keep,];
			norm = norm[keep,];
			labels = 1:length(counts[1,]);

			# Truth
			truth_vals = l2fc[keep];
			truth_names = names(l2fc); truth_names = truth_names[keep];

			if (count_type == 5) {
				Pos_Ctrl = truth_names[truth_vals < -1*log(5)/log(2)] # fold change > 5 = TP -> only increases in variance.
			} else {
				Pos_Ctrl = truth_names[abs(truth_vals) > log(5)/log(2)] # fold change > 5 = TP
			}
			Neg_Ctrl = truth_names[abs(truth_vals) < 1] # fold change < 2 = TN

			# FS
			fit <- NBumiFitModel(NBumiConvertToInteger(counts))
			fs1 <- NBumiFeatureSelectionCombinedDrop(fit) 
			fs4 <- NBumiFeatureSelectionHighVar(fit)
			hvg <- BrenneckeGetVariableGenes(norm, fdr=2, fitMeanQuantile=0.5, suppress.plot=TRUE)
			m3d <- M3DropFeatureSelection(norm, mt_method = "fdr", mt_threshold = 2, suppress.plot=TRUE)
			fs6 <- hvg$p.value; names(fs6) <- hvg$Gene;
			fs7 <- m3d$p.value; names(fs7) <- m3d$Gene;
			fs9 <- Gini_FS(norm)
			fs11 <- Monocle2_pca_FS(counts, labels, pcs=c(2,3))
			fs12 <- Cor_FS_neg(norm) # requires high mem
			fs14 <- Monocle2_pca_FS(counts, labels, pcs=c(1,2,3))


			sets = list(fs1,fs4,fs6, fs7, fs9, fs11, fs12, fs14)
			my_methods <- list()
			my_recall <- vector(length=length(sets));
			score <- vector(length=length(sets));
			for(i in 1:length(sets)) {
			        #print(i)
			        if (is.vector(sets[[i]])) {
			                tab = data.frame(Gene=names(sets[[i]]), p.value=sets[[i]])
			        } else {
			                tab = data.frame(Gene=rownames(sets[[i]]), p.value=sets[[i]][,1])
			        }
			        my_methods[[i]] = Quality_AUC(tab)
				thing <- my_methods[[i]];
				score[i] = thing$AUC
				my_recall[i] = Quality_TopRecall(tab, n=2000)
				# Aggregate score across methods, types, files
				output1_table[i, (f-1)*length(case_names)+count_type-2] = score[i]
				output2_table[i, (f-1)*length(case_names)+count_type-2] = my_recall[i]
			}
	
		}
	}
	my_FS_colours = c(Depth_col, NBVar_col, hvg_1_col, MM_col, gini_col, pca_1_col, cor_col, pca_2_col)
	my_FS_names = c("NBDrop", "NBDisp", "HVG", "M3Drop", "Gini", "PCA(2-3)", "Cor", "PCA(1-3)")

	png(paste(type[t],"_Simulation_AUC.png", sep=""), width=5, height=5, units="in", res=300) # only fixed 6 Feb
	par(mar=c(3,4,1,1))
	bar_vals <- sapply(case_names, function(c){rowMeans(output1_table[,case_set==c])})
	bar_errs <- sapply(case_names, function(c){ sqrt(rowVars(output1_table[,case_set==c])/sum(case_set==c)) })
	xes <- barplot(bar_vals, beside=TRUE, col=my_FS_colours, ylab="AUC", ylim=c(0,max(bar_vals)+0.2))
	legend("topright", my_FS_names, fill=my_FS_colours, bty="n", ncol=3)
	arrows(as.vector(xes), as.vector(bar_vals), as.vector(xes), as.vector(bar_vals)+as.vector(bar_errs)*2, angle=90, len=0.1)
	arrows(as.vector(xes), as.vector(bar_vals), as.vector(xes), as.vector(bar_vals)-as.vector(bar_errs)*2, angle=90, len=0.1)
	abline(h=0.5, col="grey35", lty=2)
	dev.off()

	png(paste(type[t],"_Simulation_Recall.png", sep=""), width=5, height=5, units="in", res=300) # Only fixed 6 Feb
	par(mar=c(3,4,1,1))
	bar_vals <- sapply(case_names, function(c){rowMeans(output2_table[,case_set==c])})
	bar_errs <- sapply(case_names, function(c){ sqrt(rowVars(output2_table[,case_set==c])/sum(case_set==c)) })
	xes <- barplot(bar_vals, beside=TRUE, col=my_FS_colours, ylab="True Positives", ylim=c(0,max(bar_vals)*1.1), main="Top 2000 Genes")
	legend("topright", my_FS_names, fill=my_FS_colours, bty="n", ncol=3)
	arrows(as.vector(xes), as.vector(bar_vals), as.vector(xes), as.vector(bar_vals)+as.vector(bar_errs)*2, angle=90, len=0.1)
	arrows(as.vector(xes), as.vector(bar_vals), as.vector(xes), as.vector(bar_vals)-as.vector(bar_errs)*2, angle=90, len=0.1)
	dev.off()
}

