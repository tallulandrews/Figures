# R3.3
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

bg__dropout_plot_base_custom <- function (expr_mat, xlim = NA, suppress.plot=FALSE) {

        gene_info = bg__calc_variables(expr_mat);

        xes = log(gene_info$s)/log(10);
        put_in_order = order(xes);
        fancy <- densCols(xes, gene_info$p, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
        colours <-  colorRampPalette(c("#BBBBBB", "#444444"))(256) #grey->grey
        dens.col = colours[dens]

        if (!suppress.plot) {
                par(fg="black")
                if (!(sum(is.na(xlim)))) {
                        plot(xes,gene_info$p, main="", ylab="", xlab="", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1), cex=0.5)
                } else {
                        plot(xes,gene_info$p, main="", ylab="", xlab="", col = dens.col,pch=16, ylim=c(0,1), cex=0.5)
                }
                title(ylab="Dropout Rate", line=2)
                title(xlab="log10(expression)", line=2)
        }
        invisible(list(p=gene_info$p, s=gene_info$s, xes=xes, data=expr_mat, order=put_in_order));
}

MM_col = "black"
SCDE_col="magenta3"
ZIFA_col="blue"

custom_plot <- function(data, file_prefix=NA, known_groups=as.character(1:length(data[1,])), is_counts=TRUE, xlim=c(-1.5,6), textpos=-0.5) {
	data_list = normalize_data(data, labels = known_groups, is.counts = is_counts)
	require("RColorBrewer")
	norm = data_list$data

	if (!is.na(file_prefix)) {
		png(paste(file_prefix,"_CustomPlot.png", sep=""), width=6, height=6, units="in", res=300);
	}
	par(xpd=T)
        par(fg="black")
	BasePlot = bg__dropout_plot_base_custom(norm);
	if(length(grep("ERCC",rownames(data_list$data), ignore.case=TRUE)) > 0) {
		erccs = grep("ERCC",rownames(data_list$data), ignore.case=TRUE)
		points(log(BasePlot$s[erccs])/log(10),BasePlot$p[erccs], col="black", cex=1.5)
	}


	textMM <- tryCatch({	
	        MM = bg__fit_MM(BasePlot$p, BasePlot$s);
		lines(BasePlot$xes[BasePlot$order],MM$predictions[BasePlot$order],lty=1,lwd=3,col=MM_col);
		paste("MM   \tK  = ",sprintf("%.2f",signif(MM$K,3)),"\t\tSSr = ",MM$SSr,"\tSAr = ",MM$SAr,"\n", sep="")
		}, error = function(cond) {
			message(cond);
			return("MM did not converge.")
		}
	)
	textSCDE <- tryCatch({	
        	SCDE = bg__fit_logistic(BasePlot$p, BasePlot$s);
		lines(BasePlot$xes[BasePlot$order],SCDE$predictions[BasePlot$order],lty=2,lwd=3,col=SCDE_col);
		paste("Logi \tB0 = ",signif(SCDE$B0,4),"\t\tSSr = ",SCDE$SSr,"\tSAr = ",SCDE$SAr,"\n",
		         "     \t\tB1 = ",signif(SCDE$B1,4),"\n", sep="")
		}, error = function(cond) {
			message(cond);
			return("Logi did not converge.")
		}
	)
	# nested tryCatch
	textZIFA <- tryCatch({	
	        ZIFA = bg__fit_ZIFA(BasePlot$p, BasePlot$s);
		lines(BasePlot$xes[BasePlot$order],ZIFA$predictions[BasePlot$order],lty=3,lwd=3,col=ZIFA_col);
		paste("ZIFA \tl  = ",signif(ZIFA$lambda,3),"\t\tSSr = ",ZIFA$SSr,"\tSAr = ",ZIFA$SAr,"\n", sep="")
		}, error = function(cond) {
			message(cond);
			return("ZIFA did not converge.")
		}
	)
	if (!is.na(file_prefix)) {
		dev.off()
	}
	return(cbind(c(MM$SSr,SCDE$SSr,ZIFA$SSr),c(MM$SAr,SCDE$SAr,ZIFA$SAr),c(MM$K,SCDE$B1,ZIFA$lambda)))
}

Depth_col = "goldenrod1"
Norm_col = "purple"

NB_plot <- function(counts, size_factor=(colSums(counts)/median(colSums(counts)))) {
        norm <- NBumiConvertToInteger(t(t(counts)/size_factor));
        fit_adjust <- NBumiFitModel(counts);
        fit_basic <- NBumiFitBasicModel(norm);
        check_adjust <- NBumiCheckFit(counts, fit_adjust, suppress.plot=TRUE)
        check_basic <- NBumiCheckFit(norm, fit_basic, suppress.plot=TRUE)
        nc = fit_adjust$vals$nc
#       plot( log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), fit_adjust$vals$djs, col="white" )
#       arrows(log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), fit_adjust$vals$djs, 
#               log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), rowSums(check_adjust$exp_ps), col="navy",
#               length=0)
#       arrows(log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), fit_basic$vals$djs, 
#               log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), rowSums(check_basic$exp_ps), col="purple", 
#               length=0)
        plot( fit_adjust$vals$tjs/nc, fit_adjust$vals$djs/nc, col="black", pch=16, xlab="Expression", log="x", ylab= "Dropout Rate", cex=0.75)
#        points( fit_adjust$vals$tjs/nc, fit_basic$vals$djs/nc, col="black", pch=16, cex=0.75)
        points( fit_adjust$vals$tjs/nc, rowSums(check_adjust$exp_ps)/nc, col=Depth_col, pch=16, cex=0.5 )
        points( fit_adjust$vals$tjs/nc, rowSums(check_basic$exp_ps)/nc, col=Norm_col, pch=16, cex=0.5 )
        err_adj <- sum(abs(rowSums(check_adjust$exp_ps)/nc-fit_adjust$vals$djs/nc))
#        err_bas <- sum(abs(rowSums(check_basic$exp_ps)/nc-fit_basic$vals$djs/nc))
        err_bas <- sum(abs(rowSums(check_basic$exp_ps)/nc-fit_adjust$vals$djs/nc))
#        legend("topright", paste(c("Depth-Adjusted\nError:", "Normalized\nError:"), round(c(err_adj, err_bas)), c("\n","\n")), col=c("goldenrod1","purple"), pch=16, bty="n")
        return(c(err_adj, err_bas))
}


#### Main-Text Ola & Blischak ####

png("Figure_FittingModels.png", width=4.5*2, height=4.5, units="in", res=300)
par(mfrow=c(1,2))
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Ola_SC.R")
Ola <- as.matrix(data[-spikes,]);
Ola_stats <- custom_plot(Ola, is_counts=TRUE, xlim=c(-1.5,6), textpos=-0.5)
title("Zero-inflation models")
legend("topright", paste(c("M3Drop\nError:", "SCDE\nError:", "ZIFA\nError:"),round(Ola_stats[,2]), c("\n","\n","\n")), lty=c(1,2,3), col=c(MM_col, SCDE_col, ZIFA_col), bty="n", lwd=3, cex=0.5)

source("../My_R_packages/M3D/R/NB_UMI.R"); require("matrixStats");
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Blishcak_UMI.R")
Blish <- counts;
blish_stats <- NBumiCompareModels(Blish)
title("Negative Binomial Models")
dev.off()


#### Full-Transcript Supplementary ####

png("SupplFigure_FullTranscript.png", width=8, height=6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,3,1))
scaling=0.5
# Shalek
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/ShalekKO_clean.RData")
Shalek_stats <- custom_plot(ShalekKO_list$data, is_counts=TRUE, xlim=c(-1.5,6), textpos=-0.5)
title(main="Shalek (FPKM)", line=0.5); mtext("A",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
Deng_stats <- custom_plot(Deng_embyro_list$data, is_counts=TRUE, xlim=c(-1.5,6), textpos=-0.5)
title(main="Deng (CPM)", line=0.5); mtext("B",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)

# Ola

# Buettner
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TeichCC_clean.RData")
T_stats = custom_plot(TeichCC_list$data[!(rownames(TeichCC_list$data) %in% pseudo),],is_counts=TRUE)
title(main="Buettner (CPM)", line=0.5); mtext("C",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(TeichCC_list)

# Pollen
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/pollen.rda")
pollen_stats = custom_plot(pollen,is_counts=FALSE)
title(main="Pollen (FPKM)", line=0.5); mtext("D",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(pollen)

# Zhong
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/zhong.rda")
ultralow = which(rowMeans(zhong) < 10^-5)
zhong = zhong[-ultralow,]
rownames(zhong) = 1:length(zhong[,1])
Z_stats = custom_plot(zhong,is_counts=FALSE)
title(main="Biase (FPKM)", line=0.5); mtext("E",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(zhong)

par(mar=c(4.5,4,1,1))
barplot(cbind(Ola_stats[,2], Shalek_stats[,2],Deng_stats[,2],T_stats[,2],pollen_stats[,2],Z_stats[,2]), beside=T, col=c(MM_col, SCDE_col, ZIFA_col), 
	names=c("Kolo", "Shalek","Deng","Buettner","Pollen", "Biase"), ylab="Absolute Error", las=3)
mtext("F",side=2, line=2, cex=2*scaling, font=2, at=5500, las=1)
legend("topleft", fill=c(MM_col, SCDE_col, ZIFA_col), c("M3Drop", "SCDE", "ZIFA"), bty="n")
dev.off()

#### UMI-Tagged Supplementary ####
png("SupplFigure_UMITag.png", width=6, height=6, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,3,1))
scaling=0.5
# Kirschner
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/kirschner.rda")
data_list = normalize_data(kirschner, labels = 1:length(kirschner[1,]), is.counts = FALSE)
kir_stats = NB_plot(data_list$data)
title(main="Klein (UMI)", line=0.5); mtext("A",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(kirschner)

# Linnarsson
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/linnarsson.rda")
data_list = normalize_data(linnarsson, labels = 1:length(linnarsson[1,]), is.counts = FALSE)
lin_stats = NB_plot(data_list$data)
title(main="Zeisel (UMI)", line=0.5); mtext("B",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(linnarsson)

# Usoskin
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/usoskin.rda")
data_list = normalize_data(usoskin, labels = 1:length(usoskin[1,]), is.counts = FALSE)
uso_stats = NB_plot(NBumiConvertToInteger(data_list$data))
title(main="Usoskin (5'Seq)", line=0.5); mtext("C",side=2, line=2, cex=2*scaling, font=2, at=1.1, las=1)

par(mar=c(5,4,1,1))
barplot(cbind(blish_stats, kir_stats, lin_stats, uso_stats), beside=T, col=c(Depth_col,Norm_col), 
	names=c("Blischak\n5' UMI", "Klein\n3' UMI", "Zeisel\n5' UMI", "Usoskin\n5' Reads"), ylab="Absolute Error", las=3)
mtext("D",side=2, line=2, cex=2*scaling, font=2, at=1500, las=1)
legend("topleft", fill=c(Depth_col, Norm_col), c("Depth-Adjusted", "Normalized"), bty="n")

dev.off()
