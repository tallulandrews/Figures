# R3.3
require("matrixStats")
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

custom_plot <- function(data, file_prefix=NA, known_groups=as.character(1:length(data[1,])), is_counts=TRUE, xlim=c(-1.5,6), textpos=-0.5, suppress.plot=FALSE) {
	data_list = normalize_data(data, labels = known_groups, is.counts = is_counts)
	require("RColorBrewer")
	norm = data_list$data

	if (!is.na(file_prefix)) {
		png(paste(file_prefix,"_CustomPlot.png", sep=""), width=6, height=6, units="in", res=300);
	}
	par(xpd=T)
        par(fg="black")
	BasePlot = bg__dropout_plot_base_custom(norm, suppress.plot=suppress.plot);
	if(length(grep("ERCC",rownames(data_list$data), ignore.case=TRUE)) > 0) {
		erccs = grep("ERCC",rownames(data_list$data), ignore.case=TRUE)
		if (!suppress.plot) {
		points(log(BasePlot$s[erccs])/log(10),BasePlot$p[erccs], col="black", cex=1.5)
		}
	}


	textMM <- tryCatch({	
	        MM = bg__fit_MM(BasePlot$p, BasePlot$s);
		if (!suppress.plot) {
		lines(BasePlot$xes[BasePlot$order],MM$predictions[BasePlot$order],lty=1,lwd=3,col=MM_col);
		}
		paste("MM   \tK  = ",sprintf("%.2f",signif(MM$K,3)),"\t\tSSr = ",MM$SSr,"\tSAr = ",MM$SAr,"\n", sep="")
		}, error = function(cond) {
			message(cond);
			return("MM did not converge.")
		}
	)
	textSCDE <- tryCatch({	
        	SCDE = bg__fit_logistic(BasePlot$p, BasePlot$s);
		if (!suppress.plot) {
		lines(BasePlot$xes[BasePlot$order],SCDE$predictions[BasePlot$order],lty=2,lwd=3,col=SCDE_col);
		}
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
		if (!suppress.plot) {
		lines(BasePlot$xes[BasePlot$order],ZIFA$predictions[BasePlot$order],lty=3,lwd=3,col=ZIFA_col);
		}
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
Norm_col = "mediumpurple1"

NB_plot <- function(counts, size_factor=(colSums(counts)/median(colSums(counts))), suppress.plot=FALSE) {
        norm <- NBumiConvertToInteger(t(t(counts)/size_factor));
        fit_adjust <- NBumiFitModel(counts);
        fit_basic <- NBumiFitBasicModel(norm);
        check_adjust <- NBumiCheckFitFS(counts, fit_adjust, suppress.plot=TRUE)
        check_basic <- NBumiCheckFitFS(norm, fit_basic, suppress.plot=TRUE)
        nc = fit_adjust$vals$nc
#       plot( log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), fit_adjust$vals$djs, col="white" )
#       arrows(log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), fit_adjust$vals$djs, 
#               log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), rowSums(check_adjust$exp_ps), col="navy",
#               length=0)
#       arrows(log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), fit_basic$vals$djs, 
#               log(fit_adjust$vals$tjs/fit_adjust$vals$nc)/log(10), rowSums(check_basic$exp_ps), col="purple", 
#               length=0)
if (!suppress.plot) {
        plot( log(fit_adjust$vals$tjs/nc)/log(10), fit_adjust$vals$djs/nc, col="black", pch=16, xlab="", ylab= "", cex=0.75)
	title( ylab= "Dropout Rate", xlab="log10(expression)", line=2)
#        points( fit_adjust$vals$tjs/nc, fit_basic$vals$djs/nc, col="black", pch=16, cex=0.75)
        points( log(fit_adjust$vals$tjs/nc)/log(10), rowSums(check_adjust$exp_ps)/nc, col=Depth_col, pch=16, cex=0.5 )
        points( log(fit_basic$vals$tjs/nc)/log(10), rowSums(check_basic$exp_ps)/nc, col=Norm_col, pch=16, cex=0.5 )
}
        err_adj <- sum(abs(rowSums(check_adjust$exp_ps)/nc-fit_adjust$vals$djs/nc))
#        err_bas <- sum(abs(rowSums(check_basic$exp_ps)/nc-fit_basic$vals$djs/nc))
        err_bas <- sum(abs(rowSums(check_basic$exp_ps)/nc-fit_adjust$vals$djs/nc))
#        legend("topright", paste(c("Depth-Adjusted\nError:", "Normalized\nError:"), round(c(err_adj, err_bas)), c("\n","\n")), col=c("goldenrod1","purple"), pch=16, bty="n")
        return(c(err_adj, err_bas))
}


#### Main-Text Ola & Blischak ####

#### Full-Transcript Supplementary ####
scaling=0.5

# Pollen
png("Poster_FullTranscript.png", width=6, height=6, units="in", res=500)
par(mar=c(4,4,3,1))
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/pollen.rda")
pollen_stats_MM = custom_plot(pollen,is_counts=FALSE)
pollen_stats_NB = NB_plot(NBumiConvertToInteger(pollen), suppress.plot=TRUE)
title(main="Pollen (FPKM)", line=0.5); #mtext("D",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(pollen)
legend("bottomleft", c("M3Drop","SCDE","ZIFA"), fill=c(MM_col, SCDE_col, ZIFA_col))
dev.off()

#### UMI-Tagged Supplementary ####
scaling=0.5

# Linnarsson
png("Poster_UMITag.png", width=6, height=6, units="in", res=500)
par(mar=c(4,4,3,1))
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/clustools-master/data/linnarsson.rda")
data_list = normalize_data(linnarsson, labels = 1:length(linnarsson[1,]), is.counts = FALSE)
lin_stats_NB = NB_plot(data_list$data)
lin_stats_MM <- custom_plot(data_list$data, is_counts=TRUE, xlim=c(-1.5,6), textpos=-0.5, suppress.plot=TRUE)
title(main="Zeisel (UMI)", line=0.5); #mtext("B",side=2, line=2, cex=2*scaling, font=2, at=1.075, las=1)
rm(linnarsson)
legend("bottomleft", c("DANB", "Base NB"), fill=c(Depth_col, Norm_col))
dev.off()
