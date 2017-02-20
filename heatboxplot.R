silhouette_boxplot <- function(x, col="grey50", ylim=c(min(x),max(x)), main="") {
	par(mar=c(0.0.5,2,0.5))
	boxplot(x, ylim=ylim, col=col, border=col, staplewex=0, outline=FALSE, xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=FALSE)

}
colorbar_plot <- function(colour_vec, horiz=FALSE) {
	if (!horiz) {
		image(rbind(1:length(colour_vec)), 
			col = colour_vec, axes = FALSE)
	} else {
		image(cbind(1:length(colour_vec)), 
			col = colour_vec, axes = FALSE)
	}
}

get_height_for_k <- function(dendro, k) {
	height_out <- -1;
        curr_k <- 1;
        dendro_list <- list(dendro)
        dendro_heights <- attr(dendro, "height")
        while( curr_k < k ){
                to_split <- which(dendro_heights == max(dendro_heights))
                to_split_dendro <- dendro_list[[to_split]]
                to_split_height <-  dendro_heights[to_split]

                children <- as.list(to_split_dendro)
                for (i in 1:length(children)) {
                        dendro_heights <- c(dendro_heights,attr(children[[i]],"height"))
                        dendro_list[[length(dendro_list)+1]] <- children[[i]]
                }
                # Remove to split
                dendro_list[to_split] <- NULL
		height_out = mean(dendro_heights[to_split], 
				max(dendro_heights[-to_split]));
                dendro_heights <- dendro_heights[-to_split]
                curr_k <- curr_k-1+length(children)
        }
	return(height_out)
}



cut_dendro_plot <- funcion(dendro, k, horiz=TRUE, full=TRUE) {
	dendro <- as.dendrogram(dendro);
	if (!full) { dendro <- cut(dendro, k=k) }
	plot(dendro, axes=FALSE, leaflab="none", horiz=horiz)

	if (full) {
		h <- get_height_for_k(dendro, k)
		abline(v=h, lty=2, col="grey40")
	}
}



epic_dendro_boxplots <- function(x, distfun = dist, hclustfun = hclust,
	k, markers, marker_col, vbar_labels, vbar_cols) {
	# Layout
	nrow = k; ncol = length(markers);
        lhei <- c(rep(box_height, times=k))
        lwid <- c(dendro_width, bar_width, rep(box_width, times=ncol))
        lmat <- matrix(c(rep(1,times=k),rep(2,times=k),3:(ncol*k)), 
			ncol=ncol, byrow=FALSE)

	# Create Dendrogram
        hcc <- hclustfun(distfun( t(x) ))
	cut_dendro_plot(hcc, k)

	# Create colour bar
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
	if (!is.factor(vbar_labels)) {
		vbar_labels <- factor(vbar_labels);
	}
	my_col <- vbar_cols[vbar_labels[colInd]]
	colorbar_plot(my_col)

	# Create marker boxes
	groups <- cutree(ddc, k=k)
	groups <- groups[colInd]
	x <- x[,colInd]
	g_order <- unique(groups)
	for (m in markers) {
		ylim = c(min(x[rownames(x) == m,]), max(x[rownames(x) == m,]))
		for (g in g_order) {
			dat <- x[rownames(x) == m, groups == g]
			silhouette_boxplot(dat, col=marker_col[markers==m], ylim=ylim)
			if (g == g_order[1]) {
				title(main=m)
			}
		}
	}
}
