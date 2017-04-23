source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Simulations/Simulations_Functions.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/DE_vs_bulk/Load_Blishcak_UMI.R")

original_data = counts

require("M3Drop")
fit <- NBumiFitModel(original_data);
Tis = fit$vals$tis
Tis_gamma = fit_gamma(Tis)
Mjs = fit$vals$tjs/fit$vals$nc
Mjs_norm = c(mean(log(Mjs)/log(10)), sd(log(Mjs)/log(10)))
mean2disp_coeffs = NBumiFitDispVsMean(fit, suppress.plot=TRUE)
min_mean = 10^-5;

mu_raw = c(0.1, 1, 10)

set.seed(1001)

n_cells = 1000000

c_depths = round(rgamma(n_cells, shape=Tis_gamma$shape, scale=Tis_gamma$scale));

mus <-(mu_raw %*% t(c_depths)/sum(c_depths))*n_cells
disp_size <- exp(log(rowMeans(mus))*mean2disp_coeffs[2]+mean2disp_coeffs[1])

base <- sapply(1:3, function(i) {sapply(mus[i,], function(m) {rnbinom(1, mu=m, size=disp_size[i])})})

truth_var = c( var(base[,1]), var(base[,2]), var(base[,3]) )
truth_drop = c( sum(base[,1]==0), sum(base[,2]==0), sum(base[,3]==0) )/n_cells

sizes = c(50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200)
n_rep = 50

outTABLElow_var = matrix(0, ncol=length(sizes), nrow=n_rep)
outTABLEmed_var = matrix(0, ncol=length(sizes), nrow=n_rep)
outTABLEhig_var = matrix(0, ncol=length(sizes), nrow=n_rep)
outTABLElow_drop = matrix(0, ncol=length(sizes), nrow=n_rep)
outTABLEmed_drop = matrix(0, ncol=length(sizes), nrow=n_rep)
outTABLEhig_drop = matrix(0, ncol=length(sizes), nrow=n_rep)

for (j in 1:length(sizes)) {
	s = sizes[j]
	for (i in 1:n_rep) {
		set = sample(1:length(base[,1]), size=s)
		diff = (colVars(base[set,])-truth_var)/truth_var
		outTABLElow_var[i,j] = diff[1]
		outTABLEmed_var[i,j] = diff[2]
		outTABLEhig_var[i,j] = diff[3]
		diff2 = (colSums(base[set,] == 0)/length(set)-truth_drop)/truth_drop
		outTABLElow_drop[i,j] = diff2[1]
                outTABLEmed_drop[i,j] = diff2[2]
                outTABLEhig_drop[i,j] = diff2[3]

	}
}

colnames(outTABLElow_var) = sizes
colnames(outTABLEmed_var) = sizes
colnames(outTABLEhig_var) = sizes
colnames(outTABLElow_drop) = sizes
colnames(outTABLEmed_drop) = sizes
colnames(outTABLEhig_drop) = sizes


#boxplot(outTABLEmed_var, ylim=c(-0.7,0.7))
#boxplot(outTABLEmed_drop, ylim=c(-0.7,0.7))

data = c(as.vector(outTABLEmed_var), as.vector(outTABLEmed_drop))
size_lab = rep(c(sizes, sizes), each=n_rep)
type_lab = rep(c("var","drop"), each=n_rep*length(sizes))

source("Colour_Scheme.R")
drop_col = Depth_col
var_col = "red3"#NBVar_col

png("Sim_Error_Drop_Var.png", width=5, height=5, units="in", res=300)
thing = boxplot(data~type_lab*size_lab, col=c(drop_col, var_col), xlab="# Cells", ylab="Sampling Error", xaxt="n")
axis(1, at=seq(from=1, to=length(sizes)*2, by=2)+0.5, labels=sizes, las=2)

legend("topright", c("Dropout Rate", "Variance"), fill=c(drop_col, var_col), bty="n")
dev.off()
