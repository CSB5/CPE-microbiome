#######################################
# Packages in use:
# vegan
# dendextend
# ggplot2
#######################################

# Function to obtain PCoA output
calc_pcoa <- function (data_raw, meta_status) {
  ff=vegdist(t(data_raw[,-1]), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
  brayMatrix=as.matrix(ff)
  
  cmds = cmdscale(ff, eig=TRUE)
  eigen = cmds$eig / sum(cmds$eig) * 100
  
  dat = (merge(cmds$points, meta_status, by.x=0, by.y="Index"))
  
  return (list("eigen"=eigen, "dat"=dat))
}

# Load data and calculate PCoA
genus_raw = read.table("g.metaphlan2.profile_merged.tsv", header = TRUE)
sample_status = read.delim("sample_status.txt")
pcoa_out = calc_pcoa(genus_raw, sample_status)
dat_ori = pcoa_out$dat
eigen = pcoa_out$eigen

# Compute hierarchical clustering
names(dat_ori)[1] <- "Index"
ori_dist_matrix = dist(dat_ori[,2:3], method = 'euclidean')
ori_hclust_avg <- hclust(ori_dist_matrix, method = 'average')
plot(ori_hclust_avg)
ori_cut_avg = cutree(ori_hclust_avg, k = 4)
plot(dat_ori$V1, dat_ori$V2, col = c("#0070C0", "#FFC000", "#C00000", "#000000")[ori_cut_avg],
     xaxt="n", yaxt="n", xlab="", ylab="")

# Plot dendrogram
dend = as.dendrogram(hclust(ori_dist_matrix, method = 'average'))
dend1 = color_branches(dend, k = 4, col = c("#000000","#FFC000","#0070C0","#C00000"),
                       groupLabels = c("IV","III","II","I"))
dend1 = color_labels(dend1, k = 4, col = c("#000000","#FFC000","#0070C0","#C00000"))
labels(dend1) = dat_ori$Index[as.numeric(labels(dend1))]
dend1 <- set(dend1, "labels_cex", 0.25)
plot(dend1)

# Add labels indicating dendrogram clusters
dat_ori$Status <- factor(dat_ori$Status , levels=c("CPE_positive", "CPE_negative","Family_member"))
dat_ori_clusts = dat_ori
dat_ori_clusts$clust = as.factor(ori_cut_avg)

# Plot PCoA
ggplot(dat_ori_clusts, aes(x=V1,y=V2))  +
  geom_density2d(color="#ABABAB") +
  geom_point(aes(col=Status, bg=Status, shape=clust), size = 3) + theme_classic() +
  scale_color_manual(values = c("indianred", "cornflowerblue", "lightgreen")) +
  scale_fill_manual(values = c("indianred", "cornflowerblue", "lightgreen")) +
  scale_shape_manual(values = c(21,22,23,24)) +
  labs(x=paste0('PCoA1 (',round(eigen[1], 1),'%)'), y=paste0('PCoA2 (',round(eigen[2], 1),'%)'))
