#####################
#### Supp Fig 2a ####
#####################
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

# dat_ori from 1a.R
dat_ori$Status <- factor(dat_ori$Status , levels=c("CPE_positive", "CPE_negative","Family_member"))

# p-values for all 3 pairwise combinations (change "Status" accordingly)
ori_subset = dat_ori[dat_ori$Status == "CPE_negative" | dat_ori$Status == "Family_member",]
wilcox.test(ori_subset$V1 ~ ori_subset$Status)

# Generate boxplot
boxplot(V1~Status,data=dat_ori, ylab="PCoA 1 value", xlab="",
        col = c("indianred","cornflowerblue","lightgreen"),ylim=c(-0.7,0.48), yaxt="n")
axis(2, at = seq(-0.5, 0.4, by = 0.25))

#####################
#### Supp Fig 2b ####
#####################
# dat_ori, genus_raw from 1a.R
ori_corrWithPC1.df <- data.frame(Genus=genus_raw$Index, Sign=integer(224), Corr=double(224))
for (i in seq(224)){
  corrVal = cor(dat_ori$V1, t(genus_raw[i,-1]))
  # Find sign of correlation between abundance and PCoA 1
  ori_corrWithPC1.df[i,2] = (corrVal > 0) * 1
  # Find magnitude of correlation between abundance and PCoA 1
  ori_corrWithPC1.df[i,3] = abs(corrVal)
}

# Sort based on magnitude of correlation
ori_corrWithPC1.df[is.na(ori_corrWithPC1.df)] = -9
ori_corrWithPC1.df = ori_corrWithPC1.df[order(ori_corrWithPC1.df$Corr, decreasing=T),]
