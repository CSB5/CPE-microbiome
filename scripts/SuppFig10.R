#####################
#### Supp Fig 10 ####
#####################
# Fig 10b
oneclass_overlap = read.delim("../data/SuppFig10_OneClass_Iso_Overlap.txt")
boxplot(overlap_proportion ~ species, data=oneclass_overlap, lwd = 2, outline=FALSE, ylim=c(0.975,1),
        col=c("#060086", "#80170E"), yaxt="n", ylab="", xlab="")
axis(2, at = seq(0.98, 1, by = 0.01))
stripchart(overlap_proportion ~ species, vertical = TRUE, data=oneclass_overlap, method = "jitter", add = TRUE,
           pch = 20, col = "#AAAAAA")

# Fig 10a
Ecoli_hist_bin_100 = read.delim("../data/SuppFig10_Ec_hist_bin100.txt")

# Multiple strains: 0506-T-V02
lastBinMod = as.numeric(Ecoli_hist_bin_100[71,5:104])
lastBinMod[100] = 12000
barplot(lastBinMod, ylim=c(0,12000), yaxt="n", col="#C00000")
axis(2, at = seq(0, 10000, by = 5000))
axis(1, at=c(0,60,120))

# Two strains: 0506-T-V03
lastBinMod = as.numeric(Ecoli_hist_bin_100[72,5:104])
lastBinMod[100] = 12000
barplot(lastBinMod, ylim=c(0,12000), yaxt="n", col="#FFC000")
axis(2, at = seq(0, 10000, by = 5000))
axis(1, at=c(0,60,120))

# One strain: 0506-T-V04
lastBinMod = as.numeric(Ecoli_hist_bin_100[73,5:104])
lastBinMod[100] = 12000
barplot(lastBinMod, ylim=c(0,12000), yaxt="n", col="#0070C0")
axis(2, at = seq(0, 10000, by = 5000))
axis(1, at=c(0,60,120))
