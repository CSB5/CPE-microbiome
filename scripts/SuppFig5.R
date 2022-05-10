#####################
#### Supp Fig 5 ####
#####################

braycurtis <- read.delim("../data/SuppFig5_BrayCurtis.txt")
braycurtis$type = factor(braycurtis$type, levels=c("fam-fam","fam-pos","fam-neg0","fam-neg2"))
boxplot(bc_dist~type, data=braycurtis, ylab="", xlab="",
        col = c("lightgreen", "indianred", "dodgerblue4", "cornflowerblue"), xaxt="n",yaxt="n",
        ylim=c(0,1.3))
axis(2, at = seq(0, 1, by = 0.5))
axis(1, at = seq(0, 4, by = 1))

# Change "type" to calculate p-values of different pairs
bcSubset = braycurtis[braycurtis$type == "fam-fam" | braycurtis$type == "fam-pos",]
wilcox.test(bcSubset$bc_dist ~ bcSubset$type, alternative = "less")
