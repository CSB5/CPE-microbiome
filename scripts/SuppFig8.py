#####################
#### Supp Fig 8 ####
#####################

ARG_abundance <- read.csv("../data/SuppFig8_ARGAbundance.csv")
sample_status = read.delim("sample_status.txt")
ARG_merge = (merge(sample_status, ARG_abundance, by.y="Index"))
# Remove "1 month post clearance" designation as that is not used in the plots
ARG_mergeNoPostOne = ARG_merge[ARG_merge$SubStatus != "CPE_negative_1",]
ARG_mergeNoPostOne$SubStatus <- factor(ARG_mergeNoPostOne$SubStatus,
          levels=c("CPE_positive","CPE_negative_0","CPE_negative_2","Family_member"))

# Supp Fig 8a
boxplot(log10(ARG)~SubStatus,data=ARG_mergeNoPostOne, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n", yaxt="n",
        ylim=c(-4,-1))
axis(2, at = seq(-4, -1, by = 1))

# Supp Fig 8b
boxplot(log10(cARG_cDB_abundance)~SubStatus,data=ARG_mergeNoPostOne, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n", yaxt="n",
        ylim=c(-7.2,-1.3))
axis(2, at = seq(-7, -2, by = 1))

# Supp Fig 8c
ARG_mergeNoPostOne$ARG_diff = ARG_mergeNoPostOne$ARG - ARG_mergeNoPostOne$cARG_cDB_abundance
boxplot(log10(ARG_diff)~SubStatus,data=ARG_mergeNoPostOne, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n", yaxt="n",
        ylim=c(-4,-1))
axis(2, at = seq(-4, -1, by = 1))

# Change "status2" to find p-values of different combinations
status1 = "CPE_positive"
#status2 = "CPE_negative_0"
#status2 = "CPE_negative_2"
#status2 = "Family_member"
argSubset = ARG_mergeNoPostOne[ARG_mergeNoPostOne$SubStatus == status1 | ARG_mergeNoPostOne$SubStatus == status2,]
wilcox.test(argSubset$ARG ~ argSubset$SubStatus)
