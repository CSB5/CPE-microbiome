#####################
#### Supp Fig 4 ####
#####################
div_noEntero_norm_Mar22 <- read.delim("../data/SuppFig4_NoEntero_Summary.txt")
sample_status = read.delim("../data/sample_status.txt")
divFinal_Mar22 = (merge(sample_status, div_noEntero_norm_Mar22, by.y="Index"))
# Remove "1 month post clearance" designation as that is not used in the plots
divFinalNoPostOne_Mar22 = divFinal_Mar22[divFinal_Mar22$SubStatus != "CPE_negative_1",]
divFinalNoPostOne_Mar22$SubStatus <- factor(divFinalNoPostOne_Mar22$SubStatus,
        levels=c("CPE_positive","CPE_negative_0","CPE_negative_2","Family_member"))

# Plot Figure 4a
boxplot(noEntero_Shannon~SubStatus,data=divFinalNoPostOne_Mar22, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))

# Plot Figure 4b
boxplot(noEntero_Count~SubStatus,data=divFinalNoPostOne_Mar22, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(0,90))
axis(2, at = seq(0, 80, by = 20))
