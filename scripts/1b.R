# Load data file
meta_diversity <- read.delim("Supplementary_file_1.txt", stringsAsFactors=TRUE)
meta_diversity$Colonization.status <- factor(meta_diversity_jan22$Colonization.status, levels = c("CPE_positive", "CPE_negative_0", "CPE_negative_1", "CPE_negative_2", "Family_member"))

# Remove the "1 month post-clearance" group for easier visualization
meta_diversity_no1 = meta_diversity[(meta_diversity$Colonization.status != "CPE_negative_1"), ]
meta_diversity_no1$Colonization.status <- factor(meta_diversity_no1$Colonization.status, levels = c("CPE_positive", "CPE_negative_0", "CPE_negative_2", "Family_member"))

# Plot
boxplot(Genus.level.Shannon.diversity~Colonization.status, data=meta_diversity_no1, ylab="", xlab="", 
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))

# Calculate Wilcoxon rank-sum test p-value (replace statusB as appropriate)
statusA = "CPE_positive"
statusB = "CPE_negative_0"
blkA = meta_diversity_no1$Genus.level.Shannon.diversity[meta_diversity_no1$Colonization.status == statusA]
blkB = meta_diversity_no1$Genus.level.Shannon.diversity[meta_diversity_no1$Colonization.status == statusB]
wilcox.test(blkA, blkB, alternative = "less")
