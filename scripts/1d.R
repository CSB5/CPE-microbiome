#######################################
# Packages in use:
# ggplot2
#######################################

stDev = read.delim("relAbunStDev.txt")
stDev$SubjectType <- factor(stDev$SubjectType, levels = c("Index", "Family"))

# Change y=EnteroAbunStDev to other column names as necessary
p <- ggplot(stDev, aes(x=SubjectType, y=EnetroAbunStDev, fill=SubjectType)) + geom_violin() + 
  geom_boxplot(width=0.05, fill="white") + theme_classic() + theme(legend.position = "none") +
  ylim(0,50) + theme(axis.text.x=element_blank())
p + scale_fill_manual(values=c("#DA70D6","lightgreen"))

# Change y=EnteroAbunStDev to other column names as necessary
wilcox.test(stDev$EnteroAbunStDev ~ stDev$SubjectType, alternative = "two.sided")
