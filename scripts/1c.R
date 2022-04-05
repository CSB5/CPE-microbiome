library(dplyr)
library(readr)
library(ggplot2)

#sed 's/UNINTEGRATED.//' ANALYSIS/HUMANN2/merge_pathway_status.tsv.lefse.ana > ANALYSIS/HUMANN2/merge_pathway_status.tsv.lefse.ana.f
dat <- read_tsv('../ANALYSIS/HUMANN2/merge_pathway_status.tsv.lefse.ana', col_names = F) %>%
    filter(X3!='-', X4>3)

#arrange(dat, group, X1, LDA) %>% 
#    mutate(name=factor(X2, levels=X2, ordered = TRUE)) -> plot.dat

                                        #Get the minus and the scalling
dat$X4 = dat$X4 - 3
dat[dat[,3]=="Negative",4] = -1 * dat[dat[,3]=="Negative",4] 



pdf(file="../PLOT_PAPER/lefse.pdf", width=15,height=4)

ggplot(dat, aes(x=X1, y=X4, fill=X3, label=X1)) + 
    geom_bar(stat='identity') + 
    geom_text(aes(y=0, x=X1), hjust=ifelse(dat$X4<0, -0, 1), nudge_y = -sign(dat$X4)*0.1) + 
    coord_flip() + 
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1), labels = c('-4','-3.5','-3/3','3.5','4'), limits = c(-1.5,1.5), name = 'LDA score') + 
    scale_fill_manual(values=c('cornflowerblue','indianred','darkgreen')) + 
    labs(x=NULL) + 
    theme_classic() + 
    facet_grid(X1~., scale='free_y',  space = "free_y") + 
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'cm'), legend.position = 'left', 
          axis.text.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0, "lines"))

dev.off()
