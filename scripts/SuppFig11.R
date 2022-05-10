#####################
#### Supp Fig 11 ####
#####################

stacked_ec_index <- read.delim("../data/SuppFig11_Ecoli_index_prop.txt")
stacked_ec_index = data.matrix(stacked_ec_index[,2:13])
row.names(stacked_ec_index) = c("GM","G","R","Y","B")
stacked_ec_index = stacked_ec_index[nrow(stacked_ec_index):1, ]
barplot(stacked_ec_index, col = c("#0070C0","#FFC000","#C00000","#BFBFBF","black"),
        ylim = c(0,1), border="white", space=0.04, font.axis=1)

stacked_kp_index <- read.delim("../data/SuppFig11_Kleb_index_prop.txt")
stacked_kp_index = data.matrix(stacked_kp_index[,2:13])
row.names(stacked_kp_index) = c("GM","G","R","Y","B")
stacked_kp_index = stacked_kp_index[nrow(stacked_kp_index):1, ]
barplot(stacked_kp_index, col = c("#0070C0","#FFC000","#C00000","#BFBFBF","black"),
        ylim = c(0,1), border="white", space=0.04, font.axis=1)

stacked_ec_fam <- read.delim("../data/SuppFig11_Ecoli_fam_prop.txt")
stacked_ec_fam = data.matrix(stacked_ec_fam[,2:13])
row.names(stacked_ec_fam) = c("GM","G","R","Y","B")
stacked_ec_fam = stacked_ec_fam[nrow(stacked_ec_fam):1, ]
barplot(stacked_ec_fam, col = c("#0070C0","#FFC000","#C00000","#BFBFBF","black"),
        ylim = c(0,1), border="white", space=0.04, font.axis=1)

stacked_kp_fam <- read.delim("../data/SuppFig11_Kleb_fam_prop.txt")
stacked_kp_fam = data.matrix(stacked_kp_fam[,2:13])
row.names(stacked_kp_fam) = c("GM","G","R","Y","B")
stacked_kp_fam = stacked_kp_fam[nrow(stacked_kp_fam):1, ]
barplot(stacked_kp_fam, col = c("#0070C0","#FFC000","#C00000","#BFBFBF","black"),
        ylim = c(0,1), border="white", space=0.04, font.axis=1)
