library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
# # # 
tree <- read.tree('C:/Users/Yatir/OneDrive/code/general/trees/review_family.nw')
rnk_phy_col_df <- read.csv('C:/Users/Yatir/OneDrive/code/general/review_family_phylogenetic_mtdta.tsv', header=TRUE, sep='\t')
fig <- ggtree(tree, layout='circular', size=.27)
fig <- fig %<+% rnk_phy_col_df
# fig1 <- fig1 + geom_tippoint(alpha=0.75, size=1.5, show.legend=TRUE)
fig <- fig + geom_tippoint(aes(colour=phyla_class), alpha=0.75, size=1.8, show.legend=TRUE)
# fig <- fig + geom_text2(aes(label=node, subset=!isTip), hjust=-.3, size=0.8) # adding nodes number in the figure
# fig <- fig + geom_tiplab(aes(colour=phyla_class), 
fig <- fig + geom_tiplab(colour='black', 
                        geom='text', fontface=3,  size=1, 
                        offset=0.01, align=FALSE, 
                        linesize=1, linetype='dotted', 
                        show.legend=FALSE) # leaf label (family/order) next to phylum/class
ggsave("C:/Users/Yatir/OneDrive/code/general/all_bacteria_tree.svg")
print(fig)
