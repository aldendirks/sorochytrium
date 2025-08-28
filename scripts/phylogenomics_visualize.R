# Load libraries
library(tidyverse)
library(treedataverse)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(cowplot)
library(aplot)
library(reshape2)
library(ggnewscale)
library(phytools)
library(ggimage)
library(castor)
library(ggrepel)
library(rsvg)
library(viridis)
library(nlme)
library(geiger)

##########################################################################################

# Consensus concatenation tree
tree <- read.tree("results/phylogenomics/concat.contree")
tree <- root_at_midpoint(tree)
tree$root.edge <- 0.01

# Concordance factors tree
tree.cf <- read.iqtree("results/phylogenomics/concord.cf.tree")
tree.cf@phylo <- root_at_midpoint(tree.cf@phylo)

# Concordance factors dataset
cf.df <- read.table("results/phylogenomics/concord.cf.stat", header=TRUE)
colnames(cf.df)[1] <- "node"
cf.df$node <- cf.df$node + 1

# Show tree node numbers
cladogram.nodes.plot <- ggtree(tree, branch.length="none") + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  geom_tiplab()
cladogram.nodes.plot
cladogram.nodes.plot2 <- ggtree(tree.cf, branch.length="none") + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  geom_tiplab()
cladogram.nodes.plot2

# Adjust tibble size for mapping onto tree
# By making it a tibble then a tree again, it expands the dataframe to include NAs. What's important is it becomes the same size as the tree, or rather same dimensions as tree tibble
tree.cf <- as_tibble(tree.cf)
tree.cf$label2 <- NA
tree.cf$label2[15:24] <- tree.cf$label[15:24]
tree.cf <- tree.cf %>% 
  separate(label2, into = c("bootstrap","gCF","sCF"), sep="/") %>%
  mutate_at(c('bootstrap', 'gCF',"sCF"), as.numeric) %>%
  mutate(combined_cf = paste(round(gCF), round(sCF), sep="/"))
tree.cf$combined_cf[tree.cf$combined_cf == "NA/NA"] <- NA
tree.cf <- as.treedata(tree.cf)

# Make pie charts
# Do not include the node column
# See: https://guangchuangyu.github.io/software/ggtree/vignettes/ggtree-inset.html
pies <- nodepie(cf.df, cols=c(2,4,6,8), alpha = 1)
pies <- lapply(pies, function(g) g + 
                 scale_fill_manual(values = c("darkgreen", "firebrick1", "pink", "grey90")))

# Error with these functions. This is apparently a fix to add clade labels. It might mess other things up?
# https://github.com/YuLab-SMU/enrichplot/issues/249
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

# Make cladogram with pie charts at nodes
cladogram.plot <- ggtree(tree, branch.length="none") + geom_tiplab(size=3, parse=TRUE) + geom_rootedge() +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["gCF"]])), size=3, nudge_y = 0.2) +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["sCF"]])), size=3, nudge_y = -0.2)+
  geom_cladelab(node=15, label="Blastocladiomycota", barsize=1, align=TRUE, offset=6) +
  geom_cladelab(node=25, label="Sanchytriomycota", barsize=1, align=TRUE, offset=6) + 
  xlim_tree(15)
cladogram.plot
cladogram.pies <- inset(cladogram.plot, pies, width=0.1, height=0.1)
cladogram.pies
ggsave("cladogram.pdf", width=10, height=5, limitsize = FALSE)

# Make cladogram with pie charts at nodes
cladogram.plot <- ggtree(tree, branch.length="none") + geom_tiplab(size=3, parse=TRUE) + geom_rootedge() +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["bootstrap"]])), size=3, nudge_y = 0.2) +
  geom_nodelab(mapping = aes(x = branch, label = tree.cf@data[["combined_cf"]]), size=3, nudge_y = -0.2)+
  geom_cladelab(node=15, label="Blastocladiomycota", barsize=1, align=TRUE, offset=6) +
  geom_cladelab(node=25, label="Sanchytriomycota", barsize=1, align=TRUE, offset=6) + 
  xlim_tree(15)
cladogram.plot
cladogram.pies <- inset(cladogram.plot, pies, width=0.1, height=0.1)
cladogram.pies
ggsave("cladogram2.pdf", width=10, height=5, limitsize = FALSE)

# Make ML tree with branch lengths
ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + xlim_tree(0.45)
ggtree(tree) + geom_text2(aes(label=node), hjust=-.3) + xlim_tree(0.45)
ggtree(tree) + geom_tiplab() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + xlim_tree(0.45)
tree.plot <- ggtree(tree) + geom_tiplab(size=3, parse=TRUE) + geom_rootedge() +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["gCF"]]), subset = round(tree.cf@data[["gCF"]]) >= 50), size=3, nudge_y = 0.3) + 
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["gCF"]]), subset = round(tree.cf@data[["gCF"]]) < 50), size=3, nudge_y = 0.3, color="firebrick1", fontface='bold') +
  geom_treescale() +
  geom_cladelab(node=15, label="Blastocladiomycota", barsize=1, align=TRUE, offset=0.7) +
  geom_cladelab(node=25, label="Sanchytriomycota", barsize=1, align=TRUE, offset=0.7) + 
  xlim_tree(2)
tree.plot
ggsave("phylogram.pdf", width=10, height=5, limitsize = FALSE)

##########################################################################################

# Consensus concatenation tree good partitions
tree <- read.tree("phylogenomics/good_partitions.treefile")
tree <- root_at_midpoint(tree)
tree$root.edge <- 0.01

# Concordance factors tree
tree.cf <- read.iqtree("phylogenomics/concord_good_partitions.cf.tree")
tree.cf@phylo <- root_at_midpoint(tree.cf@phylo)

# Concordance factors dataset
cf.df <- read.table("phylogenomics/concord_good_partitions.cf.stat", header=TRUE)
colnames(cf.df)[1] <- "node"
cf.df$node <- cf.df$node + 1

# Adjust tibble size for mapping onto tree
tree.cf <- as_tibble(tree.cf)
tree.cf$label2 <- NA
tree.cf$label2[15:24] <- tree.cf$label[15:24]
tree.cf <- tree.cf %>% 
  separate(label2, into = c("bootstrap","gCF","sCF"), sep="/") %>%
  mutate_at(c('bootstrap', 'gCF',"sCF"), as.numeric) %>%
  mutate(combined_cf = paste(round(gCF), round(sCF), sep="/"))
tree.cf$combined_cf[tree.cf$combined_cf == "NA/NA"] <- NA
tree.cf <- as.treedata(tree.cf)

# Make pie charts
# Do not include the node column
# See: https://guangchuangyu.github.io/software/ggtree/vignettes/ggtree-inset.html
pies <- nodepie(cf.df, cols=c(2,4,6,8), alpha = 1)
pies <- lapply(pies, function(g) g + 
                 scale_fill_manual(values = c("darkgreen", "firebrick1", "pink", "grey90")))

# Make cladogram with pie charts at nodes
cladogram.plot <- ggtree(tree, branch.length="none") + geom_tiplab(size=3, parse=TRUE) + geom_rootedge() +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["gCF"]])), size=3, nudge_y = 0.2) +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["sCF"]])), size=3, nudge_y = -0.2)+
  geom_cladelab(node=15, label="Blastocladiomycota", barsize=1, align=TRUE, offset=6) +
  geom_cladelab(node=25, label="Sanchytriomycota", barsize=1, align=TRUE, offset=6) + 
  xlim_tree(15)
cladogram.plot
cladogram.pies <- inset(cladogram.plot, pies, width=0.1, height=0.1)
cladogram.pies
ggsave("cladogram_good_partitions.pdf", width=10, height=5, limitsize = FALSE)

# Make cladogram with pie charts at nodes
cladogram.plot <- ggtree(tree, branch.length="none") + geom_tiplab(size=3, parse=TRUE) + geom_rootedge() +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["bootstrap"]])), size=3, nudge_y = 0.2) +
  geom_nodelab(mapping = aes(x = branch, label = tree.cf@data[["combined_cf"]]), size=3, nudge_y = -0.2)+
  geom_cladelab(node=15, label="Blastocladiomycota", barsize=1, align=TRUE, offset=6) +
  geom_cladelab(node=25, label="Sanchytriomycota", barsize=1, align=TRUE, offset=6) + 
  xlim_tree(15)
cladogram.plot
cladogram.pies <- inset(cladogram.plot, pies, width=0.1, height=0.1)
cladogram.pies
ggsave("cladogram_good_partitions2.pdf", width=10, height=5, limitsize = FALSE)

# Make ML tree with branch lengths
ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + xlim_tree(0.45)
ggtree(tree) + geom_text2(aes(label=node), hjust=-.3) + xlim_tree(0.45)
ggtree(tree) + geom_tiplab() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + xlim_tree(0.45)
tree.plot <- ggtree(tree) + geom_tiplab(size=3, parse=TRUE) + geom_rootedge() +
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["gCF"]]), subset = round(tree.cf@data[["gCF"]]) >= 50), size=3, nudge_y = 0.3) + 
  geom_nodelab(mapping = aes(x = branch, label = round(tree.cf@data[["gCF"]]), subset = round(tree.cf@data[["gCF"]]) < 50), size=3, nudge_y = 0.3, color="firebrick1", fontface='bold') +
  geom_treescale() +
  geom_cladelab(node=15, label="Blastocladiomycota", barsize=1, align=TRUE, offset=0.7) +
  geom_cladelab(node=25, label="Sanchytriomycota", barsize=1, align=TRUE, offset=0.7) + 
  xlim_tree(2)
tree.plot
ggsave("phylogram_good_partitions.pdf", width=10, height=5, limitsize = FALSE)
