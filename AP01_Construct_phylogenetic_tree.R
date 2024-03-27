library(tidyverse)
library(U.PhyloMaker)  # devtools::install_github("jinyizju/U.PhyloMaker")
library(ape)

# 读取物种清单数据
sp.list <- read.csv("orchid_data/orchid_traits.csv") %>%
  select(species,genus,family)
megatree <- read.tree("orchid_data/plant_megatree.tre")
gen.list <- read.csv("orchid_data/plant_genus_list.csv")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)


tree <- result$phylo

# 简单绘制树的图形(非必需)
plot(tree, cex = .01, show.tip.label = F)
axisPhylo()

# 保存系统树到本地
write.tree(tree, "orchid_data/Orchidaceae.tre")