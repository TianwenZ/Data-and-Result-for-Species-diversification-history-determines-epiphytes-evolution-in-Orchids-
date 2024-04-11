library(tidyverse)
library(U.PhyloMaker)  # devtools::install_github("jinyizju/U.PhyloMaker")
library(ape)

# read in the plant list
sp.list <- read.csv("orchid_data/orchid_traits.csv") %>%
  select(species,genus,family)
megatree <- read.tree("orchid_data/plant_megatree.tre")
gen.list <- read.csv("orchid_data/plant_genus_list.csv")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)


tree <- result$phylo

# simply plotting the phylogentic tree (not necessary)
plot(tree, cex = .01, show.tip.label = F)
axisPhylo()

# save it
write.tree(tree, "orchid_data/Orchidaceae.tre")