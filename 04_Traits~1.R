library(phylolm) # Note, this R package needs to be loaded before ggtree
library(ape)
library(phytools)
library(rr2)
library(tidyverse)
library(patchwork)
library(evobiR)

# read in data
traits <- read.csv("orchid_data/orchid_traits.csv") %>%
    column_to_rownames(var = "species")
# 读取发育树并进行标准化
tree <- read.tree("orchid_data/Orchidaceae.tre")
st_tree <- function(tree) {
  tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))
  return(tree)
}
st_tree <- st_tree(tree)

# name check
geiger::name.check(st_tree, traits)

# construct the models, among them, LF represents the life form, CR represents the climate region, PV indicates represents the pollination vectors,
                                    # AS represents the (pollinator) attraction strategies, 1 indicates the intercept term
LF_x_CR_CH_DR_PHY_1 <- phyloglm(terrestrial ~ temperate + subtropic + tropic + pollination_vectors + attraction_strategies + 1,
                                data = traits, phy = tree)

LF_x_PHY_1 <- phyloglm(terrestrial ~ 1,
                       data = traits, phy = tree)
LF_x_1 <- glm(terrestrial ~ 1, data = traits, family = "binomial")

PV_x_PHY_1 <- phyloglm(pollination_vectors ~ 1,
                       data = traits, phy = tree)
PV_x_1 <- glm(pollination_vectors ~ 1, data = traits, family = "binomial")

AS_x_PHY_1 <- phyloglm(attraction_strategies ~ 1,
                       data = traits, phy = tree)
AS_x_1 <- glm(attraction_strategies ~ 1, data = traits, family = "binomial")

CR_Tro_x_PHY_1 <- phyloglm(tropic ~  1,
                                data = traits, phy = tree)
CR_Tro_x_1 <- glm(tropic ~  1, data = traits, family = "binomial")

CR_Subtro_x_PHY_1 <- phyloglm(subtropic ~  1,
                                data = traits, phy = tree)
CR_Subtro_x_1 <- glm(subtropic ~  1, data = traits, family = "binomial")

CR_Tem_x_PHY_1 <- phyloglm(temperate ~  1,
                                data = traits, phy = tree)
CR_Tem_x_1 <- glm(temperate ~  1, data = traits, family = "binomial")

#  essentially evaluates whether the inclusion of additional factors in mod1 significantly improves its fit to the data compared to mod2
get_r2like <- function(mod1, mod2, df) {
  # mod1: The model
  # mod2: The model to be subtracted
  # df: Degrees of freedom for the factor, here one needs to independently judge the difference in factors between the two models.
    #For example, CR, although only two letters long, actually represents three climate regions
  # The output is a list containing model names, subtracted model name, corresponding partial coefficients, log-likelihood difference, and p-values.
  w <- c(as.character(substitute(mod1)),
         as.character(substitute(mod2)),
         as.character(R2_lik(mod1, mod2)),
         as.character(logLik(mod1)[[1]] - logLik(mod2)[[1]]),
         as.character(pchisq(2 * (logLik(mod1)[[1]] - logLik(mod2)[[1]]),
                             df = df, lower.tail = F))
  )
  return(w)
}

df_fullphy <- data.frame(get_r2like(LF_x_PHY_1, LF_x_1, df = 1),
                         get_r2like(PV_x_PHY_1, PV_x_1, df = 1),
                         get_r2like(AS_x_PHY_1, AS_x_1, df = 1),
                         get_r2like(CR_Tro_x_PHY_1, CR_Tro_x_1, df = 1),
                         get_r2like(CR_Subtro_x_PHY_1, CR_Subtro_x_1, df = 1),
                         get_r2like(CR_Tem_x_PHY_1, CR_Tem_x_1, df = 1),
                         row.names = c("full_mod", "reduce_mod", "R2_lik", "delta_log_lik", "p")) %>%
  t() %>%
  as_tibble() %>%
  mutate(R2_lik = parse_number(R2_lik),
         delta_log_lik = parse_number(delta_log_lik),
         p = parse_number(p)) %>%
  mutate_if(is.numeric, round, 3)

print(df_fullphy)

write.csv(df_fullphy, "orchid_result/sf_table2_phylo_p_.csv", row.names = FALSE)




#calculate the k value, Pagel’s λ and alpha values of all traits
# k value
terrestrial <- setNames(traits$terrestrial, rownames(traits))
tropic <- setNames(traits$tropic, rownames(traits))
subtropic <- setNames(traits$subtropic, rownames(traits))
temperate <- setNames(traits$temperate, rownames(traits))
pollination_vectors <- setNames(traits$pollination_vectors, rownames(traits))
attraction_strategies <- setNames(traits$attraction_strategies, rownames(traits))

terrestrial_K <- phylosig(st_tree, terrestrial, method = "K", test = T, nsim = 999)
tropic_K <- phylosig(st_tree, tropic, method = "K", test = T, nsim = 999)
subtropic_K <- phylosig(st_tree, subtropic, method = "K", test = T, nsim = 999)
temperate_K <- phylosig(st_tree, temperate, method = "K", test = T, nsim = 999)
pollination_vectors_K <- phylosig(st_tree, pollination_vectors, method = "K", test = T, nsim = 999)
attraction_strategies_K <- phylosig(st_tree, attraction_strategies, method = "K", test = T, nsim = 999)

print(terrestrial_K)
print(tropic_K)
print(subtropic_K)
print(temperate_K)
print(attraction_strategies_K)
print(pollination_vectors_K)


# Pagel’s λ

terrestrial_lambda <- phylolm(terrestrial~1,data=traits,phy=st_tree,model="lambda")
tropic_lambda <- phylolm(tropic~1,data=traits,phy=st_tree,model="lambda")
subtropic_lambda <- phylolm(subtropic~1,data=traits,phy=st_tree,model="lambda")
temperate_lambda <- phylolm(temperate~1,data=traits,phy=st_tree,model="lambda")
pollination_vectors_lambda <- phylolm(pollination_vectors~1,data=traits,phy=st_tree,model="lambda")
attraction_strategies_lambda <- phylolm(attraction_strategies~1,data=traits,phy=st_tree,model="lambda")

print(terrestrial_lambda)
print(tropic_lambda)
print(subtropic_lambda)
print(temperate_lambda)
print(attraction_strategies_lambda)
print(pollination_vectors_lambda)


# alpha
alphas <- NULL
all_label <- c("terrestrial", "tropic", "subtropic", "subtropic", "temperate", "pollination_vectors", "attraction_strategies")

for (i in 1:length(all_label)) {
  f <- paste(all_label[i], "~", "1")
  model <- phyloglm(f, data = traits_filtered, phy = tree_filtered)
  alpha <- model$alpha
  alphas <- c(alphas, alpha)
}
df_alpha <- tibble(name = all_label, alpha = alphas)


