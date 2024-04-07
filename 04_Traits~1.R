library(phylolm) # 注意，这一行需要先在ggtree前载入
library(ape)
library(phytools)
library(rr2)
library(tidyverse)
library(patchwork)
library(evobiR)

# 读取数据
traits <- read.csv("orchid_data/orchid_traits.csv") %>%
    column_to_rownames(var = "species")
# 读取发育树并进行标准化
tree <- read.tree("orchid_data/Orchidaceae.tre")
st_tree <- function(tree) {
  tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))
  return(tree)
}
st_tree <- st_tree(tree)

# 命名检查
geiger::name.check(st_tree, traits)

# 构建模型两两之间的比较，其中LF表示因变量生长形式，CR表示气候分布区，CH表示有无传粉媒介，DR表示是否为欺骗性传粉，1表示截距项
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

# 计算模型的差异
get_r2like <- function(mod1, mod2, df) {
  # - mod1: 模型
  # - mod2: 需要减去的模型
  # - df: 因子的自由度，这里需要自定判断两个模型之间因子的差异量，例如本文中CR虽然名字只有两个字母，但其实表示了3个气候区
  # 输出结果是一个包含模型名、减去模型名、对应偏方系数、逻辑似然差值和p值的列表
  w <- c(as.character(substitute(mod1)),
         as.character(substitute(mod2)),
         as.character(R2_lik(mod1, mod2)),
         as.character(logLik(mod1)[[1]] - logLik(mod2)[[1]]),
         as.character(pchisq(2 * (logLik(mod1)[[1]] - logLik(mod2)[[1]]),
                             df = df, lower.tail = F))
  )
  return(w)
}

# 输出结果
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




#传统系统发育信号值计算
# 计算5个形状的k值
terrestrial <- setNames(traits$terrestrial, rownames(traits))
tropic <- setNames(traits$tropic, rownames(traits))
subtropic <- setNames(traits$subtropic, rownames(traits))
temperate <- setNames(traits$temperate, rownames(traits))
pollination_vectors <- setNames(traits$pollination_vectors, rownames(traits))
attraction_strategies <- setNames(traits$attraction_strategies, rownames(traits))

# 按照 st_tree$tip.label 的顺序重新排列整个 traits 数据框
traits_reordered <- traits[match(st_tree$tip.label, rownames(traits)), ]
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


#计算5个性状的Pagel’s λ值

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


# 计算5个性状的α值
alphas <- NULL
all_label <- c("terrestrial", "tropic", "subtropic", "subtropic", "temperate", "pollination_vectors", "attraction_strategies")

for (i in 1:length(all_label)) {
  f <- paste(all_label[i], "~", "1")
  model <- phyloglm(f, data = traits_filtered, phy = tree_filtered)
  alpha <- model$alpha
  alphas <- c(alphas, alpha)
}

# 创建结果数据框
df_alpha <- tibble(name = all_label, alpha = alphas)


