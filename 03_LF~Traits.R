library(phylolm) # 注意，这一行需要先在ggtree前载入
library(ape)
library(phytools)
library(rr2)
library(tidyverse)
library(patchwork)

# read in data
traits <- read_csv("orchid_data/orchid_traits.csv") %>%
  mutate(species_name = species) %>%
  column_to_rownames(var = "species_name")

# read in the phylogentic tree and standardize the tree
tree <- read.tree("orchid_data/Orchidaceae.tre")
st_tree <- function(tree) {
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))
  return(tree)
}
tree <- st_tree(tree)

# name check
geiger::name.check(tree, traits)

# construct the models, among them, LF represents the life form, CR represents the climate region, PV indicates represents the pollination vectors,
                                  # AS represents the (pollinator) attraction strategies, 1 indicates the intercept term
LF_x_CR_PV_AS_PHY_1 <- phyloglm(terrestrial ~ temperate + subtropic + tropic + pollination_vectors + attraction_strategies + 1,
                                data = traits, phy = tree)
LF_x_CR_PV_AS_1 <- glm(terrestrial ~ temperate + subtropic + tropic + pollination_vectors + attraction_strategies + 1,
                       data = traits, family = "binomial")
LF_x_CR_PV_PHY_1 <- phyloglm(terrestrial ~ temperate + subtropic + tropic + pollination_vectors + 1,
                             data = traits, phy = tree)
LF_x_PV_AS_PHY_1 <- phyloglm(terrestrial ~ pollination_vectors + attraction_strategies + 1,
                             data = traits, phy = tree)
LF_x_CR_AS_PHY_1 <- phyloglm(terrestrial ~ temperate + subtropic + tropic + attraction_strategies + 1,
                             data = traits, phy = tree)
LF_x_CR_PV_1 <- glm(terrestrial ~ temperate + subtropic + tropic + pollination_vectors + 1,
                    data = traits, family = "binomial")
LF_x_CR_AS_1 <- glm(terrestrial ~ temperate + subtropic + tropic + attraction_strategies + 1,
                    data = traits, family = "binomial")
LF_x_CR_PHY_1 <- phyloglm(terrestrial ~ temperate + subtropic + tropic + 1,
                          data = traits, phy = tree)
LF_x_PV_AS_1 <- glm(terrestrial ~ pollination_vectors + attraction_strategies + 1,
                    data = traits, family = "binomial")
LF_x_PV_PHY_1 <- phyloglm(terrestrial ~ pollination_vectors + 1,
                          data = traits, phy = tree)
LF_x_AS_PHY_1 <- phyloglm(terrestrial ~ attraction_strategies + 1,
                          data = traits, phy = tree)
LF_x_CR_1 <- glm(terrestrial ~ temperate + subtropic + tropic + 1,
                 data = traits, family = "binomial")
LF_x_PV_1 <- glm(terrestrial ~ pollination_vectors + 1,
                 data = traits, family = "binomial")
LF_x_AS_1 <- glm(terrestrial ~ attraction_strategies + 1,
                 data = traits, family = "binomial")
LF_x_PHY_1 <- phyloglm(terrestrial ~ 1,
                       data = traits, phy = tree)
LF_x_1 <- glm(terrestrial ~ 1, data = traits, family = "binomial")

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

df_fullphy <- data.frame(get_r2like(LF_x_CR_PV_AS_PHY_1, LF_x_1, df = 6),
                         get_r2like(LF_x_CR_PV_AS_PHY_1, LF_x_CR_PV_AS_1, df = 1),
                         get_r2like(LF_x_CR_PV_AS_PHY_1, LF_x_CR_PV_PHY_1, df = 1),
                         get_r2like(LF_x_CR_PV_AS_PHY_1, LF_x_CR_AS_PHY_1, df = 1),
                         get_r2like(LF_x_CR_PV_AS_PHY_1, LF_x_PV_AS_PHY_1, df = 3),
                         get_r2like(LF_x_CR_PV_AS_1, LF_x_1, df = 5),
                         get_r2like(LF_x_CR_PV_AS_1, LF_x_CR_PV_1, df = 1),
                         get_r2like(LF_x_CR_PV_AS_1, LF_x_CR_AS_1, df = 1),
                         get_r2like(LF_x_CR_PV_AS_1, LF_x_PV_AS_1, df = 3),
                         get_r2like(LF_x_CR_PV_PHY_1, LF_x_CR_PV_1, df = 1),
                         get_r2like(LF_x_CR_PV_PHY_1, LF_x_PV_PHY_1, df = 3),
                         get_r2like(LF_x_CR_PV_PHY_1, LF_x_CR_PHY_1, df = 1),
                         get_r2like(LF_x_CR_AS_PHY_1, LF_x_CR_AS_1, df = 1),
                         get_r2like(LF_x_CR_AS_PHY_1, LF_x_CR_PHY_1, df = 1),
                         get_r2like(LF_x_CR_AS_PHY_1, LF_x_AS_PHY_1, df = 1),
                         get_r2like(LF_x_PV_AS_PHY_1, LF_x_PV_PHY_1, df = 1),
                         get_r2like(LF_x_PV_AS_PHY_1, LF_x_AS_PHY_1, df = 1),
                         get_r2like(LF_x_CR_PV_1, LF_x_CR_1, df = 1),
                         get_r2like(LF_x_CR_PV_1, LF_x_PV_1, df = 3),
                         get_r2like(LF_x_CR_AS_1, LF_x_CR_1, df = 1),
                         get_r2like(LF_x_CR_AS_1, LF_x_AS_1, df = 3),
                         get_r2like(LF_x_CR_PHY_1, LF_x_CR_1, df = 1),
                         get_r2like(LF_x_CR_PHY_1, LF_x_PHY_1, df = 3),
                         get_r2like(LF_x_PV_AS_1, LF_x_PV_1, df = 1),
                         get_r2like(LF_x_PV_AS_1, LF_x_AS_1, df = 1),
                         get_r2like(LF_x_PV_PHY_1, LF_x_PV_1, df = 1),
                         get_r2like(LF_x_PV_PHY_1, LF_x_PHY_1, df = 1),
                         get_r2like(LF_x_AS_PHY_1, LF_x_AS_1, df = 1),
                         get_r2like(LF_x_AS_PHY_1, LF_x_PHY_1, df = 1),
                         get_r2like(LF_x_CR_1, LF_x_1, df = 3),
                         get_r2like(LF_x_PV_1, LF_x_1, df = 1),
                         get_r2like(LF_x_AS_1, LF_x_1, df = 1),
                         get_r2like(LF_x_PHY_1, LF_x_1, df = 1),
                         row.names = c("full_mod", "reduce_mod", "R2_lik", "delta_log_lik", "p")) %>%
  t() %>%
  as_tibble() %>%
  mutate(R2_lik = parse_number(R2_lik),
         delta_log_lik = parse_number(delta_log_lik),
         p = parse_number(p)) %>%
  mutate_if(is.numeric, round, 3)


# writing out the result
write.csv(df_fullphy, "orchid_result/fig2_&_sp_table1.csv", row.names = FALSE)

# plotting
# full model
subset_full <- df_fullphy[2:5, ]
p_full <- ggplot(subset_full, aes(x = reduce_mod, y = R2_lik)) +
  geom_col(fill = "#F0BB62") +
  # The order of the x-axis in the bar chart is arranged from largest to smallest according to R2_lik
  scale_x_discrete(limits = subset_full$reduce_mod[order(subset_full$R2_lik, decreasing = T)],
                   labels = c("Phylogeny", "Attraction Strategies", "Climate", "Pollination Vectors")) +
  geom_text(aes(label = paste0(round(R2_lik*100,1),"%")), vjust = -0.5, family = "serif") +
  labs(x = "remain", y = expression(R[lik]^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 14,family = "serif")) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))
ggsave("orchid_result/fig2_a.svg",p_full, width = 4, height = 4)

# model with out phylogeny conservation
subset_rmphy <- df_fullphy[7:9, ]
p_rmphy <- ggplot(subset_rmphy, aes(x = reduce_mod, y = R2_lik)) +
  geom_col(fill = "#F0BB62") +
  # 柱形图的 x 轴顺序按照 R2_lik 从大到小排列
  scale_x_discrete(limits = subset_rmphy$reduce_mod[order(subset_rmphy$R2_lik, decreasing = T)],
                   labels = c("Climate", "Attraction Strategies", "Pollination Vectors")) +
  # 添加标签,使用百分比
  geom_text(aes(label = paste0(round(R2_lik*100,1),"%")), vjust = -0.5, family = "serif") +
  labs(x = "remain", y = expression(R[lik]^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 14,family = "serif")) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))
ggsave("orchid_result/spf_rmphy.svg",p_rmphy, width = 4, height = 4)


# combine the plot
p_merge <- p + p_rmphy +
  plot_layout(guides = 'collect') +
  plot_annotation(title = "", tag_levels = "A")
  #theme(aspect.ratio = 1)
# writing out the plot
ggsave("orchid_result/compare_full_rmphy_merge.svg", p_merge, width = 9, height = 4)

#  clustered bar plot
merge <- left_join(subset_full, subset_rmphy, by = "remain")
merge_rmornot <- subset(merge, select = c(remain, R2_lik.x, R2_lik.y))
p_duidie <-ggplot(merge_rmornot, aes(x = remain)) +
  geom_col(aes(y = R2_lik.x, fill = "R2_lik.x"), position = "dodge", stat = "identity") +
  geom_col(aes(y = R2_lik.y, fill = "R2_lik.y"), position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("R2_lik.x" = "blue", "R2_lik.y" = "red")) +
  labs(x = "Factors",
       y = expression(R[lik]^2))

ggsave("orchid_result/compare_full_rmphy_duidie.svg", p_duidie, width = 9, height = 4)
