## 01 Traits distribution
library(tidyverse)
library(scales)
library(patchwork)
library(dplyr)

# read in orchid traits data
data <- read.csv("orchid_data/orchid_traits.csv")
# data summary
data %>%
    group_by(terrestrial) %>%
    summarise(counts = n())
# climate region distribution
(p_climate <- data %>%
  mutate(terrestrial = ifelse(terrestrial == 0, "Epiphytic", "Terrestrial")) %>%
  group_by(terrestrial) %>%
  summarise(across(c(tropic, temperate, subtropic), sum)) %>%
  pivot_longer(cols = c(tropic, temperate, subtropic),
               values_to = "counts",
               names_to = "Climate") %>%
  ggplot(mapping = aes(x = factor(Climate, levels = c("tropic","subtropic", "temperate")),
                       y = counts,
                       fill = terrestrial)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = counts), position = position_fill(vjust = 0.5),family = "serif") +
  labs(x = "Climate", y = "", fill = "Life Form") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("") +
  theme(text = element_text(size = 10,family = "serif")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
)

# ggsave("orchid_result/climate_region.png", p_climate, width = 6, height = 4)

# pollination vectors distribution
(p_pollination_vectors <- data %>%
  mutate(terrestrial = ifelse(terrestrial == 0, "Epiphytic", "Terrestrial"),
         pollination_vectors = ifelse(pollination_vectors == 1, "Biotic", "Abiotic")) %>%
  group_by(terrestrial,pollination_vectors) %>%
  summarize(counts = n()) %>%
  ggplot(mapping = aes(x = factor(pollination_vectors, levels = c("Biotic", "Abiotic")),
                       y = counts,
                       fill = terrestrial)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = counts), position = position_fill(vjust = 0.5),family = "serif") +
  labs(x = "Pollination Vectors", y = "", fill = "Life Form") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("") +
  theme(text = element_text(size = 10,family = "serif")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
)
# ggsave("orchid_result/pollination_vectors.png", p_pollination_vectors, width = 6, height = 4)

# pollinatior attraction strategies distribution
(p_attraction_strategies <- data %>%
  mutate(terrestrial = ifelse(terrestrial == 0, "Epiphytic", "Terrestrial"),
         attraction_strategies = ifelse(attraction_strategies == 1, "Deception", "Reward")) %>%
  group_by(terrestrial,attraction_strategies) %>%
  summarize(counts = n()) %>%
  na.omit() %>% #因为这个数据集中有缺失值，所以要去掉
  ggplot(mapping = aes(x = factor(attraction_strategies, levels = c("Deception", "Reward")),
                       y = counts,
                       fill = terrestrial)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = counts), position = position_fill(vjust = 0.5),family = "serif") +
  labs(x = "Attraction Strategies", y = "", fill = "Life Form") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("") +
  theme(text = element_text(size = 10,family = "serif")) +
  theme(axis.text.x = element_text(angle = 10,hjust = 1))
)
# ggsave("orchid_result/attraction_strategies.png", p_attraction_strategies, width = 6, height = 4)

p_merge <- p_climate / p_pollination_vectors / p_attraction_strategies +
  plot_layout(guides = 'collect') +
  plot_annotation(title = "", tag_levels = "A")
# write out the figure
ggsave("orchid_result/fig1_abc_distribution_of_lf_traits.svg", p_merge, width = 10, height = 25)

## construct the phylogenetic tree
library(tidyverse)
library(U.PhyloMaker)  # devtools::install_github("jinyizju/U.PhyloMaker")
library(ape)

# read in the species list
sp.list <- read.csv("orchid_data/orchid_traits.csv") %>%
  select(species,genus,family)
megatree <- read.tree("orchid_data/plant_megatree.tre")
gen.list <- read.csv("orchid_data/plant_genus_list.csv")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)


tree <- result$phylo

# visualize the tree(not nesswary)
plot(tree, cex = .01, show.tip.label = F)
axisPhylo()

# save the phylogenetic tree
write.tree(tree, "orchid_data/Orchidaceae.tre")

## 02 traits' correlation

## kendall correlation coefficient & p-value

library(readr)
library(corrplot)
library(ggplot2)

data <- read_csv("orchid_data/orchid_traits.csv")
selected_vars <- data[, c('terrestrial', 'epiphytic', 'tropic', 'subtropic', 'temperate', 'pollination_vectors', 'attraction_strategies')]
colnames(selected_vars) <- c("Terrestrial", "Epiphytic", 'Tropic', 'Subtropic', 'Temperate', "Pollination Vectors", "Attraction Strategies")

#Initialize matrices to store the Kendall correlation results and p-values
cor_matrix <- matrix(nrow = ncol(selected_vars), ncol = ncol(selected_vars))
p_matrix <- matrix(nrow = ncol(selected_vars), ncol = ncol(selected_vars))
rownames(cor_matrix) <- colnames(selected_vars)
colnames(cor_matrix) <- colnames(selected_vars)
rownames(p_matrix) <- colnames(selected_vars)
colnames(p_matrix) <- colnames(selected_vars)

# Calculate Kendall correlation and p-values for each pair of variables, and fill the matrices
for (i in 1:ncol(selected_vars)) {
  for (j in i:ncol(selected_vars)) {
    if (i == j) {
      cor_matrix[i, j] <- 1  # The correlation of a variable with itself is 1
      p_matrix[i, j] <- NA  # P-value is not applicable for the same variable comparison
    } else {
      test_result <- cor.test(selected_vars[[i]], selected_vars[[j]], method = "kendall", exact = FALSE)
      cor_matrix[i, j] <- test_result$estimate
      cor_matrix[j, i] <- test_result$estimate  # Fill both the upper and lower triangle of the matrix
      p_matrix[i, j] <- test_result$p.value
      p_matrix[j, i] <- test_result$p.value  # Same for p-values
    }
  }
}


cor_matrix <- round(cor_matrix, digits = 3)
p_matrix <- round(p_matrix, digits = 3)

print(cor_matrix)
print(p_matrix)  # Also print the p-value matrix

corrplot(cor_matrix, method = "circle")
svg("orchid_result/fig1e_correlations_kendall.svg", width = 9.6, height = 6.4)
corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 1.2)
dev.off()

# The relative contribution of orchids factors to variation of life form
library(phylolm)
library(ape)
library(phytools)
library(rr2)
library(tidyverse)
library(patchwork)

# read in orchid data
traits <- read_csv("orchid_data/orchid_traits.csv") %>%
  mutate(species_name = species) %>%
  column_to_rownames(var = "species_name")

# read in the phylogeny and standardize the tree
tree <- read.tree("orchid_data/Orchidaceae.tre")
st_tree <- function(tree) {
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))
  return(tree)
}
tree <- st_tree(tree)

# check the latin name of the species
geiger::name.check(tree, traits)

# Build models for pairwise comparisons, where "LF" represents the dependent variable of growth form, "CR" represents the climatic distribution area,
# "CH" indicates the presence or absence of pollination vectors, "DR" denotes deceptive pollination, and "1" represents the intercept term.
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

# Construct a function to calculate the differences between models
get_r2like <- function(mod1, mod2, df) {
  # - mod1: full model
  # - mod2: The model that needs to be subtracted
  # - df:
  # Degrees of freedom of the factor, where the difference in factors between two models is assessed. For instance, in this paper, "CR" actually represents three climatic zones
  # The output is a list that includes the model name, subtracted model name, corresponding regression coefficient, log-likelihood difference, and p-value
  w <- c(as.character(substitute(mod1)),
         as.character(substitute(mod2)),
         as.character(R2_lik(mod1, mod2)),
         as.character(logLik(mod1)[[1]] - logLik(mod2)[[1]]),
         as.character(pchisq(2 * (logLik(mod1)[[1]] - logLik(mod2)[[1]]),
                             df = df, lower.tail = F))
  )
  return(w)
}

# run the relative contribution of factors to variation of life form
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


# output the result
write.csv(df_fullphy, "orchid_result/fig2_&_sp_table1.csv", row.names = FALSE)

# visualization
# full model consider the phylogenetic conservation
subset_full <- df_fullphy[2:5, ]
p_full <- ggplot(subset_full, aes(x = reduce_mod, y = R2_lik)) +
  geom_col(fill = "#F0BB62") +
  scale_x_discrete(limits = subset_full$reduce_mod[order(subset_full$R2_lik, decreasing = T)],
                   labels = c("Phylogeny", "Attraction Strategies", "Climate", "Pollination Vectors")) +
  geom_text(aes(label = paste0(round(R2_lik*100,1),"%")), vjust = -0.5, family = "serif") +
  labs(x = "remain", y = expression(R[lik]^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 14,family = "serif")) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))
ggsave("orchid_result/fig2_a.svg",p_full, width = 4, height = 4)

# model without the phylogenetic conservation
subset_rmphy <- df_fullphy[7:9, ]
p_rmphy <- ggplot(subset_rmphy, aes(x = reduce_mod, y = R2_lik)) +
  geom_col(fill = "#F0BB62") +
  scale_x_discrete(limits = subset_rmphy$reduce_mod[order(subset_rmphy$R2_lik, decreasing = T)],
                   labels = c("Climate", "Attraction Strategies", "Pollination Vectors")) +
  geom_text(aes(label = paste0(round(R2_lik*100,1),"%")), vjust = -0.5, family = "serif") +
  labs(x = "remain", y = expression(R[lik]^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 14,family = "serif")) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))
ggsave("orchid_result/spf_rmphy.svg",p_rmphy, width = 4, height = 4)

# merge two pictures together
p_merge <- p + p_rmphy +
  plot_layout(guides = 'collect') +
  plot_annotation(title = "", tag_levels = "A")
  #theme(aspect.ratio = 1)
# write out
ggsave("orchid_result/compare_full_rmphy_merge.svg", p_merge, width = 9, height = 4)

merge <- left_join(subset_full, subset_rmphy, by = "remain")
merge_rmornot <- subset(merge, select = c(remain, R2_lik.x, R2_lik.y))
p_duidie <-ggplot(merge_rmornot, aes(x = remain)) +
  geom_col(aes(y = R2_lik.x, fill = "R2_lik.x"), position = "dodge", stat = "identity") +
  geom_col(aes(y = R2_lik.y, fill = "R2_lik.y"), position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("R2_lik.x" = "blue", "R2_lik.y" = "red")) +
  labs(x = "Factors",
       y = expression(R[lik]^2))

ggsave("orchid_result/compare_full_rmphy_duidie.svg", p_duidie, width = 9, height = 4)

## 04 Traits‘ phylogenetic conservatism

library(phylolm)
library(ape)
library(phytools)
library(rr2)
library(tidyverse)
library(patchwork)

# read in orchid data
traits <- read.csv("orchid_data/orchid_traits.csv") %>%
    column_to_rownames(var = "species")

# read in the phylogeny and standardize the tree
tree <- read.tree("orchid_data/Orchidaceae.tre")
st_tree <- function(tree) {
  tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))
  return(tree)
}
st_tree <- st_tree(tree)

# check the latin name of the species
geiger::name.check(st_tree, traits)

# Build models for pairwise comparisons, where "LF" represents the dependent variable of growth form, "CR" represents the climatic distribution area,
# "CH" indicates the presence or absence of pollination vectors, "DR" denotes deceptive pollination, and "1" represents the intercept term.
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

# Construct a function to calculate the differences between models
get_r2like <- function(mod1, mod2, df) {
  w <- c(as.character(substitute(mod1)),
         as.character(substitute(mod2)),
         as.character(R2_lik(mod1, mod2)),
         as.character(logLik(mod1)[[1]] - logLik(mod2)[[1]]),
         as.character(pchisq(2 * (logLik(mod1)[[1]] - logLik(mod2)[[1]]),
                             df = df, lower.tail = F))
  )
  return(w)
}

# # run the phylogenetic conservatism of all traits
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