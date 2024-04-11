library(tidyverse)
library(scales)
library(patchwork)
library(dplyr)

# 读取数据
data <- read.csv("orchid_data/orchid_traits_without_deceit.csv")
data1 <- read.csv("orchid_data/Is_Deceit.csv")
# 增加是否为欺骗者信息
data <- left_join(data, data1[, c("species", "attraction_strategies")], by = "species")
# 将增加了是否为欺骗者信息的数据集导出
write.csv(data, "orchid_data/orchid_traits.csv", row.names = FALSE)

# 统计信息
data %>%
    group_by(terrestrial) %>%
    summarise(counts = n())
# 气候分布区
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

# 是否具有传粉媒介
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

# 是否为欺骗性传粉
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
# 绘制svg格式的结果
ggsave("orchid_result/fig1_abc_distribution_of_lf_traits.svg", p_merge, width = 10, height = 25)