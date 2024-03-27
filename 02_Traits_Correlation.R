########kendall & p值

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

#Calculate Kendall correlation and p-values for each pair of variables, and fill the matrices
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
