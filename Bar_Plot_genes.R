library(readxl)
library(dplyr)
library(ggplot2)

# Read data
genes <- read_excel("C:/Users/Anna Graf/OneDrive/Desktop/Master/Masterarbeit/RNAi_screen_hits_and_function.xlsx")

# Summarize by function
function_summary <- genes %>%
  group_by(gene_function) %>%
  summarise(
    Count = n(),
    Genes = paste(sort(name), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(Count), gene_function)

function_summary <- function_summary %>%
  mutate(
    gene_function = factor(gene_function,
                           levels = rev(gene_function))
  )

# Plot
ggplot(function_summary,
       aes(x = Count,
           y = gene_function,
           fill = Count)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = Genes),
    hjust = 0,
    nudge_x = 0.2,
    size = 3
  ) +
  scale_fill_gradient(
    low = "#F8BBD0",
    high = "#C2185B"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.6))) +
  labs(
    x = "Number of genes",
    y = "Gene function",
    title = "Gene functions of Hits"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )
