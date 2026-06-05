## Make cartoon with simulated eQTL data for Figure 1a

library(tidyverse)
library(ggbeeswarm)

set.seed(47)

# simulate gene expression data for each genotype
AA <- rnorm(10, mean = 15)
Aa <- rnorm(10, mean = 10)
aa <- rnorm(10, mean = 5)

# combine into one table
genex <- tibble(AA, Aa, aa) %>% 
  pivot_longer(cols = everything(), names_to = "genotype", values_to = "expression") %>% 
  mutate(genotype = fct_inorder(genotype))

# make plot
ggplot(genex, aes(x = genotype, y = expression)) +
  geom_beeswarm(color = "#4685c5", cex = 6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black") +
  labs(y = "Gene expression") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(filename = "other_analyses/cartoons_fig1a/fig1a_eqtl.pdf", width = 2, height = 1.5)
