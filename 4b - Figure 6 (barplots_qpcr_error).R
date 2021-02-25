#################
## Config
#################
library(tidyverse)
library(magrittr)
library(forcats)
library(here)

setwd(here())

################
## Default vars
################
gene_names <- c(
  "PAX7"
  ,"RUNX3"
  ,"ARNT2"
  #,"MEF2C"
  ,"GLI3"
  ,"PBX3"
  ,"CREB3L1"
)

diseases <- c(
  "Ewing"
  ,"Hepatoblastoma"
  ,"Neuroblastoma"
  ,"Meduloblastoma"
)

##########################
## Pivoting qPCR results
##########################
results_wide <- read_tsv("obs_wide.tsv")

results <- results_wide %>% 
  pivot_longer(paste0("n", 1:3)) %>%
  group_by(gene, disease, cell) %>%
  summarise(
    mean = mean(value)
    ,se   = sd(value)/sqrt(n())
  ) %>%
  ungroup

results %<>% filter(gene %in% gene_names)

results %<>% arrange(disease, cell)

results %<>% mutate(
  gene    = fct_relevel(gene, gene_names)
  ,disease = fct_relevel(disease, diseases)
  ,cell    = fct_reorder(paste0(cell, " - ", disease), as.numeric(disease))
)

#####################
## Plotting
#####################
ggplot(results, aes(x = cell, y = mean)) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se)
    ,stat  = "identity"
    ,color = "#333333"
  ) +
  geom_bar(
    stat  = "identity"
    ,fill  = "white"
    ,color = "#333333"
  ) + 
  facet_wrap(gene ~ ., ncol = 2, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Disease") +
  ylab("Relative mRNA level")

ggsave("qpcr_se.pdf", width = 0.95 * 5, height = 8, useDingbats = F)