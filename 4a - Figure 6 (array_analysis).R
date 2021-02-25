#################
## Config
#################
library(hgu133plus2.db)
library(tidyverse)
library(magrittr)
library(GEOquery)
library(EnvStats)
library(broom)
library(affy)
library(here)

overwrite_rds <- FALSE

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

geo_disease_mapping <- c(
  "GSE34620" = "Ewing"
  ,"GSE75271" = "Hepatoblastoma"
  ,"GSE16476" = "Neuroblastoma"
  ,"GSE14827" = "Osteosarcoma"
  ,"GSE29684" = "Retinoblastoma"
  ,"GSE66533" = "Rhabdomyosarcoma"
  ,"GSE53224" = "Wilms' tumor"
  ,"GSE3526"  = "Normal tissue"
  #,"GSE108089"= "ped_no_ewing"
)

geo_ids <- names(geo_disease_mapping)


####################################
## Downloading and extracting files
####################################
for(geo_id in geo_ids) {
  # Downloading
  geo_file <- paste0("downloaded_data/",geo_id,"_RAW.tar")
  if(!file.exists(geo_file)){
    download.file(
      url      = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=",geo_id,"&format=file")
      ,destfile = geo_file
      ,mode     = "wb"
    )
  }
  # Extracting
  if(!dir.exists(paste0("downloaded_data/",geo_id))){
    untar(geo_file, exdir = paste0("downloaded_data/",geo_id))
  }
}


##################
## Probe mapping
##################
probe_mapping <- AnnotationDbi::select(
   x       = hgu133plus2.db
  ,keys    = gene_names
  ,columns = "PROBEID"
  ,keytype = "ALIAS"
)


##########################
## Retrieving expression
##########################
cel_files <- list.files("downloaded_data", pattern = "CEL", recursive = T, full.names = T)

# If file does not exist, process and save. Else, just load it.
if(!file.exists("expression.rds") | overwrite_rds){
  rawdata <- ReadAffy(filenames = cel_files, compress = TRUE)
  
  expression_set <- rma(rawdata)
  
  # Releasing memory
  rm(rawdata)
  gc()
  
  expression <- exprs(expression_set)
  
  expression <- expression[probe_mapping[,"probe_id"],]
  
  # Releasing memory
  rm(expression_set)
  gc()
  
  saveRDS(expression, file = "expression.rds")
} else {
  expression <- readRDS("expression.rds")
}


####################################
## GEO <-> CEL <-> Disease mapping
####################################
geo_cel_mapping <- data.frame(
   cel_file = basename(cel_files)
  ,geo_id   = basename(dirname(cel_files))
)

geo_cel_mapping %<>% mutate(disease = geo_disease_mapping[geo_id])

cel_disease_mapping <- geo_cel_mapping %$% set_names(disease, cel_file)


############
## Pivoting
############
expression_long <- as.data.frame(expression)

expression_long %<>% mutate(probe_id = rownames(.))

expression_long %<>% pivot_longer(-probe_id, "cel_file")

expression_long %<>% mutate(disease = cel_disease_mapping[cel_file])


#########################
## Finding best probes
#########################
score <- function(x) geoSD(x)^(1/geoMean(x))

probe_scores <- expression_long %>%
  filter(disease == "Ewing") %>%
  group_by(probe_id) %>%
  summarise(value = score(value)) %>%
  inner_join(probe_mapping, by = c("probe_id" = "PROBEID"))

best_probes <- probe_scores %>%
  group_by(ALIAS) %>%
  mutate(best = value == max(value)) %>%
  select(probe_id, best, gene_symbol = ALIAS)

expression_long %<>% inner_join(best_probes)


#####################
## Plotting
#####################
# Plotting order
expression_long %<>% mutate(
   gene_symbol = fct_relevel(gene_symbol, gene_names)
  ,disease     = fct_relevel(disease, geo_disease_mapping)
)

# Selecting best probes
expression_long_best_probes <- filter(expression_long, best == TRUE)

# Base plot elements
base_plot <- ggplot() +
  geom_boxplot(aes(x = disease, y = value)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Disease") +
  ylab("Normalized probe expression")

# All probes
base_plot %+%
  expression_long +
  facet_wrap(gene_symbol ~ paste0(probe_id, ifelse(best, " (best)", "")))

# Saving PDF
ggsave("normalized_expression_all_probes.pdf", width = 5, height = 8, useDingbats = F)

# Best probes only
base_plot %+%
  expression_long_best_probes +
  facet_wrap(gene_symbol ~ ., ncol = 2)

# Saving PDF
ggsave("normalized_expression_best_probes.pdf", width = 5, height = 8, useDingbats = F)


########################
## Statistical testing
########################
best_probes_summary <- expression_long_best_probes %>%
  group_by(disease, gene_symbol) %>%
  summarise(
    count  = n(),
    mean   = mean(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    IQR    = IQR(value, na.rm = TRUE)
  )

# For each gene, pairwise test its expression between all diseases
pairwise_tests_by_gene <- expression_long_best_probes %>%
  group_by(gene_symbol) %>%
  summarise(
    test = pairwise.wilcox.test(value, disease, p.adjust.method = "bonf") %>% list
  ) %>%
  mutate(test = map(test, tidy)) %>%
  unnest(test)

# Print diseases in which gene expression
# is statistically significant for ALL pairwise combinations (p < 0.01)
significant_tests <- pairwise_tests_by_gene %>%
  pivot_longer(starts_with("group"), values_to = "disease") %>%
  group_by(gene_symbol, disease) %>%
  summarise(significant = all(p.value < 0.01))

significant_tests %>% filter(significant == TRUE)
