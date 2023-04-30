library(janitor)
library(tidyverse)

cnaDF <- read.table(
  file = "resources/SCNA.gz",
  header = T,
  sep = "\t",
  comment.char = "",
  quote = "",
  stringsAsFactors = F
) %>%
  clean_names()

cnaLongDF <- cnaDF  %>% 
  select(
    gene_id = id,
    starts_with("cpt"),
    starts_with("x")
  ) %>% 
  pivot_longer(
    cols = starts_with(c("x", "cpt")),
    names_to = "tumor",
    values_to = "cna_abundance"
  )

rnaDF <- read.table(
  file = "resources/mRNA.gz",
  header = T,
  sep = "\t",
  comment.char = "",
  quote = "",
  stringsAsFactors = F
) %>%
  clean_names()

proteinDF <- read.table(
  file = "resources/proteins.gz",
  header = T,
  sep = "\t",
  comment.char = "",
  quote = "",
  stringsAsFactors = F
) %>%
  clean_names()

proteinLongDF <- proteinDF  %>% 
  select(
    protein_id = id,
    gene_id = gene_symbol,
    starts_with("cpt"),
    starts_with("x")
  ) %>% 
  pivot_longer(
    cols = starts_with(c("x", "cpt")),
    names_to = "tumor",
    values_to = "protein_abundance"
  )

rnaLongDF <- rnaDF  %>% 
  select(
    gene_id = id,
    gene_description = description,
    starts_with("cpt"),
    starts_with("x")
  ) %>% 
  pivot_longer(
    cols = starts_with(c("x", "cpt")),
    names_to = "tumor",
    values_to = "mRNA_abundance"
  )

cnaRnaProteinDF <- cnaLongDF %>% 
  full_join(
    rnaLongDF,
    by = c("gene_id", "tumor")
  ) %>% 
  full_join(
    proteinLongDF,
    by = c("gene_id", "tumor")
  )

cnaRnaProteinCorrelationDF <- cnaRnaProteinDF  %>% 
  group_by(
    protein_id
  ) %>% 
  summarize(
    r_cna_rna = cor(
      x = cna_abundance, 
      y = mRNA_abundance, 
      use = "pairwise.complete.obs", 
      method = "spearman"
    ),
    r_cna_protein = cor(
      x = cna_abundance, 
      y = protein_abundance, 
      use = "pairwise.complete.obs", 
      method = "spearman"
    )
  )

correlationDF <- cnaRnaProteinCorrelationDF  %>%
  pivot_longer(
    cols = c("r_cna_rna", "r_cna_protein"),
    names_to = "comparison",
    values_to = "r"
  )  %>% 
  mutate(
    comparison = factor(comparison, levels = c("r_cna_rna", "r_cna_protein"))
  ) %>% 
  filter(
    !is.na(r)
  )

write.table(
  x = correlationDF,
  file = gzfile("resources/correlation.gz"),
  col.names = T,
  row.names = F,
  sep = "\t"
)

