library(janitor)
library(tidyverse)
library(mclust)

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

cnaRnaProteinCorrelationDF <- cnaRnaProteinCorrelationDF %>% 
  filter(
    !is.na(r_cna_rna) & !is.na(r_cna_protein)
  ) %>%
  mutate(
    attenuation_coefficient = r_cna_protein - r_cna_rna
  )

gmm <- densityMclust(cnaRnaProteinCorrelationDF$attenuation_coefficient)

summary(gmm, parameters = TRUE)

cnaRnaProteinCorrelationDF$attenuation_p <- pnorm(
  q = cnaRnaProteinCorrelationDF$attenuation_coefficient,
  mean = gmm$parameters$mean[3],
  sd = sqrt(gmm$parameters$variance$sigmasq[3])
)

cnaRnaProteinCorrelationDF$attenuation_category <- factor(
  x = ifelse(cnaRnaProteinCorrelationDF$attenuation_p <= 0.05, "Attenuated", "Background"),
  levels = c("Attenuated", "Background")
)

write.table(
  x = cnaRnaProteinCorrelationDF,
  file = gzfile("resources/attenuation.gz"),
  col.names = T,
  row.names = F,
  sep = "\t"
)


