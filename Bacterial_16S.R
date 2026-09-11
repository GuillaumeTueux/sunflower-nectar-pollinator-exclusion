# 16S bacterial community: phyloseq import, contamination filtering
# (extraction blanks + PCR negative controls), rarefaction check, genus-level
# composition (Acinetobacter / chloroplast focus), and beta-diversity
# (PERMANOVA, UniFrac PCoA) per cultivar.

library(phyloseq)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(iNEXT)
library(tibble)
library(vegan)
library(car)
library(emmeans)
library(purrr)
library(scales)
library(forcats)
library(ggtext)
library(patchwork)

pkgs <- c("phyloseq", "ape", "dplyr", "tidyr", "ggplot2", "iNEXT", "tibble",
          "vegan", "car", "emmeans", "purrr", "scales", "forcats")
pkg_versions <- data.frame(
  package = pkgs,
  version = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)),
  stringsAsFactors = FALSE
)
pkg_versions

options(scipen = 999)
set.seed(64121)

theme_natcomm <- theme_bw(base_size = 35) +
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 25, face = "italic"),
    axis.text.y = element_text(size = 34),
    legend.title = element_text(size = 35),
    legend.text = element_text(size = 34),
    strip.text = element_text(size = 36, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(4, 4, 4, 4, "mm")
  )

italicize_genus <- function(x, exclude = c("Chloroplast", "Mitochondria", "Other", "Unknown Genus")) {
  dplyr::if_else(x %in% exclude, x, paste0("*", x, "*"))
}

# --- Import data -------------------------------------------------------
# otu_table.txt / metadata.txt / taxonomy.tsv / tree.nwk: QIIME2 exports
data_otu <- read.table("data/otu_table.txt", header = TRUE, sep = "\t",
                       fileEncoding = "latin1", check.names = FALSE, row.names = 1)

data_grp <- read.table("data/metadata.txt", header = TRUE, sep = "\t",
                       stringsAsFactors = TRUE, row.names = 1)
data_grp$Trial <- factor(substr(rownames(data_grp), 1, 6))

data_taxo <- read.table("data/taxonomy.tsv", header = TRUE, sep = "\t", quote = "",
                        comment.char = "", fill = TRUE, row.names = 1, check.names = FALSE)

MyTree <- read.tree("data/tree.nwk")

# --- Build the phyloseq object ---------------------------------------------
data_taxo_separated <- data_taxo %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", remove = FALSE, fill = "right", extra = "drop") %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), ~ sub("^.*__", "", .))) %>%
  filter(Kingdom == "Bacteria") %>%
  select(-Taxon)

physeq <- phyloseq(
  otu_table(as.matrix(data_otu), taxa_are_rows = TRUE),
  tax_table(as.matrix(data_taxo_separated)),
  sample_data(data_grp, errorIfNULL = TRUE),
  MyTree
)

physeq
sample_names(physeq)

data_grp %>%
  filter(genotype %in% c("CELESTO", "IDILLIC")) %>%
  count(genotype, Trial, modality)

# extract OTU/taxonomy table for a given sample (reads > 0 only)
extract_otu_tax <- function(sid) {
  otu_counts <- as.integer(otu_table(physeq)[, sid])
  cbind(
    data.frame(OTU = taxa_names(physeq), reads = otu_counts),
    as.data.frame(tax_table(physeq))
  ) %>%
    as_tibble() %>%
    filter(reads > 0) %>%
    arrange(desc(reads)) %>%
    mutate(sample = sid) %>%
    select(sample, everything())
}

# --- Mock community: check recovery, then remove ---------------------------
mock_samples <- c("Mock_2025", "Mock_2024", "Mock_2023")
stopifnot(all(mock_samples %in% sample_names(physeq)))

otu_tax_mocks <- bind_rows(lapply(mock_samples, extract_otu_tax))
print(otu_tax_mocks, n = Inf)

physeq <- prune_samples(!(sample_names(physeq) %in% mock_samples), physeq)
physeq_raw <- physeq

# --- Contamination filtering -------------------------------------------------
water_samples <- c("Water_2024", "Water_2025", "Water_2025_2")
negative_samples <- c("24MO11_T1", "24MO11_T2", "25ZM05_T1", "25ZM05_T2")
stopifnot(all(c(water_samples, negative_samples) %in% sample_names(physeq)))

# OTUs exceeding 5% relative abundance in any extraction blank are removed,
# except chloroplast OTUs (co-amplified plant plastidial 16S, not contamination)
otu_mat <- as(otu_table(physeq), "matrix")
if (!taxa_are_rows(physeq)) otu_mat <- t(otu_mat)

water_mat <- otu_mat[, water_samples, drop = FALSE]
water_tot <- colSums(water_mat)
water_prop <- sweep(water_mat, 2, water_tot, FUN = "/")
water_prop[, water_tot == 0] <- 0

taxa_to_remove <- rownames(water_prop)[apply(water_prop, 1, function(x) any(x > 0.05))]

chloro_otus <- rownames(tax_table(physeq))[tax_table(physeq)[, "Order"] == "Chloroplast"]
taxa_to_remove <- setdiff(taxa_to_remove, chloro_otus)

cat("OTUs removed after chloroplast exclusion:\n")
print(tax_table(physeq)[taxa_to_remove, ])

physeq <- prune_taxa(!(taxa_names(physeq) %in% taxa_to_remove), physeq)
physeq <- prune_samples(!(sample_names(physeq) %in% water_samples), physeq)

# inspect PCR negative controls before removing them
otu_tax_negatives <- bind_rows(lapply(negative_samples, extract_otu_tax))
print(otu_tax_negatives, n = Inf)

physeq <- prune_samples(!(sample_names(physeq) %in% negative_samples), physeq)

# minimum read/OTU thresholds
physeq <- prune_samples(sample_sums(physeq) >= 1000, physeq)
physeq <- prune_taxa(taxa_sums(physeq) >= 10, physeq)

reads_per_sample <- sample_sums(physeq)
otus_per_sample <- estimate_richness(physeq, measures = "Observed")$Observed

summary(reads_per_sample)
summary(otus_per_sample)
mean(reads_per_sample)
mean(otus_per_sample)

data.frame(SampleID = sample_names(physeq), Reads = reads_per_sample, OTUs = otus_per_sample)

n_samples_genotype_modality <- data.frame(sample_data(physeq)) %>%
  mutate(genotype = factor(genotype), modality = factor(modality)) %>%
  count(genotype, modality, name = "n_samples") %>%
  arrange(genotype, modality)
n_samples_genotype_modality

# relative-abundance version of the filtered dataset, reused throughout
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
sample_names(physeq_rel) <- sub("_16S$", "", sample_names(physeq_rel))

# --- Rarefaction curves ------------------------------------------------------
ps_rare <- prune_samples(sample_sums(physeq) > 0, physeq)
mat <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) mat <- t(mat)

final_depth <- sample_sums(ps_rare)[rownames(mat)]
step <- 250
df_list <- vector("list", nrow(mat))
names(df_list) <- rownames(mat)

for (s in rownames(mat)) {
  x <- mat[s, ]
  maxd <- sum(x)
  d <- as.numeric(final_depth[s])
  if (is.na(d) || d < 2 || maxd < 2) next
  d <- min(d, maxd)
  depths <- unique(pmax(1, c(seq(1, d, by = step), d)))
  r <- suppressWarnings(vegan::rarefy(x, sample = depths))
  df_list[[s]] <- data.frame(sample = s, depth = depths, richness = as.numeric(r), final = d)
}

rare_df <- bind_rows(df_list) %>% filter(!is.na(richness))

ggplot(rare_df, aes(depth, richness, group = sample)) +
  geom_line(linewidth = 0.4, alpha = 0.4) +
  geom_point(data = rare_df %>% group_by(sample) %>% slice_max(depth, n = 1, with_ties = FALSE),
             aes(depth, richness), size = 1.2, alpha = 0.7) +
  labs(x = "Reads", y = "Expected richness") +
  theme_bw() + theme_natcomm

# --- Figure S3, panel A: genus-level composition per sample ----------------
df16S <- psmelt(physeq_rel)

dat16S <- df16S %>%
  mutate(Genus = if_else(is.na(Genus) | Genus == "", "Unknown Genus", as.character(Genus))) %>%
  group_by(Sample, modality, genotype, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample, modality, genotype) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup()

focal_genera <- c("Acinetobacter", "Chloroplast")

dat16S_plot <- dat16S %>%
  mutate(Genus = if_else(Genus %in% focal_genera, Genus, "Other")) %>%
  group_by(Sample, modality, genotype, Genus) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  mutate(
    Genus = factor(Genus, levels = c("Acinetobacter", "Chloroplast", "Other")),
    modality = factor(as.character(modality),
                      levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access"))
  )

genus_labels_A <- setNames(italicize_genus(levels(dat16S_plot$Genus)), levels(dat16S_plot$Genus))

barplot16S <- ggplot(dat16S_plot, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Acinetobacter" = "#E69F00", "Chloroplast" = "#4DAF4A", "Other" = "#BDBDBD"),
                    labels = genus_labels_A, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0),
                     name = "Relative abundance (%)") +
  labs(x = "Samples", fill = "Genus") +
  facet_wrap(genotype ~ modality, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 25, margin = margin(t = 20)),
        axis.title.y = element_text(size = 25, margin = margin(r = 25)),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 28),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.title = element_text(size = 30),
        legend.text = ggtext::element_markdown(size = 28),
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(0.8, "lines"),
        strip.text.x = element_text(size = 9, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"),
        plot.margin = margin(5, 5, 5, 5, "mm"))
barplot16S

top5_16S <- dat16S %>%
  group_by(modality, genotype, Genus) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  group_by(modality, genotype) %>%
  slice_max(order_by = MeanRelAbund, n = 5) %>%
  arrange(modality, genotype, desc(MeanRelAbund)) %>%
  mutate(MeanRelAbund = round(MeanRelAbund * 100, 1))
print(top5_16S, n = Inf)

dat16S %>%
  filter(modality == "No_Visitors", genotype == "CELESTO", Genus == "Chloroplast") %>%
  summarise(n = n(), min = min(RelAbund) * 100, max = max(RelAbund) * 100, mean = mean(RelAbund) * 100)

# --- Figure S3, panel B: genus composition with chloroplast subtracted -----
chloro_otus_final <- rownames(tax_table(physeq))[tax_table(physeq)[, "Order"] == "Chloroplast"]
physeq_noChloro <- prune_taxa(!(taxa_names(physeq) %in% chloro_otus_final), physeq)
stopifnot(all(sample_sums(physeq_noChloro) > 0))

physeq_noChloro_rel <- transform_sample_counts(physeq_noChloro, function(x) x / sum(x))
sample_names(physeq_noChloro_rel) <- sub("_16S$", "", sample_names(physeq_noChloro_rel))

df16S_noChloro <- psmelt(physeq_noChloro_rel) %>%
  mutate(Genus = if_else(is.na(Genus) | Genus == "", "Unknown Genus", as.character(Genus)))

dat16S_genus <- df16S_noChloro %>%
  group_by(Sample, modality, genotype, Genus) %>%
  summarise(RelAbund = sum(Abundance), .groups = "drop")

top_genera <- dat16S_genus %>%
  group_by(Genus) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  slice_max(order_by = MeanRelAbund, n = 10) %>%
  pull(Genus)

dat16S_genus_plot <- dat16S_genus %>%
  mutate(Genus = if_else(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(Sample, modality, genotype, Genus) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  mutate(
    Genus = factor(Genus, levels = c(top_genera, "Other")),
    modality = factor(as.character(modality),
                      levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access"))
  )

genus_palette <- c(scales::hue_pal()(length(top_genera)), "grey70")
names(genus_palette) <- c(top_genera, "Other")
genus_labels_B <- setNames(italicize_genus(levels(dat16S_genus_plot$Genus)), levels(dat16S_genus_plot$Genus))

barplot16S_genus <- ggplot(dat16S_genus_plot, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = genus_palette, labels = genus_labels_B, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0),
                     name = "Relative abundance (%)") +
  labs(x = "Samples", fill = "Genus") +
  facet_wrap(genotype ~ modality, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 25, margin = margin(t = 20)),
        axis.title.y = element_text(size = 25, margin = margin(r = 25)),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 28),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.title = element_text(size = 30),
        legend.text = ggtext::element_markdown(size = 28),
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(0.8, "lines"),
        strip.text.x = element_text(size = 9, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"),
        plot.margin = margin(5, 5, 5, 5, "mm"))
barplot16S_genus

FigureS3 <- barplot16S / barplot16S_genus &
  theme(plot.tag = element_text(face = "italic"))
FigureS3

# --- Acinetobacter and chloroplast: relative abundance by treatment --------
df_target <- psmelt(physeq_rel) %>%
  mutate(Genus = as.character(Genus), Order = as.character(Order),
         modality = as.character(modality), genotype = as.character(genotype)) %>%
  mutate(target = case_when(
    Order == "Chloroplast" ~ "Chloroplast",
    Genus == "Acinetobacter" ~ "Acinetobacter",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(target)) %>%
  group_by(Sample, modality, genotype, target) %>%
  summarise(pct = sum(Abundance) * 100, .groups = "drop")

all_samples <- psmelt(physeq_rel) %>% distinct(Sample, modality, genotype)

df_target_full <- expand_grid(all_samples, target = c("Chloroplast", "Acinetobacter")) %>%
  left_join(df_target, by = c("Sample", "modality", "genotype", "target")) %>%
  mutate(
    pct = replace_na(pct, 0),
    modality = gsub("_", " ", modality),
    modality = factor(modality, levels = c("No Visitors", "All Visitors", "Day Visitors", "Night Visitors")),
    target = factor(target, levels = c("Chloroplast", "Acinetobacter"))
  )

# mean relative abundance per treatment x cultivar (cited in Results)
prop_means <- df_target_full %>%
  group_by(modality, genotype, target) %>%
  summarise(mean_pct = mean(pct), sd_pct = sd(pct), n_samples = n(), .groups = "drop") %>%
  arrange(target, modality, genotype)
prop_means

df_summary <- df_target_full %>%
  group_by(modality, genotype, target) %>%
  summarise(mean_pct = mean(pct), se_pct = sd(pct) / sqrt(n()), .groups = "drop") %>%
  mutate(
    modality = factor(modality,
                      levels = c("No Visitors", "All Visitors", "Day Visitors", "Night Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
    target = factor(target, levels = c("Chloroplast", "Acinetobacter"))
  )

df_target_full <- df_target_full %>%
  mutate(modality = factor(modality,
                           levels = c("No Visitors", "All Visitors", "Day Visitors", "Night Visitors"),
                           labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")))

pd <- position_dodge(width = 0.8)

p <- ggplot(df_summary, aes(x = modality, y = mean_pct, fill = target)) +
  geom_col(position = pd, width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
                position = pd, width = 0.25, linewidth = 0.6) +
  geom_point(data = df_target_full, aes(x = modality, y = pct, group = target),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             shape = 21, size = 2.5, alpha = 0.85, fill = "white", color = "black", inherit.aes = FALSE) +
  facet_wrap(~ genotype) +
  scale_fill_manual(values = c("Chloroplast" = "#4DAF4A", "Acinetobacter" = "#E69F00"),
                    labels = c("Chloroplast", expression(italic("Acinetobacter")))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 100)) +
  labs(x = NULL, y = "Relative abundance (%)", fill = NULL) +
  theme_natcomm +
  theme(axis.text.x = element_text(size = 25, face = "plain", angle = 45, hjust = 1),
        legend.position = "right")
p

# --- Acinetobacter OTU identity checks -----------------------------------
tax_table(physeq) %>% as.data.frame() %>% filter(Genus == "Acinetobacter")

psmelt(physeq_rel) %>%
  filter(Genus == "Acinetobacter") %>%
  group_by(OTU, modality, genotype) %>%
  summarise(mean_relabund = mean(Abundance), .groups = "drop") %>%
  arrange(modality, genotype, desc(mean_relabund)) %>%
  print(n = Inf)

psmelt(physeq_rel) %>%
  filter(Order == "Chloroplast") %>%
  group_by(OTU, modality, genotype) %>%
  summarise(mean_relabund = mean(Abundance), .groups = "drop") %>%
  arrange(modality, genotype, desc(mean_relabund)) %>%
  print(n = Inf)

# share of total reads carried by the two dominant OTUs
psmelt(physeq_rel) %>%
  mutate(dominant = OTU %in% c("77415788a2999e145f594d8f3738cd02e1e334f8",
                               "90c39638c0a1a1b6a2ca253c711811238359c694")) %>%
  group_by(Sample, dominant) %>%
  summarise(RelAbund = sum(Abundance), .groups = "drop") %>%
  group_by(dominant) %>%
  summarise(mean = mean(RelAbund), .groups = "drop")

psmelt(physeq_rel) %>%
  mutate(group = case_when(
    OTU == "77415788a2999e145f594d8f3738cd02e1e334f8" ~ "Acinetobacter_dominant",
    Genus == "Acinetobacter" ~ "Acinetobacter_other",
    OTU == "90c39638c0a1a1b6a2ca253c711811238359c694" ~ "Chloroplast_dominant",
    Order == "Chloroplast" ~ "Chloroplast_other",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(total_relabund = sum(Abundance), .groups = "drop") %>%
  mutate(prop_within_group = case_when(
    grepl("Acinetobacter", group) ~ total_relabund / sum(total_relabund[grepl("Acinetobacter", group)]),
    grepl("Chloroplast", group) ~ total_relabund / sum(total_relabund[grepl("Chloroplast", group)])
  ))

# dominant Acinetobacter OTU in Celesto (continuous + daytime access)
dominant_acineto_otu <- psmelt(physeq_rel) %>%
  filter(Genus == "Acinetobacter", genotype == "CELESTO", modality %in% c("All_Visitors", "Day_Visitors")) %>%
  group_by(OTU) %>%
  summarise(mean_relabund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_relabund)) %>%
  slice(1) %>%
  pull(OTU)

dominant_acineto_otu  # expected: "77415788a2999e145f594d8f3738cd02e1e334f8"
tax_table(physeq)[dominant_acineto_otu, ]

# same OTU in Idillic: detection rate and abundance range when detected
ac1_detection_by_modality <- psmelt(physeq_rel) %>%
  filter(OTU == dominant_acineto_otu, genotype == "IDILLIC", modality %in% c("All_Visitors", "Day_Visitors")) %>%
  group_by(modality) %>%
  summarise(n_detected = sum(Abundance > 0), n_total = n(), .groups = "drop")

ac1_detection_combined <- psmelt(physeq_rel) %>%
  filter(OTU == dominant_acineto_otu, genotype == "IDILLIC", modality %in% c("All_Visitors", "Day_Visitors")) %>%
  summarise(n_detected = sum(Abundance > 0), n_total = n())

ac1_detection_by_modality
ac1_detection_combined

ac2_range_detected_only <- psmelt(physeq_rel) %>%
  filter(OTU == dominant_acineto_otu, genotype == "IDILLIC",
         modality %in% c("All_Visitors", "Day_Visitors"), Abundance > 0) %>%
  group_by(modality) %>%
  summarise(n_detected = n(), min_relabund = min(Abundance) * 100,
            max_relabund = max(Abundance) * 100, mean_relabund = mean(Abundance) * 100, .groups = "drop")
print(ac2_range_detected_only, n = Inf, digits = 6)

# --- Beta-diversity, per cultivar -------------------------------------------
physeq_counts <- prune_samples(sample_sums(physeq) > 0, physeq)
physeq_counts_rel <- transform_sample_counts(physeq_counts, function(x) x / sum(x))

sample_df_rel <- as(sample_data(physeq_counts_rel), "data.frame")
sample_df_counts <- as(sample_data(physeq_counts), "data.frame")

for (geno in c("IDILLIC", "CELESTO")) {
  cat("\n\n##################################################\n")
  cat("# Genotype:", geno, "\n")
  cat("##################################################\n")
  
  ps_g <- subset_samples(physeq_counts, genotype == geno)
  ps_g <- prune_samples(sample_sums(ps_g) > 0, ps_g)
  ps_g_rel <- transform_sample_counts(ps_g, function(x) x / sum(x))
  
  sdf_rel <- as(sample_data(ps_g_rel), "data.frame")
  sdf_counts <- as(sample_data(ps_g), "data.frame")
  
  d_bray <- phyloseq::distance(ps_g_rel, "bray")
  
  otu_g <- as(otu_table(ps_g), "matrix")
  if (taxa_are_rows(ps_g)) otu_g <- t(otu_g)
  d_jacc <- vegdist(otu_g, method = "jaccard", binary = TRUE)
  
  d_uuni <- phyloseq::UniFrac(ps_g, weighted = FALSE, normalized = TRUE)
  d_wuni <- phyloseq::UniFrac(ps_g, weighted = TRUE, normalized = TRUE)
  
  cat("\n--- Betadisper ---\n")
  for (nm in c("bray", "jacc", "uuni", "wuni")) {
    d <- get(paste0("d_", nm))
    bd <- betadisper(d, sdf_counts$modality)
    pt <- permutest(bd, permutations = 9999)
    cat(nm, "~ modality: p =", pt$tab$`Pr(>F)`[1], "\n")
  }
  
  cat("\n--- PERMANOVA ---\n")
  for (nm in c("bray", "jacc", "uuni", "wuni")) {
    d <- get(paste0("d_", nm))
    sdf <- if (nm == "bray") sdf_rel else sdf_counts
    cat(nm, ":\n")
    print(adonis2(d ~ modality + Trial, data = sdf, permutations = 99999, by = "margin"))
  }
}

# --- PCoA, unweighted UniFrac, per cultivar ---------------------------------
pcoa_list <- list()

for (geno in c("CELESTO", "IDILLIC")) {
  ps_g <- subset_samples(physeq_counts, genotype == geno)
  ps_g <- prune_samples(sample_sums(ps_g) > 0, ps_g)
  
  d_uuni <- phyloseq::UniFrac(ps_g, weighted = FALSE, normalized = TRUE)
  pco <- cmdscale(d_uuni, eig = TRUE, k = 2)
  var_exp <- round(pco$eig / sum(pco$eig[pco$eig > 0]) * 100, 1)
  
  df_pco <- data.frame(PCoA1 = pco$points[, 1], PCoA2 = pco$points[, 2], Sample = rownames(pco$points)) %>%
    left_join(data.frame(sample_data(ps_g)) %>% rownames_to_column("Sample"), by = "Sample") %>%
    mutate(genotype = geno,
           ax1_lab = paste0("PCoA1 (", var_exp[1], "%)"),
           ax2_lab = paste0("PCoA2 (", var_exp[2], "%)"))
  
  pcoa_list[[geno]] <- df_pco
}

df_pcoa_uuni <- bind_rows(pcoa_list) %>%
  mutate(
    panel = factor(genotype, levels = c("CELESTO", "IDILLIC")),
    modality = factor(modality,
                      levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
    Trial = factor(Trial, levels = c("24MO11", "25ZM05"), labels = c("2024", "2025"))
  )

axis_labels <- df_pcoa_uuni %>% distinct(panel, ax1_lab, ax2_lab)

modality_colours <- c(
  "No access" = "#000000", "Continuous access" = "#CC79A7",
  "Daytime access" = "#F5C400", "Nighttime access" = "#0072B2"
)

ggplot(df_pcoa_uuni, aes(x = PCoA1, y = PCoA2, colour = modality)) +
  geom_point(aes(shape = Trial), size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = modality), type = "t", level = 0.95, linewidth = 0.6) +
  scale_colour_manual(values = modality_colours) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~ genotype, scales = "free") +
  labs(x = "PCoA1", y = "PCoA2", colour = "Treatment", shape = "Trial Year") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme_bw(base_size = 35) +
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 34),
    legend.title = element_text(size = 35),
    legend.text = element_text(size = 34),
    strip.text = element_text(size = 36, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(4, 4, 4, 4, "mm")
  )

# --- Chloroplast presence check, objects used in the PERMANOVA loop --------
sum(tax_table(physeq_counts)[, "Order"] == "Chloroplast", na.rm = TRUE)

for (geno in c("IDILLIC", "CELESTO")) {
  ps_g <- subset_samples(physeq_counts, genotype == geno)
  ps_g <- prune_samples(sample_sums(ps_g) > 0, ps_g)
  n_chloro <- sum(tax_table(ps_g)[, "Order"] == "Chloroplast", na.rm = TRUE)
  cat(geno, ": n OTU Chloroplast =", n_chloro, "\n")
}