# ITS fungal metabarcoding: import and contamination filtering, fungal alpha
# diversity (Hill numbers) across pollinator access treatments, beta
# diversity (PERMANOVA/UniFrac) and community composition, and camera-trap
# visit rates as continuous predictors of microbiome structure.
#
# Filtering strategy: water/extraction blanks at 5% global threshold,
# sampling blanks ("tem") at 5% per trial, Malasseziomycetes exclusion,
# OTUs with <10 total reads and samples with <1000 reads discarded.

library(phyloseq)
library(ape)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)
library(car)
library(emmeans)
library(multcomp)
library(iNEXT)
library(ggh4x)
library(ggtext)
library(patchwork)
library(pairwiseAdonis)
library(purrr)

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

pal_modality <- c(
  "No access" = "#000000", "Continuous access" = "#CC79A7",
  "Daytime access" = "#F5C400", "Nighttime access" = "#0072B2"
)

# =============================================================================
# 1. Import and filtering
# =============================================================================

# --- Import data -------------------------------------------------------
# otu_table_merged.txt / taxonomy_merged.tsv: post Metschnikowia-merge
#   feature table and taxonomy (see merge_metschnikowia_otus.py)
# metadata.txt: sample metadata, including "modality" (sampling-blank rows
#   coded "tem") and "genotype"
# tree.nwk: representative-sequence tree, pruned to the samples present
data_otu <- read.table("data/otu_table_merged.txt", header = TRUE, sep = "\t",
                       fileEncoding = "latin1", check.names = FALSE, row.names = 1)

data_grp <- read.table("data/metadata.txt", header = TRUE, sep = "\t",
                       stringsAsFactors = TRUE, row.names = 1)
data_grp$Trial <- factor(substr(rownames(data_grp), 1, 6))

data_taxo <- read.table("data/taxonomy_merged.tsv", header = TRUE, fill = TRUE, row.names = 1)

MyTree <- read.tree("data/tree.nwk")
tree_merged <- drop.tip(MyTree, setdiff(MyTree$tip.label, rownames(data_otu)))

data_taxo_separated <- data_taxo %>%
  separate(ID, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", remove = FALSE) %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), ~ sub("^.*__", "", .))) %>%
  select(-ID) %>%
  filter(Kingdom == "Fungi")

physeq <- phyloseq(
  otu_table(as.matrix(data_otu), taxa_are_rows = TRUE),
  tax_table(as.matrix(data_taxo_separated)),
  sample_data(data_grp, errorIfNULL = TRUE),
  tree_merged
)

physeq
sample_names(physeq)

data_grp %>%
  filter(genotype %in% c("SY_CELESTO", "ES_IDILLIC")) %>%
  count(genotype, Trial, modality)

extract_otu_tax <- function(ps, sid) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  tibble(sample = sid, OTU = rownames(otu), reads = as.integer(otu[, sid])) %>%
    filter(reads > 0) %>%
    left_join(as.data.frame(as.matrix(tax_table(ps))) %>% rownames_to_column("OTU"), by = "OTU") %>%
    arrange(desc(reads))
}

mock_samples <- c("Mock_2025", "Mock_2024", "Mock_2023")
stopifnot(all(mock_samples %in% sample_names(physeq)))
otu_tax_mocks <- bind_rows(lapply(mock_samples, \(s) extract_otu_tax(physeq, s)))

physeq <- prune_samples(!(sample_names(physeq) %in% mock_samples), physeq)

water_samples <- c("Water_2024", "Water_2025")
stopifnot(all(water_samples %in% sample_names(physeq)))
otu_tax_waters <- bind_rows(lapply(water_samples, \(s) extract_otu_tax(physeq, s)))

# --- Water/extraction blank filter: 5% global, OTU removed dataset-wide ----
otu_mat <- as(otu_table(physeq), "matrix")
if (!taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
thr_rel <- 0.05

water_in_ps <- intersect(water_samples, sample_names(physeq))
water_counts <- otu_mat[, water_in_ps, drop = FALSE]
water_libs <- colSums(water_counts)
water_libs[water_libs == 0] <- 1L
water_prop <- sweep(water_counts, 2, water_libs, "/")
in_water <- apply(water_prop, 1, function(x) any(x > thr_rel))

removed_otus_water <- names(in_water[in_water])
cat(length(removed_otus_water), "OTUs removed by the water filter (global)\n")

physeq_pw <- prune_taxa(!in_water, physeq)
physeq_pw <- prune_samples(!(sample_names(physeq_pw) %in% water_in_ps), physeq_pw)
cat("After water filter — taxa:", ntaxa(physeq_pw), "| samples:", nsamples(physeq_pw), "\n")

# --- Sampling blank filter: 5% per trial, OTU zeroed only within that trial's samples ---
sd_pw <- data.frame(sample_data(physeq_pw))
sd_pw$SampleID <- rownames(sd_pw)

neg_all <- sd_pw$SampleID[tolower(as.character(sd_pw$modality)) == "tem"]
cat("\nSampling blanks detected:", length(neg_all), "\n")

otu_pw <- as(otu_table(physeq_pw), "matrix")
if (!taxa_are_rows(physeq_pw)) otu_pw <- t(otu_pw)

remove_mat <- matrix(FALSE, nrow = nrow(otu_pw), ncol = ncol(otu_pw), dimnames = dimnames(otu_pw))

for (tr in unique(as.character(sd_pw$Trial))) {
  samples_tr <- sd_pw$SampleID[as.character(sd_pw$Trial) == tr]
  neg_tr <- intersect(samples_tr, neg_all)
  ech_tr <- setdiff(samples_tr, neg_tr)
  if (length(neg_tr) == 0) next
  
  neg_counts <- otu_pw[, neg_tr, drop = FALSE]
  libs <- colSums(neg_counts)
  libs[libs == 0] <- 1L
  neg_prop <- sweep(neg_counts, 2, libs, "/")
  in_neg_tr <- apply(neg_prop, 1, function(x) any(x > thr_rel))
  
  if (any(in_neg_tr) && length(ech_tr) > 0) remove_mat[in_neg_tr, ech_tr] <- TRUE
  cat("Trial", tr, ":", length(neg_tr), "blanks,", length(ech_tr), "samples,",
      sum(in_neg_tr), "OTUs removed locally\n")
}

otu_pw_filtered <- otu_pw
otu_pw_filtered[remove_mat] <- 0L
otu_back <- if (taxa_are_rows(physeq_pw)) otu_pw_filtered else t(otu_pw_filtered)
otu_table(physeq_pw) <- otu_table(otu_back, taxa_are_rows = taxa_are_rows(physeq_pw))

physeq_tmp <- prune_samples(!(sample_names(physeq_pw) %in% neg_all), physeq_pw)

# --- Malasseziomycetes exclusion and minimum-read thresholds ----------------
physeq_final <- prune_taxa(taxa_sums(physeq_tmp) > 0, physeq_tmp)
cat("\nSamples after blank filters (no controls):", nsamples(physeq_final), "\n")
cat("Taxa after blank filters:", ntaxa(physeq_final), "\n")

physeq_final <- subset_taxa(physeq_final, is.na(Class) | Class != "Malasseziomycetes")

otu_tot <- taxa_sums(physeq_final)
keep_otu <- otu_tot >= 10
physeq_final <- prune_taxa(keep_otu, physeq_final)

repeat {
  n_samp_before <- nsamples(physeq_final)
  n_otu_before <- ntaxa(physeq_final)
  physeq_final <- prune_samples(sample_sums(physeq_final) > 1000, physeq_final)
  physeq_final <- prune_taxa(taxa_sums(physeq_final) >= 10, physeq_final)
  n_samp_after <- nsamples(physeq_final)
  n_otu_after <- ntaxa(physeq_final)
  cat("Samples:", n_samp_before, "->", n_samp_after, " | OTUs:", n_otu_before, "->", n_otu_after, "\n")
  if (n_samp_after == n_samp_before && n_otu_after == n_otu_before) break
}

cat("\nFinal — samples:", nsamples(physeq_final), "| OTUs:", ntaxa(physeq_final), "\n")

sd_final <- data.frame(sample_data(physeq_final)) %>%
  mutate(genotype = factor(genotype), modality = factor(modality))

n_samples_genotype_modality <- sd_final %>%
  count(genotype, modality, name = "n_samples") %>%
  arrange(genotype, modality)
n_samples_genotype_modality

otu25 <- as(otu_table(physeq_final), "matrix")
if (!taxa_are_rows(physeq_final)) otu25 <- t(otu25)

meta25 <- data.frame(sample_data(physeq_final))
meta25$SampleID <- rownames(meta25)

# --- Rarefaction curves ------------------------------------------------------
ps_rare <- prune_samples(sample_sums(physeq_final) > 0, physeq_final)
mat_rare <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) mat_rare <- t(mat_rare)
final_depth <- sample_sums(ps_rare)[rownames(mat_rare)]

step <- 250
df_list <- vector("list", nrow(mat_rare))
names(df_list) <- rownames(mat_rare)

for (s in rownames(mat_rare)) {
  x <- mat_rare[s, ]
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
  labs(x = "Subsampled reads (post-filter depth)", y = "Expected richness") +
  theme_bw() + theme_natcomm

# =============================================================================
# 2. Alpha diversity
# =============================================================================

hill25 <- iNEXT::estimateD(otu25, q = 0:2, datatype = "abundance", base = "coverage", level = 0.95) |>
  arrange(Assemblage, Order.q, desc(Method)) |>
  group_by(Assemblage, Order.q) |>
  slice(1) |>
  ungroup() |>
  pivot_wider(id_cols = Assemblage, names_from = Order.q, values_from = qD, names_prefix = "ENSq") |>
  rename(SampleID = Assemblage) |>
  mutate(Richness_corr = ENSq0, ENS1 = ENSq1, ENS2 = ENSq2) |>
  select(SampleID, Richness_corr, ENS1, ENS2)

data_alpha25 <- meta25 |>
  left_join(hill25, by = "SampleID") |>
  mutate(Trial = factor(Trial), modality = factor(modality), genotype = factor(genotype))

all_alpha <- data_alpha25 %>%
  mutate(
    modality = factor(modality, levels = c("All_Visitors", "Day_Visitors", "Night_Visitors", "No_Visitors")),
    Trial = factor(Trial), genotype = factor(genotype),
    Richness_sqrt = sqrt(Richness_corr),
    ENS1_log = log(ENS1), ENS2_log = log(ENS2)
  )

# --- Treatment x genotype ANOVA, per metric, with Tukey pairwise contrasts ---
aov_rich <- aov(Richness_sqrt ~ modality * genotype + Trial, data = all_alpha)
summary(aov_rich)
shapiro.test(residuals(aov_rich))
car::leveneTest(residuals(aov_rich) ~ all_alpha$modality)
emm_rich <- emmeans(aov_rich, ~ modality | genotype)
pairs(emm_rich, adjust = "tukey")
cld_rich <- multcomp::cld(emm_rich, adjust = "tukey", Letters = letters)

aov_ens1 <- aov(ENS1_log ~ modality * genotype + Trial, data = all_alpha)
summary(aov_ens1)
shapiro.test(residuals(aov_ens1))
car::leveneTest(residuals(aov_ens1) ~ all_alpha$modality)
emm_ens1 <- emmeans(aov_ens1, ~ modality | genotype)
pairs(emm_ens1, adjust = "tukey")
cld_ens1 <- multcomp::cld(emm_ens1, adjust = "tukey", Letters = letters)

aov_ens2 <- aov(ENS2_log ~ modality * genotype + Trial, data = all_alpha)
summary(aov_ens2)
shapiro.test(residuals(aov_ens2))
car::leveneTest(residuals(aov_ens2) ~ all_alpha$modality)
emm_ens2 <- emmeans(aov_ens2, ~ modality | genotype)
pairs(emm_ens2, adjust = "tukey")
cld_ens2 <- multcomp::cld(emm_ens2, adjust = "tukey", Letters = letters)

# --- Figure: alpha diversity by modality, faceted by cultivar x index -------
plot_df <- data_alpha25 %>%
  filter(!is.na(modality), !is.na(genotype), !is.na(Trial)) %>%
  select(genotype, Trial, modality, Richness_corr, ENS1, ENS2) %>%
  pivot_longer(cols = c(Richness_corr, ENS1, ENS2), names_to = "Indice", values_to = "Diversite") %>%
  mutate(
    Indice = factor(Indice, levels = c("Richness_corr", "ENS1", "ENS2"),
                    labels = c("Richness (q=0)", "ENS1 (q=1)", "ENS2 (q=2)")),
    modality = factor(modality, levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
    genotype = factor(genotype, levels = c("SY_CELESTO", "ES_IDILLIC"), labels = c("CELESTO", "IDILLIC")),
    Trial = factor(Trial, levels = c("24MO11", "25ZM05"), labels = c("2024", "2025"))
  )

y_lims <- plot_df %>% group_by(Indice) %>% summarise(ymax = max(Diversite, na.rm = TRUE) * 1.20, .groups = "drop")
ymax_rich <- y_lims$ymax[y_lims$Indice == "Richness (q=0)"]
ymax_ens1 <- y_lims$ymax[y_lims$Indice == "ENS1 (q=1)"]
ymax_ens2 <- y_lims$ymax[y_lims$Indice == "ENS2 (q=2)"]

letters_df <- bind_rows(
  as.data.frame(cld_rich) %>% mutate(Indice = "Richness (q=0)"),
  as.data.frame(cld_ens1) %>% mutate(Indice = "ENS1 (q=1)"),
  as.data.frame(cld_ens2) %>% mutate(Indice = "ENS2 (q=2)")
) %>%
  mutate(
    .group = trimws(.group),
    modality = factor(modality, levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
    genotype = factor(genotype, levels = c("SY_CELESTO", "ES_IDILLIC"), labels = c("CELESTO", "IDILLIC")),
    Indice = factor(Indice, levels = c("Richness (q=0)", "ENS1 (q=1)", "ENS2 (q=2)"))
  )

label_pos <- plot_df %>%
  group_by(Indice, genotype, modality) %>%
  summarise(y_max = max(Diversite, na.rm = TRUE), .groups = "drop")

letters_plot <- letters_df %>%
  left_join(label_pos, by = c("Indice", "genotype", "modality")) %>%
  left_join(y_lims, by = "Indice") %>%
  group_by(Indice, genotype) %>%
  mutate(n_groups = n_distinct(.group)) %>%
  ungroup() %>%
  mutate(y_label = y_max + 0.10 * ymax)

letters_indiv <- letters_plot %>% filter(n_groups > 1)

ns_summary <- letters_plot %>%
  filter(n_groups == 1) %>%
  group_by(Indice, genotype, ymax) %>%
  summarise(y_label = max(y_label, na.rm = TRUE), .groups = "drop") %>%
  mutate(x_pos = 2.5)

p_alpha <- ggplot(plot_df, aes(x = modality, y = Diversite, fill = modality)) +
  geom_boxplot(outlier.shape = NA, width = 0.68, alpha = 0.85, colour = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Trial), fill = "white", colour = "black", width = 0.16,
              size = 2.5, alpha = 0.9, stroke = 0.5) +
  geom_text(data = letters_indiv, aes(x = modality, y = y_label, label = .group),
            inherit.aes = FALSE, size = 8, fontface = "bold") +
  geom_text(data = ns_summary, aes(x = x_pos, y = y_label, label = "ns"),
            inherit.aes = FALSE, size = 8, fontface = "italic") +
  facet_grid2(genotype ~ Indice, scales = "free_y", independent = "y") +
  facetted_pos_scales(y = list(
    Indice == "Richness (q=0)" ~ scale_y_continuous(limits = c(0, ymax_rich), expand = expansion(mult = c(0, 0))),
    Indice == "ENS1 (q=1)" ~ scale_y_continuous(limits = c(0, ymax_ens1), expand = expansion(mult = c(0, 0))),
    Indice == "ENS2 (q=2)" ~ scale_y_continuous(limits = c(0, ymax_ens2), expand = expansion(mult = c(0, 0)))
  )) +
  scale_fill_manual(values = pal_modality, guide = "none") +
  scale_shape_manual(values = c("2024" = 21, "2025" = 24), name = "Trial") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Fungal alpha diversity") +
  theme_bw(base_size = 35) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 30),
    strip.text.x = ggtext::element_markdown(size = 26, face = "bold"),
    strip.text.y = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    strip.clip = "off",
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(t = 4, r = 4, b = 15, l = 4, unit = "mm"),
    panel.spacing.x = unit(15, "mm"),
    panel.spacing.y = unit(8, "mm")
  )
print(p_alpha)

# --- No-access baseline: cultivar effect alone, no treatment term ------------
alpha_base <- data_alpha25 %>%
  filter(modality == "No_Visitors", !is.na(genotype), !is.na(Trial)) %>%
  mutate(
    genotype = factor(genotype, levels = c("SY_CELESTO", "ES_IDILLIC")),
    Trial = factor(Trial),
    Richness_sqrt = sqrt(Richness_corr),
    ENS1_log = log(ENS1), ENS2_log = log(ENS2)
  )

table(alpha_base$genotype, alpha_base$Trial)

aov_rich_base <- aov(Richness_sqrt ~ genotype + Trial, data = alpha_base)
summary(aov_rich_base)
shapiro.test(residuals(aov_rich_base))
car::leveneTest(residuals(aov_rich_base) ~ alpha_base$genotype)
emmeans(aov_rich_base, ~ genotype)

aov_ens1_base <- aov(ENS1_log ~ genotype + Trial, data = alpha_base)
summary(aov_ens1_base)
shapiro.test(residuals(aov_ens1_base))
car::leveneTest(residuals(aov_ens1_base) ~ alpha_base$genotype)
emmeans(aov_ens1_base, ~ genotype)

aov_ens2_base <- aov(ENS2_log ~ genotype + Trial, data = alpha_base)
summary(aov_ens2_base)
shapiro.test(residuals(aov_ens2_base))
car::leveneTest(residuals(aov_ens2_base) ~ alpha_base$genotype)
emmeans(aov_ens2_base, ~ genotype)

pal_geno <- c("CELESTO" = "#1B9E77", "IDILLIC" = "#7570B3")

plot_base <- alpha_base %>%
  select(genotype, Trial, Richness_corr, ENS1, ENS2) %>%
  pivot_longer(c(Richness_corr, ENS1, ENS2), names_to = "Indice", values_to = "Diversite") %>%
  mutate(
    Indice = factor(Indice, levels = c("Richness_corr", "ENS1", "ENS2"),
                    labels = c("Richness (q=0)", "ENS1 (q=1)", "ENS2 (q=2)")),
    genotype = factor(genotype, levels = c("SY_CELESTO", "ES_IDILLIC"), labels = c("CELESTO", "IDILLIC")),
    Trial = factor(Trial, levels = c("24MO11", "25ZM05"), labels = c("2024", "2025"))
  )

ggplot(plot_base, aes(genotype, Diversite, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, colour = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Trial), colour = "black", width = 0.08, size = 2.5, alpha = 0.9, stroke = 0.5) +
  facet_wrap(~ Indice, scales = "free_y") +
  scale_fill_manual(values = pal_geno, guide = "none") +
  scale_shape_manual(values = c("2024" = 16, "2025" = 17), name = "Trial") +
  labs(x = NULL, y = "Fungal alpha diversity") +
  theme_bw(base_size = 35) +
  theme(axis.title.y = element_text(size = 30), axis.text = element_text(size = 28),
        strip.text = element_text(size = 28, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        axis.line = element_line(linewidth = 0.7), axis.ticks = element_line(linewidth = 0.7),
        legend.text = element_text(size = 24), legend.title = element_text(size = 26))

# =============================================================================
# 3. Beta diversity and community composition
# =============================================================================

physeq_counts <- prune_samples(sample_sums(physeq_final) > 0, physeq_final)

ps_i <- prune_samples(sample_data(physeq_counts)$genotype == "ES_IDILLIC", physeq_counts)
ps_i <- prune_samples(sample_sums(ps_i) > 0, ps_i)
ps_i <- prune_taxa(taxa_sums(ps_i) > 0, ps_i)
ps_i_rel <- transform_sample_counts(ps_i, function(x) x / sum(x))
sdf_i <- as(sample_data(ps_i), "data.frame"); sdf_i$SampleID <- rownames(sdf_i)
sdf_i_rel <- as(sample_data(ps_i_rel), "data.frame")
otu_i <- { m <- as(otu_table(ps_i), "matrix"); if (taxa_are_rows(ps_i)) t(m) else m }

ps_c <- prune_samples(sample_data(physeq_counts)$genotype == "SY_CELESTO", physeq_counts)
ps_c <- prune_samples(sample_sums(ps_c) > 0, ps_c)
ps_c <- prune_taxa(taxa_sums(ps_c) > 0, ps_c)
ps_c_rel <- transform_sample_counts(ps_c, function(x) x / sum(x))
sdf_c <- as(sample_data(ps_c), "data.frame"); sdf_c$SampleID <- rownames(sdf_c)
sdf_c_rel <- as(sample_data(ps_c_rel), "data.frame")
otu_c <- { m <- as(otu_table(ps_c), "matrix"); if (taxa_are_rows(ps_c)) t(m) else m }

d_bray_i <- phyloseq::distance(ps_i_rel, "bray")
d_jacc_i <- vegdist(otu_i, method = "jaccard", binary = TRUE)
d_uuni_i <- phyloseq::UniFrac(ps_i, weighted = FALSE, normalized = TRUE)
d_wuni_i <- phyloseq::UniFrac(ps_i, weighted = TRUE, normalized = TRUE)

d_bray_c <- phyloseq::distance(ps_c_rel, "bray")
d_jacc_c <- vegdist(otu_c, method = "jaccard", binary = TRUE)
d_uuni_c <- phyloseq::UniFrac(ps_c, weighted = FALSE, normalized = TRUE)
d_wuni_c <- phyloseq::UniFrac(ps_c, weighted = TRUE, normalized = TRUE)

for (nm in c("bray_i", "jacc_i", "uuni_i", "wuni_i", "bray_c", "jacc_c", "uuni_c", "wuni_c")) {
  d <- get(paste0("d_", nm))
  grp <- if (grepl("_i$", nm)) sdf_i$modality else sdf_c$modality
  print(permutest(betadisper(d, grp), permutations = 9999))
}

adonis2(d_bray_i ~ modality + Trial, data = sdf_i_rel, permutations = 99999, by = "margin")
adonis2(d_jacc_i ~ modality + Trial, data = sdf_i, permutations = 99999, by = "margin")
adonis2(d_uuni_i ~ modality + Trial, data = sdf_i, permutations = 99999, by = "margin")
adonis2(d_wuni_i ~ modality + Trial, data = sdf_i, permutations = 99999, by = "margin")

adonis2(d_bray_c ~ modality + Trial, data = sdf_c_rel, permutations = 99999, by = "margin")
adonis2(d_jacc_c ~ modality + Trial, data = sdf_c, permutations = 99999, by = "margin")
adonis2(d_uuni_c ~ modality + Trial, data = sdf_c, permutations = 99999, by = "margin")
adonis2(d_wuni_c ~ modality + Trial, data = sdf_c, permutations = 99999, by = "margin")

# --- Pairwise comparisons, BH-adjusted across the 6 pairs within each metric ---
# pairwise.adonis2()'s own p.adjust.m argument does not adjust the printed
# Pr(>F); p.adjust() is applied explicitly here to get correct q-values.
pairwise_beta_pq <- function(dist_obj, data_df, genotype_label, metric_label) {
  res <- pairwise.adonis2(dist_obj ~ modality, data = data_df, permutations = 9999)
  res$parent_call <- NULL
  p_raw <- sapply(res, function(x) x$`Pr(>F)`[1])
  data.frame(Genotype = genotype_label, Metric = metric_label,
             Comparison = names(p_raw), p = as.numeric(p_raw), row.names = NULL)
}

pairwise_beta_results <- bind_rows(
  pairwise_beta_pq(d_bray_c, sdf_c_rel, "CELESTO", "Bray-Curtis"),
  pairwise_beta_pq(d_jacc_c, sdf_c, "CELESTO", "Jaccard"),
  pairwise_beta_pq(d_uuni_c, sdf_c, "CELESTO", "Unweighted UniFrac"),
  pairwise_beta_pq(d_wuni_c, sdf_c, "CELESTO", "Weighted UniFrac"),
  pairwise_beta_pq(d_bray_i, sdf_i_rel, "IDILLIC", "Bray-Curtis"),
  pairwise_beta_pq(d_jacc_i, sdf_i, "IDILLIC", "Jaccard"),
  pairwise_beta_pq(d_uuni_i, sdf_i, "IDILLIC", "Unweighted UniFrac"),
  pairwise_beta_pq(d_wuni_i, sdf_i, "IDILLIC", "Weighted UniFrac")
) %>%
  group_by(Genotype, Metric) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(Genotype, Metric, Comparison)

print(pairwise_beta_results, n = Inf)
write.csv(pairwise_beta_results, "output/pairwise_beta_p_q_ITS.csv", row.names = FALSE)

# --- PCoA figures, per cultivar x metric --------------------------------------
make_pcoa_plot <- function(dist_obj, ps_obj, meta_df, title_label) {
  ord <- ordinate(ps_obj, method = "PCoA", distance = dist_obj)
  eig_pos <- ord$values$Eigenvalues[ord$values$Eigenvalues > 0]
  pct <- round(eig_pos / sum(eig_pos) * 100, 1)
  
  scores <- as.data.frame(ord$vectors[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$SampleID <- rownames(scores)
  
  plot_data <- scores %>%
    left_join(meta_df[, c("SampleID", "modality", "Trial")], by = "SampleID") %>%
    mutate(
      modality = factor(modality, levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                        labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
      Trial = factor(Trial, levels = c("24MO11", "25ZM05"), labels = c("2024", "2025"))
    )
  
  ggplot(plot_data, aes(x = PC1, y = PC2, colour = modality)) +
    stat_ellipse(level = 0.95, linewidth = 0.6, linetype = "dashed", show.legend = FALSE) +
    geom_point(aes(shape = Trial, fill = modality), size = 4, alpha = 0.85, stroke = 0.5) +
    scale_colour_manual(values = pal_modality, guide = "none") +
    scale_fill_manual(values = pal_modality, name = "Treatment",
                      guide = guide_legend(override.aes = list(shape = 21, size = 5, stroke = 0.5))) +
    scale_shape_manual(values = c("2024" = 21, "2025" = 24), name = "Trial") +
    labs(title = title_label, x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
    theme_bw(base_size = 35) +
    theme(plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 30), axis.text = element_text(size = 28),
          panel.border = element_rect(colour = "black", linewidth = 0.8),
          axis.line = element_line(linewidth = 0.7), axis.ticks = element_line(linewidth = 0.7),
          plot.margin = unit(rep(4, 4), "mm"),
          legend.text = element_text(size = 24), legend.title = element_text(size = 26))
}

subt <- function(p, title) p + labs(title = NULL, subtitle = title) +
  theme(plot.subtitle = element_text(size = 24, hjust = 0.5, colour = "grey40", face = "italic"))

p_bray_i4 <- make_pcoa_plot(d_bray_i, ps_i, sdf_i, "IDILLIC") %>% subt("Bray-Curtis")
p_jacc_i4 <- subt(make_pcoa_plot(d_jacc_i, ps_i, sdf_i, "IDILLIC"), "Jaccard")
p_uuni_i4 <- subt(make_pcoa_plot(d_uuni_i, ps_i, sdf_i, "IDILLIC"), "Unweighted UniFrac")
p_wuni_i4 <- subt(make_pcoa_plot(d_wuni_i, ps_i, sdf_i, "IDILLIC"), "Weighted UniFrac")

p_bray_c4 <- make_pcoa_plot(d_bray_c, ps_c, sdf_c, "CELESTO") %>% subt("Bray-Curtis")
p_jacc_c4 <- subt(make_pcoa_plot(d_jacc_c, ps_c, sdf_c, "CELESTO"), "Jaccard")
p_uuni_c4 <- subt(make_pcoa_plot(d_uuni_c, ps_c, sdf_c, "CELESTO"), "Unweighted UniFrac")
p_wuni_c4 <- subt(make_pcoa_plot(d_wuni_c, ps_c, sdf_c, "CELESTO"), "Weighted UniFrac")

p_pcoa_all <- (p_bray_c4 | p_jacc_c4 | p_uuni_c4 | p_wuni_c4) /
  (p_bray_i4 | p_jacc_i4 | p_uuni_i4 | p_wuni_i4) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave("output/pcoa_supp_all.png", p_pcoa_all, width = 28, height = 14, dpi = 300, bg = "white")

# --- Dominant taxa by cultivar x modality --------------------------------------
tax_final <- as.data.frame(as.matrix(tax_table(physeq_final)))
tax_final$OTU <- rownames(tax_final)

otu_mat_final <- as(otu_table(physeq_final), "matrix")
if (!taxa_are_rows(physeq_final)) otu_mat_final <- t(otu_mat_final)

meta_final <- as(sample_data(physeq_final), "data.frame")
meta_final$SampleID <- rownames(meta_final)
meta_final$total_reads <- sample_sums(physeq_final)[meta_final$SampleID]

otu_long <- as.data.frame(t(otu_mat_final)) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "OTU", values_to = "reads") %>%
  left_join(meta_final[, c("SampleID", "genotype", "modality", "total_reads")], by = "SampleID") %>%
  mutate(prop = reads / total_reads)

top_otus <- otu_long %>%
  group_by(genotype, modality, OTU) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  left_join(tax_final[, c("OTU", "Genus", "Species", "Class", "Order")], by = "OTU") %>%
  arrange(genotype, modality, desc(mean_prop)) %>%
  group_by(genotype, modality) %>%
  slice_head(n = 10)
print(top_otus, n = 80)

genus_summary <- otu_long %>%
  left_join(tax_final[, c("OTU", "Genus")], by = "OTU") %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", "Unassigned", Genus)) %>%
  group_by(SampleID, genotype, modality, total_reads, Genus) %>%
  summarise(reads_genus = sum(reads), .groups = "drop") %>%
  mutate(prop_genus = reads_genus / total_reads) %>%
  group_by(genotype, modality, Genus) %>%
  summarise(mean_prop = round(mean(prop_genus) * 100, 2),
            median_prop = round(median(prop_genus) * 100, 2),
            pct_present = round(mean(reads_genus > 0) * 100, 1), .groups = "drop") %>%
  arrange(genotype, modality, desc(mean_prop))
print(genus_summary, n = 10)

# --- Metschnikowia detection frequency and relative abundance -----------------
focal_combined <- otu_long %>%
  left_join(tax_final[, c("OTU", "Genus")], by = "OTU") %>%
  filter(Genus == "Metschnikowia") %>%
  group_by(SampleID, genotype, modality, Genus) %>%
  summarise(prop = sum(prop), .groups = "drop") %>%
  mutate(
    modality = factor(modality, levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
    genotype = factor(genotype, levels = c("SY_CELESTO", "ES_IDILLIC"), labels = c("CELESTO", "IDILLIC"))
  )

prev_summary <- focal_combined %>%
  group_by(genotype, modality) %>%
  summarise(n_present = sum(prop > 0), n_total = n(), pct = n_present / n_total * 100, .groups = "drop")

abund_mets_global <- focal_combined %>%
  group_by(genotype, modality) %>%
  summarise(mean_prop = mean(prop), se_prop = sd(prop) / sqrt(n()), .groups = "drop")

lab_detection <- "*Metschnikowia* detection frequency (%)"
lab_abundance <- "*Metschnikowia* relative abundance (%)"

prev_plot <- prev_summary %>% select(genotype, modality, value = pct) %>% mutate(metric = lab_detection)
abund_plot <- abund_mets_global %>% select(genotype, modality, value = mean_prop) %>%
  mutate(value = value * 100, metric = lab_abundance)
abund_se <- abund_mets_global %>% select(genotype, modality, se_prop) %>% mutate(se_prop = se_prop * 100)

combined_plot <- bind_rows(prev_plot, abund_plot) %>%
  mutate(metric = factor(metric, levels = c(lab_detection, lab_abundance)),
         genotype = factor(genotype, levels = c("CELESTO", "IDILLIC")),
         modality = factor(modality, levels = names(pal_modality))) %>%
  left_join(abund_se, by = c("genotype", "modality"))

ggplot(combined_plot, aes(x = modality, y = value, fill = modality)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.4) +
  geom_errorbar(data = combined_plot %>% filter(metric == lab_abundance),
                aes(ymin = pmax(value - se_prop, 0), ymax = value + se_prop), width = 0.18, linewidth = 0.8) +
  geom_point(data = focal_combined %>% mutate(value = prop * 100, metric = factor(lab_abundance, levels = c(lab_detection, lab_abundance))),
             aes(x = modality, y = value), position = position_jitter(width = 0.12),
             size = 2.2, shape = 21, fill = "white", colour = "black", stroke = 0.5, inherit.aes = FALSE) +
  facet_grid(genotype ~ metric) +
  scale_fill_manual(values = pal_modality, guide = "none") +
  scale_y_continuous(limits = c(-2, 100), expand = expansion(mult = c(0, 0.05)), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 35) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 30),
    strip.text.x = ggtext::element_markdown(size = 26, face = "bold"),
    strip.text.y = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90"), strip.clip = "off",
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7), axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(t = 4, r = 4, b = 15, l = 4, unit = "mm"),
    panel.spacing.x = unit(30, "mm"), panel.spacing.y = unit(8, "mm")
  )

# --- Top 15 fungal families, stacked barplot (supplementary) ------------------
ps_fam <- tax_glom(physeq_final, taxrank = "Family", NArm = FALSE)
ps_fam_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

otu_fam <- as(otu_table(ps_fam_rel), "matrix")
if (taxa_are_rows(ps_fam_rel)) otu_fam <- t(otu_fam)

tax_fam <- as.data.frame(as.matrix(tax_table(ps_fam_rel)))
tax_fam$OTU <- rownames(tax_fam)
tax_fam$Family_label <- ifelse(is.na(tax_fam$Family) | tax_fam$Family == "", "Unassigned", tax_fam$Family)
tax_fam$Family_label <- gsub("_fam_Incertae_sedis", " (unclassified)", tax_fam$Family_label)
tax_fam$Family_label <- gsub("_", " ", tax_fam$Family_label)

meta_fam <- as(sample_data(ps_fam_rel), "data.frame")
meta_fam$SampleID <- rownames(meta_fam)

fam_long <- as.data.frame(otu_fam) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "OTU", values_to = "prop") %>%
  left_join(tax_fam[, c("OTU", "Family_label")], by = "OTU") %>%
  left_join(meta_fam[, c("SampleID", "genotype", "modality")], by = "SampleID")

top15_fam <- fam_long %>%
  group_by(Family_label) %>%
  summarise(total = sum(prop), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(Family_label)

fam_summary_15 <- fam_long %>%
  mutate(Family_group = ifelse(Family_label %in% top15_fam, Family_label, "Other")) %>%
  group_by(SampleID, genotype, modality, Family_group) %>%
  summarise(prop = sum(prop), .groups = "drop") %>%
  group_by(genotype, modality, Family_group) %>%
  summarise(mean_prop = mean(prop) * 100, .groups = "drop") %>%
  mutate(
    modality = factor(modality, levels = c("No_Visitors", "All_Visitors", "Day_Visitors", "Night_Visitors"),
                      labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
    genotype = factor(genotype, levels = c("SY_CELESTO", "ES_IDILLIC"), labels = c("CELESTO", "IDILLIC")),
    Family_group = factor(Family_group, levels = rev(c(top15_fam, "Other")))
  )

pal_15 <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    "#44AA99", "#882255", "#117733", "#DDCC77", "#332288", "#AA4499", "#88CCEE", "#661100"),
  top15_fam
)
pal_15[["Unassigned"]] <- "#7B2D8B"
pal_15[["Other"]] <- "grey80"

ggplot(fam_summary_15, aes(x = modality, y = mean_prop, fill = Family_group)) +
  geom_bar(stat = "identity", colour = NA, width = 0.8) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_manual(values = pal_15, name = "Family",
                    guide = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(colour = NA))) +
  labs(x = NULL, y = "Mean relative abundance (%)") +
  theme_bw(base_size = 22) +
  theme(
    axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 18, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 18), strip.text = element_text(size = 22, face = "bold"),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7), axis.ticks = element_line(linewidth = 0.7),
    plot.margin = unit(rep(4, 4), "mm"),
    legend.text = element_text(size = 18), legend.title = element_text(size = 20), legend.key.size = unit(6, "mm")
  )

# =============================================================================
# 4. Camera-trap visit rates as predictors of microbiome structure
# =============================================================================

# --- Import data -------------------------------------------------------
# plant_means.csv: per-plant mean visit rates (see pollinator_visits.R)
camera <- read.csv("data/plant_means.csv", stringsAsFactors = FALSE)

sample_data(physeq_final) <- sample_data(
  data.frame(sample_data(physeq_final)) %>%
    rownames_to_column("SampleID") %>%
    left_join(camera, by = "PLANT_ID") %>%
    column_to_rownames("SampleID")
)

physeq_cam <- subset_samples(physeq_final, modality %in% c("All_Visitors", "Day_Visitors", "Night_Visitors") & !is.na(mean_bee))
physeq_cam <- prune_taxa(taxa_sums(physeq_cam) > 0, physeq_cam)
sample_data(physeq_cam)$mean_day <- sample_data(physeq_cam)$mean_bee + sample_data(physeq_cam)$mean_bumblebee

cat("physeq_cam:", nsamples(physeq_cam), "samples,", ntaxa(physeq_cam), "OTUs\n")

alpha_cam <- data_alpha25 %>%
  inner_join(data.frame(sample_data(physeq_cam)) %>% select(SampleID, mean_bee, mean_bumblebee, mean_moth), by = "SampleID") %>%
  mutate(mean_day = mean_bee + mean_bumblebee)

# --- Alpha diversity ~ visit rate, day and night windows, per cultivar -------
prep_alpha_subset <- function(geno, mods) {
  alpha_cam %>%
    filter(genotype == geno, modality %in% mods) %>%
    mutate(Richness_sqrt = sqrt(Richness_corr), ENS1_log = log(ENS1), ENS2_log = log(ENS2))
}

day_all_cel <- prep_alpha_subset("SY_CELESTO", c("All_Visitors", "Day_Visitors"))
day_all_idil <- prep_alpha_subset("ES_IDILLIC", c("All_Visitors", "Day_Visitors"))
night_all_cel <- prep_alpha_subset("SY_CELESTO", c("All_Visitors", "Night_Visitors"))
night_all_idil <- prep_alpha_subset("ES_IDILLIC", c("All_Visitors", "Night_Visitors"))

fit_cel_rich <- lm(Richness_sqrt ~ modality + Trial + mean_day, data = day_all_cel)
fit_cel_ens1 <- lm(ENS1_log ~ modality + Trial + mean_day, data = day_all_cel)
fit_cel_ens2 <- lm(ENS2_log ~ modality + Trial + mean_day, data = day_all_cel)
summary(fit_cel_rich); shapiro.test(residuals(fit_cel_rich)); ncvTest(fit_cel_rich)
summary(fit_cel_ens1); shapiro.test(residuals(fit_cel_ens1)); ncvTest(fit_cel_ens1)
summary(fit_cel_ens2); shapiro.test(residuals(fit_cel_ens2)); ncvTest(fit_cel_ens2)

fit_idil_rich <- lm(Richness_sqrt ~ modality + mean_day + Trial, data = day_all_idil)
fit_idil_ens1 <- lm(ENS1_log ~ modality + mean_day + Trial, data = day_all_idil)
fit_idil_ens2 <- lm(ENS2_log ~ modality + mean_day + Trial, data = day_all_idil)
summary(fit_idil_rich); shapiro.test(residuals(fit_idil_rich)); ncvTest(fit_idil_rich)
summary(fit_idil_ens1); shapiro.test(residuals(fit_idil_ens1)); ncvTest(fit_idil_ens1)
summary(fit_idil_ens2); shapiro.test(residuals(fit_idil_ens2)); ncvTest(fit_idil_ens2)

fit_cel_night_rich <- lm(Richness_sqrt ~ modality + mean_moth + Trial, data = night_all_cel)
fit_cel_night_ens1 <- lm(ENS1_log ~ modality + mean_moth + Trial, data = night_all_cel)
fit_cel_night_ens2 <- lm(ENS2_log ~ modality + mean_moth + Trial, data = night_all_cel)
summary(fit_cel_night_rich); shapiro.test(residuals(fit_cel_night_rich)); ncvTest(fit_cel_night_rich)
summary(fit_cel_night_ens1); shapiro.test(residuals(fit_cel_night_ens1)); ncvTest(fit_cel_night_ens1)
summary(fit_cel_night_ens2); shapiro.test(residuals(fit_cel_night_ens2)); ncvTest(fit_cel_night_ens2)

fit_idil_night_rich <- lm(Richness_sqrt ~ modality + mean_moth + Trial, data = night_all_idil)
fit_idil_night_ens1 <- lm(ENS1_log ~ modality + mean_moth + Trial, data = night_all_idil)
fit_idil_night_ens2 <- lm(ENS2_log ~ modality + mean_moth + Trial, data = night_all_idil)
summary(fit_idil_night_rich); shapiro.test(residuals(fit_idil_night_rich)); ncvTest(fit_idil_night_rich)
summary(fit_idil_night_ens1); shapiro.test(residuals(fit_idil_night_ens1)); ncvTest(fit_idil_night_ens1)
summary(fit_idil_night_ens2); shapiro.test(residuals(fit_idil_night_ens2)); ncvTest(fit_idil_night_ens2)

# --- Predicted average marginal effect of visit rate, with t-based CI -------
compute_avg_pred <- function(data, response_var, focal_var, covariates = c("modality", "Trial"), n_grid = 100) {
  data <- droplevels(data)
  fm <- reformulate(c(covariates, focal_var), response = response_var)
  fit <- lm(fm, data = data)
  x_grid <- seq(min(data[[focal_var]], na.rm = TRUE), max(data[[focal_var]], na.rm = TRUE), length.out = n_grid)
  t_crit <- qt(0.975, df.residual(fit))
  
  purrr::map_dfr(x_grid, function(x_val) {
    newdata <- data
    newdata[[focal_var]] <- x_val
    X <- model.matrix(delete.response(terms(fit)), newdata)
    pred_i <- as.numeric(X %*% coef(fit))
    cov_pred <- X %*% vcov(fit) %*% t(X)
    data.frame(x = x_val, y_hat = mean(pred_i), se = sqrt(sum(cov_pred)) / nrow(data))
  }) %>%
    mutate(ymin = y_hat - t_crit * se, ymax = y_hat + t_crit * se)
}

build_pred_df <- function(subset_data, response_vars, focal_var, genotype_label) {
  purrr::map_dfr(response_vars, function(resp) {
    df <- compute_avg_pred(subset_data, resp, focal_var)
    df$ENS_metric <- resp
    df
  }) %>% mutate(genotype = genotype_label)
}

build_raw_df <- function(subset_data, response_vars, focal_var, genotype_label) {
  subset_data <- droplevels(subset_data)
  subset_data %>%
    select(all_of(c(focal_var, response_vars))) %>%
    rename(x = all_of(focal_var)) %>%
    pivot_longer(cols = all_of(response_vars), names_to = "ENS_metric", values_to = "y_raw") %>%
    mutate(genotype = genotype_label)
}

metric_levels <- c("Richness_corr", "ENS1", "ENS2")
metric_labels <- c("Richness_corr" = "Richness", "ENS1" = "ENS1", "ENS2" = "ENS2")

diurnal_pred <- bind_rows(
  build_pred_df(day_all_cel, metric_levels, "mean_day", "CELESTO"),
  build_pred_df(day_all_idil, metric_levels, "mean_day", "IDILLIC")
) %>% mutate(genotype = factor(genotype, levels = c("CELESTO", "IDILLIC")),
             ENS_metric = factor(ENS_metric, levels = metric_levels, labels = metric_labels))

diurnal_raw <- bind_rows(
  build_raw_df(day_all_cel, metric_levels, "mean_day", "CELESTO"),
  build_raw_df(day_all_idil, metric_levels, "mean_day", "IDILLIC")
) %>% mutate(genotype = factor(genotype, levels = c("CELESTO", "IDILLIC")),
             ENS_metric = factor(ENS_metric, levels = metric_levels, labels = metric_labels))

moth_pred <- bind_rows(
  build_pred_df(night_all_cel, metric_levels, "mean_moth", "CELESTO"),
  build_pred_df(night_all_idil, metric_levels, "mean_moth", "IDILLIC")
) %>% mutate(genotype = factor(genotype, levels = c("CELESTO", "IDILLIC")),
             ENS_metric = factor(ENS_metric, levels = metric_levels, labels = metric_labels))

moth_raw <- bind_rows(
  build_raw_df(night_all_cel, metric_levels, "mean_moth", "CELESTO"),
  build_raw_df(night_all_idil, metric_levels, "mean_moth", "IDILLIC")
) %>% mutate(genotype = factor(genotype, levels = c("CELESTO", "IDILLIC")),
             ENS_metric = factor(ENS_metric, levels = metric_levels, labels = metric_labels))

make_pred_plot <- function(pred_data, raw_data, x_label) {
  ggplot(pred_data, aes(x = x, y = y_hat)) +
    geom_line(color = "black", linewidth = 1.5) +
    geom_point(data = raw_data, aes(x = x, y = y_raw), color = "black", alpha = 1, size = 2) +
    facet_grid2(genotype ~ ENS_metric, scales = "free", independent = "y") +
    labs(x = x_label, y = "Fungal alpha diversity") +
    theme_natcomm
}

make_pred_plot(diurnal_pred, diurnal_raw, "Bee and bumblebee visitation (detections/photo)")
make_pred_plot(moth_pred, moth_raw, "Moth visitation (detections/photo)")

# --- Beta diversity ~ visit rate (continuous), Trial as covariate -----------
build_sub_ps <- function(ps, geno, mods) {
  sdf <- data.frame(sample_data(ps))
  keep <- sdf$genotype == geno & sdf$modality %in% mods
  sub_ps <- prune_samples(rownames(sdf)[keep], ps)
  prune_taxa(taxa_sums(sub_ps) > 0, sub_ps)
}

ps_day_all_cel <- build_sub_ps(physeq_cam, "SY_CELESTO", c("All_Visitors", "Day_Visitors"))
ps_day_all_idil <- build_sub_ps(physeq_cam, "ES_IDILLIC", c("All_Visitors", "Day_Visitors"))
ps_night_all_cel <- build_sub_ps(physeq_cam, "SY_CELESTO", c("All_Visitors", "Night_Visitors"))
ps_night_all_idil <- build_sub_ps(physeq_cam, "ES_IDILLIC", c("All_Visitors", "Night_Visitors"))

sample_data(ps_night_all_cel)$mean_moth <- as.numeric(sample_data(ps_night_all_cel)$mean_moth)
sample_data(ps_night_all_idil)$mean_moth <- as.numeric(sample_data(ps_night_all_idil)$mean_moth)

test_beta_covariate <- function(ps_obj, covariate, dist_method) {
  ps_rel <- transform_sample_counts(ps_obj, function(x) x / sum(x))
  d <- phyloseq::distance(ps_rel, dist_method)
  sdf <- as(sample_data(ps_obj), "data.frame")
  adonis2(as.formula(paste0("d ~ ", covariate, " + Trial")), data = sdf, permutations = 99999, by = "margin")
}

test_beta_covariate_jaccard <- function(ps_obj, covariate) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) otu <- t(otu)
  d <- vegdist(t(otu), method = "jaccard", binary = TRUE)
  sdf <- as(sample_data(ps_obj), "data.frame")
  adonis2(as.formula(paste0("d ~ ", covariate, " + Trial")), data = sdf, permutations = 99999, by = "margin")
}

test_beta_covariate_unifrac <- function(ps_obj, covariate, weighted) {
  d <- phyloseq::UniFrac(ps_obj, weighted = weighted, normalized = TRUE)
  sdf <- as(sample_data(ps_obj), "data.frame")
  adonis2(as.formula(paste0("d ~ ", covariate, " + Trial")), data = sdf, permutations = 99999, by = "margin")
}

# extract R2/F/P for the covariate term from an adonis2 table without relying
# on broom (broom does not reliably parse adonis2's column names)
extract_adonis_term <- function(res, term) {
  tab <- as.data.frame(res)
  row <- tab[term, ]
  data.frame(R2 = row$R2, F = row$F, `Pr(>F)` = row$`Pr(>F)`, check.names = FALSE)
}

beta_covariate_summary <- bind_rows(
  extract_adonis_term(test_beta_covariate(ps_day_all_cel, "mean_day", "bray"), "mean_day") %>% mutate(Genotype = "CELESTO", Window = "Diurnal", Metric = "Bray-Curtis"),
  extract_adonis_term(test_beta_covariate_jaccard(ps_day_all_cel, "mean_day"), "mean_day") %>% mutate(Genotype = "CELESTO", Window = "Diurnal", Metric = "Jaccard"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_day_all_cel, "mean_day", weighted = FALSE), "mean_day") %>% mutate(Genotype = "CELESTO", Window = "Diurnal", Metric = "Unweighted UniFrac"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_day_all_cel, "mean_day", weighted = TRUE), "mean_day") %>% mutate(Genotype = "CELESTO", Window = "Diurnal", Metric = "Weighted UniFrac"),
  extract_adonis_term(test_beta_covariate(ps_day_all_idil, "mean_day", "bray"), "mean_day") %>% mutate(Genotype = "IDILLIC", Window = "Diurnal", Metric = "Bray-Curtis"),
  extract_adonis_term(test_beta_covariate_jaccard(ps_day_all_idil, "mean_day"), "mean_day") %>% mutate(Genotype = "IDILLIC", Window = "Diurnal", Metric = "Jaccard"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_day_all_idil, "mean_day", weighted = FALSE), "mean_day") %>% mutate(Genotype = "IDILLIC", Window = "Diurnal", Metric = "Unweighted UniFrac"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_day_all_idil, "mean_day", weighted = TRUE), "mean_day") %>% mutate(Genotype = "IDILLIC", Window = "Diurnal", Metric = "Weighted UniFrac"),
  extract_adonis_term(test_beta_covariate(ps_night_all_cel, "mean_moth", "bray"), "mean_moth") %>% mutate(Genotype = "CELESTO", Window = "Nocturnal", Metric = "Bray-Curtis"),
  extract_adonis_term(test_beta_covariate_jaccard(ps_night_all_cel, "mean_moth"), "mean_moth") %>% mutate(Genotype = "CELESTO", Window = "Nocturnal", Metric = "Jaccard"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_night_all_cel, "mean_moth", weighted = FALSE), "mean_moth") %>% mutate(Genotype = "CELESTO", Window = "Nocturnal", Metric = "Unweighted UniFrac"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_night_all_cel, "mean_moth", weighted = TRUE), "mean_moth") %>% mutate(Genotype = "CELESTO", Window = "Nocturnal", Metric = "Weighted UniFrac"),
  extract_adonis_term(test_beta_covariate(ps_night_all_idil, "mean_moth", "bray"), "mean_moth") %>% mutate(Genotype = "IDILLIC", Window = "Nocturnal", Metric = "Bray-Curtis"),
  extract_adonis_term(test_beta_covariate_jaccard(ps_night_all_idil, "mean_moth"), "mean_moth") %>% mutate(Genotype = "IDILLIC", Window = "Nocturnal", Metric = "Jaccard"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_night_all_idil, "mean_moth", weighted = FALSE), "mean_moth") %>% mutate(Genotype = "IDILLIC", Window = "Nocturnal", Metric = "Unweighted UniFrac"),
  extract_adonis_term(test_beta_covariate_unifrac(ps_night_all_idil, "mean_moth", weighted = TRUE), "mean_moth") %>% mutate(Genotype = "IDILLIC", Window = "Nocturnal", Metric = "Weighted UniFrac")
) %>%
  select(Genotype, Window, Metric, R2, F, `Pr(>F)`)

print(beta_covariate_summary, n = Inf)