# Floral VOC (GC-MS) analysis: import and log-transform, outlier removal,
# per-compound treatment effect within each cultivar (BH-adjusted), Sidak
# pairwise letters, heatmap of significant compounds, and beta-diversity
# (Bray-Curtis PERMANOVA) on the full volatile blend.

library(readr)
library(dplyr)
library(tibble)
library(purrr)
library(emmeans)
library(ggplot2)
library(tidyr)
library(multcomp)
library(multcompView)
library(pheatmap)
library(vegan)
library(patchwork)

set.seed(64621)

# --- Import data -------------------------------------------------------
# Data_VOC.csv: compounds as rows, samples as columns, plus compound
#   metadata columns (Metabolite_ID, Retention_time, Name, FORMULA, Library,
#   CAS, Forward, Reverse, Avg, RI_expermental, RI_literature, Diff_RI,
#   mass, max_sample_name, max_value)
# metadata_VOC.csv: one row per sample
#   Nom_Ech, Identifier, nr, genotype, rep, Treatment_modality
raw <- read_delim(
  "data/Data_VOC.csv",
  delim = ";",
  locale = locale(decimal_mark = "."),
  na = "",                        # 4 compounds are literally named "NA"
  col_types = cols(.default = "c")
)

metadata <- read_delim("data/metadata_VOC.csv", delim = ";", locale = locale(decimal_mark = "."))

# --- Compound metadata ---------------------------------------------------
compound_info_cols <- c("Metabolite_ID", "Retention_time", "Name", "FORMULA",
                        "Library", "CAS", "Forward", "Reverse", "Avg",
                        "RI_expermental", "RI_literature", "Diff_RI", "mass",
                        "max_sample_name", "max_value")

data_metabolite <- raw %>%
  dplyr::select(all_of(compound_info_cols)) %>%
  mutate(across(c(Retention_time, Forward, Reverse, Avg, RI_expermental,
                  RI_literature, Diff_RI, mass, max_value), as.numeric))

metabolite_names <- setNames(data_metabolite$Name, paste0("VOC_", data_metabolite$Metabolite_ID))
retention_time <- setNames(data_metabolite$Retention_time, paste0("VOC_", data_metabolite$Metabolite_ID))

# Annotate a VOC_XXXX code with its library name, or "Unknown (RT ... min)"
# if unnamed. The literal "NA" name (4 compounds) is not treated as missing.
full_annotation <- function(compound_codes) {
  name <- trimws(metabolite_names[compound_codes])
  rt <- retention_time[compound_codes]
  ifelse(is.na(name), paste0("Unknown (RT ", round(rt, 1), " min)"), name)
}

# --- Abundance matrix: samples as rows, compounds as columns -------------
sample_cols <- setdiff(names(raw), compound_info_cols)

data <- raw %>%
  dplyr::select(Metabolite_ID, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  column_to_rownames("Metabolite_ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Nom_Ech")

colnames(data)[-1] <- paste0("VOC_", colnames(data)[-1])

data_full <- metadata %>%
  dplyr::select(Nom_Ech, Identifier, nr, genotype, rep, Treatment_modality) %>%
  left_join(data, by = "Nom_Ech")

stopifnot(n_distinct(metadata$Nom_Ech) == nrow(metadata))
stopifnot(n_distinct(data$Nom_Ech) == nrow(data))
stopifnot(nrow(data_full) == nrow(data))
stopifnot(sum(is.na(data_full$genotype)) == 0)

# replace "-" with "_" in modality codes for compatibility with car::Anova/emmeans
data_full <- data_full %>%
  mutate(Treatment_modality = gsub("-", "_", Treatment_modality))

voc_cols <- grep("^VOC_", names(data_full), value = TRUE)

data_full_log <- data_full %>%
  mutate(across(all_of(voc_cols), ~ log10(. + 1))) %>%
  filter(genotype != "QC")

n_distinct(data_full_log$Nom_Ech)
length(voc_cols)
table(data_full_log$genotype)
table(data_full_log$genotype, data_full_log$Treatment_modality)

# --- Outlier screening (PCA on all compounds) -----------------------------
run_pca_dist_modalite <- function(df) {
  mat <- df %>% dplyr::select(all_of(voc_cols)) %>% as.matrix()
  pca <- prcomp(mat, scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$Nom_Ech <- df$Nom_Ech
  scores$genotype <- unique(df$genotype)
  scores$Treatment_modality <- df$Treatment_modality
  scores$dist_genotype <- mahalanobis(scores[, c("PC1", "PC2")],
                                      center = colMeans(scores[, c("PC1", "PC2")]),
                                      cov = cov(scores[, c("PC1", "PC2")]))
  scores %>%
    group_by(Treatment_modality) %>%
    mutate(dist_modalite = mahalanobis(cbind(PC1, PC2),
                                       center = colMeans(cbind(PC1, PC2)),
                                       cov = cov(cbind(PC1, PC2)))) %>%
    ungroup()
}

scores_tous <- bind_rows(
  run_pca_dist_modalite(filter(data_full_log, genotype == "CEL")),
  run_pca_dist_modalite(filter(data_full_log, genotype == "IDI"))
)

scores_tous %>%
  group_by(genotype) %>%
  arrange(desc(dist_modalite)) %>%
  dplyr::select(genotype, Nom_Ech, Treatment_modality, dist_genotype, dist_modalite) %>%
  slice_head(n = 5) %>%
  ungroup()

mat_global <- data_full_log %>% dplyr::select(all_of(voc_cols)) %>% as.matrix()
pca_global <- prcomp(mat_global, scale. = TRUE)

scores_global <- as.data.frame(pca_global$x[, 1:2])
colnames(scores_global) <- c("PC1", "PC2")
scores_global$Nom_Ech <- data_full_log$Nom_Ech
scores_global$genotype <- data_full_log$genotype
scores_global$Treatment_modality <- data_full_log$Treatment_modality
scores_global$est_retire <- scores_global$Nom_Ech == "20260519_helex_inra_florets_series_57"

var_exp <- round(100 * summary(pca_global)$importance[2, 1:2], 1)

ggplot(scores_global, aes(PC1, PC2, color = genotype, shape = Treatment_modality)) +
  geom_point(data = filter(scores_global, !est_retire), size = 3, alpha = 0.85) +
  geom_point(data = filter(scores_global, est_retire), size = 5, color = "black") +
  geom_text(data = filter(scores_global, est_retire), aes(x = PC1, y = PC2, label = "series_57 (outlier)"),
            vjust = -1.3, color = "black", fontface = "bold", size = 4, inherit.aes = FALSE) +
  scale_color_manual(values = c(CEL = "#D55E00", IDI = "#0072B2")) +
  labs(title = "Global PCA -- PC1_PC2",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       color = "Genotype", shape = "Modality") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# --- Models per compound, per cultivar ------------------------------------
data_full_log_CEL <- filter(data_full_log, genotype == "CEL",
                            Nom_Ech != "20260519_helex_inra_florets_series_57")
data_full_log_IDI <- filter(data_full_log, genotype == "IDI")

data_modele <- bind_rows(data_full_log_CEL, data_full_log_IDI)

data_long <- data_modele %>%
  dplyr::select(Nom_Ech, genotype, Treatment_modality, all_of(voc_cols)) %>%
  pivot_longer(all_of(voc_cols), names_to = "Metabolite_col", values_to = "valeur")

ggplot(data_long, aes(valeur, fill = genotype)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c(CEL = "#D55E00", IDI = "#0072B2")) +
  labs(title = "Pooled distribution of VOC values (log10(abundance + 1))",
       x = "log10(abundance + 1)", y = "Number of values (all compounds pooled)",
       fill = "Genotype") +
  theme_minimal(base_size = 13)

models_CEL <- map(voc_cols, ~ lm(reformulate("Treatment_modality", response = .x), data = data_full_log_CEL))
names(models_CEL) <- voc_cols
models_IDI <- map(voc_cols, ~ lm(reformulate("Treatment_modality", response = .x), data = data_full_log_IDI))
names(models_IDI) <- voc_cols

shapiro_CEL <- map_dbl(models_CEL, ~ shapiro.test(residuals(.x))$p.value)
shapiro_IDI <- map_dbl(models_IDI, ~ shapiro.test(residuals(.x))$p.value)

# --- Omnibus treatment effect, per compound x cultivar --------------------
anova_omnibus_CEL <- map2_dfr(models_CEL, names(models_CEL), function(mod, nom) {
  broom::tidy(anova(mod)) %>% filter(term == "Treatment_modality") %>% mutate(Metabolite_col = nom)
})
anova_omnibus_IDI <- map2_dfr(models_IDI, names(models_IDI), function(mod, nom) {
  broom::tidy(anova(mod)) %>% filter(term == "Treatment_modality") %>% mutate(Metabolite_col = nom)
})

anova_omnibus <- bind_rows(
  anova_omnibus_CEL %>% mutate(genotype = "CEL", shapiro_p = shapiro_CEL[Metabolite_col]),
  anova_omnibus_IDI %>% mutate(genotype = "IDI", shapiro_p = shapiro_IDI[Metabolite_col])
) %>%
  group_by(genotype) %>%
  mutate(q_value = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_CEL <- anova_omnibus %>% filter(genotype == "CEL", q_value < 0.05) %>% pull(Metabolite_col)
sig_IDI <- anova_omnibus %>% filter(genotype == "IDI", q_value < 0.05) %>% pull(Metabolite_col)

length(sig_CEL)
length(sig_IDI)

anova_omnibus %>%
  dplyr::select(Metabolite_col, genotype, p.value, q_value) %>%
  filter(Metabolite_col %in% union(sig_CEL, sig_IDI)) %>%
  arrange(genotype, Metabolite_col) %>%
  print(n = Inf)

# --- Sidak pairwise comparisons and compact letter groupings --------------
get_cld <- function(mod, nom, geno) {
  emm <- emmeans(mod, ~ Treatment_modality)
  cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
  as.data.frame(cld_out) %>% mutate(Metabolite_col = nom, genotype = geno)
}

cld_CEL <- map_dfr(sig_CEL, ~ get_cld(models_CEL[[.x]], .x, "CEL"))
cld_IDI <- map_dfr(sig_IDI, ~ get_cld(models_IDI[[.x]], .x, "IDI"))
cld_all <- bind_rows(cld_CEL, cld_IDI)

cld_all %>%
  dplyr::select(Metabolite_col, genotype, Treatment_modality, emmean, .group) %>%
  arrange(genotype, Metabolite_col, Treatment_modality)

# --- Compounds shared between cultivars -----------------------------------
shared_sig <- intersect(sig_CEL, sig_IDI)
shared_sig
length(shared_sig)
round(100 * length(shared_sig) / length(union(sig_CEL, sig_IDI)), 1)
round(100 * length(shared_sig) / length(sig_CEL), 1)
round(100 * length(shared_sig) / length(sig_IDI), 1)
round(100 * length(shared_sig) / length(voc_cols), 1)

# --- Supplementary table ---------------------------------------------------
supp_table_voc_wide <- cld_all %>%
  left_join(
    anova_omnibus %>% dplyr::select(Metabolite_col, genotype, model_P = p.value, model_q = q_value),
    by = c("Metabolite_col", "genotype")
  ) %>%
  mutate(
    Cultivar = dplyr::recode(genotype, CEL = "CELESTO", IDI = "IDILLIC"),
    Modality = dplyr::recode(Treatment_modality,
                             N_O = "No access", O = "Continuous access",
                             O_D = "Daytime access", O_N = "Nighttime access"),
    Group = trimws(.group)
  ) %>%
  dplyr::select(Compound = Metabolite_col, Cultivar, Modality, Group, model_P, model_q) %>%
  pivot_wider(names_from = Modality, values_from = Group) %>%
  dplyr::select(Compound, Cultivar,
                `No access`, `Continuous access`, `Daytime access`, `Nighttime access`,
                `Model P` = model_P, `Model q` = model_q) %>%
  mutate(
    `Model P` = signif(`Model P`, 3),
    `Model q` = signif(`Model q`, 3)
  ) %>%
  arrange(Cultivar, Compound)

supp_table_voc_wide <- supp_table_voc_wide %>%
  mutate(Annotation = full_annotation(Compound), .after = Compound)

stopifnot(nrow(supp_table_voc_wide) == length(sig_CEL) + length(sig_IDI))

write_excel_csv(supp_table_voc_wide, "output/supp_table_voc_contrasts.csv")

# --- Heatmap of significant compounds --------------------------------------
mat_heatmap <- data_modele %>%
  dplyr::select(all_of(voc_cols)) %>%
  as.matrix() %>%
  t()
colnames(mat_heatmap) <- data_modele$Nom_Ech

annotation_col <- data.frame(
  Genotype = data_modele$genotype,
  Modality = factor(data_modele$Treatment_modality,
                    levels = c("N_O", "O", "O_D", "O_N"),
                    labels = c("No access", "Continuous access", "Daytime access", "Nighttime access")),
  row.names = data_modele$Nom_Ech
)

annotation_colors <- list(
  Genotype = c(CEL = "#D55E00", IDI = "#009E73"),
  Modality = c("No access" = "#000000", "Continuous access" = "#CC79A7",
               "Daytime access" = "#F5C400", "Nighttime access" = "#0072B2")
)

ordre <- data_modele %>%
  mutate(Modality = factor(Treatment_modality, levels = c("N_O", "O", "O_D", "O_N"))) %>%
  arrange(genotype, Modality, rep) %>%
  pull(Nom_Ech)

mat_heatmap <- mat_heatmap[, ordre]
annotation_col <- annotation_col[ordre, , drop = FALSE]

mat_cel <- mat_heatmap[sig_CEL, annotation_col$Genotype == "CEL"]
mat_idi <- mat_heatmap[sig_IDI, annotation_col$Genotype == "IDI"]

annotation_col_cel <- annotation_col[colnames(mat_cel), "Modality", drop = FALSE]
annotation_col_idi <- annotation_col[colnames(mat_idi), "Modality", drop = FALSE]
colnames(annotation_col_cel) <- "Treatment"
colnames(annotation_col_idi) <- "Treatment"

annotation_colors_mod <- list(Treatment = annotation_colors$Modality)

height_cel <- max(4, 1.2 + length(sig_CEL) * 0.22)
height_idi <- max(4, 1.2 + length(sig_IDI) * 0.22)

label_width_in <- function(labels, chars_per_inch = 11) max(nchar(labels)) / chars_per_inch
mat_width_in <- function(n_col, col_width_in = 0.22) n_col * col_width_in

width_cel <- mat_width_in(ncol(mat_cel)) + label_width_in(rownames(mat_cel)) + 2.8
width_idi <- mat_width_in(ncol(mat_idi)) + label_width_in(rownames(mat_idi)) + 2.8

pheatmap(mat_cel, scale = "row", annotation_col = annotation_col_cel,
         annotation_colors = annotation_colors_mod, show_rownames = TRUE, show_colnames = FALSE,
         cluster_cols = FALSE, clustering_distance_rows = "correlation", main = "CELESTO",
         fontsize = 12, fontsize_row = 9, border_color = "black",
         filename = "output/heatmap_celesto.png", width = width_cel, height = height_cel, res = 300)

pheatmap(mat_idi, scale = "row", annotation_col = annotation_col_idi,
         annotation_colors = annotation_colors_mod, show_rownames = TRUE, show_colnames = FALSE,
         cluster_cols = FALSE, clustering_distance_rows = "correlation", main = "IDILLIC",
         fontsize = 12, fontsize_row = 9, border_color = "black",
         filename = "output/heatmap_idillic.png", width = width_idi, height = height_idi, res = 300)

# --- Volatile blend beta-diversity: Bray-Curtis PERMANOVA, per cultivar ---
data_bc <- data_full %>%
  filter(genotype != "QC", Nom_Ech != "20260519_helex_inra_florets_series_57") %>%
  mutate(Treatment_modality = gsub("-", "_", Treatment_modality))

stopifnot(all(!is.na(data_bc[voc_cols])))
stopifnot(all(data_bc[voc_cols] >= 0))

mat_bc <- data_bc %>% dplyr::select(all_of(voc_cols)) %>% as.matrix()
rownames(mat_bc) <- data_bc$Nom_Ech
mat_rel <- decostand(mat_bc, method = "total")

d_all <- vegdist(mat_rel, method = "bray")  # all samples, used for the genotype-only test below

idx_CEL <- data_bc$genotype == "CEL"
meta_CEL <- data_bc[idx_CEL, ]
d_CEL <- vegdist(mat_rel[idx_CEL, ], method = "bray")

permutest(betadisper(d_CEL, meta_CEL$Treatment_modality), permutations = 9999)
adonis2(d_CEL ~ Treatment_modality, data = meta_CEL, permutations = 9999)

idx_IDI <- data_bc$genotype == "IDI"
meta_IDI <- data_bc[idx_IDI, ]
d_IDI <- vegdist(mat_rel[idx_IDI, ], method = "bray")

permutest(betadisper(d_IDI, meta_IDI$Treatment_modality), permutations = 9999)
adonis2(d_IDI ~ Treatment_modality, data = meta_IDI, permutations = 9999)

# --- PCoA, Bray-Curtis ------------------------------------------------------
modality_labels <- c(N_O = "No access", O = "Continuous access",
                     O_D = "Daytime access", O_N = "Nighttime access")
modality_colors <- c("No access" = "#000000", "Continuous access" = "#CC79A7",
                     "Daytime access" = "#F5C400", "Nighttime access" = "#0072B2")

build_pcoa_df <- function(d, meta) {
  pcoa <- cmdscale(d, k = 2, eig = TRUE)
  var_exp <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)
  df <- as.data.frame(pcoa$points)
  colnames(df) <- c("PCo1", "PCo2")
  df$Treatment_modality <- factor(modality_labels[meta$Treatment_modality], levels = modality_labels)
  list(df = df, var_exp = var_exp)
}

pcoa_CEL <- build_pcoa_df(d_CEL, meta_CEL)
pcoa_IDI <- build_pcoa_df(d_IDI, meta_IDI)

p_CEL <- ggplot(pcoa_CEL$df, aes(PCo1, PCo2, color = Treatment_modality)) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = modality_colors) +
  labs(title = "CELESTO", x = paste0("PCo1 (", pcoa_CEL$var_exp[1], "%)"),
       y = paste0("PCo2 (", pcoa_CEL$var_exp[2], "%)"), color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_IDI <- ggplot(pcoa_IDI$df, aes(PCo1, PCo2, color = Treatment_modality)) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = modality_colors) +
  labs(title = "IDILLIC", x = paste0("PCo1 (", pcoa_IDI$var_exp[1], "%)"),
       y = paste0("PCo2 (", pcoa_IDI$var_exp[2], "%)"), color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_CEL + p_IDI + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# --- PERMANOVA, genotype effect only ----------------------------------------
permutest(betadisper(d_all, data_bc$genotype), permutations = 9999)
adonis2(d_all ~ genotype, data = data_bc, permutations = 9999)