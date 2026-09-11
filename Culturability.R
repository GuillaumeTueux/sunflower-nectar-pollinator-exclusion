# Nectar culturability: growth scores across pollinator access treatments,
# per cultivar. Wilcoxon pairwise comparisons on per-sample means, BH-adjusted.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(multcompView)

set.seed(64121)

# --- Import data -------------------------------------------------------
# growth_scores.csv: one row per plate
#   Genotype, Modality, Sample_Name, Trial, Counts (growth score, 0-4 scale)
data <- read.csv2("data/growth_scores.csv", header = TRUE)

data <- data %>%
  mutate(Modality = factor(
    dplyr::recode(as.character(Modality),
                  "No_Visitors"    = "No access",
                  "All_Visitors"   = "Continuous access",
                  "Day_Visitors"   = "Daytime access",
                  "Night_Visitors" = "Nighttime access"),
    levels = c("No access", "Continuous access", "Daytime access", "Nighttime access")
  ))

modality_colours <- c(
  "No access"         = "#000000",
  "Continuous access" = "#CC79A7",
  "Daytime access"    = "#F5C400",
  "Nighttime access"  = "#0072B2"
)

# one row per sample: mean growth score across the four culture media
data_mean <- data %>%
  group_by(Genotype, Modality, Sample_Name, Trial) %>%
  summarise(mean_Counts = mean(Counts, na.rm = TRUE), .groups = "drop")

data_mean %>% dplyr::count(Genotype, Modality)

data_celesto <- data_mean %>% filter(Genotype == "CELESTO")
data_idillic <- data_mean %>% filter(Genotype == "IDILLIC")

# --- Pairwise Wilcoxon tests --------------------------------------------
modality_pairs <- combn(
  c("No access", "Continuous access", "Daytime access", "Nighttime access"),
  2, simplify = FALSE
)

wilcox_pairwise <- function(df, pair, var) {
  g1 <- df %>% filter(Modality == pair[1]) %>% pull({{ var }})
  g2 <- df %>% filter(Modality == pair[2]) %>% pull({{ var }})
  wt <- wilcox.test(g1, g2, exact = FALSE)
  data.frame(
    group1 = pair[1], group2 = pair[2],
    n_group1 = length(g1), n_group2 = length(g2),
    W = as.numeric(wt$statistic), p.value = wt$p.value
  )
}

res_c_Counts <- bind_rows(lapply(modality_pairs, wilcox_pairwise, df = data_celesto, var = mean_Counts))
res_c_Counts$p.adj <- p.adjust(res_c_Counts$p.value, method = "BH")

res_i_Counts <- bind_rows(lapply(modality_pairs, wilcox_pairwise, df = data_idillic, var = mean_Counts))
res_i_Counts$p.adj <- p.adjust(res_i_Counts$p.value, method = "BH")

print(res_c_Counts)
print(res_i_Counts)

# --- Compact letter display ---------------------------------------------
get_cld <- function(res_wilcox, var_name, geno) {
  pvec <- setNames(
    res_wilcox$p.adj,
    paste(gsub(" ", "_", res_wilcox$group1), gsub(" ", "_", res_wilcox$group2), sep = "-")
  )
  cld <- multcompLetters(pvec, threshold = 0.05)$Letters
  data.frame(
    Modality = gsub("_", " ", names(cld)),
    letter = cld, Genotype = geno, Variable = var_name, row.names = NULL
  )
}

cld_c_Counts <- get_cld(res_c_Counts, "Growth score", "CELESTO")
cld_i_Counts <- get_cld(res_i_Counts, "Growth score", "IDILLIC")

cld_all <- bind_rows(cld_c_Counts, cld_i_Counts)
print(cld_all)

# --- Supplementary table --------------------------------------------------
group_means <- data_mean %>%
  group_by(Genotype, Modality) %>%
  summarise(
    mean_val = mean(mean_Counts, na.rm = TRUE),
    sd_val = sd(mean_Counts, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Variable = "Growth score", mean_fmt = sprintf("%.2f ± %.2f", mean_val, sd_val)) %>%
  dplyr::select(Genotype, Variable, Modality, mean_fmt)

supp_table <- bind_rows(
  res_c_Counts %>% mutate(Genotype = "CELESTO", Variable = "Growth score"),
  res_i_Counts %>% mutate(Genotype = "IDILLIC", Variable = "Growth score")
) %>%
  mutate(
    sig = ifelse(p.adj < 0.001, "***", ifelse(p.adj < 0.01, "**", ifelse(p.adj < 0.05, "*", "ns"))),
    p.value_fmt = ifelse(p.value < 0.001, "< 0.001", formatC(p.value, format = "f", digits = 3)),
    p.adj_fmt = ifelse(p.adj < 0.001, "< 0.001", formatC(p.adj, format = "f", digits = 3)),
    W = round(W, 1)
  ) %>%
  left_join(group_means, by = c("Genotype", "Variable", "group1" = "Modality")) %>%
  dplyr::rename(Mean_1 = mean_fmt) %>%
  left_join(group_means, by = c("Genotype", "Variable", "group2" = "Modality")) %>%
  dplyr::rename(Mean_2 = mean_fmt) %>%
  dplyr::select(Genotype, Variable, group1, n_group1, Mean_1, group2, n_group2, Mean_2,
                W, p.value_fmt, p.adj_fmt, sig) %>%
  dplyr::rename(
    Comparison_1 = group1, n1 = n_group1,
    Comparison_2 = group2, n2 = n_group2,
    `p-value` = p.value_fmt, `p.adj (BH)` = p.adj_fmt, Significance = sig
  )

stopifnot(!anyNA(supp_table$Mean_1), !anyNA(supp_table$Mean_2))

print(supp_table)
write.csv(supp_table, "output/supp_table_culturability_wilcoxon.csv", row.names = FALSE)

# --- Figure: growth score per plate, by cultivar and treatment -----------
summary_count <- data_mean %>%
  group_by(Genotype, Modality) %>%
  summarise(mean_count = mean(mean_Counts), se_count = sd(mean_Counts) / sqrt(n()), .groups = "drop")

ymax_count <- data_mean %>%
  group_by(Genotype, Modality) %>%
  summarise(top = max(mean_Counts, na.rm = TRUE), .groups = "drop")

cld_count <- cld_all %>%
  left_join(summary_count, by = c("Genotype", "Modality")) %>%
  left_join(ymax_count, by = c("Genotype", "Modality")) %>%
  group_by(Genotype) %>%
  mutate(y_pos = pmax(mean_count + se_count, top) + 0.08 * max(top)) %>%
  ungroup()

p_count <- ggplot(summary_count, aes(x = Modality, y = mean_count, fill = Modality)) +
  geom_col(width = 0.68, alpha = 0.85, colour = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(mean_count - se_count, 0), ymax = mean_count + se_count),
                width = 0.18, linewidth = 0.8) +
  geom_jitter(data = data_mean, aes(x = Modality, y = mean_Counts),
              width = 0.08, size = 2.5, alpha = 0.9, shape = 21,
              fill = "white", colour = "black", stroke = 0.5, inherit.aes = FALSE) +
  geom_text(data = cld_count, aes(x = Modality, y = y_pos, label = letter),
            vjust = -0.4, size = 10, inherit.aes = FALSE) +
  scale_fill_manual(values = modality_colours) +
  facet_wrap(~ Genotype) +
  labs(x = NULL, y = "Growth score per plate") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
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

print(p_count)

# --- Ceiling effect: Celesto under diurnal access -------------------------
ceiling_plate <- data %>%
  filter(Genotype == "CELESTO", Modality %in% c("Daytime access", "Continuous access")) %>%
  group_by(Modality) %>%
  summarise(n_plates = n(), pct_at_4 = mean(Counts == 4, na.rm = TRUE) * 100, .groups = "drop")
print(ceiling_plate)

ceiling_sample <- data_mean %>%
  filter(Genotype == "CELESTO", Modality %in% c("Daytime access", "Continuous access")) %>%
  group_by(Modality) %>%
  summarise(
    n = n(), n_at_ceiling = sum(mean_Counts == 4),
    samples_below = paste(round(mean_Counts[mean_Counts < 4], 2), collapse = ", "),
    .groups = "drop"
  )
print(ceiling_sample)

# --- Growth score relative to no-access, Celesto ---------------------------
means_celesto <- data_mean %>%
  filter(Genotype == "CELESTO") %>%
  group_by(Modality) %>%
  summarise(mean_count = mean(mean_Counts, na.rm = TRUE), n = n(), .groups = "drop")

no_access_mean <- means_celesto$mean_count[means_celesto$Modality == "No access"]

means_celesto <- means_celesto %>%
  mutate(pct_vs_no_access = round((mean_count - no_access_mean) / no_access_mean * 100, 1))
print(means_celesto)

q_continuous_vs_no <- res_c_Counts %>%
  filter((group1 == "No access" & group2 == "Continuous access") |
           (group1 == "Continuous access" & group2 == "No access")) %>%
  pull(p.adj)

q_daytime_vs_no <- res_c_Counts %>%
  filter((group1 == "No access" & group2 == "Daytime access") |
           (group1 == "Daytime access" & group2 == "No access")) %>%
  pull(p.adj)

cat("Continuous vs No access: q =", round(q_continuous_vs_no, 4), "\n")
cat("Daytime vs No access: q =", round(q_daytime_vs_no, 4), "\n")