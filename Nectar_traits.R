# Nectar volume, sugar mass and sugar composition (sucrose, glucose, fructose)
# for the Celesto and Idillic sunflower cultivars.

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

set.seed(64121)

# --- Import data -------------------------------------------------------
# nectar_volume_mass.csv: one row per capillary
#   Genotype, Plant_Rep, Sample_Name, Capillar_ID, Sampling_Date,
#   Volume_Nectar_ul, mass_Sugar_ug, Sowing
data_nectar <- read.table(
  "data/nectar_volume_mass.csv",
  header = FALSE, sep = ";", dec = ".", skip = 1,
  stringsAsFactors = FALSE, fileEncoding = "latin1"
)
colnames(data_nectar) <- c(
  "Genotype", "Plant_Rep", "Sample_Name", "Capillar_ID",
  "Sampling_Date", "Volume_Nectar_ul", "mass_Sugar_ug", "Sowing"
)

# nectar_sugar_composition.xlsx: one row per technical replicate
#   Sample, Genotype, Sucrose, Fructose, Glucose (g/L, enzymatic assay)
df_sugar <- read_xlsx("data/nectar_sugar_composition.xlsx")

# --- Plot themes ---------------------------------------------------------
theme_nectar_traits <- theme_bw(base_size = 35) +
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 34),
    axis.text.y = element_text(size = 34),
    legend.title = element_text(size = 34),
    legend.text = element_text(size = 30),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(size = 36, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(4, 4, 4, 4, "mm")
  )

theme_sugar_composition <- theme_bw(base_size = 35) +
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 34),
    axis.text.y = element_text(size = 34),
    legend.title = element_text(size = 35),
    legend.text = element_text(size = 34),
    strip.text = element_text(size = 36, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(4, 4, 4, 4, "mm")
  )

# ===========================================================================
# Nectar volume and sugar mass
# ===========================================================================

data_nectar <- data_nectar %>%
  mutate(
    Sampling_Date = as.Date(Sampling_Date, format = "%d/%m/%Y"),
    Volume_Nectar_ul = as.numeric(Volume_Nectar_ul),
    mass_Sugar_ug = as.numeric(mass_Sugar_ug),
    Genotype = as.factor(Genotype),
    Sowing = as.factor(Sowing)
  )

# plant-level means, averaged across capillaries
data_nectar_plant <- data_nectar %>%
  group_by(Genotype, Sowing, Plant_Rep, Sample_Name) %>%
  summarise(
    Volume_mean_ul = mean(Volume_Nectar_ul, na.rm = TRUE),
    Sugar_mean_ug = mean(mass_Sugar_ug, na.rm = TRUE),
    n_capillaries = n(),
    .groups = "drop"
  )

wilcox.test(Volume_mean_ul ~ Genotype, data = data_nectar_plant, exact = FALSE)
wilcox.test(Sugar_mean_ug ~ Genotype, data = data_nectar_plant, exact = FALSE)
wilcox.test(Volume_mean_ul ~ Sowing, data = data_nectar_plant, exact = FALSE)
wilcox.test(Sugar_mean_ug ~ Sowing, data = data_nectar_plant, exact = FALSE)

data_nectar_plant %>%
  group_by(Genotype) %>%
  summarise(
    mean_vol = mean(Volume_mean_ul), sd_vol = sd(Volume_mean_ul),
    mean_sug = mean(Sugar_mean_ug), sd_sug = sd(Sugar_mean_ug)
  )

pal_genotype <- c("CELESTO" = "grey70", "IDILLIC" = "grey30")
shape_sowing <- c("April" = 16, "June" = 17)

df_panel_nectar <- data_nectar_plant %>%
  mutate(
    Genotype = factor(Genotype, levels = c("CELESTO", "IDILLIC")),
    Sowing = factor(Sowing, levels = c("April", "June"))
  ) %>%
  pivot_longer(
    cols = c(Volume_mean_ul, Sugar_mean_ug),
    names_to = "Trait", values_to = "Value"
  ) %>%
  mutate(
    Trait = factor(
      Trait,
      levels = c("Volume_mean_ul", "Sugar_mean_ug"),
      labels = c("Nectar volume (µL/floret)", "Sugar mass (µg/floret)")
    )
  )

ggplot(df_panel_nectar, aes(x = Genotype, y = Value)) +
  geom_boxplot(
    aes(fill = Genotype), width = 0.6, linewidth = 0.6, color = "black",
    outlier.shape = NA, alpha = 0.7
  ) +
  geom_jitter(aes(shape = Sowing), color = "black", width = 0.15, size = 5, alpha = 0.9) +
  facet_wrap(~ Trait, nrow = 2, scales = "free_y", strip.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  scale_fill_manual(values = pal_genotype, guide = "none") +
  scale_shape_manual(values = shape_sowing, name = "Sowing") +
  labs(x = NULL, y = NULL) +
  theme_classic() + theme_nectar_traits +
  theme(strip.placement = "outside")

# ===========================================================================
# Sugar composition (sucrose vs hexose)
# ===========================================================================

df_prop <- df_sugar %>%
  mutate(
    Sucrose_clean = pmax(Sucrose, 0),
    Fructose_clean = pmax(Fructose, 0),
    Glucose_clean = pmax(Glucose, 0),
    Total_sugar = Sucrose_clean + Fructose_clean + Glucose_clean
  ) %>%
  filter(Total_sugar > 0) %>%
  mutate(
    Sucrose_prop = Sucrose_clean / Total_sugar,
    Fructose_prop = Fructose_clean / Total_sugar,
    Glucose_prop = Glucose_clean / Total_sugar
  )

df_mean <- df_prop %>%
  pivot_longer(
    cols = c(Sucrose_prop, Fructose_prop, Glucose_prop),
    names_to = "Sugar", values_to = "Prop"
  ) %>%
  mutate(Sugar = dplyr::recode(
    Sugar, Sucrose_prop = "Sucrose", Fructose_prop = "Fructose", Glucose_prop = "Glucose"
  )) %>%
  group_by(Sample, Genotype, Sugar) %>%
  summarise(Prop = mean(Prop, na.rm = TRUE), .groups = "drop")

# hexose = glucose + fructose
df_mean_hex <- df_mean %>%
  pivot_wider(names_from = Sugar, values_from = Prop) %>%
  mutate(Hexose = coalesce(Glucose, 0) + coalesce(Fructose, 0)) %>%
  select(Sample, Genotype, Sucrose, Hexose) %>%
  pivot_longer(cols = c(Sucrose, Hexose), names_to = "Sugar", values_to = "Prop") %>%
  mutate(Sugar = factor(Sugar, levels = c("Sucrose", "Hexose")), Prop_pct = 100 * Prop)

df_bar_se <- df_mean_hex %>%
  group_by(Genotype, Sugar) %>%
  summarise(
    mean_pct = 100 * mean(Prop, na.rm = TRUE),
    sd_pct = 100 * sd(Prop, na.rm = TRUE),
    n = sum(!is.na(Prop)),
    se_pct = sd_pct / sqrt(n),
    .groups = "drop"
  )

ggplot(df_bar_se, aes(x = Genotype, y = mean_pct, fill = Sugar)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(
    data = df_mean_hex,
    aes(x = Genotype, y = Prop_pct, group = Sugar),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    shape = 21, size = 4, fill = "white", colour = "black", inherit.aes = FALSE
  ) +
  scale_fill_grey(start = 0.2, end = 0.8) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.24))) +
  labs(x = NULL, y = "Relative sugar proportion (%)", fill = "Sugar") +
  theme_sugar_composition

kruskal.test(Prop ~ Genotype, data = filter(df_mean_hex, Sugar == "Hexose"))