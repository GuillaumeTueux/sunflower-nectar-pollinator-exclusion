# Pollinator visitation from camera-trap detections: cleaning, day/night
# windows for the continuous-access treatment, per-plant visit rates, and
# consistency checks between continuous access and the single-window
# treatments (daytime vs continuous, nighttime vs continuous).

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(patchwork)

# --- Import data -------------------------------------------------------
# pollinator_detections_raw.csv: one row per photo (";"-separated)
#   plant_id, capture_date, image_name,
#   predpolli_3cls_v_09_25_Bee, predpolli_3cls_v_09_25_Bumblebee,
#   predpolli_3cls_v_09_25_Moth_Butterfly
data_raw <- read.csv("data/pollinator_detections_raw.csv", header = TRUE, sep = ";")

# metadata_2024.xlsx / metadata_2025.xlsx: one row per plant
#   PLANT_ID, Genotype, Modality
metadata_2024 <- readxl::read_excel("data/metadata_2024.xlsx", sheet = 1)
metadata_2025 <- readxl::read_excel("data/metadata_2025.xlsx", sheet = 1)
metadata_all <- bind_rows(metadata_2024, metadata_2025)

# --- Clean timestamps and known artefacts --------------------------------
data_raw <- data_raw %>%
  mutate(
    capture_date = as.POSIXct(capture_date, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"),
    Date = as.Date(capture_date),
    time_rounded = sprintf("%02d:00", hour(capture_date)),
    Time = as.POSIXct(paste(Date, time_rounded), format = "%Y-%m-%d %H:%M", tz = "Europe/Paris")
  ) %>%
  filter(!year(Date) %in% c(2000, 2016))  # camera date-stamp artefacts

# plants excluded due to camera malfunction (flash failure, wrong date stamp, misalignment)
plants_to_remove <- c(
  "24MO11.0604.14", "25ZM05.10.3.06", "25ZM05.10.3.17",
  "25ZM05.11.4.02", "25ZM05.11.5.09", "25ZM05.11.6.02"
)

data_raw <- data_raw %>% filter(!plant_id %in% plants_to_remove)

# --- Merge metadata and format -------------------------------------------
data <- data_raw %>%
  left_join(metadata_all %>% dplyr::select(PLANT_ID, Genotype, Modality),
            by = c("plant_id" = "PLANT_ID")) %>%
  rename(
    PLANT_ID = plant_id,
    Bee = predpolli_3cls_v_09_25_Bee,
    Bumblebee = predpolli_3cls_v_09_25_Bumblebee,
    Moth.Butterfly = predpolli_3cls_v_09_25_Moth_Butterfly
  ) %>%
  mutate(capture_dt = as.POSIXct(capture_date, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"))

# fall back to the timestamp embedded in the image file name when capture_date is missing
data <- data %>%
  mutate(capture_dt = if_else(
    is.na(capture_dt),
    as.POSIXct(sub(".*_([0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{2}-[0-9]{2}-[0-9]{2}).*", "\\1", image_name),
               format = "%Y-%m-%d-%H-%M-%S", tz = "Europe/Paris"),
    capture_dt
  ))

# --- Day/night window for continuous-access plants ------------------------
# 22:00 cutoff, matching the daytime/nighttime treatment split
anchors <- data %>%
  filter(Modality == "All_Visitors") %>%
  group_by(PLANT_ID) %>%
  summarise(cutoff_night = as.POSIXct(paste(as.Date(min(capture_dt)), "22:00:00"),
                                      format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"))

# --- Per-plant mean visit rates -------------------------------------------
plant_means <- data %>%
  left_join(anchors, by = "PLANT_ID") %>%
  mutate(period = if_else(Modality == "All_Visitors",
                          if_else(capture_dt < cutoff_night, "day", "night"),
                          NA_character_)) %>%
  group_by(PLANT_ID, Genotype, Modality) %>%
  summarise(
    mean_bee       = mean(Bee[Modality != "All_Visitors" | period == "day"]),
    mean_bumblebee = mean(Bumblebee[Modality != "All_Visitors" | period == "day"]),
    mean_moth      = mean(Moth.Butterfly[Modality != "All_Visitors" | period == "night"]),
    .groups = "drop"
  )

write.csv(plant_means, "output/plant_means.csv", row.names = FALSE)

plant_means %>%
  group_by(Genotype, Modality) %>%
  summarise(
    n = n(),
    mean_bee_g = mean(mean_bee), max_bee_g = max(mean_bee),
    mean_bumblebee_g = mean(mean_bumblebee), max_bumblebee_g = max(mean_bumblebee),
    mean_moth_g = mean(mean_moth), max_moth_g = max(mean_moth),
    .groups = "drop"
  ) %>%
  arrange(Modality, Genotype) %>%
  print(n = Inf, width = Inf)

# --- Figure: visit rate per photo, by cultivar and treatment --------------
plant_summary <- plant_means %>%
  pivot_longer(cols = c(mean_bee, mean_bumblebee, mean_moth),
               names_to = "Pollinator", values_to = "mean_per_photo") %>%
  mutate(Pollinator = dplyr::recode(Pollinator,
                                    mean_bee = "Bee", mean_bumblebee = "Bumblebee",
                                    mean_moth = "Moth.Butterfly"))

make_panel <- function(df, genotype_filter) {
  ggplot(df %>% filter(Genotype == genotype_filter),
         aes(x = Modality, y = mean_per_photo, fill = Modality)) +
    geom_boxplot(outlier.shape = NA, width = 0.65) +
    geom_jitter(width = 0.15, alpha = 1, size = 1.5, color = "black") +
    facet_grid(Pollinator ~ ., scales = "free_y",
               labeller = labeller(Pollinator = c(
                 "Moth.Butterfly" = "Moth", "Bee" = "Bee", "Bumblebee" = "Bumblebee"
               ))) +
    scale_fill_manual(values = c(
      "Day_Visitors" = "#F5C400", "Night_Visitors" = "#0072B2", "All_Visitors" = "#CC79A7"
    )) +
    scale_x_discrete(labels = c(
      "Day_Visitors" = "Daytime access", "Night_Visitors" = "Nighttime access",
      "All_Visitors" = "Continuous access"
    )) +
    labs(x = NULL, y = "Mean visits per photo", title = genotype_filter) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 30),
      axis.title = element_text(size = 36),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      plot.title = element_text(size = 32, face = "bold", hjust = 0),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      plot.margin = unit(rep(4, 4), "mm")
    )
}

p_celesto <- make_panel(plant_summary, "CELESTO")
p_idillic <- make_panel(plant_summary, "IDILLIC")

p_celesto / p_idillic

ggsave("output/pollinator_visits_per_photo.png", width = 14, height = 20, units = "in", dpi = 300)

# --- Supplementary table: plants/photos/detections per cultivar x treatment x trial ---
supp_table <- data %>%
  mutate(Trial = substr(PLANT_ID, 1, 6)) %>%
  left_join(anchors, by = "PLANT_ID") %>%
  mutate(period = if_else(Modality == "All_Visitors",
                          if_else(capture_dt < cutoff_night, "day", "night"),
                          NA_character_)) %>%
  group_by(Genotype, Modality, Trial) %>%
  summarise(
    n_plants = n_distinct(PLANT_ID),
    total_photos = n(),
    total_bee = sum(Bee),
    total_bumblebee = sum(Bumblebee),
    total_moth = sum(Moth.Butterfly),
    total_photos_night = sum(period == "night", na.rm = TRUE),
    total_photos_day = sum(period == "day", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Genotype, Modality, Trial)

write.csv(supp_table, "output/supp_table_visits.csv", row.names = FALSE)

# --- Consistency between continuous access and single-window treatments ---
compare_modality <- function(df, guild_var, mod_a, mod_b, genotype_filter) {
  d <- df %>% filter(Genotype == genotype_filter, Modality %in% c(mod_a, mod_b))
  test <- wilcox.test(d[[guild_var]] ~ d$Modality)
  data.frame(
    Genotype = genotype_filter, Guild = guild_var,
    Comparison = paste(mod_a, "vs", mod_b),
    W = test$statistic, p_value = test$p.value
  )
}

results_similarity <- bind_rows(
  compare_modality(plant_means, "mean_bee", "All_Visitors", "Day_Visitors", "CELESTO"),
  compare_modality(plant_means, "mean_bee", "All_Visitors", "Day_Visitors", "IDILLIC"),
  compare_modality(plant_means, "mean_bumblebee", "All_Visitors", "Day_Visitors", "CELESTO"),
  compare_modality(plant_means, "mean_bumblebee", "All_Visitors", "Day_Visitors", "IDILLIC"),
  compare_modality(plant_means, "mean_moth", "All_Visitors", "Night_Visitors", "CELESTO"),
  compare_modality(plant_means, "mean_moth", "All_Visitors", "Night_Visitors", "IDILLIC")
) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

print(results_similarity, row.names = FALSE)