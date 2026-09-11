# Metschnikowia OTU merging based on a FastTree2 (GTR+GAMMA) phylogeny of
# OTU sequences plus UNITE reference sequences. The merging threshold is
# calibrated per tree as the maximum intraspecific patristic distance among
# M. reukaufii references, and OTU clades below that threshold are merged
# into single features via greedy, non-overlapping node selection.

library(ape)
library(dplyr)
library(tidyr)
library(phytools)

# --- Import data -------------------------------------------------------
# FastTree_GTR_gamma.tree: FastTree2 output, OTU sequences (SHA1-hash
#   headers) + UNITE reference sequences for Metschnikowia and related
#   genera (Clavispora, Kodamaea)
raw <- readLines("data/FastTree_GTR_gamma.tree")

# --- Step 1: load and root the tree ---------------------------------------
full_string <- paste(raw, collapse = "")
full_string <- sub("^\xEF\xBB\xBF", "", full_string)          # strip BOM
full_string <- gsub("\\|", "_", full_string)
full_string <- gsub("_k__Fungi[^:,)]+", "", full_string)       # strip taxonomy suffix from tip labels
full_string <- trimws(full_string)

tree <- read.tree(text = full_string)
cat("Rooted before correction:", is.rooted(tree), "\n")

tree <- midpoint.root(tree)
cat("Rooted after correction:", is.rooted(tree), "\n")
cat("Total tips:", length(tree$tip.label), "\n")

# --- Step 2: separate OTU tips from reference tips -------------------------
otu_tips <- grep("^[0-9a-f]{40}$", tree$tip.label, value = TRUE)  # SHA1 hash = OTU
ref_tips <- setdiff(tree$tip.label, otu_tips)

cat("OTUs:", length(otu_tips), "\n")
cat("References:", length(ref_tips), "\n")

ref_species <- data.frame(
  tip = ref_tips,
  species = sub("^([A-Z][a-z]+_[a-z]+)_.*", "\\1", ref_tips),
  stringsAsFactors = FALSE
)

# --- Step 3: patristic distance matrix --------------------------------------
dist_mat <- cophenetic.phylo(tree)

# --- Step 4: calibrate the merging threshold --------------------------------
# Expand outward from the MRCA of the M. reukaufii references until the
# first non-reukaufii named species is encountered; the threshold is the
# max pairwise OTU distance at the node just before that point, plus a
# negligible epsilon to avoid ties at the boundary.

all_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
node_descendants <- lapply(all_nodes, function(n) extract.clade(tree, n)$tip.label)
names(node_descendants) <- all_nodes

node_info <- data.frame(
  node = all_nodes,
  n_tips = sapply(node_descendants, length),
  n_otus = sapply(node_descendants, function(d) sum(d %in% otu_tips)),
  stringsAsFactors = FALSE
)

get_ancestors <- function(tree, node) {
  root <- length(tree$tip.label) + 1
  anc <- c()
  current <- node
  while (current != root) {
    parent_row <- which(tree$edge[, 2] == current)
    if (length(parent_row) == 0) break
    parent <- tree$edge[parent_row, 1]
    anc <- c(anc, parent)
    current <- parent
  }
  anc
}

get_ref_species_in_node <- function(node) {
  refs_in <- intersect(node_descendants[[as.character(node)]], ref_tips)
  sp <- unique(ref_species$species[ref_species$tip %in% refs_in])
  sp[!grepl("_sp$", sp)]
}

reukaufii_refs <- ref_species$tip[ref_species$species == "Metschnikowia_reukaufii"]
anchor_node <- getMRCA(tree, reukaufii_refs)
cat("Anchor node (MRCA of reukaufii references):", anchor_node, "\n\n")

path <- c(anchor_node, get_ancestors(tree, anchor_node))

threshold <- NA
prev_max_dist <- NA

for (n in path) {
  otus_in <- intersect(node_descendants[[as.character(n)]], otu_tips)
  n_otu <- length(otus_in)
  max_d <- if (n_otu >= 2) max(dist_mat[otus_in, otus_in]) else NA
  species_here <- get_ref_species_in_node(n)
  other_species <- setdiff(species_here, "Metschnikowia_reukaufii")
  
  cat(sprintf("Node %d | n_otus=%d | max_dist=%s | other species: %s\n",
              n, n_otu,
              ifelse(is.na(max_d), "NA", sprintf("%.5f", max_d)),
              ifelse(length(other_species) == 0, "none", paste(other_species, collapse = ", "))))
  
  if (length(other_species) > 0) {
    threshold <- prev_max_dist + 1e-9
    cat(">>> Other species detected here. Threshold retained (previous node + epsilon):", threshold, "<<<\n")
    break
  }
  prev_max_dist <- max_d
}

# --- Step 5: identify OTU clades to merge -----------------------------------
cat("\nThreshold applied:", threshold, "\n")

candidates <- node_info %>%
  filter(n_otus >= 2)

candidates$max_dist_otus <- sapply(candidates$node, function(n) {
  d <- node_descendants[[as.character(n)]]
  otus_in <- intersect(d, otu_tips)
  if (length(otus_in) < 2) return(Inf)
  max(dist_mat[otus_in, otus_in])
})

candidates <- candidates %>%
  filter(max_dist_otus < threshold) %>%
  arrange(desc(n_otus))

# --- Step 6: greedy, non-overlapping selection ------------------------------
assigned_otus <- character(0)
keep <- logical(nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  all_desc <- node_descendants[[as.character(candidates$node[i])]]
  otus_in <- intersect(all_desc, otu_tips)
  if (!any(otus_in %in% assigned_otus)) {
    keep[i] <- TRUE
    assigned_otus <- c(assigned_otus, otus_in)
  }
}

final_groups <- candidates[keep, ]
cat("\n", nrow(final_groups), "groups identified,", length(assigned_otus), "OTUs merged\n")

# --- Step 7: annotate each group with its nearest reference species --------

find_nearest_ref <- function(otu_members, ref_tips, dist_mat, ref_species) {
  sub_dist <- dist_mat[otu_members, ref_tips, drop = FALSE]
  mean_dist_per_tip <- colMeans(sub_dist)
  best_ref_tip <- names(which.min(mean_dist_per_tip))
  best_dist <- mean_dist_per_tip[best_ref_tip]
  best_species <- ref_species$species[ref_species$tip == best_ref_tip]
  list(ref_tip = best_ref_tip, species = best_species, distance = best_dist)
}

# Annotation is always based on the true minimal MRCA of the group's OTUs,
# never on the node retained by the greedy selection (step 6), which can be
# an inflated ancestor spanning unrelated references.
annotate_otu_group <- function(otus_in) {
  minimal_node <- if (length(otus_in) >= 2) getMRCA(tree, otus_in) else NA
  minimal_desc <- if (!is.na(minimal_node)) extract.clade(tree, minimal_node)$tip.label else otus_in
  refs_in <- intersect(minimal_desc, ref_tips)
  
  if (length(refs_in) > 0) {
    ref_species_in <- unique(ref_species$species[ref_species$tip %in% refs_in])
    ref_species_in <- ref_species_in[!grepl("_sp$", ref_species_in)]  # drop uninformative "Metschnikowia_sp"
    ref_label <- if (length(ref_species_in) > 0) paste(ref_species_in, collapse = " / ") else "no_named_ref"
    dist_to_ref <- 0
  } else {
    info <- find_nearest_ref(otus_in, ref_tips, dist_mat, ref_species)
    ref_label <- ifelse(info$distance < threshold, info$species, "no_ref")
    dist_to_ref <- info$distance
  }
  
  list(minimal_node = minimal_node, refs_in = refs_in, ref_label = ref_label, dist_to_ref = dist_to_ref)
}

warning_log <- list()

group_annotations <- lapply(seq_len(nrow(final_groups)), function(i) {
  all_desc <- node_descendants[[as.character(final_groups$node[i])]]
  otus_in <- intersect(all_desc, otu_tips)
  
  ann <- withCallingHandlers(
    annotate_otu_group(otus_in),
    warning = function(w) {
      warning_log[[length(warning_log) + 1]] <<- data.frame(
        group_id = i, n_otus = length(otus_in),
        otus_in = paste(substr(otus_in, 1, 8), collapse = ","),
        warning_msg = conditionMessage(w),
        stringsAsFactors = FALSE
      )
      invokeRestart("muffleWarning")
    }
  )
  
  data.frame(
    group_id = i, node_used = final_groups$node[i], minimal_node = ann$minimal_node,
    node_inflated = final_groups$node[i] != ann$minimal_node,
    n_otus = length(otus_in), max_dist = final_groups$max_dist_otus[i],
    ref_in_clade = length(ann$refs_in) > 0, ref_nearby = ann$ref_label,
    dist_to_ref = ann$dist_to_ref, stringsAsFactors = FALSE
  )
})

warning_df <- bind_rows(warning_log)
cat(nrow(warning_df), "warnings captured\n")
print(warning_df)

group_summary <- bind_rows(group_annotations)
cat("\n--- Group summary (minimal node, corrected) ---\n")
print(group_summary)

# --- Step 8: final merge table -----------------------------------------------
merge_list <- lapply(seq_len(nrow(final_groups)), function(i) {
  all_desc <- node_descendants[[as.character(final_groups$node[i])]]
  otus_in <- intersect(all_desc, otu_tips)
  ann <- annotate_otu_group(otus_in)
  
  data.frame(
    otu_original = otus_in, otu_merged = otus_in[1], group_id = i,
    ref_nearby = ann$ref_label, dist_to_ref = ann$dist_to_ref,
    stringsAsFactors = FALSE
  )
})

merge_table <- bind_rows(merge_list)

singleton_otus <- setdiff(otu_tips, assigned_otus)

singleton_annotations <- lapply(singleton_otus, function(otu) {
  ann <- annotate_otu_group(otu)
  data.frame(
    otu_original = otu, otu_merged = otu, group_id = NA,
    ref_nearby = ann$ref_label, dist_to_ref = ann$dist_to_ref,
    stringsAsFactors = FALSE
  )
})

singletons <- bind_rows(singleton_annotations)
merge_table <- bind_rows(merge_table, singletons)

cat("\n=== FINAL SUMMARY ===\n")
cat("Initial OTUs:", length(otu_tips), "\n")
cat("Final OTUs:", length(unique(merge_table$otu_merged)), "\n")
cat("  - groups:", nrow(final_groups), "containing", length(assigned_otus), "OTUs\n")
cat("  - singletons:", length(singleton_otus), "\n")

# --- Step 9: circular tree figure (supplementary) ---------------------------
# ggtree conflicts with phyloseq if it is still loaded in the session
try(detach("package:phyloseq", unload = TRUE), silent = TRUE)
try(unloadNamespace("phyloseq"), silent = TRUE)

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)

pal_manual <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
  "#A65628", "#F781BF", "#00CED1", "#FFD700", "#1B9E77",
  "#7570B3", "#D95F02", "#E6AB02", "#666666", "#A6761D",
  "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462"
)

tip_meta <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(
    merge_table %>% dplyr::select(label = otu_original, group_id, ref_nearby),
    by = "label"
  ) %>%
  mutate(
    ref_nearby_short = case_when(
      grepl("pulcherrima|fructicola|kunwiensis|picachoensis|leonuri|sinensis|rubicola|citriensis|chrysoperlae|pimensis|_sp$", ref_nearby) ~ "Metschnikowia pulcherrima complex",
      grepl("chrysomelidarum|rancensis|vanudenii", ref_nearby) ~ "Metschnikowia spp. (unresolved clade)",
      TRUE ~ gsub("_", " ", ref_nearby)
    ),
    type = case_when(
      label %in% ref_tips ~ "reference",
      !is.na(group_id) ~ paste0("Cluster ", group_id, " — ", ref_nearby_short),
      TRUE ~ "singleton"
    ),
    display = ifelse(
      label %in% ref_tips,
      gsub("_", " ", sub("^([A-Z][a-z]+_[a-z]+).*", "\\1", label)),
      ""
    )
  )

cluster_levels_full <- tip_meta %>%
  filter(!label %in% ref_tips, !is.na(group_id)) %>%
  pull(type) %>% unique() %>%
  .[order(as.integer(sub("Cluster ([0-9]+).*", "\\1", .)))]

type_levels <- c(cluster_levels_full, "singleton", "reference")
type_cols <- c(
  setNames(pal_manual[seq_len(length(cluster_levels_full))], cluster_levels_full),
  "singleton" = "grey80",
  "reference" = "#1a1a1a"
)

tip_meta$type <- factor(tip_meta$type, levels = type_levels)

p <- ggtree(tree, layout = "circular", branch.length = "none",
            linewidth = 0.2, color = "grey40") %<+% tip_meta

p <- p +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(y = label, fill = type), width = 1.5, offset = 0.02) +
  scale_fill_manual(
    values = type_cols, name = "Assigned species", na.value = "grey90",
    guide = guide_legend(
      override.aes = list(size = 10),
      title.theme = element_text(size = 28, face = "bold"),
      label.theme = element_text(size = 25, face = "italic")
    )
  ) +
  theme(legend.position = "right", legend.key.size = unit(1, "cm"))

p <- p + geom_tiplab(
  aes(label = display), size = 4.5, offset = 1.75, hjust = 0,
  fontface = "italic", color = "grey30"
)

ggsave("output/metsch_tree.tiff", plot = p, width = 30, height = 30, units = "in", compression = "lzw")
ggsave("output/metsch_tree.png", plot = p, width = 30, height = 30, units = "in")

# --- Step 10: export --------------------------------------------------------
qiime_merge <- merge_table %>%
  dplyr::select(feature_id = otu_original, merged_id = otu_merged)

write.table(qiime_merge, file = "output/metschnikowia_merge_map.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(merge_table, file = "output/metschnikowia_merge_detail.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)