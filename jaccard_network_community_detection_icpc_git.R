# Goal: Build Jaccard-based ICPC code/medication networks per age group and relate cluster burden to SES and region.
# Why: Network communities provide an interpretable view of co-occurring patterns that standard rates can miss.
# Depends on: Source assets and upstream cohort artifacts.
# Upstream cohort artifacts in 'pathname'.
#    - data.RDS: cohort with ICPC/topics content and stratifiers.
#    - data_lookup.RDS: Patient-year denominators used in rates/models.
# Output: In-memory network objects, cluster assignments, model summaries, and figures (no explicit file write).
#    - n/a: This script currently runs clustering/modelling workflows without persisting dedicated output files.

# Load libraries ----
library(haven)
library(dplyr)
library(skimr)
library(tidyverse)
library(tidymodels)
library(tidygraph)
library(purrr)
library(future)
library(furrr)
library(workflows)
library(broom.helpers)
library(lubridate)
library(forcats)
library(ggplot2)
library(data.table)
library(readxl)
library(duckdb)
library(UpSetR)
library(ComplexUpset)
library(bit64)
library(vctrs)
library(viridis)
library(igraph)
library(ggraph)
library(vegan)
library(ggforce)
library(IsingFit)
library(patchwork)
library(broom)
library(MASS)
library(lmtest)
library(sandwich)
library(gtsummary)
library(cluster)
library(writexl)

conflicted::conflict_prefer("select", "dplyr") # Prefer dplyr's functions
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(igraph::degree)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(bit64::`:`)

# LOAD DATA ----
journals_ses_country_ages_days <- readRDS("path/data.RDS")              # Final dataset with gp contacts
lookup_patient_years_final     <- readRDS("path/data_lookup.rds")       # Final population

# Patient-years per age_group x SES x region combination
patient_years_strata <- lookup_patient_years_final |>
  group_by(age_group, region_latest, ses_category_total) |>
  summarise(
    total_patient_years = sum(patient_years, na.rm = TRUE),
    n_patients          = n(),
    n_distinct_patients = n_distinct(patientID), # To check if this gives the same output as n()
    .groups = "drop")

testthat::expect_equal(patient_years_strata$n_patients, patient_years_strata$n_distinct_patients)
testthat::expect_true(all(patient_years_strata$total_patient_years > 0))
testthat::expect_equal(nrow(patient_years_strata |>
                              group_by(age_group, region_latest, ses_category_total) |>
                              summarise(n = n(), .groups = "drop") |>
                              filter(n > 1)), 0L)

# Colours of SES and country categories for visualisation
country_colours <- c("Netherlands" = "#73D055FF", "Europe" = "#B8DE29FF", "Outside Europe" = "#FDE725FF")
ses_colours     <- c("Low" = "#1F968BFF", "Middle" = "#29AF7FFF", "High" = "#55C667FF")

# Set seed for reproducibility
set.seed(1301)

icpc_binary_df <- journals_ses_country_ages_days |>
  select(patientID, age_group, icpc_clean, icpc_group, ses_category_total, region_latest, cpi_id_crypt, gp_note_date) |>
  mutate(icpc_num = as.integer(str_extract(icpc_clean, "[0-90]+"))) |>
  filter(icpc_num < 30 | icpc_num > 69) |>                          # Remove process codes (30-69)
  filter(icpc_num != 0) |>
  filter(!str_starts(icpc_clean, "Y"))

# Check if there are no ICPC codes starting with Y (= male problems)
testthat::expect_equal(icpc_binary_df |> filter(str_starts(icpc_clean, "Y")) |> nrow(), 0L)

contact_df <- icpc_binary_df |>
  select(cpi_id_crypt, patientID, age_group, ses_category_total, region_latest, gp_note_date)

check <- icpc_binary_df |>
  count(icpc_num)

testthat::expect_true(all(check$icpc_num %in% c(1:29, 69:99)))

# Data preparation ----
## Build patient ICPC incidence matrix (binary) ----
# Matrix for 18-25 years
icpc_matrix_18_25 <- icpc_binary_df |>
  filter(age_group == "18-25 years") |>
  distinct(patientID, icpc_clean) |>
  mutate(value = 1) |>
  pivot_wider(
    names_from  = icpc_clean,
    values_from = value,
    values_fill = 0)

testthat::expect_true(all(as.matrix(select(icpc_matrix_18_25, -patientID)) %in% c(0, 1)))
testthat::expect_equal(ncol(icpc_matrix_18_25) - 1L, n_distinct(filter(icpc_binary_df, age_group == "18-25 years")$icpc_clean))
testthat::expect_equal(nrow(icpc_matrix_18_25), n_distinct(filter(icpc_binary_df, age_group == "18-25 years")$patientID))

# Matrix for 26-35 years
icpc_matrix_26_35 <- icpc_binary_df |>
  filter(age_group == "26-35 years") |>
  distinct(patientID, icpc_clean) |>
  mutate(value = 1) |>
  pivot_wider(
    names_from  = icpc_clean,
    values_from = value,
    values_fill = 0)

testthat::expect_true(all(as.matrix(select(icpc_matrix_26_35, -patientID)) %in% c(0, 1)))
testthat::expect_equal(ncol(icpc_matrix_26_35) - 1L, n_distinct(filter(icpc_binary_df, age_group == "26-35 years")$icpc_clean))
testthat::expect_equal(nrow(icpc_matrix_26_35), n_distinct(filter(icpc_binary_df, age_group == "26-35 years")$patientID))

# Matrix for 36-45 years
icpc_matrix_36_45 <- icpc_binary_df |>
  filter(age_group == "36-45 years") |>
  distinct(patientID, icpc_clean) |>
  mutate(value = 1) |>
  pivot_wider(
    names_from  = icpc_clean,
    values_from = value,
    values_fill = 0)

testthat::expect_true(all(as.matrix(select(icpc_matrix_36_45, -patientID)) %in% c(0, 1)))
testthat::expect_equal(ncol(icpc_matrix_36_45) - 1L, n_distinct(filter(icpc_binary_df, age_group == "36-45 years")$icpc_clean))
testthat::expect_equal(nrow(icpc_matrix_36_45), n_distinct(filter(icpc_binary_df, age_group == "36-45 years")$patientID))

# Matrix for 46-55 years
icpc_matrix_46_55 <- icpc_binary_df |>
  filter(age_group == "46-55 years") |>
  distinct(patientID, icpc_clean) |>
  mutate(value = 1) |>
  pivot_wider(
    names_from  = icpc_clean,
    values_from = value,
    values_fill = 0)

testthat::expect_true(all(as.matrix(select(icpc_matrix_46_55, -patientID)) %in% c(0, 1)))
testthat::expect_equal(ncol(icpc_matrix_46_55) - 1L, n_distinct(filter(icpc_binary_df, age_group == "46-55 years")$icpc_clean))
testthat::expect_equal(nrow(icpc_matrix_46_55), n_distinct(filter(icpc_binary_df, age_group == "46-55 years")$patientID))

# =============================================================================
# Pipeline per age group:
#   prevalence filter -> pairwise Jaccard similarity between codes ->
#   edge-thresholded weighted undirected graph -> Louvain community detection
#   (resolution 1.5) -> node/edge/coordinate output
# =============================================================================

# =============================================================================
# I. ICPC CODES CLUSTERS
# =============================================================================

# A. Helper function: Jaccard co-occurrence ----
# Takes a binary patient x code matrix (rows = patients, columns = ICPC codes),
# and age group label for plot titles, and an edge weight threshold below which
# edges are removed from the network to to reduce noise.
build_jaccard_network <- function(m,
                                  age_label,
                                  edge_threshold = 0.075,
                                  min_code_prev  = 0.01,
                                  resolution     = 1.5) {
  
  message("Processing: ", age_label)
  
  # Remove patient ID column
  m <- m |>
    select(-patientID) |>
    as.matrix()
  
  # Remove rare codes below prevalence threshold
  code_prev <- colMeans(m)
  m <- m[, code_prev >= min_code_prev]
  
  message(" Codes retained after prevalence filter: ", ncol(m))
  
  # Step 1: Jaccard similarity between ICPC codes (codes as rows via t(m)) ----
  jac <- 1 - as.matrix(vegdist(t(m),                 # vegdist computes pairwise Jaccard distances between codes.
                               method = "jaccard",   # Subtracting from 1 converts distance to similarity: 1 = identical co-occurrence, 0 = no overlap.
                               binary = TRUE))        # ... based on binary co-occurrence across patients
  
  diag(jac) <- 0                                      # The diagonal (self-similarity) is set to 0 to prevent self-loops.
  
  jac[jac < edge_threshold] <- 0                      # Apply the minimum similarity threshold
  
  # Step 2: build weighted undirected graph: remove isolated nodes ----
  g <- graph_from_adjacency_matrix(
    jac,
    mode     = "undirected",
    weighted = TRUE,
    diag     = FALSE)
  
  g <- delete_vertices(g, degree(g) == 0)
  
  n_before_isolation <- vcount(g)
  g <- delete_vertices(g, degree(g) == 0)
  message(" Isolated nodes removed: ", n_before_isolation - vcount(g),
          " (", round(1 - vcount(g) / n_before_isolation, 2), " of prevalence-filtered codes)")
  
  # Step 3: Louvain community detection on edge-weighted graph ----
  set.seed(1301)
  clust_louvain <- cluster_louvain(g,
                                   weights    = E(g)$weight,   # Use edge weights for community detection
                                   resolution = resolution)    # Controls the number of communities: try values between 1.2 and 2.0 (1.5 is appropriate)
  V(g)$cluster <- membership(clust_louvain)
  
  g <- delete_vertices(g, is.na(V(g)$cluster)) # Remove clusters = NA
  
  # Standardise cluster levels: reorder by size (largest = community 1)
  community_sizes <- sort(table(V(g)$cluster), decreasing = TRUE)
  new_labels      <- setNames(seq_along(community_sizes),
                              names(community_sizes))
  V(g)$cluster    <- new_labels[as.character(V(g)$cluster)]
  
  # Quality check: fragmentation + modularity
  cluster_sizes <- table(V(g)$cluster)
  message(" n_clusters: ", length(cluster_sizes),
          " | largest cluster share: ", round(max(cluster_sizes) / sum(cluster_sizes), 2),
          " | n_singletons: ", sum(cluster_sizes == 1))
  message(" modularity (resolution = ", resolution, "): ",
          round(igraph::modularity(g, V(g)$cluster, weights = E(g)$weight), 3))
  
  # Step 4: node metadata: numeric ICPC part and ICPC type for node shape ----
  V(g)$icpc_numeric <- as.numeric(str_extract(V(g)$name, "\\d+"))
  V(g)$type <- case_when(
    V(g)$icpc_numeric < 30 ~ "Symptom diagnosis",
    V(g)$icpc_numeric > 69 ~ "Disease diagnosis",
    .default               = "Other")
  
  g_tbl <- as_tbl_graph(g) |>
    activate(edges) |>
    mutate(
      from_cluster  = .N()$cluster[from],
      to_cluster    = .N()$cluster[to],
      intra_cluster = from_cluster == to_cluster)
  
  # Step 5: Network PLOT: hulls per community, edge weighted by Jaccard similarity ----
  # nodes coloured by community and shaped by ICPC type
  p <- ggraph(g_tbl, layout = "graphopt",
              spring.length = 4,
              charge = 0.05) +
    
    # Make an area around the nodes of the community
    geom_mark_hull(
      aes(x      = x,
          y      = y,
          group  = as.factor(cluster),
          fill   = as.factor(cluster),
          colour = NA),
      concavity = 5,
      alpha     = 0.1,
      expand    = unit(4, "mm")) +
    geom_edge_arc(
      aes(width = weight),
      alpha    = 0.2,
      strength = 0.4) +
    geom_node_point(
      aes(colour = as.factor(cluster),
          shape  = type),
      size = 9) +
    # Annotation of the ICPC code in the node
    geom_node_text(
      aes(label = name),
      repel    = FALSE,
      size     = 3,
      colour   = "white",
      fontface = "bold") +      # , max.overlaps = 20
    scale_colour_viridis_d(option = "D", guide = "none", drop = FALSE) +
    scale_fill_viridis_d(option = "D", guide = "none", drop = FALSE) +
    scale_edge_width(range = c(0.2, 2), guide = "none") +
    scale_shape_manual(
      values = c("Symptom diagnosis" = 16,
                 "Disease diagnosis" = 15,
                 Other = 3),
      name = "ICPC type") +
    labs(title = age_label) +
    theme_graph(base_size = 11) +
    theme(
      legend.position  = "bottom",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA))
  
  # Graph faceted per cluster # New
  p_facet <- ggraph(g_tbl,
                    layout = "graphopt",
                    spring.length = 4,
                    charge        = 0.05) +
    geom_edge_link(
      aes(width = weight, filter = intra_cluster),
      alpha       = 0.3,
      show.legend = FALSE) +
    geom_node_point(
      aes(colour = as.factor(cluster), shape = type),
      size = 10,
      show.legend = FALSE) +
    geom_node_text(
      aes(label = name),
      repel    = FALSE,
      size     = 4,
      colour   = "white",
      fontface = "bold") +
    facet_nodes(~ cluster,     # One panel per cluster
                ncol = 1) +
    scale_colour_viridis_d(option = "D", guide = "none", drop = FALSE) +
    scale_fill_viridis_d(option = "D", guide = "none", drop = FALSE) +
    scale_edge_width(range = c(0.2, 2), name = "Jaccard similarity") +
    scale_shape_manual(
      values = c("Symptom diagnosis" = 16,
                 "Disease diagnosis" = 15,
                 Other               = 3),
      name = "ICPC type") +
    coord_cartesian(clip = "off") +
    labs(title = "Cluster") + # , caption = "Node proximity reflects Jaccard co-occurrence similarity"
    theme_graph(base_size = 12) +
    theme(
      plot.title            = element_text(size = 12),
      plot.caption          = element_text(hjust = 0),
      plot.caption.position = "plot",
      strip.text            = element_text(size = 12),
      strip.text.x          = element_blank(),
      strip.text.y          = element_blank(),
      strip.background      = element_blank(),
      legend.position       = "none",
      panel.spacing         = unit(1, "lines")) # Increase vertical space between panels
  
  list(
    graph       = g,
    tbl_graph   = g_tbl,
    communities = clust_louvain,
    plot        = p,
    plot_facet  = p_facet,
    age_group   = age_label)}

## Named list of binary patient x code matrix per age group ----
matrices <- list(
  "18-25 years" = icpc_matrix_18_25,
  "26-35 years" = icpc_matrix_26_35,
  "36-45 years" = icpc_matrix_36_45,
  "46-55 years" = icpc_matrix_46_55)

## Run the function for all age groups ----
# Default parameter settings (edge_threshold = 0.075, min_code_prev = 0.01, resolution = 1.5)
# FINAL (2026-08):
results_network <- matrices |>
  imap(\(m, ag) build_jaccard_network(m, age_label = ag))

# Console output of the final run:
# Processing: 18-25 years
# Codes retained after prevalence filter: 171
# Isolated nodes removed: 0 (0 of prevalence-filtered codes)
# n_clusters: 10 | largest cluster share: 0.25 | n_singletons: 0
# modularity (resolution = 1.5): 0.126
# Processing: 26-35 years
# Codes retained after prevalence filter: 194
# Isolated nodes removed: 0 (0 of prevalence-filtered codes)
# n_clusters: 11 | largest cluster share: 0.18 | n_singletons: 0
# modularity (resolution = 1.5): 0.071
# Processing: 36-45 years
# Codes retained after prevalence filter: 200
# Isolated nodes removed: 0 (0 of prevalence-filtered codes)
# n_clusters: 12 | largest cluster share: 0.18 | n_singletons: 0
# modularity (resolution = 1.5): 0.179
# Processing: 46-55 years
# Codes retained after prevalence filter: 215
# Isolated nodes removed: 0 (0 of prevalence-filtered codes)
# n_clusters: 12 | largest cluster share: 0.17 | n_singletons: 1
# modularity (resolution = 1.5): 0.094

# Test other parameters
# results_network <- matrices |>
#   imap(\(m, ag) build_jaccard_network(
#     m,
#     age_label      = ag,
#     edge_threshold = 0.085,
#     min_code_prev  = 0.01,
#     resolution     = 1.5))

# Node level membership and ICPC type
cluster_tables <- results_network |>
  map(\(res) {
    res$tbl_graph |>
      activate(nodes) |>
      as_tibble() |>
      transmute(
        icpc_code    = name,
        icpc_numeric = icpc_numeric,
        icpc_type    = type,
        cluster      = cluster) |>
      arrange(cluster, icpc_code)})

# Cluster-level summary table per age group
cluster_summary_tables <- results_network |>
  map(\(res) {
    node_df <- res$tbl_graph |> activate(nodes) |> as_tibble()
    edge_df <- res$tbl_graph |> activate(edges) |> as_tibble()
    
    intra_stats <- edge_df |>
      filter(intra_cluster) |>
      group_by(cluster = from_cluster) |>
      summarise(
        n_intra_edges      = n(),
        mean_jaccard_intra = mean(weight),
        .groups = "drop")
    
    node_df |>
      count(cluster, name = "n_codes") |>
      mutate(pct_of_nodes = round(100 * n_codes / sum(n_codes), 1)) |>
      left_join(intra_stats, by = "cluster") |>
      mutate(across(c(n_intra_edges, mean_jaccard_intra), ~ replace_na(.x, 0))) |>
      arrange(cluster)})

# Output of the cluster edges
edge_tables <- results_network |>
  map(\(res) {
    node_lookup <- res$tbl_graph |> activate(nodes) |> as_tibble() |> mutate(node_id = row_number())
    
    res$tbl_graph |> activate(edges) |> as_tibble() |>
      left_join(node_lookup |> select(node_id, from_code = name), by = c("from" = "node_id")) |>
      left_join(node_lookup |> select(node_id, to_code = name), by = c("to" = "node_id")) |>
      filter(intra_cluster) |>
      transmute(cluster = from_cluster, from_code, to_code, weight)})

# Node layout coordinates (x, y)
node_coords <- results_network |>
  imap(\(res, ag) {
    set.seed(1301)
    layout_tbl <- create_layout(res$tbl_graph, layout = "graphopt",
                                spring.length = 4,
                                charge = 0.05)
    layout_tbl |>
      as_tibble() |>
      select(name, x, y, cluster, type) |>
      mutate(age_group = ag)}) |>
  bind_rows()

# Alternative parameter settings explored (not used for the final results):
# results_network <- matrices |> imap(\(m, ag) build_jaccard_network(
#   m,
#   age_label      = ag,
#   edge_threshold = 0.1,
#   min_code_prev  = 0.01,
#   resolution     = 1.3))
# Gives best modularity scores, but clusters are not very well interpretable

# results_network <- matrices |> imap(\(m, ag) build_jaccard_network(
#   m,
#   age_label      = ag,
#   edge_threshold = 0.09,
#   min_code_prev  = 0.01,
#   resolution     = 1.2))

# Print all combined plots ----
walk(results_network, ~ print(.x$plot))
