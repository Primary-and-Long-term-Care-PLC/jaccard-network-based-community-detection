# Jaccard network-based community detection of GP-registered morbidity

Analysis code accompanying the manuscript on the **clinical burden among women in Dutch general-practice (GP) care across the life course**.

This repository contains the R code used to construct Jaccard co-occurrence networks of ICPC diagnosis codes and to detect communities (clusters) of frequently co-occurring codes, separately for each age group, using the Louvain algorithm. It also contains the downstream code that exports node, edge and layout tables so the network figures can be reproduced outside the secure analysis environment.

> **Manuscript:** *W.E. Berends-Hoekstra, S. Weiland, J.J. Hamstra,I.E. Sommer, T.C. Olde Hartman, J.W.M. Muris, H.M.M. Vos, A. de Jonge, M.H. Blanker, L.L. Peters. More Than the Sum of its Parts: A Dutch Population-Based Study of Women’s Clinical Burden Across the Life Course and Sociodemographic Background, [journal - in preparation], [year - in preparation]*
> **Correspondence:** W.E. Berends-Hoekstra (WEB)

## Repository contents

| File | Purpose |
|------|---------|
| `jaccard_network_community_detection_icpc.R` |  Loads packages and data, derives patient-year denominators, and builds the binary patient × ICPC-code matrices per age group, builds the Jaccard co-occurrence network, runs Louvain community detection, computes diagnostics, and exports node/edge/coordinate tables. |

## Method summary

For each age group:

1. **Binary patient × code matrix** — one row per woman, one column per ICPC code, 1 if the code occurred at least once.
2. **Jaccard similarity** between codes — `1 - vegan::vegdist(t(m), method = "jaccard", binary = TRUE)`; i.e. the number of patients with both codes divided by the number with at least one.
3. **Prevalence filter** — codes present in < 1% of patients are removed (unstable similarities).
4. **Edge threshold** — edges with a Jaccard similarity below the selected threshold are removed (0.075 for ICPC; a single threshold is applied across age groups for comparability).
5. **Community detection** — Louvain (`igraph::cluster_louvain`, `resolution = 1.5`) on the edge-weighted, undirected graph; communities are renumbered by size.
6. **Diagnostics** — number of communities, largest-cluster share, number of singletons, and modularity.

A fixed seed (`set.seed(1301)`) is set for reproducibility of the clustering and the graph layout.

## Requirements

- **R ≥ 4.1** (uses the native pipe `|>` and the `\(x)` lambda syntax).
- Key packages: `tidyverse`, `igraph`, `tidygraph`, `ggraph`, `ggforce`, `vegan`, `viridis`, `writexl`. The full list of `library()` calls is at the top of `00_preamble_load_prepare_icpc.R`.

```r
install.packages(c("tidyverse", "igraph", "tidygraph", "ggraph",
                   "ggforce", "vegan", "viridis", "writexl", "testthat"))
```

## Data availability

The analyses were run inside the **secure remote-access (DRE) environment of Statistics Netherlands (CBS)**. The underlying individual-level GP and register data **cannot be shared** and are not included in this repository. Access to CBS microdata is available to authorised researchers under the conditions of CBS (https://www.cbs.nl/microdata).

The absolute file paths in the scripts (`H:/UMCG/...`) are those of the analysis environment and are **not** accessible externally; adapt them to run the code on your own data.

## Reproducibility notes

- No individual-level data or output files are included.

## How to cite

If you use this code, please cite the manuscript and this repository. 

> Berends-Hoekstra, W.E., *et al.* (2027).

## License
MIT License

---
Primary-and-Long-term-Care-PLC · University Medical Center Groningen (UMCG)
