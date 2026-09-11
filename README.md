# Analysis code for: "Same visitors, different outcomes: floral phenotype gates pollinator-vectored nectar microbial establishment"

Analysis code for a two-year field exclusion experiment testing how diurnal
and nocturnal pollinator guilds shape nectar microbiome assembly in two
sunflower cultivars (Celesto, Idillic) contrasting in nectar volume and
sugar composition.

## Data

Raw sequencing data: ENA study [PRJEB121466](https://www.ebi.ac.uk/ena/browser/view/PRJEB121466)
Nectar traits, visitation counts, culturability scores, reference sequences: [Recherche Data Gouv, doi:10.57745/LGKSOT](https://doi.org/10.57745/LGKSOT)

Scripts expect input files under a local `data/` folder and write outputs to
`output/`; neither is included in this repository. File names and expected
columns are documented in a header comment at the top of each script.

## Scripts

| Script | Produces | Manuscript figures/tables |
|---|---|---|
| `metschnikowia_tree_merge.R` | Phylogeny-informed OTU merging for *Metschnikowia* | Figure S1 |
| `merge_metschnikowia_otus.py` | Merged feature table / taxonomy / sequences from the R output | — |
| `its_microbiome_analysis.R` | ITS import, filtering, alpha/beta diversity, camera-trap regressions | Figure 4A–C, S4–S6 |
| `bacterial_16S.R` | 16S import, filtering, composition, beta diversity | Figure 3, S3 |
| `voc_analysis.R` | Floral VOC analysis | Figure 5, S7, Table S8 |
| `nectar_traits.R` | Nectar volume, sugar mass, sugar composition | Figure 1B |
| `pollinator_visits.R` | Camera-trap visit rates, day/night windows | Figure 2, Table S2 |
| `culturability.R` | Culturable microbial load, growth scores | Figure 4D, Table S5 |

Execution order: `metschnikowia_tree_merge.R` → `merge_metschnikowia_otus.py`
→ `its_microbiome_analysis.R`. The other scripts are independent of each
other and of this chain, except `its_microbiome_analysis.R`, which also
reads the output of `pollinator_visits.R` (`plant_means.csv`) for the
camera-trap section.

## Upstream steps not included here

- OTU clustering (VSEARCH, 97%), chimera removal (UChime), curation (mumu)
  and taxonomic assignment (QIIME2): custom Python scripts at
  http://forge.inrae.fr/aurelien.carlier/attracthol
- Sequence alignment (MAFFT v7.505), trimming (Trimal, >95% gaps) and
  phylogeny inference (FastTree2, `-gtr -gamma`) upstream of
  `metschnikowia_tree_merge.R`
- Pollinator detection/classification model (YOLO11x): forge.inrae.fr/astr/public/pollicrop
- GC-MS raw data processing (MetAlign, MSClust)
- Acinetobacter strain identification: EzBioCloud 16S-based ID (web service)

## Software

R v4.5.1. Key packages: phyloseq v1.52.0, vegan v2.7.1, iNEXT v3.0.1,
pairwiseAdonis v0.4.1, ape v5.8.1, phytools v2.5.2, ggtree v3.16.3,
car v3.1.3, emmeans v2.0.0, multcomp v1.4.29, ggplot2 v4.0.3, ggh4x v0.3.1,
patchwork v1.3.2, pheatmap v1.0.13. Full list and versions in each script's
header / the manuscript's Methods.

## Citation

[DOI / citation once published]
