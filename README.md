# Correlated evolution of thermal tolerance and egg size links resistance and recovery to extreme heat in amphibians



## Description

This repository contains the metadata, raw data, phylogenetic information, and R scripts necessary to reproduce the analyses and figures for the manuscript: "Pascucci, G., Tejedo, M., Gvoždík, L. Correlated evolution of thermal tolerance and egg size links resistance and recovery to extreme heat in amphibians".

This material can be freely shared, adapted, or reprocessed (including for commercial use) as long as appropriate credit is given to the authors.
License details: [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

**You can cite this repository as:**
Pascucci, G. (2026). Data from: Correlated evolution of thermal tolerance and egg size links resistance and recovery to extreme heat in amphibians (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.23039857

**Repository contact:** Gioele Pascucci ([gioelepascucci@live.it](mailto:gioelepascucci@live.it))

## Data collection & processing summary

Primary physiological and life-history traits data were compiled across 83 amphibian species from published literature and established database sources. Thermal tolerance values were standardized for acclimation temperature, duration, and physiological endpoints using `heatTolerance.R`. Spatial range maps and occurrence records were retrieved from GBIF (gbif.org, 2026) and IUCN (2025), then integrated with elevation raster from WorldClim 2.1 (Fick & Hijmans, 2017), using `speciesRange.R` to compute elevation ranges. All statistical analyses were conducted within a phylogenetically informed framework (e.g, Phylogenetic Generalized Least Squares) in R v4.5.3 (R Core Team, 2026; complete package version details are documented in `sessionR_info.txt`). For full details on variable definitions and methodologies, refer to the main article and the metadata sheets within each data file provided in this repository.

---

## Directory structure

```text

├── sessionR_info.txt       : R version and package versions used for the analysis.
                              

/data
├── datasetTraits.xlsx      : Compiled dataset containing the 83 species of | Sheet-1: "Trait_FINAL"
│                             amphibians, including raw trait values,       | Sheet-2: "Raw_values"
│                             references and metadata.                      | Sheet-3: "Metadata"
│   
├── Jetz&Pyron2018.nex      : A subset of the amphibian phylogeny from      
│                             Jetz & Pyron (2018).                          
│   
└── literatureAnalysis.xlsx : Compiled results from the literature          | Sheet-1: "Literature_RESULT"
                              analysis.                                     | Sheet-2: "Metadata"
                                                                              

/figures
└── Output directory intended for the figures generated for the main 
    manuscript and the Supplementary Information (Text S1).

/script
├── heatTolerance.R : Standardizes CTmax values from different studies to 
│                     ensure consistent acclimation temperature, acclimation 
│                     time, and physiological endpoints.
│   
├── speciesRange.R  : Integrates spatial data from multiple sources (IUCN 
│                     range maps, GBIF occurrences, and WorldClim Digital 
│                     Elevation Models) to calculate species distribution.
│   
└── main.R          : Reproduce all primary analyses of the study, 
                      including PGLS models, plot generation, literature 
                      analysis summaries, diagnostic tests, and sensitivity 
                      analyses.

```
> **Note:** Additional external data must be downloaded into the `/data` folder before running the scripts (see instructions below).

---

## Setup and analysis workflow

To successfully replicate the study, please follow the steps in the exact order below.

### Step 1: Download required external data

Before running any scripts, you must download the following external datasets and place them directly into the `/data` directory.
<p></p>

1. **AmphiTherm dataset** (Pottier et al., 2025)
   * Download the version v1.0.1 `.zip` file from [Zenodo](https://doi.org/10.5281/zenodo.16751427).
<p></p>

2. **IUCN Amphibian range maps** (IUCN, 2025; accessed August 2026)
   * Download the shapefile `.zip` for all amphibians from the [IUCN Red List](https://www.iucnredlist.org/resources/spatial-data-download#amphibians).
<p></p>

3. **WorldClim elevation raster** (Fick & Hijmans, 2017)
   * Download the version 2.1, elevation data (30-seconds resolution) from [WorldClim](https://www.worldclim.org/data/worldclim21.html).
<p></p>

### Step 2: Execute scripts

Once all data files are located in the `/data` folder, run the R scripts in the following sequence:
<p></p>

1. **`heatTolerance.R`**
   * **Action:** Standardizes the thermal tolerance data.
<p></p>

2. **`speciesRange.R`**
   * **Action:** Calculates the spatial distribution and altitudinal ranges for the species.
<p></p>

3. **`main.R`**
   * **Action:** Executes the core phylogenetic and statistical analyses and generates the final manuscript figures.

---
