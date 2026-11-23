# Correlated evolution of thermal tolerance and egg size links resistance and recovery to extreme heat in amphibians



## Description

This repository contains the raw data, phylogenetic information, and R scripts necessary to reproduce the analyses and figures of the manuscript "Correlated evolution of thermal tolerance and egg size links resistance and recovery to extreme heat in amphibians" (Pascucci, G., Tejedo, M., Gvoždík, L.).

This material can be freely shared, adapted, or reprocessed (including for commercial use) as long as appropriate credit is given to the authors.
License details: [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

**You can cite this repository as:**
Pascucci, G., Tejedo, M., & Gvoždík, L. (2026). Data from: Correlated evolution of thermal tolerance and egg size links resistance and recovery to extreme heat in amphibians (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.23039857

**Repository contact:** Gioele Pascucci ([gioelepascucci@live.it](mailto:gioelepascucci@live.it))

## Abstract
We tested 83 species of amphibians (Anura, 43; Caudata, 40), representing 19 familes. Using phylogenetic comparative analyses, across multiple evolutionary models, we investigated correlations between critical thermal maximum (CTmax) and three life-history traits, clutch size, egg size, and age at maturity, while controlling for geographic gradients and body size (snout-vent length).

## Data collection & processing

Primary physiological and life-history traits data were compiled from published literature and established databases.
Thermal tolerance values were standardized for acclimation temperature, duration, and physiological endpoints using `heatTolerance.R`.
Species range maps and occurrence records were retrieved from IUCN (2025) and GBIF (gbif.org, 2026), then integrated with elevation raster from WorldClim 2.1 (Fick & Hijmans, 2017), using `speciesRange.R`.
Phylogenetic comparative analyses (e.g., PGLS and multivariate PGLS) based on the Jetz and Pyron (2018) phylogeny were implemented in R v4.5.3 (R Core Team, 2026). The complete analyses workflow and R package environment are provided in `main.R` and `sessionR_info.txt`, respectively.

## Directory structure

```text

├── sessionR_info.txt       : R version and package versions used for the analyses.
                              

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
                      including PGLS models, plot generation, diagnostic tests, and sensitivity 
                      analyses.

```
> **Note:** Additional external data must be downloaded into the `/data` folder before running the scripts (see instructions below).

---

## Setup and analyses workflow

To successfully replicate the study, please follow the steps in the exact order below.

### Step 1: Download required external data

Before running any scripts, you must download the following external datasets and place them directly into the `/data` directory (keep .zip files compressed).
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
   * **Action:** Standardize the thermal tolerance data.
<p></p>

2. **`speciesRange.R`**
   * **Action:** Calculate the spatial distribution and elevation ranges for the species.
<p></p>

3. **`main.R`**
   * **Action:** Execute the core phylogenetic and statistical analyses and produce the manuscript figures.

---