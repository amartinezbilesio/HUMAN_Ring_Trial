# HUMAN Ring Trial

Reproducible R/Quarto workflow for processing mass spectrometry data from the HUMAN Ring Trial and generating spectral libraries and `compDb` objects.

## Overview

This repository contains the code and metadata required to reproduce the spectral library generation workflow used for the HUMAN Ring Trial.

The workflow is based on the RforMassSpectrometry ecosystem and is designed to operate on mzML input files provided separately from this repository.

## Repository structure

```text
HUMAN_Ring_Trial/
├── 1_preprocessing/
│   ├── afekta/
│   │      ├── HE/
│   │      └── Others/
│   ├── cembio/HE/
│   ├── hmgu/HE/
│   ├── icl/HE/
│   ├── NAPS_info.xlsx
│   ├── generic_preprocessing.qmd
│   ├── setup.R
│   └── standards.xlsx
│
├── 2_annotation/sirius_all_ms2/
│
├── 3_library_generation/
│   ├── build_spectral_library_provenance.R
│   └── spectral_library_provenance.xlsx
│
├── .gitignore
├── .Rprofile
├── DESCRIPTION.md
├── README.md
├── README_data.md
├── README_workflow.md
└── renv.lock
```

## Requirements

The workflow requires:

* R
* RStudio or another R environment capable of running Quarto
* Quarto
* the R/Bioconductor packages specified by the project
* access to the required HUMAN Ring Trial mzML files

## Installation

Clone this repository and open the `.Rproj` file.

Restore the R environment with:

```r
renv::restore()
```

This installs the package versions recorded in `renv.lock`.

## Input data

The mzML files are not stored in this GitHub repository.

Users must obtain the required mzML data separately.

See:

```text
README_data.md
```

for the expected data structure and information about data availability.

## Configure the data location

Open the Quarto document ring_trial_library_building.qmd and set:

```r
data_dir <- "PATH/TO/HUMAN_RING_TRIAL_DATA"
```

For example:

```r
data_dir <- "D:/HUMAN_Ring_Trial_Data"
```

## Provenance

The expected mzML files and their associated metadata are documented in:

```text
3_library_generation/spectral_library_provenance.xlsx
```

The provenance table contains the information required by the workflow to associate each sample with its laboratory, file name, and polarity.

## Running the workflow

After configuring the input data location, open the relevant Quarto document and render it.

The workflow will:

1. read the provenance metadata;
2. locate the corresponding mzML files;
3. process the mass spectrometry data;
4. generate the required spectral-library objects;
5. generate the corresponding `compDb` output.

The exact outputs depend on the workflow version and configuration.

## Reproducibility

Package versions are recorded using `renv`.

To recreate the project environment:

```r
renv::restore()
```

The repository intentionally separates code and metadata from the large mzML input files.

## Data availability

See:

```text
README_data.md
```

## License

[MIT License]
