# HUMAN Ring Trial

## Input data

The HUMAN Ring Trial mzML files are not included in this GitHub repository because of their size and data-distribution considerations.

## Required data

To run the workflow, users must have access to the mzML files used by the analysis.

The location of the top-level data directory must be specified locally in:

```text
data_dir <- "PATH/TO/METASCI_DATA"
```

For example:

```r
data_dir <- "D:/HUMAN_Ring_Trial_Data"
```

The expected directory structure is:

```text
data_dir/
├── AFEKTA/
│
├── CEMBIO/
│
├── HMGU/
│
├── ICL/
│
└── AUTH/

```

## File provenance

The expected mzML filenames and associated metadata are defined in:

```text
project_dir <- "PATH/TO/PROJECT_DIR"

library_dir <-
    file.path(project_dir,
              "3_library_generation")

provenance_file <-
    file.path(
        library_dir,
        "spectral_library_provenance.xlsx"
    )
```

The workflow uses this provenance table to locate the corresponding mzML files and associate them with the relevant laboratory and polarity.

## Data availability

If the mzML files are publicly available, the corresponding repository, accession numbers, or download instructions should be provided here.

If access to the data is restricted, users should obtain access from the data provider before attempting to reproduce the workflow.

```text
LAB             MTBLS_ID              MTBLS_assay_name              MTBLS_public_status             MTBLS_access_link
AFEKTA          MTBLS14343		                                      TRUE    	                      https://www.ebi.ac.uk/metabolights/MTBLS14343
CEMBIO          REQ20251003213571		                                FALSE
HMGU            MTBLS15694		                                      FALSE	                          https://www.ebi.ac.uk/metabolights/MTBLS15694
ICL             MTBLS14274		                                      FALSE	                          https://www.ebi.ac.uk/metabolights/MTBLS14274
AUTH		                                                            FALSE	
```
