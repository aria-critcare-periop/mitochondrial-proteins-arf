# Mitochondrial proteins and inflammatory phenotypes in acute respiratory failure

Analysis code for the study **Plasma mitochondrial-associated proteins overlap with the hyperinflammatory phenotype in acute respiratory failure**.

This repository contains six R scripts and five Mplus input files.

## Data and software

Download the Olink NPX and clinical data from [Zenodo](https://doi.org/10.5281/zenodo.7076472), and the MitoCarta3.0 inventory. Place these files in a `data/` folder:

- `MGH_COVID_OLINK_NPX.txt`
- `MGH_COVID_Clinical_Info.txt`
- `Human.MitoCarta3.0.xls`

The analyses require R, Mplus, and the R packages used in the scripts.

## Run order

Run the files from the repository folder:

1. Run `01_prepare_cohort.R`. This creates `mgh_lca_input.dat`.
2. Run `lca_k1.inp` through `lca_k5.inp` in Mplus. Keep the resulting `.out` files and `class_k2.dat` in the repository folder.
3. Run `02_inflammatory_phenotype.R`.
4. Run `03_mitochondrial_cluster.R`.
5. Run `04_differential_abundance.R`, `05_figure1.R`, and `06_cluster_diagnostics.R`.
