# Our task

**Analysis of T-cells immunological landscape and its association with graft-versus-host disease**

Transplantation of allogeneic hematopoietic stem cells allows achieving a biological cure in patients with leukemia. But the main reason for the decrease of disease-free survival and quality of life is chronic graft-versus-host disease (chronic GVHD), which develops in more than half of patients. Chronic GVHD remains a clinical diagnosis but attempts are being made worldwide to find diagnostic markers of this complication. The pathogenesis of this condition is not fully understood, but the role of T-lymphocytes has been proven.

A large dataset of immunological data is presented (analysis of the subpopulation composition of T-lymphocytes - 162 subpopulations) for 70 patients after stem cells transplantation. Chronic GVHD developed in 28 of these 70 patients.

**Purpose**: application of modern data processing and visualization technologies for the prediction and characterization of chronic GVHD.

[More](https://https://bioinf.institute/hack/teams "BioHack Team page")

Task 1: Identify landscape of T cells and their combination in patients after allo-HSCT.
Task 2: Using biostatistics and machine learning methods identify the relationship between T cells and their combination in patients after allo-HSCT and development of cGVHD.
Task 3: Visualize using t-SNE, PCA and etc.  T cells populations in patients who develop cGVHD and who not. 

# Dataset preprocessing

## Potential duplicates 

We have dropped potential duplicates for which 90% of the data or higher had Spearman rho values different by 0.05 or lower. 8 pairs different by markers of maturity or activity, but not by the identified cell type, were found. We have dropped the populations that had less information about the population.
Here are the potentially duplicated cell types: 
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/potential_dupes_tight.png "heatmap")

## Median deviation cutoff

We used the analog of dispersion from the median-centered statistics paradigm: mad of the cell population per sample median + 1. For day 90 and day 180 of blood sampling, we have plotted the median deviation values and identified samples with the highest and lowest variability. After the intersection, 24 cell populations were discarded.

Median deviation for day_90:

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/median_deviation_day_90.png.png "median_deviation_day_90.png")

Median deviation for day_180:

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/median_deviation_day_180.png.png "median_deviation_day_180.png")
