# Our task

**Analysis of T-cells immunological landscape and its association with graft-versus-host disease**

Mentor: Mikhail Drokov, MD

Transplantation of allogeneic hematopoietic stem cells allows achieving a biological cure in patients with leukemia. But the main reason for the decrease of disease-free survival and quality of life is chronic graft-versus-host disease (chronic GVHD), which develops in more than half of patients. Chronic GVHD remains a clinical diagnosis but attempts are being made worldwide to find diagnostic markers of this complication. The pathogenesis of this condition is not fully understood, but the role of T-lymphocytes has been proven.

A large dataset of immunological data is presented (analysis of the subpopulation composition of T-lymphocytes - 162 subpopulations) for 70 patients after stem cells transplantation. Chronic GVHD developed in 28 of these 70 patients.

**Purpose**: application of modern data processing and visualization technologies for the prediction and characterization of chronic GVHD.

[More](https://bioinf.institute/hack/teams "BioHack Team page")

1. Identify landscape of T cells and their combination in patients after allo-HSCT.
   
2. Using biostatistics and machine learning methods identify the relationship between T cells and their combination in patients after allo-HSCT and development of cGVHD.
   
3. Visualize using t-SNE, PCA and etc.  T cells populations in patients who develop cGVHD and who not. 

## Presentation

The project was presented at [BioHack Hackathon 2023](https://bioinf.institute/hack/teams "BioHack_2023")

Presentation video: [YouTube](https://www.youtube.com/live/E-VtzpypL6M?feature=share&t=5301)

[PDF of the slides](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/T_cells_cGVHD_Hackathon_2023_final.pdf "BioHack_2023")

_**The project won 1st Prize.**_


# Dataset preprocessing

We performed inital dataset preprocessing and worked with blood samples collected on day 90 for patients who did not develop cGVHD yet. The metrics for cell population was concentration of identified T cell population concentration per ml of blood. Each population had information about its maturity and activity marker. Overall, there were 66 patients. We received only aggregated knowledge about dataset clinical characteristics except day and status of the onset. Mostly, patients were women (60%) and had either AML (56%) or cll (42.6%). 91% of patients received reduced-intensity conditioning, and 80% received similar prophylaxy. The last patient dat follow-up was considered the day 275 since the blood collection.

## Potential duplicates 

We have dropped potential duplicates for which 90% of the data or higher had Spearman rho values different by 0.05 or lower. 8 pairs different by markers of maturity or activity, but not by the identified cell type, were found. We have dropped the populations that had less information about the population.
Here are the potentially duplicated cell types: 
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/potential_dupes_tight.png "heatmap")

## Median deviation cutoff

We used the analog of dispersion from the median-centered statistics paradigm: mad of the cell population per sample median + 1. For day 90 and day 180 of blood sampling, we have plotted the median deviation values and identified samples with the highest and lowest variability. After the intersection, 24 cell populations were discarded.

Median deviation for day_90 (cutoffs <0.5; >4):

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/median_deviation_day_90.png "median_deviation_day_90")

Median deviation for day_180 (cutoffs <0.5; >1.6):

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/median_deviation_day_180.png "median_deviation_day_180")


## Log-scaling

We have identified outliers in raw non-scaled data and thus log2-scaled the data.

Before:
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/eda_distplot_non_scaled.png "eda_distplot_non_scaled")

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/eda_pca_plot_non_scaled.png "pca_distplot_non_scaled")

After:
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/eda_distplot_log_scaled.png "eda_distplot_log_scaled")

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/eda_pca_plot_log_scaled.png "pca_distplot_log_scaled")

## Feature engineering

To work with the important cGVHD-predicting populations, we have created new features using ssGSEA score approach, common to gene expression analysis. Rank-based approach helped us to overcome the issue with the data: for same cell type and maturity state one does not know if detected markers intersect or not (as we've established with duplicate detection). Rank sums helps us to decrease each value contribution to the value. 
We calculated ssGSEA score for all "active" state T cells (HLA-DR+, CD226+, TIGIT-), "suppressed" (TIGIT+, PD1+, CD35+, CD226-), all CD8+ cells and CD4+ cells.

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/ssgsea.jpg "ssgsea")

# Time-to-event analysis

One of the approaches for the task was to perform time-to-event approach.

## KME
First, we were searching for dependencies using Kaplan Meier approach. Using log-rank test we have identified populations with a predictive tendency, for example, CD8+ HLA-DR+ cells and Th1 TM were associated with reduced risk. We used 3 cutoffs (Q1, Q2, Q3) and FDR was applied inside every cutoff. The majority of populations (out of 132 cell types) were not influencing the survival. 

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/km1.png "km1")

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/km2.png "km2")

## CPH

### Assumptions
We have also used the Cox proportional hazards regression univariate models. We checked the CPH assumptions. If the check failed, the population was categorized and the model was checked again. If the model failed the Shoenfild residuals tests, such populations were dropped. 

### Permutations tests was baselines
For each feature, 10 permutations rounds were performed and the mean permutated value coeffiecient, as well as Harrel's C-index, indicating model predicting power were tracked.

### CPH results
Here are the populations that had a significant log-rank test results and for which C-index was higher than permutated value C-index:

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/forest_plot.png "forest_plot")

Higher CD4+ TE maturity PD1- TIGIT+ levels associated with associated with cGVHD onset. Higher levels of CD8+ cells with different stages of maturity associated with with lack of cGVHD. Also, several engineered features had a higher association with the target.

# Clustering

## NMF clustering

Using sparse NMF approach, we have clustered populations and patients with 2 clusters. It was the best metric when the matrix prediction was tested 10 times using cophenetic correlation and dispersion as metrics of interests for number of clusters under 5 (considering low amount of patients).
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/nmf_ranks_40.png "nmf_ranks")

The resulting clusters were as follows (100 iterations):
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/nmf_select_clustering_heatmap.png "heatmap")

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/nmf_select_clustering_heatmap_event_patch_plot.png "heatmap")

The populations that were plotted were selected by feature fraction (max/sum by row) derived from the resulting matrices. Distribution of cGVHD events was significant (Fisher test, FDR=0.045 for each cluster).

## Leiden clustering

Another approach was to use graph-based clustering. 
![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/leiden.jpg "leiden")

## Clustering Results

Both clusterings identified a predominantly CD4 T cell naive cluster with worse outcome (cGVHD) and a CD8 T cell mature cluster with a better outcome (no cGVHD).

# Classification approach
To detect non-linear dependencies, we wanted to try other ML approaches.

## Baseline
We have selected features using 500 iteration of BorutaShap. To assess baseline, we have used pycaret tool and identified the scores. Low Dummy classifier scores says that the task is solvable. 

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/pycaret.jpg "pycaret")

## Classifier building
We have tried to make a classifier (XGBoost and RF) and optimize it using optuna (1000 rounds). While training scores were high, models faild the validation dataset estimation. Possibly due to very low amount of samples.

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/ROC_PR_curves.png "ROC_PR_curves")

## Shap
Modelling help us to assess possible interesting non-linear dependencies between features:

![alt text](https://github.com/onion-42/cGVHD_T_cell_populations_BioHackathon_2023/blob/main/plots/shap.png "shap")

# Results

CD8+ T-cells (any maturity, predominantly active) are associated with better outcomes after allo-HSCT transplantation.
CD4+ T-cells (low effector function) are associated with chronic GVHD onset.
