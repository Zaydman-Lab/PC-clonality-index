# PC-clonality-score
Code repository for PC clonality index package to accompany manuscript by Azimi et al

## Introduction

The serum free light chain (sFLC) assay is a screening test for monoclonal gammopathy. The interpretation of this assay is done by taking a ratio of kappa to lambda free light chains (sFLC-ratio); however, interpretation of this calculation is confounded by renal impairment, a condition that skews the sFLC-ratio towards kappa. To this end, we used principal component analysis (PCA) to separate "normal" variation due to renal function (defined by principal component 1) from "abnormal" variation due to clonality (defined by principal component 2), and used deviation from principal component 1 to define a new metric and reference interval for sFLC interpretation. This package allows users to apply and evaluate the PCA-based approach for sFLC interpretation using their own local data. 

## Getting Started
### Dependencies
- Python 3.10+
- Optparse
- Matplotlib
- Seaborn
- NumPy
- Pandas
- SciPy

### Installation
To run, clone the PC-clonality-score Github repository:

```git clone https://github.com/Zaydman-Lab/PC-clonality-score.git```

## Usage
PC-clonality-score can be used to accomplish 3 main functions:
1. Using local data to generate local PC2-based reference interval and equation for calculating PC2 metric.
2. Evaluating performance (diagnostic sensitivity and specificity) of locally-derived OR WashU-derived PC2 metric and interval on local data.
3. Calculating PC2 metric and assigning as normal/abnormal using locally-derived OR WashU-derived PC2 metric equation and interval.

### Function 1: Defining local PC2-based reference interval and equation
To define a local PC2-based reference interval and equation, perform the following steps using a local non-MG cohort:
1. Create an CSV file with X rows and 2 columns.
- Each row should contain a unique sFLC sample. 
- The left and right columns should contain kappa and lambda sFLC values, respectively. 
- Place the CSV file into the "Data" directory of the repository. 
2. Run pc_clonality.py with the "-n" flag followed by the filename of CSV file created in Step 1. For example:

```python pc_clonality.py -n "your_non_MG_cohort_filename.csv"```
3. Check output in the "Output" folder of the repository
- "pc_vars.csv" will contain the PC2 interval bounds ('PCA_RI_low' and 'PCA_RI_high') as well as all of the variables needed to calculate the PC2 metric for new cases. 
- The equation to calculate the PC2 metric for new cases takes the following form: A*((log(kappa)-kappa_mean)/kappa_std) + B((log(lambda)-lambda_mean)/lambda_std)

###
