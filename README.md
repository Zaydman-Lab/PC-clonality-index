# PC-clonality-score
Code repository for PC clonality index package to accompany manuscript by Azimi et al.

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
PC-clonality-score can be used to accomplish 3 main use cases:
1. Using local data to generate local PC2-based reference interval and equation for calculating PC2 metric.
2. Evaluating performance (diagnostic sensitivity and specificity) of locally-derived OR WashU-derived PC2 metric and interval on local data.
3. Calculating PC2 metric and assigning as normal/abnormal using locally-derived OR WashU-derived PC2 metric equation and interval.

The steps for performing each of these use cases is described below.

### Use Case 1: Defining local PC2-based reference interval and equation
To define a local PC2-based reference interval and equation, perform the following steps:
1. Define a local "non-MG" cohort.
2. Create a CSV file with X rows and 2 columns.
- Each row should contain a unique sFLC sample from the non-MG cohort. 
- The left and right columns should contain kappa and lambda sFLC values, respectively. 
- Place the CSV file into the "Data" directory of the repository. 
3. Run pc_clonality.py with the "-n" flag followed by the filename of CSV file created in Step 1. For example:

```python pc_clonality.py -n "your_non_MG_cohort_filename.csv"```
4. Inspect the output in the "Output" folder of the repository
- "pc_vars.csv" will contain the PC2 interval bounds ('PCA_RI_low' and 'PCA_RI_high') as well as all of the variables needed to calculate the PC2 metric for new cases. 
- The equation to calculate the PC2 metric for new cases takes the following form: A*((log(kappa)-kappa_mean)/kappa_std) + B((log(lambda)-lambda_mean)/lambda_std)

### Use Case 2: Evaluating diagnostic performance (sensitivity and specificity) of PC2 metric and interval
The diagnostic performance of the PC2 metric and interval can be evaluated using either 1) a locally-defined PC2 metric equation and reference interval, or 2) the PC2 metric and interval derived using the WashU non-MG cohort. To evaluate diagnostic performance using either of these methods, perform the following steps:
1. Define a local "MG" cohort.
2. Create a CSV file with X rows and 2 columns.
- Each row should contain a unique sFLC sample from the MG cohort. 
- The left and right columns should contain kappa and lambda sFLC values, respectively. 
- Place the CSV file into the "Data" directory of the repository. 
3. To use evaluate diagnostic performance using a locally-defined PC2 metric and interval, run:

```python pc_clonality.py -n "your_non_MG_cohort_filename.csv" -m "your_MG_cohort_filename.csv"```

To use the WashU-derived PC2 metric and interval, simply run the above while omitting the -n flag:

```python pc_clonality.py -m "your_MG_cohort_filename.csv"```

4. Inspect the output in the "Output" folder of the repository:
- "performance.csv" will contain the sensitivity and specificity for the manufacturer's sFLC-ratio-based interval ("sFLC_Sp" and "sFLC_Se") and PC2-based metric and interval ("PC_Sp", "PC_Se"). 

### Use Case 3: Calculating the PC2 metric for new cases
PC-clonality-index provides a function for calculating the PC2 metric for new cases based on their kappa and lambda sFLC results. Similar to use case #2, this use case can be performed using the PC2 metric and interval derived from either 1) a locally-defined PC2 metric equation and reference interval, or 2) the PC2 metric and interval derived using the WashU non-MG cohort. To calculate the PC2 metric for new cases using either of  these methods, perform the following steps:
1. Create a CSV file with X rows and 2 columns.
- Each row should contain a unique sFLC sample. 
- The left and right columns should contain kappa and lambda sFLC values, respectively. 
- Place the CSV file into the "Data" directory of the repository. 
2. To calculate the PC2 metric using a locally-defined PC2 metric and interval, run:
```python pc_clonality.py -n "your_non_MG_cohort_filename.csv" -c "your_cases_filename.csv"```

To use the WashU-derived PC2 metric and interval, simply run the above while omitting the -n flag:

```python pc_clonality.py -m "your_MG_cohort_filename.csv"```

3. Inspect the output in the "Output" folder of the repository:
- "cases_pc2.csv" will contain four columns containing the kappa values, lambda values, PC2 metric, and an abnormal flag (0 for normal, 1 for abnormal) for each case, respectively. 

### Combining use cases
Use cases 1, 2, and 3 can be combined into one step. For example, to define a local PC2-based interval, evaluate the diagnostic performance, and calculate the PC2 metric for new cases, run the following:

```python pc_clonality.py -n "your_non_MG_cohort_filename.csv" -m "your_MG_cohort_filename.csv" -c "your_cases_filename.csv"```

### Defining custom interval bounds
By default, PC-clonality-index will calculate the PC2-metric-based reference interval using a 95% (2.5-97.5%ile) diagnostic interval. These bounds can be customized by providing "-l" and -u" flags. For example, to define a PC2 reference interval using a 90% diagnostic interval, run the following:

```python pc_clonality.py -n "your_non_MG_cohort_filename.csv" -l 5 -u 95```

## Contact
For any questions, comments, or suggestions, please feel free to contact [Mark Zaydman](mailto:zaydmanm@wustl.edu) or [Vahid Azimi](mailto:a.vahid@wustl.edu).
