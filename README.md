# PC-clonality-index
Code repository for PC clonality index package to accompany manuscript by Azimi et al.

## Introduction

This package can be used to perform three use cases: 1) Derive an equation for calculating the PC2 clonality index and a corresponding reference interval using serum free light chain assay results from a non-MG cohort, as described in Azimi et al. (reference), 2) Calculate the PC2 clonality index for a set of serum free light chain results, and 3) Calculate the sensitivity and specificity of the PC2 clonality index-based reference interval and the manufacturer's sFLC-ratio-based interval using non-MG and MG cohorts. 

The behavior of the program will be determined by the files and input flags you give it. We have provided a template file called "WashU.p"; you will use this file if you want to use your own local reference cohort for your own institutional data. The default behavior if "non_mg.csv" is not provided is to use the WashU data that is described in Azimi et al. (1). The non-MG cohort will ideally consist of an sFLC result from at least 120 patients (2) that do not have monoclonal gammopathy and have varying degrees of renal function. "MG.csv" will be used if you want to compare the sensitivity of the PC2 clonality index reference interval with the manufacturer's sFLC-ratio-based reference interval (specificity will be defined by the interval percentage specified [default 95%]). For additional details on non-MG and MG cohort definitions, please refer to Azimi et al. (reference). "cases.csv" is a file you can upload to calculate the PC clonality index for any given set of serum free light chain results. 

<img width="867" alt="Screen Shot 2022-11-30 at 4 40 29 PM" src="https://user-images.githubusercontent.com/10732984/204924499-e5ffab22-ebb5-4664-aa74-0aefa95b1484.png">


## Getting Started
### Installation
To run, clone the PC-clonality-score Github repository to a local directory:

```git clone https://github.com/Zaydman-Lab/PC-clonality-score.git```
### Environment and Dependencies
1. [Download Conda] (https://conda.io/projects/conda/en/latest/user-guide/install/download.html) and follow [installation instructions] (https://conda.io/projects/conda/en/latest/user-guide/install/index.html#).
2. Create the conda environment by running the .yml file:
```conda env create -f environment.yml```
This will create the Conda environment and install the following dependencies:
- Python 3
- Optparse
- Matplotlib
- Seaborn
- NumPy
- Pandas
- SciPy

3. Activate the Conda environment:

```conda activate sflc```

4. Run the code (as described below).

## Usage
### Case 1: Defining local PC2-based reference interval and equation
#### Inputs
- Optional: non_mg.csv
#### Usage
```python pc_clonality.py -n path_to_non_mg.csv```
#### Outputs
- cases_1.png: .png image of plot of non-MG cohort with manufacturer's sFLC-ratio-based and PC2-based reference intervals superimposed
- vars.csv: contains the PC2 clonality index-based interval ('PCA_RI_low' and 'PCA_RI_high') as well as the variables needed to calculate the PC2 clonality index for new cases. The equation takes the form of:

```math
A*(\frac{ln(kappa)-{B}}{C}) + D*(\frac{ln(lambda)-{E}}{F})
```

Where 'ln' represents the natural logarithm.

Description of variables (WashU-specific values)

A: (0.71)

B: mean of log(kappa) values in non-MG cohort (1.31)

C: standard deviation of log(kappa) values in non-MG cohort (0.87)

D: (-0.71)

E: mean of log(lambda) values in non-MG cohort (0.88)

F: standard deviation of log(lambda) values in non-MG cohort (0.68)

The WashU-specific PCA_RI_low and PCA_RI_high are -0.542 and 0.584, respectively.

### Case 2: Calculating PC clonality index for new cases
#### Inputs
- Required: cases.csv
- Optional: non_MG.csv
#### Usage
```python pc_clonality.py -n path_to_non_mg.csv -c path_to_cases.csv```
#### Output
- pc2_cases.csv: contains four columns containing the kappa values, lambda values, PC2 clonality index, and an abnormal flag (0 for normal, 1 for abnormal) for each case.
- cases_2.png: .png image of plot of non-MG cohort and MG cohort with manufacturer's sFLC-ratio-based and PC2-based reference intervals superimposed

### Case 3: Evaluating Sensitivity and Specificity
#### Inputs
- Required: mg.csv
- Optional: non_mg.csv
#### Usage
```python pc_clonality.py -n path_to_non_mg.csv -m path_to_mg.csv```
#### Output
- "performance.csv" will contain the sensitivity and specificity for the manufacturer's sFLC-ratio-based interval ("sFLC_Sp" and "sFLC_Se") and PC2 clonality index-based reference interval ("PC_Sp", "PC_Se"). 
- cases_3.png: .png image of plot of non-MG cohort and cases with manufacturer's sFLC-ratio-based and PC2-based reference intervals superimposed

### Defining custom interval bounds
By default, PC-clonality-index will calculate the PC2-metric-based reference interval using a 95% (2.5-97.5%ile) diagnostic interval. These bounds can be customized by providing "-l" and "-u" flags. For example, to define the PC2 metric equation and a corresponding 90% reference interval, run the following:

```python pc_clonality.py -n path_to_non_mg.csv -l 5 -u 95```

### Combining use cases
Use cases 1, 2, and 3 can be combined into a single step. For example, to define a local PC2-based interval, evaluate the diagnostic performance, and calculate the PC2 metric for new cases, run the following:

```python pc_clonality.py -n path_to_non_mg.csv -m path_to_mg.csv -c path_to_cases.csv```

## References
1. Azimi V, Slade M, Fiala M, Fortier J, Stockeryl-Goldstein K, Frater JL, Brestoff JR, Jackups R, Zaydman MA. A novel metric improves the accuracy and robustness to renal function of serum free light chain assay interpretation. 
2. Horowitz, GL. Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory, 3rd Edition. Clinical and Laboratory Standards Institute, 2008.

## Contact
For any questions, comments, or suggestions, please feel free to contact Mark Zaydman via this repository or [email](mailto:zaydmanm@wustl.edu) or Vahid Azimi via [email](mailto:a.vahid@wustl.edu).
