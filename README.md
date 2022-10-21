# PC-clonality-score
Code repository for PC clonality index package to accompany manuscript by Azimi et al

## Introduction

The serum free light chain (sFLC) assay is a screening test for monoclonal gammopathy. The interpretation of this assay is done by taking a ratio of kappa to lambda free light chains (sFLC-ratio); however, interpretation of this calculation is confounded by renal impairment, a condition that skews the sFLC-ratio towards kappa. To this end, we used principal component analysis (PCA) to separate "normal" variation due to renal function (defined by principal component 1) from "abnormal" variation due to clonality (defined by principal component 2), and used deviation from principal component 1 to define a new metric and reference interval for sFLC interpretation. This package allows users to apply and evaluate the PCA-based approach for sFLC interpretation using their own local data. 

## Getting Started
###Dependencies
- Python 3.10+
- Optparse
- Matplotlib
- Seaborn
- NumPy
- Pandas
- SciPy

## Installation
To run, clone the PC-clonality-score Github repository
'''git clone https://github.com/Zaydman-Lab/PC-clonality-score.git'''
