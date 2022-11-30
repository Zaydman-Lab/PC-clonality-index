"""
pc_clonality.py
Created by Vahid Azimi
11/29/22

Purpose: main point of entry into PC-clonality-index package

This script will execute one of the three cases described in the .README based on the user input

Functions:


generate_embedding(nonmg_path,lb,ub,output_path):
    '''Returns tranform object and embeddings'''
"""

#%% imports
from optparse import OptionParser
import pandas as pd 
import matplotlib.pyplot as plt
import pickle
from typing import Tuple, Callable
import visualize
import transforms
import performance
import numpy as np


#%% helper functions
def df2array(df: pd.DataFrame)->np.array:
	"""Returns np array of kappa and labmda values derived from input dataframe"""
	X = np.column_stack((df['kappa'].to_numpy(), df['lambda'].to_numpy())) 
	return(X)
#%%
def derive_interval(nonmg: pd.DataFrame,lb: float,ub: float)->Tuple[np.array,list[float,float],Callable[[np.array,float],np.array],Callable[[np.array,float],np.array]]:
	"""Derive PC2-based refence interval using WU nonmg cohort or user input nonmg cohort"""
	X_nonmg = df2array(nonmg)
	L_nonmg = transforms.log_transform(X_nonmg) #apply log transform to X_nonmg
	z_transform = transforms.make_ztransform(L_nonmg) #compute z transform of L_nonmg
	Z_nonmg = z_transform(L_nonmg) #apply z transform to L_nonmg
	pc_transform, U, S, Vh = transforms.make_pctransform(Z_nonmg,return_decomposition=True) #compute pc transform of Z_nonmg
	pc2 = pc_transform(Z_nonmg)[:,1] #pc2 projections for nonmg cohort
	pc2_RI = [np.percentile(pc2,lb),np.percentile(pc2,ub)] #compute reference interval for pc2
	Vh_raw=transforms.log_transform(z_transform(Vh,inverse=True),inverse=True) #project right singular vectors from z space into raw sFLC space
	equation_parameters = {
		'A': Vh_raw[0,0], #First right singular vector in raw space
		'B': np.mean(transforms.log_transform(X_nonmg)[:,0]), #Mean of log transformed kappa values for non_mg cohort
		'C': np.std(transforms.log_transform(X_nonmg)[:,0]), #Standard deviation of log transformed kappa values for non_mg cohort
		'D': Vh_raw[0,1], #Second right singular vector in raw space                      
		'E': np.mean(transforms.log_transform(X_nonmg)[:,1]), #Mean of log transformed lambda values for non_mg cohort
		'F': np.std(transforms.log_transform(X_nonmg)[:,1]) #Standard deviation of log transformed lambda values for non_mg cohort
	} 
	visualize.plot_interval(X_nonmg, pc2_RI, z_transform, pc_transform)
	with open('./Output/pc_clonality_index.txt','w') as f:
		f.write('PC clonality index equation parameters\n')
		f.write('--------------------------------------\n')
		for key in equation_parameters.keys():
			f.write(f"\t%s = %.2f\n" % (key,equation_parameters[key]))
		f.write('\n\nPC clonality index reference interval\n')
		f.write('--------------------------------------\n')
		f.write(f"%s %%ile = %.2f\n" % (lb,pc2_RI[0]))
		f.write(f"%s %%ile = %.2f" % (ub,pc2_RI[1]))		
	return(X_nonmg,pc2_RI, equation_parameters, z_transform, pc_transform)

#%%
def evaluate_interval(nonmg: pd.DataFrame, mg: pd.DataFrame, lb: float, ub: float):
	"""Evaluate PC2-based interval for MG diagnosis"""
	X_mg=df2array(nonmg)
	X_nonmg, pc2_RI, _, z_transform, pc_transform = derive_interval(nonmg,lb,ub)
	performance = {}
	performance.update(performance.SeSp_sFLCR(X_nonmg,X_mg,0.26,1.65)) #evaluate performance of manufacturer's sFLC-ratio-based interval
	performance.update(performance.SeSp_PCA(X_nonmg, X_mg, pc2_RI, pc_transform, z_transform)) #evaluate performance of pc2-based interval
	pd.DataFrame(data=performance, index=['Measure']).T.to_csv('./Output/performance.csv')
	visualize.evaluation(X_nonmg,X_mg, pc2_RI, z_transform, pc_transform)

#%%
def apply_interval(nonmg: pd.DataFrame,cases: pd.DataFrame, lb: float, ub: float):
	"""Apply PC2-based interval to new patient data"""
	_, pc2_RI, _, z_transform, pc_transform = derive_interval(nonmg,lb,ub)
	X_cases = df2array(cases)
	L_cases = transforms.log_transform(X_cases) #apply log transform to X_cases
	Z_cases = z_transform(L_cases) #apply z transform to L_cases
	cases['pc2'] = pc_transform(Z_cases)[:,1] #pc2 projections for nonmg cohort
	cases['abnormal?'] = ((cases['pc2'] > pc2_RI[0]) & (cases['pc2'] < pc2_RI[1])) #add abnormality flag per PC2-based interval
	cases.to_csv('./Output/annotated_cases')
	visualize.cases(cases,pc2_RI, z_transform, pc_transform)



#%%
def main():
	parser = OptionParser()
	parser.add_option(
		"-n", "--nonmg",
		dest = "nonmg_fpath",
		help = "[optional] path to non-MG cohort csv file, default='.Data/WashU.p'"
	)
	parser.add_option(
		"-m", "--mg",
		dest = "mg_fpath",
		help = "[optional] path to MG cohort csv file, default=None"
	)   
	parser.add_option(
		"-c", "--cases",
		dest = "cases_fpath",
		help = "[optional] path to csv file for cases, default=None"
	)  
	parser.add_option(
		"-l", "--lower",
		dest = "user_lb",
		help = "[optional] % lower bound for reference interval, default=2.5"
	)  
	parser.add_option(
		"-u", "--upper",
		dest = "user_ub",
		help = "[optional] % upper bound for reference interval, default=97.5"
	)    
	(options, args) = parser.parse_args()

	if options.user_lb and options.user_ub:
		lb=float(options.user_lb)
		ub=float(options.user_ub)
	else:
		lb=2.5
		ub=97.5			   
		    
	if options.nonmg_fpath:
		nonmg = pd.read_csv(options.nonmg_fpath)
	else:
		with open('./Data/WashU.p', 'rb') as file:
			nonmg=pickle.load(file)
	print(derive_interval(nonmg,lb,ub))
			   
	if options.mg_fpath:
		mg = pd.read_csv(options.mg_fpath)
		evaluate_interval(nonmg,mg,lb,ub)

	if options.cases_fpath:
		cases = pd.read_csv(options.cases_fpath)
		apply_interval(nonmg,cases,lb,ub)

if __name__ == '__main__':
    main()
# %%
