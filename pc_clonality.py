"""
pc_clonality.py
Created by Vahid Azimi and Mark Zaydman
11/29/22

Purpose: main point of entry into PC-clonality-index package

This script will enable the user to derive a pc2-based equation and interval using either a user defined or
the WashU reference nonmg cohort, to evaluate the sensitivity and specificity of the pc2-based
interval for mg diagnosis using either a user defined or the Wash U mg cohorts, and to apply 
the pc2 based interval a new set of cases. The output behavior is controlled by a set of 
optional command line arguments.

Functions:
	parse_args()->Tuple[str]:
		'''Returns parsed arguments'''
	df2array(df: pd.DataFrame)->np.array:
		'''Returns np array of kappa and labmda values derived from input dataframe'''
	derive_interval(nonmg: pd.DataFrame,lb: float,ub: float)->Tuple[np.array,list[float,float],Callable[[np.array,bool],np.array],Callable[[np.array,bool],np.array]]:
		'''Derive PC2-based equation and refence interval using WU nonmg cohort or user input nonmg cohort'''
	evaluate_interval(nonmg: pd.DataFrame, mg: pd.DataFrame, lb: float, ub: float)->None:
		'''Evaluate PC2-based interval for MG diagnosis'''
	def apply_interval(nonmg: pd.DataFrame,cases: pd.DataFrame, lb: float, ub: float):
		'''Apply PC2-based interval to new patient data'''
	main(nonmg_fpath: str = None,mg_fpath: str = None,cases_fpath: str = None,lb: float = 2.5,ub: float = 97.5)->None:
		'''Generates output per user specified options'''

"""

#%% imports
from optparse import OptionParser
import pandas as pd 
import pickle
from typing import Tuple, Callable
import visualize
import transforms
import evaluate
import numpy as np

def parse_args()->Tuple[str]:
	"""Returns parsed arguments"""
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
	return(options)	

def df2array(df: pd.DataFrame)->np.array:
	"""Returns np array of kappa and labmda values derived from input dataframe"""
	X = np.column_stack((df['kappa'].to_numpy(), df['lambda'].to_numpy())) 
	return(X)

def derive_interval(nonmg: pd.DataFrame,lb: float,ub: float)->Tuple[np.array,list[float,float],Callable[[np.array,bool],np.array],Callable[[np.array,bool],np.array]]:
	"""Derive PC2-based equation and refence interval using WU nonmg cohort or user input nonmg cohort"""
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
	visualize.plot_sflc(X_nonmg, pc2_RI, z_transform, pc_transform, parameters=equation_parameters)
	with open('./Output/pc_clonality_index.txt','w') as f:
		f.write('PC clonality index equation parameters\n')
		f.write('--------------------------------------\n')
		for key in equation_parameters.keys():
			f.write(f"\t%s = %.2f\n" % (key,equation_parameters[key]))
		f.write('\n\nPC clonality index reference interval\n')
		f.write('--------------------------------------\n')
		f.write(f"\t%s %%ile = %.2f\n" % (lb,pc2_RI[0]))
		f.write(f"\t%s %%ile = %.2f" % (ub,pc2_RI[1]))		
	return(X_nonmg,pc2_RI, equation_parameters, z_transform, pc_transform)

def evaluate_interval(nonmg: pd.DataFrame, mg: pd.DataFrame, lb: float, ub: float)->None:
	"""Evaluate PC2-based interval for MG diagnosis"""
	X_mg=df2array(mg)
	X_nonmg, pc2_RI, _, z_transform, pc_transform = derive_interval(nonmg,lb,ub)
	performances=[]
	performances.append(pd.DataFrame.from_dict(evaluate.SeSp_sFLCR(X_nonmg,X_mg,0.26,1.65),orient='index')) #evaluate manufacturer's sFLC-ratio based interval
	performances.append(pd.DataFrame.from_dict(evaluate.SeSp_PCA(X_nonmg, X_mg, pc2_RI, pc_transform, z_transform),orient='index')) #evaluate pc2 interval (aka pc clonality index)
	performance=pd.concat(performances,axis=0)
	with open('./Output/performance.txt','w') as outfile:
		performance.to_string(outfile)
	visualize.plot_sflc(X_nonmg, pc2_RI, z_transform, pc_transform, X_mg=X_mg, performance=performance)

def apply_interval(nonmg: pd.DataFrame,cases: pd.DataFrame, lb: float, ub: float):
	"""Apply PC2-based interval to new patient data"""
	X_nonmg, pc2_RI, _, z_transform, pc_transform = derive_interval(nonmg,lb,ub)
	X_cases = df2array(cases)
	L_cases = transforms.log_transform(X_cases) #apply log transform to X_cases
	Z_cases = z_transform(L_cases) #apply z transform to L_cases
	cases['pc2'] = pc_transform(Z_cases)[:,1] #pc2 projections for nonmg cohort
	cases['abnormal?'] = ((cases['pc2'] > pc2_RI[0]) & (cases['pc2'] < pc2_RI[1])) #add abnormality column per PC2-based interval
	cases.to_csv('./Output/annotated_cases.csv')
	visualize.plot_sflc(X_nonmg,pc2_RI, z_transform, pc_transform, X_cases=X_cases)

def main(nonmg_fpath: str = None,mg_fpath: str = None,cases_fpath: str = None,lb: float = 2.5,ub: float = 97.5)->None:
	"""Generates output per user specified options"""
	if nonmg_fpath:
		nonmg = pd.read_csv(options.nonmg_fpath)
	else:
		with open('./Data/WashU.p', 'rb') as file:
			nonmg=pickle.load(file)

	#generate and evaluate pc2-based interval			   
	if mg_fpath:
		mg = pd.read_csv(options.mg_fpath)
		evaluate_interval(nonmg,mg,lb,ub)

	#generate and apply pc2-based interval
	if cases_fpath:
		cases = pd.read_csv(options.cases_fpath)
		apply_interval(nonmg,cases,lb,ub)

	#only generate pc2-based interval and equation
	if mg_fpath==None and cases_fpath==None:
		derive_interval(nonmg,lb,ub)



if __name__ == '__main__':
	options=parse_args()

	if options.user_lb and options.user_ub:
		lb=float(options.user_lb)
		ub=float(options.user_ub)
	else:
		lb=2.5
		ub=97.5	

	main(
		nonmg_fpath=options.nonmg_fpath,
		mg_fpath=options.mg_fpath,
		cases_fpath=options.cases_fpath,
		lb=lb,
		ub=ub
	)
