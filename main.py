"""
pc_clonality.py
Created by Vahid Azimi
11/29/22

Purpose: main point of entry into PC-clonality-index package

This script will execute one of the three cases described in the .README based on the user input

Functions:


generate_embedding(nonmg_path,lb,ub,output_path):
    '''Returns tranfform object and embeddings'''
"""

#%%
from optparse import OptionParser
import os
import generate
import pandas as pd 
import evaluate
import embed
import matplotlib.pyplot as plt
import pickle
from typing import Tuple
import visualize



#%%

def df2array(df: pd.DataFrame)->np.array:
	"""Returns np array of kappa and labmda values derived from input dataframe"""
	X = np.column_stack((df['kappa'].to_numpy(), df['lambda'].to_numpy())) 
	return(X)


def case_1(nonmg: pd.DataFrame,lb: float,ub: float)->Tuple[pd.DataFrame]:
	"""Returns PC2-based reference interval, dictionary of PC2 equation parameters, and functions for Z and pc transforms of non-MG cohort"""
	X_nonmg = df2array(nonmg)
	L_nonmg = transforms.log_transform(X_nonmg) #apply log transform to X_nonmg
	Z_transform = transforms.make_ztransform(L_nonmg) #compute z transform of L_nonmg
	Z_nonmg = Z_tranform(log_nonmg) #apply z transform to L_nonmg
	pc_transform, U, S, Vh = transforms.make_pctransform(Z_nonmg,return_decomposition=True) #compute pc transform of Z_nonmg
	pc2 = pc_transform(Z_nonmg)[:,1] #pc2 projections for nonmg cohort
	pc2_RI = [np.percentile(pc2,lb),np.percentile(pc2,ub) #compute reference interval for pc2
	Vh_raw=log_transform(z_transform(Vh,inverse=True),inverse=True) #project right singular vectors from z space into raw sFLC space
	equation_parameters = {'A': Vh_raw[0,0], #First right singular vector in raw space
		'B': np.mean(log_transform(X_normal)[:,0]), #Mean of log transformed kappa values for non_mg cohort
		'C': np.std(log_transform(X_normal)[:,0]), #Standard deviation of log transformed kappa values for non_mg cohort
		'D': Vh_raw[0,1], #Second right singular vector in raw space                      
		'E': np.mean(log_transform(X_normal)[:,1]), #Mean of log transformed lambda values for non_mg cohort
		'F': np.std(log_transform(X_normal)[:,1])} #Standard deviation of log transformed lambda values for non_mg cohort
	visualize.plot_case1(X_nonmg, pc2_RI, z_transform, pc_transform)
	return(X_nonmg,pc2_RI, equation_parameters, z_transform, pc_transform)

def case_2(non_mg: pd.DataFrame, mg: pd.DataFrame, lb: float, ub: float)->dict:
	"""Returns dictionary of performance metrics for manufacturers sFLC-ratio-based and PC2-based reference intervals"""
	X_mg=df2array(non_mg)
	X_nonmg, pc2_RI, equation_parameters, z_transform, pc_transform = case_1(non_mg,lb,ub)
        performance = {}
	performance.update(evaluate.SeSp_sFLCR(X_nonmg,X_mg,0.26,1.65) #evaluate performance of manufacturer's sFLC-ratio-based interval
        performance.update(evaluate.SeSp_PCA(X_nonmg, X_mg, pc2_RI, pc_transform, z_transform) #evaluate performance of pc2-based interval
        return(performance)

#%%
def case3(nonmg: pd.DataFrame,cases: pd.DataFrame)->pd.DataFrame:
	"""Returns PC2 embeddings and normality flags for input cases based on nonmg data model"""
	X_nonmg, pc2_RI, equation_parameters, z_transform, pc_transform = case_1(non_mg,lb,ub)	
	X_cases = df2array(cases)
	L_cases = transforms.log_transform(X_cases) #apply log transform to X_cases
	Z_cases = Z_tranform(L_cases) #apply z transform to L_cases
	cases['pc2'] = pc_transform(Z_cases)[:,1] #pc2 projections for nonmg cohort
	cases['abnormal?'] = ((df['pc2'] > pc2_RI[0]) & (df['pc2'] < pc2_RI[1]))		   
    	return(cases)

#%%
def main():

    parser = OptionParser()
    parser.add_option("-n", "--nonmg",
                    dest = "nonmg_fpath",
                    help = "[optional] path to non-MG cohort csv file, default='.Data/WashU.p'")
    parser.add_option("-m", "--mg",
                    dest = "mg_fpath",
                    help = "[optional] path to MG cohort csv file, default=None")   
    parser.add_option("-c", "--cases",
                    dest = "cases_fpath",
                    help = "[optional] path to csv file for cases, default=None")  
    parser.add_option("-l", "--lower",
                    dest = "user_lb",
                    help = "[optional] % lower bound for reference interval, default=2.5")  
    parser.add_option("-u", "--upper",
                    dest = "user_ub",
                    help = "[optional] % upper bound for reference interval, default=97.5")    
    (options, args) = parser.parse_args()

    if options.nonmg_fpath:
        nonmg = pd.read_csv(options.nonmg_fpath)
    else:
        with open('./Data/WashU.p', 'rb') as file:
            nonmg=pickle.load(file)
    
    if options.mg_fpath:
        mg = pd.read_csv(options.mg_fpath)
    else:
        mg = None
        
    if options.cases_fpath:
        cases = pd.read_csv(options.cases_fpath)
    else:
        cases = None

    if options.user_lb and options.user_ub:
        lb=float(options.user_lb)
        ub=float(options.user_ub)
    else:
        lb=2.5
        ub=97.5      

    PCA_RI, pc_dict, X_normal, z_transform, pc_transform = generate_embedding(nonmg_path,lb,ub,output_path)
    pc_dict.update({'PCA_RI_low':PCA_RI[0], 'PCA_RI_high':PCA_RI[1]})

    if mg_path!=None:
        performance_dict=evaluate_ri(X_normal, PCA_RI, pc_transform, z_transform, mg_path)
        performance_dict=pd.DataFrame(data=performance_dict, index=['Measure']).T
        performance_dict.to_csv(output_path+'/performance.csv')
    
    if options.cases_fn!=None:
        cases_path = options.cases_fn
        df_cases = embed_cases(cases_path,z_transform, pc_transform, PCA_RI)
        df_cases.to_csv(output_path+'/cases_pc2.csv',index=False)

    df_pc_dict=pd.DataFrame(data=pc_dict, index=['var']).T
    df_pc_dict.to_csv(output_path+'/pc_vars.csv')

if __name__ == '__main__':
    main()
# %%
