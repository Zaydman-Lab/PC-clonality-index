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



#%%

def case_1(nonmg: pd.DataFrame,lb: float,ub: float)->Tuple[pd.DataFrame]:
    """Returns embeddings of non"""
	X_nonmg = np.column_stack((df['kappa'].to_numpy(), df['lambda'].to_numpy())) #convert to np array
	L_nonmg = transforms.log_transform(X_nonmg) #apply log transform to X_nonmg
	Z_transform = transforms.make_ztransform(L_nonmg) #compute z transform of L_nonmg
	Z_nonmg = Z_tranform(log_nonmg) #apply z transform to L_nonmg
	pc_transform, U, S, Vh = transforms.make_pctransform(Z_nonmg,return_decomposition=True) #compute z transform of Z_nonmg
	pc2 = pc_transform(Z_nonmg)[:,1] #pc2 projections for nonmg cohort
	pc2_RI = [np.percentile(pc2,lb),np.percentile(pc2,ub) #computer reference interval for pc2
	Vh_raw=log_transform(z_transform(Vh,inverse=True),inverse=True) #project right singular vectors from z space into raw sFLC space
	equation_dict = {'A': Vh_raw[0,0], #First right singular vector in raw space
		'B': np.mean(log_transform(X_normal)[:,0]), #Mean of log transformed kappa values for non_mg cohort
		'C': np.std(log_transform(X_normal)[:,0]), #Standard deviation of log transformed kappa values for non_mg cohort
		'D': Vh_raw[0,1], #Second right singular vector in raw space                      
		'D': np.mean(log_transform(X_normal)[:,1]), #Mean of log transformed lambda values for non_mg cohort
		'F': np.std(log_transform(X_normal)[:,1])} #Standard deviation of log transformed lambda values for non_mg cohort
	generate.plot_nonmg(X_nonmg, pc2_RI, z_transform, pc_transform)
	return(pc2_RI, equation_dict, X_nonmg, z_transform, pc_transform)

#%%
def evaluate_ri(X_normal, RI, pc_transform, z_transform, mg_path=None):
    '''input: X_normal, RI, pc_transform, z_transform, mg_path
    output: SeFLC, SpFLC, SePC, SpPC'''

    lb=0.26
    ub=1.65
    df_abnormal=generate.load_func(mg_path)
    X_abnormal=generate.create_nparray(df_abnormal)
    if X_normal.shape[1]==2 and X_abnormal.shape[1]==2:
        performance_dict = evaluate.SeSp_sFLCR(X_normal,lb,ub,X_abnormal)
        PC_performance_dict = evaluate.SeSp_PCA(X_normal,RI,pc_transform,z_transform,X_abnormal)
        performance_dict.update(PC_performance_dict)
        return(performance_dict)

#%%
def embed_cases(cases_path,z_transform, pc_transform, PCA_RI):
    ## optional. Requires a different .csv file as input. Will output a new file with
    ## PC2 scores and a plot of where cases fall with respect to Katzmann and PC2 ref.
    ## ranges
    df_cases=embed.generate_embedding(cases_path,z_transform, pc_transform, PCA_RI)
    return(df_cases)

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
