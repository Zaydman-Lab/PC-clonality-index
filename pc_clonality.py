
#%%
from optparse import OptionParser
import os
import generate
import pandas as pd 
import evaluate
import embed
import matplotlib.pyplot as plt
import pickle

#%%
def generate_embedding(nonmg_path,lb,ub,output_path):
    '''input: nonmg_path, lb, ub
    output: PCA_RI, equation_dict, X_normal, z_transform, pc_transform'''


    df_nonmg = generate.load_func(nonmg_path)
    X_normal = generate.create_nparray(df_nonmg)
    z_transform, pc_transform = generate.generate_transforms(X_normal)
    pc2_col_norm = generate.calc_pc2(X_normal, z_transform, pc_transform)
    PCA_RI = generate.pc2_interval(lb, ub, pc2_col_norm)
    equation_dict=generate.equation(df_nonmg, X_normal, z_transform, pc_transform)
    generate.plot_cohorts(X_normal, PCA_RI, z_transform, pc_transform, output_path)
    return(PCA_RI, equation_dict, X_normal, z_transform, pc_transform)

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
    data_path=os.getcwd()+'/Data'
    output_path=os.getcwd()+'/Output'
    parser = OptionParser()

    parser.add_option("-n", "--nonmg",
                    dest = "nonmg_fn",
                    help = "non-MG cohort csv filename")
    parser.add_option("-m", "--mg",
                    dest = "mg_fn",
                    help = "MG cohort csv filename")   
    parser.add_option("-c", "--cases",
                    dest = "cases_fn",
                    help = "cases for embedding csv filename")  
    parser.add_option("-l", "--lower",
                    dest = "user_lb",
                    help = "user-defined lower bound")  
    parser.add_option("-u", "--upper",
                    dest = "user_ub",
                    help = "user-defined upper bound")    

    (options, args) = parser.parse_args()

    if options.nonmg_fn!=None:
        nonmg_path = options.nonmg_fn
    else:
        nonmg_path = './Data/WashU.p'
    
    if options.mg_fn!=None:
        mg_path = options.mg_fn
    else:
        mg_path = None

    lb=2.5
    ub=97.5
    if options.user_lb!=None and options.user_ub!=None:
        lb=float(options.user_lb)
        ub=float(options.user_ub)

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
