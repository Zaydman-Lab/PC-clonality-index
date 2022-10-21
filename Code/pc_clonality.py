
#%%
from optparse import OptionParser
import os
import generate
import pandas as pd 
import evaluate
import embed
import matplotlib.pyplot as plt

#%%
def generate_embedding(nonmg_path):
    ## if user passes non-MG cohort csv, will create embedding,
    ## return RI, equation, and plot based on inputted data
    ## if user does not pass non-MG cohort, embedding, 
    ## RI, equation and plot will be based on WU data


    df_nonmg = generate.load_func(nonmg_path)
    X_normal = generate.create_nparray(df_nonmg)
    z_transform, pc_transform = generate.generate_transforms(X_normal)
    pc2_col_norm = generate.calc_pc2(X_normal, z_transform, pc_transform)
    PCA_RI = generate.pc2_interval(2.5, 97.5, pc2_col_norm)
    equation_dict=generate.equation(df_nonmg, X_normal, z_transform, pc_transform)
    generate.plot_cohorts(X_normal, PCA_RI, z_transform, pc_transform)
    return(PCA_RI, equation_dict, X_normal, z_transform, pc_transform)

#%%
def evaluate_ri(X_normal, RI, pc_transform, z_transform, mg_path=None):
    ## optional. If MG cohort is not passed, only specificity will be calculated.
    ## if both MG and non-MG cohorts passed, both sens. and spec. will be calculated.
    ## if either cohort contains a creatinine column, performance will be broken out
    ## by eGFR

    ## is eGFR is not included
    lb=0.26
    ub=1.65
    df_abnormal=generate.load_func(mg_path)
    X_abnormal=generate.create_nparray(df_abnormal)
    if X_normal.shape[1]==2 and X_abnormal.shape[1]==2:
        SeFLC, SpFLC = evaluate.SeSp_sFLCR(X_normal,lb,ub,X_abnormal)
        SePC, SpPC = evaluate.SeSp_PCA(X_normal,RI,pc_transform,z_transform,X_abnormal)
        return(SeFLC, SpFLC, SePC, SpPC)

#%%
def embed_cases(cases_path,z_transform, pc_transform, PCA_RI):
    ## optional. Requires a different .csv file as input. Will output a new file with
    ## PC2 scores and a plot of where cases fall with respect to Katzmann and PC2 ref.
    ## ranges
    df_cases=embed.generate_embedding(cases_path,z_transform, pc_transform, PCA_RI)
    return(df_cases)

#%%
def main():
    path = os.getcwd()
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

    (options, args) = parser.parse_args()

    ## if a non-MG cohort is passes, use that to generate embeddings
    ## else, use the WashU cohort
    if options.nonmg_fn!=None:
        nonmg_path = path+'/'+options.nonmg_fn+'.csv'
    else:
        nonmg_path = path+'/non_MG_wu.csv'
    
    if options.mg_fn!=None:
        mg_path = path+'/'+options.mg_fn+'.csv'
    else:
        mg_path = None

    PCA_RI, equation_dict, X_normal, z_transform, pc_transform = generate_embedding(nonmg_path)
    equation_dict.update({'PCA_RI_low':PCA_RI[0], 'PCA_RI_high':PCA_RI[1]})
    if mg_path!=None:
        SeFLC, SpFLC, SePC, SpPC = evaluate_ri(X_normal, PCA_RI, pc_transform, z_transform, mg_path)
        equation_dict.update({'SpFLC': SpFLC, 'SeFLC': SeFLC, 'SpPC': SpPC, 'SePC': SePC})
    
    if options.cases_fn!=None:
        cases_path = path+'/'+options.cases_fn+'.csv'
        df_cases = embed_cases(cases_path,z_transform, pc_transform, PCA_RI)
        df_cases.to_csv('cases_pc2.to_csv',index=False)

    df_equation=pd.DataFrame(data=equation_dict, index=['var']).T
    df_equation.to_csv('output.csv')

if __name__ == '__main__':
    main()
# %%
