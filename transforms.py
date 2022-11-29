"""
transforms.py

Created by Vahid Azimi

Purpose: Helper functions for generating transforms used in PC clonality index


"""




import numpy as np
import pandas as pd
import scipy.linalg

def log_transform(INPUT,inverse=False):
    OUTPUT=np.copy(INPUT)
    if inverse:
        OUTPUT=np.exp(INPUT)
    else:
        OUTPUT=np.log(INPUT)
    return(OUTPUT)

def make_ztransform(X):
    def z_transform(INPUT,inverse=False):
        OUTPUT=np.copy(INPUT)
        if inverse:
            for col in range(INPUT.shape[1]):
                OUTPUT[:,col]=INPUT[:,col]*np.std(X[:,col])+np.mean(X[:,col])
        else:
            for col in range(INPUT.shape[1]):
                OUTPUT[:,col]=(INPUT[:,col]-np.mean(X[:,col]))/np.std(X[:,col])
        return(OUTPUT)
    return(z_transform)

def make_pctransform(X):
    U,S,Vh=scipy.linalg.svd(X,full_matrices=False)
    def pc_transform(INPUT,inverse=False, return_decomposition=False):
        if inverse:
            OUTPUT=INPUT@Vh
        else:
            OUTPUT=np.copy(INPUT)
            for row in range(INPUT.shape[0]):
                OUTPUT[row,:]=INPUT[row,:]@Vh.T
        if return_decomposition:
            return(OUTPUT, U, S, Vh)
        else:
            return(OUTPUT)
    return(pc_transform)

def find_RI(X,lb,ub):
    RI=[np.percentile(X,lb),
        np.percentile(X,ub)]
    return(RI)



#%% convert dataframes into numpy array
def create_nparray(df):

	X = np.column_stack((df['kappa'].to_numpy(), df['lambda'].to_numpy()))
	return(X)

#%% perform transforms using non-MG cohort
def generate_transforms(X_normal):
	z_transform = make_ztransform(log_transform(X_normal))
	pc_transform = make_pctransform(z_transform(log_transform(X_normal)))

	return (z_transform, pc_transform)

def calc_pc2(X,z_transform,pc_transform):

	pc2_col=pc_transform(z_transform(log_transform(X)))[:,1]
	return(pc2_col)

#%% determine PC2 reference interval
def pc2_interval(lb, ub, pc2_col_norm):

	PCA_RI=find_RI(pc2_col_norm,lb,ub)
	return(PCA_RI)

def equation(df_nonmg, X_normal, z_transform, pc_transform):
	kappa_mean = np.mean(log_transform(X_normal)[:,0])
	lambda_mean = np.mean(log_transform(X_normal)[:,1])
	kappa_std = np.std(log_transform(X_normal)[:,0])
	lambda_std = np.std(log_transform(X_normal)[:,1])
	pc, U, S, Vh = pc_transform(z_transform(log_transform(X_normal)),return_decomposition=True)
	Vh = log_transform(z_transform(Vh,inverse=True),inverse=True)
	equation_dict = {'A': Vh[0,0], 'B': Vh[0,1], 'kappa_mean': kappa_mean, 'kappa_std': kappa_std, 'lambda_mean': lambda_mean, 'lambda_std': lambda_std}
	return(equation_dict)
	
#%%    
def plot_cohorts(X_normal, PCA_RI, z_transform, pc_transform, output_path):

    Katz_lb=0.26
    Katz_ub=1.65
    fig,ax=plt.subplots()
    fig.patch.set_facecolor('xkcd:white')

    fig.set_size_inches(4,4)
    Katz_LB=np.column_stack((np.linspace(min(X_normal[:,0]),max(X_normal[:,0]),num=50),1/Katz_lb*np.linspace(min(X_normal[:,0]),max(X_normal[:,0]),num=50)))
    Katz_UB=np.column_stack((np.linspace(min(X_normal[:,0]),max(X_normal[:,0]),num=50),1/Katz_ub*np.linspace(min(X_normal[:,0]),max(X_normal[:,0]),num=50)))
    PCA_LB=np.column_stack((np.linspace(-3,7,num=50),np.ones(50)*PCA_RI[0]))
    PCA_UB=np.column_stack((np.linspace(-3,7,num=50),np.ones(50)*PCA_RI[1]))
    sns.scatterplot(x=X_normal[:,0],y=X_normal[:,1])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$[\kappa]\quad$(g/dL)')
    ax.set_ylabel(r'$[\lambda]\quad$(g/dL)')

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    xlims=ax.get_xlim()
    ylims=ax.get_ylim()

    ax.plot(Katz_LB[:,0],Katz_LB[:,1],color='k', linestyle=':')
    ax.plot(Katz_UB[:,0],Katz_UB[:,1],color='k', linestyle=':')
    LB=log_transform(z_transform(pc_transform(PCA_LB,inverse=True),inverse=True),inverse=True)
    UB=log_transform(z_transform(pc_transform(PCA_UB,inverse=True),inverse=True),inverse=True)
    if isinstance(LB,np.ndarray):
        ax.plot(LB[:,0],LB[:,1],'--k')
    if isinstance(UB,np.ndarray):
        ax.plot(UB[:,0],UB[:,1],'--k')
    line_katz = Line2D([0,1],[0,1],linestyle=':', color='black')
    line_pc2 = Line2D([0,1],[0,1],linestyle='--', color='black')
    ax.legend([line_katz, line_pc2],['Katzmann', 'PC2'])
    plt.savefig(output_path+'/pc_space.png')
