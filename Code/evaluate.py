import numpy as np
import generate

def calc_ratio(X):
    return(X[:,0]/X[:,1])

def calc_pc(X):
    pass

def SeSp_sFLCR(X_normal,lb,ub, X_abnormal=None):
    Sp=sum(np.logical_and(calc_ratio(X_normal)>lb,calc_ratio(X_normal)<ub)/len(calc_ratio(X_normal)))
    Se=None
    if X_abnormal.shape[1]==2:
        Se=sum(~np.logical_and(calc_ratio(X_abnormal)>lb,calc_ratio(X_abnormal)<ub)/len(calc_ratio(X_abnormal)))
    return(Se,Sp)

def SeSp_PCA(X_normal,RI, pc_transform, z_transform,X_abnormal=None):
    Sp_pc2=sum(np.logical_and(pc_transform(z_transform(generate.log_transform(X_normal)))[:,1]>RI[0],
            pc_transform(z_transform(generate.log_transform(X_normal)))[:,1]<RI[1]))/len(pc_transform(z_transform(generate.log_transform(X_normal)))[:,1])
    Se=None
    if X_abnormal.shape[1]==2:
        Se_pc2=sum(~np.logical_and(pc_transform(z_transform(generate.log_transform(X_abnormal)))[:,1]>RI[0],
            pc_transform(z_transform(generate.log_transform(X_abnormal)))[:,1]<RI[1]))/len(pc_transform(z_transform(generate.log_transform(X_abnormal)))[:,1])
    return(Se_pc2,Sp_pc2)


