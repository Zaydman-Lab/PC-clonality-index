import numpy as np

def calc_ratio(X):
    return(X[:,0]/X[:,1])

def SeSp_sFLCR(X_normal,lb,ub, X_abnormal=None):
    Sp=sum(np.logical_and(calc_ratio(X_normal)>lb,calc_ratio(X_normal)<ub)/len(calc_ratio(X_normal)))
    Se=None
    if X_abnormal.shape[1]==2:
        Se=sum(~np.logical_and(calc_ratio(X_abnormal)>lb,calc_ratio(X_abnormal)<ub)/len(calc_ratio(X_abnormal)))
    sFLC_performance_dict = {'sFLC_Sp':Sp, 'sFLC_Se':Se}
    return(sFLC_performance_dict)

def SeSp_PCA(X_normal,RI, pc_transform, z_transform,X_abnormal=None):
    Sp_pc2=sum(np.logical_and(pc_transform(z_transform(generate.log_transform(X_normal)))[:,1]>RI[0],
            pc_transform(z_transform(generate.log_transform(X_normal)))[:,1]<RI[1]))/len(pc_transform(z_transform(generate.log_transform(X_normal)))[:,1])
    Se_pc2=None
    if X_abnormal.shape[1]==2:
        Se_pc2=sum(~np.logical_and(pc_transform(z_transform(generate.log_transform(X_abnormal)))[:,1]>RI[0],
            pc_transform(z_transform(generate.log_transform(X_abnormal)))[:,1]<RI[1]))/len(pc_transform(z_transform(generate.log_transform(X_abnormal)))[:,1])
    PC_performance_dict = {'PC_Sp':Sp_pc2, 'PC_Se':Se_pc2}
    return(PC_performance_dict)


