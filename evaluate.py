import numpy as np
import transforms

def calc_ratio(X):
    return(X[:,0]/X[:,1])

def SeSp_sFLCR(X_nonmg,X_mg,lb,ub):
    print(np.logical_and(calc_ratio(X_nonmg)>lb,calc_ratio(X_nonmg)<ub).shape)
    Sp=sum(np.logical_and(calc_ratio(X_nonmg)>lb,calc_ratio(X_nonmg)<ub)/len(calc_ratio(X_nonmg)))
    Se=sum(~np.logical_and(calc_ratio(X_mg)>lb,calc_ratio(X_mg)<ub)/len(calc_ratio(X_mg)))
    sFLC_performance_dict = {'sFLC_Sp':Sp, 'sFLC_Se':Se}
    return(sFLC_performance_dict)

def SeSp_PCA(X_nonmg,X_mg,pc2_RI, pc_transform, z_transform):
    Sp_pc2=sum(np.logical_and(pc_transform(z_transform(transforms.log_transform(X_nonmg)))[:,1]>pc2_RI[0],
            pc_transform(z_transform(transforms.log_transform(X_nonmg)))[:,1]<pc2_RI[1]))/len(pc_transform(z_transform(transforms.log_transform(X_nonmg)))[:,1])
    Se_pc2=None
    Se_pc2=sum(~np.logical_and(pc_transform(z_transform(transforms.log_transform(X_mg)))[:,1]>pc2_RI[0],
        pc_transform(z_transform(transforms.log_transform(X_mg)))[:,1]<pc2_RI[1]))/len(pc_transform(z_transform(transforms.log_transform(X_mg)))[:,1])
    PC_performance_dict = {'PC_Sp':Sp_pc2, 'PC_Se':Se_pc2}
    return(PC_performance_dict)


