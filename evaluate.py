import numpy as np
import transforms

def calc_ratio(X):
    return(X[:,0]/X[:,1])

def SeSp_sFLCR(X_nonmg,X_mg,lb,ub):
    Sp=sum(np.logical_and(calc_ratio(X_nonmg)>lb,calc_ratio(X_nonmg)<ub)/len(calc_ratio(X_nonmg)))
    Se=sum(~np.logical_and(calc_ratio(X_mg)>lb,calc_ratio(X_mg)<ub)/len(calc_ratio(X_mg)))
    performance = {"Manufacturer's sFLC-ratio interval":{'Sensitivity':'%.2f'%Se,'Specificity':'%.2f'%Sp}}
    return(performance)

def SeSp_PCA(X_nonmg,X_mg,pc2_RI, pc_transform, z_transform):
    Sp=sum(np.logical_and(pc_transform(z_transform(transforms.log_transform(X_nonmg)))[:,1]>pc2_RI[0],
            pc_transform(z_transform(transforms.log_transform(X_nonmg)))[:,1]<pc2_RI[1]))/len(pc_transform(z_transform(transforms.log_transform(X_nonmg)))[:,1])
    Se=sum(~np.logical_and(pc_transform(z_transform(transforms.log_transform(X_mg)))[:,1]>pc2_RI[0],
        pc_transform(z_transform(transforms.log_transform(X_mg)))[:,1]<pc2_RI[1]))/len(pc_transform(z_transform(transforms.log_transform(X_mg)))[:,1])
    performance = {"PC clonality index interval":{'Sensitivity':'%.2f'%Se, 'Specificity':'%.2f'%Sp}}
    return(performance)


