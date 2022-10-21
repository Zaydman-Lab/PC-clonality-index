import generate

def generate_embedding(cases_path,z_transform,pc_transform,PCA_RI):

    df=generate.load_func(cases_path)
    X=generate.create_nparray(df)
    pc2_col = generate.calc_pc2(X,z_transform,pc_transform)
    df['pc2']=pc2_col
    filt_norm = (df['pc2'] < PCA_RI[0]) | (df['pc2'] > PCA_RI[1])
    df['abnormal?'] = 0
    df.loc[filt_norm, 'abnormal?'] = 1
    return(df)