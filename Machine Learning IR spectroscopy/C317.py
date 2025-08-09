import numpy as np
import pandas as pd
import os
import sklearn.decomposition
import sklearn.model_selection
import sklearn.neighbors

def Interpolate(dataframe):
    index_list=list(dataframe.index)
    Min,Max = min(index_list),max(index_list)
    minimum_index, maximum_index = round(Min,0),round(Max,0)
    ints=list(range(400,4001))
    NaNlist=[np.NaN for i in ints]
    df=pd.DataFrame(data=NaNlist,columns=[dataframe.columns[0]],index=ints)
    merged=pd.concat([dataframe,df])
    merged=merged.sort_index()
    df2=merged.interpolate()
    new_data=[df2.loc[int] for int in ints]
    df3=pd.DataFrame(data=new_data,columns=[dataframe.columns[0]],index=ints)
    return df3

def Normalise(dataframe):
    def numerical_integral(x_values,y_values):
        summation=0
        for i in range(1,len(x_values)):
            summation+=0.5*(y_values[i]+y_values[i-1])*(x_values[i]-x_values[i-1])
        return summation
    dataframe[dataframe.columns[0]]=dataframe[dataframe.columns[0]]/((numerical_integral(dataframe.index,dataframe.values)))
    return dataframe

def Narrow(dataframe):
    new_indices=list(range(630,881))
    new_data=[dataframe.iloc[int-400] for int in new_indices]
    dataframe=pd.DataFrame(data=new_data, columns=[dataframe.columns[0]],index=new_indices)
    return dataframe

def perform_pca(dataframe,n):
    PCA=sklearn.decomposition.PCA(n)
    x=PCA.fit_transform(dataframe.T)
    df=pd.DataFrame(data=x.T, columns=dataframe.columns)
    return df

def load_spectra(n,m): #slicing added :-4/5 to reomve the .txt from column headers as well as number #m=0 for full name, 1 for smaller
    filenames=[file.name for file in os.scandir('Data')]
    IRData=[Narrow(Normalise(Interpolate(pd.read_csv(f"data/{file}",skiprows=4,delimiter="\s+",names=[file[:-(4+2*m)]],index_col=0)))) for file in filenames]
    x=IRData.pop(0)
    while IRData!=[]:
        x=pd.concat([x,IRData.pop(0)],axis=1)
    if n==0:
        return x
    else:
        return perform_pca(x,n)
def load_new_spectra(n,m): #slicing added :-4/5 to reomve the .txt from column headers as well as number #m=0 for full name, 1 for smaller
    filenames=[file.name for file in os.scandir('new_data') if file[-4:]=='.txt']
    IRData=[Narrow(Normalise(Interpolate(pd.read_csv(f"new_data/{file}",skiprows=4,delimiter="\s+",names=[file[:-(4+2*m)]],index_col=0)))) for file in filenames]
    x=IRData.pop(0)
    while IRData!=[]:
        x=pd.concat([x,IRData.pop(0)],axis=1)
    if n==0:
        return x
    else:
        return perform_pca(x,n)

def MachineLearn(dataframe,t_size,k):
    UniqueNames=list(set(dataframe.columns))
    train_samples, test_samples = sklearn.model_selection.train_test_split(UniqueNames,test_size=t_size)
    train_labels=[]
    test_labels=[]
    for name in train_samples:
        train_labels+=5*[name[0]]
    for name in test_samples:
        test_labels+=5*[name[0]]
    OBJ=sklearn.neighbors.KNeighborsClassifier(k)
    OBJ.fit(dataframe[train_samples].T,train_labels)
    return OBJ.score(dataframe[test_samples].T,test_labels)