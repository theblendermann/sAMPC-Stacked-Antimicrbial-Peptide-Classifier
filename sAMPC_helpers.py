#Function file for the helper functions of sAMPC
import os
import numpy as np
import pandas as pd
import joblib#type: ignore
from Bio import SeqIO
#from typing import Optional
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler

def fasta_to_tsv(prodigal_out: str): # Function to read fasta files and convert them into DataFrames
    proteins_df=pd.DataFrame() #Create an empty DataFrame
    sequences: list[str]=[] #Empty List to store sequnce data
    sequence_names: list[str]=[] #Empty list to store sequence headers
    

    for record in SeqIO.parse(prodigal_out,'fasta'): #Parses input file, has to be in .fasta format
        fixed_sequences=str(record.seq).lstrip('M').rstrip('*') #Removes the starting methionine and ending stop codon
        if len(fixed_sequences)<=255 and len(fixed_sequences)>=5: #Only keeps protein sequences which are less than 255 amino acids and more than 5 amino acids in length
            sequence_names.append(f'{record.id}')
            sequences.append(fixed_sequences) 
    
    proteins_df['Names']=sequence_names #Names column is populated with protein names
    proteins_df['SEQ']=sequences #SEQ column is populated with protein sequences

    return proteins_df

def excel_reader(file_path: str, sheet_name: str):# Function to read excel files (.xlsx files)
    features=pd.ExcelFile(file_path) #Takes in file path
    features_sheet_name={}# Creates a blank dictionary to store sheets found in the excel file
    for features_sheet in features.sheet_names:
        features_sheet_name[features_sheet]=features.parse(features_sheet)
    features_df=features_sheet_name[f'{sheet_name}'] #Loads a sheet called "Features", it is the default sheet that the feature extraction functions create.

    return features_df

def data_scaler(data: pd.DataFrame, scaler_type: str, scaler_save_path: str, pnratio: str)-> tuple[np.ndarray,...]| np.ndarray:

    if scaler_type=='MinMax':
        scaler=MinMaxScaler(feature_range=(-1,1))
    elif scaler_type=='StandardScaler':
        scaler=StandardScaler()
    elif scaler_type=='RobustScaler':
        scaler=RobustScaler()
    else:
        print('Choose between RobustScaler, MinMaxScaler, and StandardScaler')

    training_df=scaler.fit_transform(data)
    joblib.dump(scaler,f'{scaler_save_path}/{scaler_type}_{pnratio}.pkl')

    return training_df, scaler #Will return the ndarray after fitting the scaler and the scaler it self [0]-> pd.Dataframe, [1]-> scaler

def stack_helper(model_parent_path: str,X_data: pd.DataFrame, model_type: str, feature_set: str): #Helper training function for the stacking training function
    ensemble_predictions=pd.DataFrame()
    for model in sorted(os.listdir(model_parent_path)):
        if model.endswith('.joblib'):
            model_path=os.path.join(model_parent_path, model)
            if model_type=='NN':
                model_scaler=joblib.load(model_path)
                X_data_scaled=model_scaler['scaler'].transform(X_data)
                ensemble_predictions[f'{feature_set} {model}']=model_scaler['model'].predict_proba(X_data_scaled)[:,1] #For NNs, Model and Scaler are stored in the same file, model_scaler['scaler'] will call the scaler, and model_scaler['model'] will call the model
            else:
                classifier=joblib.load(model_path)
                ensemble_predictions[f'{feature_set} {model}']=classifier.predict_proba(X_data)[:,1]
    #ensemble_predictions[f'{type}_{feature_set}_avg']=ensemble_predictions.mean(axis=1)
    return ensemble_predictions
