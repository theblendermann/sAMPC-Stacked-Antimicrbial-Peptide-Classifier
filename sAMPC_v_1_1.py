#Function file for the Stacked Antimicrobial Peptide classifier sAMPC
import joblib
import pandas as pd
from sAMPC_helpers import excel_reader
from sAMPC_helpers import data_scaler
from sAMPC_helpers import stack_helper



def stack_predictions(RF_MACREL_path: str, RF_AMPeP_path: str, NN_MACREL_path: str, NN_AMPeP_path: str, sAMPeP_path: str, 
                      sMACREL_path: str, MACREL_features: str, AMPeP_features: str, test: bool, scaler_math: str, scaler_loc: str):

    stack_df_AMPeP=pd.DataFrame()#make a df to load MetaModel Features
    stack_df_MACREL=pd.DataFrame()#make a df to load MetaModel Features
    stack_out=pd.DataFrame()#make a df to store MetaModel predictions
    
    data_MACREL=excel_reader(MACREL_features, sheet_name='Features')
    data_AMPeP=excel_reader(AMPeP_features, sheet_name='Features')

    if test==True:
        X_MACREL = data_MACREL.drop(columns=['Names', 'SEQ', 'Class'])
        X_MACREL_scaled= data_scaler(data=X_MACREL, use_prescaled=True, scaler_type=scaler_math, scaler_save_path=scaler_loc)
        Y_MACREL = data_MACREL['Class']
        X_AMPeP = data_AMPeP.drop(columns=['Names', 'SEQ', 'Class'])
        X_AMPeP_scaled= data_scaler(data=X_AMPeP, use_prescaled=True, scaler_type=scaler_math, scaler_save_path=scaler_loc)
        Y_AMPeP = data_AMPeP['Class']

        stack_out=data_MACREL[['Names', 'SEQ', 'Class']].copy()
    elif test==False:
        X_MACREL = data_MACREL.drop(columns=['Names', 'SEQ'])
        X_MACREL_scaled= data_scaler(data=X_MACREL, use_prescaled=True, scaler_type=scaler_math, scaler_save_path=scaler_loc)
        X_AMPeP = data_AMPeP.drop(columns=['Names', 'SEQ']) #prep for training
        X_AMPeP_scaled= data_scaler(data=X_AMPeP, use_prescaled=True, scaler_type=scaler_math, scaler_save_path=scaler_loc)
        stack_out=data_MACREL[['Names','SEQ']].copy()
    
    rf_ensemble_macrel=stack_helper(model_parent_path=RF_MACREL_path,X_data=X_MACREL,model_type='RF',feature_set='MACREL')
    nn_ensemble_macrel=stack_helper(model_parent_path=NN_MACREL_path,X_data=X_MACREL_scaled,type='NN',feature_set='MACREL')

    rf_ensemble_ampep=stack_helper(model_parent_path=RF_AMPeP_path,X_data=X_AMPeP,model_type='RF',feature_set='AMPeP')
    nn_ensemble_ampep=stack_helper(model_parent_path=NN_AMPeP_path,X_data=X_AMPeP_scaled,model_type='NN',feature_set='AMPeP')
    
    if test==True:
        stack_df_MACREL=pd.concat([X_MACREL.reset_index(drop=True),rf_ensemble_macrel.reset_index(drop=True), nn_ensemble_macrel.reset_index(drop=True),Y_MACREL.reset_index(drop=True)],axis=1)
        stack_df_AMPeP=pd.concat([X_AMPeP.reset_index(drop=True), rf_ensemble_ampep.reset_index(drop=True), nn_ensemble_ampep.reset_index(drop=True),Y_AMPeP.reset_index(drop=True)],axis=1)

        sMACREL_X=stack_df_MACREL.drop(columns='Class')
        sAMPeP_X=stack_df_AMPeP.drop(columns='Class')
    elif test==False:
        stack_df_MACREL=pd.concat([X_MACREL.reset_index(drop=True), rf_ensemble_macrel.reset_index(drop=True), nn_ensemble_macrel.reset_index(drop=True)],axis=1)
        stack_df_AMPeP=pd.concat([X_AMPeP.reset_index(drop=True), rf_ensemble_ampep.reset_index(drop=True), nn_ensemble_ampep.reset_index(drop=True)],axis=1)

        sMACREL_X=stack_df_MACREL.copy()
        sAMPeP_X=stack_df_AMPeP.copy()

    sAMPeP=joblib.load(sAMPeP_path)
    stack_out['sAMPeP_proba']=sAMPeP.predict_proba(sAMPeP_X)[:,1]
    stack_out['sAMPeP_prediction']=sAMPeP.predict(sAMPeP_X)
    sMACREL=joblib.load(sMACREL_path)
    stack_out['sMACREL_proba']=sMACREL.predict_proba(sMACREL_X)[:,1]
    stack_out['sMACREL_prediction']=sMACREL.predict(sMACREL_X)

    stack_out['Predictions']=stack_out[['sAMPeP_proba','sMACREL_proba']].mean(axis=1)
    
   
    return stack_out