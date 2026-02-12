import pandas as pd
from itertools import combinations
from collections import defaultdict
from sAMPC_helpers import excel_reader
from sklearn.metrics import cohen_kappa_score


def model_group_combinations(predictions_path: str, sheet_name: str, save_path: str):
    model_prediction=pd.DataFrame()
    model_prediction=excel_reader(file_path=predictions_path, sheet_name=sheet_name)
    column_names: list[str]=model_prediction.columns.tolist()
    model_groups=pd.DataFrame(list(combinations(column_names,2)), columns=['ColA', 'ColB'])
    model_groups.to_csv(f"{save_path}/Grouping_2.csv", index=False, header=False)

def fliess_kappa(df):
    temporary_df = pd.DataFrame()
    total_samples,total_raters = df.shape

    predictions_amp = (df>=0.5).sum(axis=1)
    predictions_namp = (df<0.5).sum(axis=1)
    temporary_df['AMP'] = predictions_amp
    temporary_df['NAMP'] = predictions_namp

    p_amp = ((temporary_df['AMP'].sum(axis=0))/(total_samples*total_raters))**2 #Summation of all the values in the col, division by total_raters*total_samples and then taking a square
    p_namp = ((temporary_df['NAMP'].sum(axis=0))/(total_samples*total_raters))**2
    random_guess = p_amp + p_namp

    extent_agreement = ((temporary_df[['AMP','NAMP']]**2).sum(axis=1)-total_raters)/(total_raters*(total_raters-1))
    extent_agreement_mean = extent_agreement.mean()

    kappa = (extent_agreement_mean - random_guess)/(1-random_guess)

    return kappa

def model_agreement(predictions_path: str, sheet_name: str, combination_path: str):
    
    model_data=excel_reader(file_path=predictions_path, sheet_name=sheet_name)
    model_groups=pd.read_csv(combination_path)
    model_pair_scores= defaultdict(dict)
    for index, row in model_groups.iterrows():
        kappa=cohen_kappa_score(model_data[row['ColA']], model_data[row['ColB']])
        model_pair=[row["ColA"], row["ColB"]]
        model_pair_scores[(row["ColA"], row["ColB"])]={
            'cohens_kappa_score':kappa,
            'model_pair_list': model_pair
        }
    
    best_pair_key: tuple =max(model_pair_scores, key= lambda x: model_pair_scores[x]['cohens_kappa_score'])
    pair: list[str]=model_pair_scores[best_pair_key]['model_pair_list']

    combination_df=model_data[pair]
    checking_list=(model_data.drop(columns=pair).drop(columns=['Class', 'SEQ','Names'])).columns.tolist()
    for i in range (10): #still working on it ik
        group_scores= defaultdict(dict)
        for model in checking_list:
            temp_df=pd.concat([combination_df, model_data[model]], axis=1)
            group=tuple(temp_df.columns.tolist())
            fliess_kappa_score=fliess_kappa()
            group_scores[group]={
                'fliess_kappa': fliess_kappa_score,
                'model_to_add': model
                }
        
        best_group_key: tuple = max(group_scores, key = lambda x: group_scores[x]['fliess_kappa'])
        pair.append(group_scores[best_group_key]['model_to_add'])

        combination_df = model_data[pair]
        checking_list.remove(group_scores[best_group_key]['model_to_add'])
    
    return combination_df.columns.tolist()
