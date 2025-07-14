#Function file for the sMACREL classifier, allows the user to train models on their own data or use pre-trained models
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import math
import peptides
import joblib
import os

def fasta_to_tsv(prodigal_out):
    proteins_df=pd.DataFrame()
    sequences=[]
    sequence_names=[]
    

    for record in SeqIO.parse(prodigal_out,'fasta'):
        if len(record.seq)<=255 and len(record.seq)>=5:
            sequence_names.append(f'{record.id}')
            sequences.append((str(record.seq).lstrip('M').rstrip('*')))
    
    proteins_df['Names']=sequence_names
    proteins_df['SEQ']=sequences

    return proteins_df

def excel_reader(file_path):
    features=pd.ExcelFile(file_path)
    features_sheet_name={}
    for features_sheet in features.sheet_names:
        features_sheet_name[features_sheet]=features.parse(features_sheet)
    features_df=features_sheet_name['Features']

    return features_df

#Feature extraction <still a work in progress>
def macrel_feature_extract(protein_data,save_path):
    
    Dataset_df=fasta_to_tsv(prodigal_out=protein_data)
    pep_seqs=Dataset_df['SEQ']
    
    CTDD_descriptors={'solvent accessibility':{'buried':{'A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'},
                                          'exposed':{'P', 'K', 'Q', 'E', 'N', 'D'},
                                          'intermediate':{'M', 'P', 'S', 'T', 'H', 'Y'}},
                      'free energy of transfer':{'fetclass1':{'I','L','V','W','A','G','T','M'},
                                                 'fetclass2':{'F','Y','S','Q','C','N'},
                                                 'fetclass3':{'P','H','K','E','D','R'}}
                      }
    global_decriptors={'tiny aar':{'A','C','G','S','T'},
                       'small amino acid residues':{'A','B','C','D','G','N','P','S','T','V'},
                       'aliphatic amino acid residues':{'A','I','L','V'},
                       'aromatic amino acid residues':{'F','H','W','Y'},
                       'nonpolar amino acid residues':{'A','C','F','G','I','L','M','P','V','W','Y'},
                       'polar amino acid residues':{'D','E','H','K','N','Q','R','S','T','Z'},
                       'charged amino acid residues':{'B','D','E','H','K','R','Z'},
                       'basic amino acid residues':{'H','K','R'},
                       'acidic amino acid residues':{'B','D','E','Z'}
                       }
    
    boman_index=[]
    aliphaticindex=[]
    instabilityindex=[]
    peptide_charge=[]
    isoelectricpoint=[]
    hydrophobicity=[]
    hydrophobicmoment=[]
    
    all_first_occurence_scores = {}
    all_global_scores={}
    for descriptor in CTDD_descriptors: #Populate the dictionary with empty lists and the names of all the classes per descriptor
        for class_name in CTDD_descriptors[descriptor]:
            key= f'{descriptor} {class_name}'
            all_first_occurence_scores[key]=[]
    
    for descriptor in global_decriptors:
        all_global_scores[f'{descriptor}']=[]
    
    for seq in pep_seqs:
        extracted_features={}
        extracted_global_features={}
        first_occurence_scores={}
        global_scores={}
        
        for descriptor in CTDD_descriptors.keys():
            extracted_features[descriptor]={}
            for class_name in CTDD_descriptors[descriptor].keys():
                class_values=CTDD_descriptors[descriptor][class_name]
                position_CTDD_classes=[]
                for i,aa in enumerate(seq):
                    if aa in class_values:
                        position_CTDD_classes.append(i+1)
                if position_CTDD_classes:
                    extracted_features[descriptor][class_name]=position_CTDD_classes
                else:
                    extracted_features[descriptor][class_name]=[]
                
        
        for descriptor in extracted_features.keys():
            first_occurence_scores[descriptor]={}
            for class_name in extracted_features[descriptor].keys():
                if len(extracted_features[descriptor][class_name])!=0:
                    first=min(extracted_features[descriptor][class_name])
                    first_score=(first/len(seq))*100
                    first_occurence_scores[descriptor][class_name]=first_score
                else:
                    first=0
                    first_score=(first/len(seq))*100
                    first_occurence_scores[descriptor][class_name]=first_score
        
        
        for descriptor in global_decriptors.keys(): #Getting postions
            extracted_global_features[descriptor]={}
            descriptor_values=global_decriptors[descriptor]
            position_global_descriptors=[]
            for i,aa in enumerate(seq):
                if aa in descriptor_values:
                    position_global_descriptors.append(i+1)
            if position_global_descriptors:
                extracted_global_features[descriptor]=position_global_descriptors
            else:
                extracted_global_features[descriptor]=[]     
        
        
        for descriptor in extracted_global_features.keys(): #The actual math for feature extraction
            global_scores[descriptor]={}
            if len(extracted_global_features[descriptor])!=0:
                global_score=len(extracted_global_features[descriptor])/len(seq)
                global_scores[descriptor]=global_score
            else:
                global_score=0
                global_scores[descriptor]=global_score
        
                
        #Appending the scores to a list to use for adding to the Dataframe        
        for descriptor in first_occurence_scores.keys():
            for class_name in first_occurence_scores[descriptor].keys():
                first_sc=first_occurence_scores[descriptor][class_name]
                all_first_occurence_scores[f'{descriptor} {class_name}'].append(first_sc)
        for descriptor in global_scores.keys():
                global_sc=global_scores[descriptor]
                all_global_scores[descriptor].append(global_sc)
                
        peptide_charge.append(ProteinAnalysis(seq).charge_at_pH(pH=7.0))
        instabilityindex.append(peptides.Peptide(seq).instability_index())
        hydrophobicity.append(peptides.Peptide(seq).hydrophobicity())
        hydrophobicmoment.append(peptides.Peptide(seq).hydrophobic_moment(angle=100,window=11))
        boman_index.append(peptides.Peptide(seq).boman())
        aliphaticindex.append(peptides.Peptide(seq).aliphatic_index())
        isoelectricpoint.append(ProteinAnalysis(seq).isoelectric_point())
        
    for column, values in all_first_occurence_scores.items():
        Dataset_df[column] = values
    for column, values in all_global_scores.items():
        Dataset_df[column] = values
    
    Dataset_df['peptide charge']=peptide_charge
    Dataset_df['isoelectric point']=isoelectricpoint
    Dataset_df['instability index']=instabilityindex
    Dataset_df['boman index']=boman_index
    Dataset_df['hydrophobicity']=hydrophobicity
    
    Dataset_df['aliphatic index']=aliphaticindex
    Dataset_df['hydrophobic moment']=hydrophobicmoment
    
    with pd.ExcelWriter(f"{save_path}/Extracted_MACREL_Features.xlsx", engine='openpyxl', mode='w') as feature_writer:
        Dataset_df.to_excel(feature_writer,index=False,sheet_name='Features')


def ampep_feature_extract(protein_data,save_path):
    
    Dataset_df=fasta_to_tsv(prodigal_out=protein_data)
    pep_seqs=Dataset_df['SEQ']

    descriptors={'hydrophobicity':{'polar':{'R', 'K', 'E', 'D', 'Q', 'N'},
                                   'neutral':{'G', 'A', 'S', 'T', 'P', 'H', 'Y'},
                                   'hydrophobic':{'C', 'L', 'V', 'I', 'M', 'F', 'W'}},
                 'Normalized van der Waals volume':{'vdwclass1':{'G', 'A', 'S', 'T', 'P', 'D'},
                                         'vdwclass2':{'N', 'V', 'E', 'Q', 'I', 'L'},
                                         'vdwclass3':{'M', 'H', 'K', 'F', 'R', 'Y', 'W'}},
                 'polarity':{'pclass1':{'L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'},
                             'pclass2':{'P', 'A', 'T', 'G', 'S'},
                             'pclass3':{'H', 'Q', 'R', 'K', 'N', 'E', 'D'}},
                 'polarizability':{'plclass1':{'G', 'A', 'S', 'D', 'T'},
                                   'plclass2':{'C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'},
                                   'plclass3':{'K', 'M', 'H', 'F', 'R', 'Y', 'W'}},
                 'charge':{'positive':{'K', 'R'}, 
                           'neutral':{'A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'},
                           'negative':{'D', 'E'}},
                 'secondary structure':{'helix':{'E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'},
                                        'strand':{'V', 'I', 'Y', 'C', 'W', 'F', 'T'},
                                        'coil':{'G', 'N', 'P', 'S', 'D'}},
                 'solvent accessibility':{'buried':{'A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'},
                                          'exposed':{'P', 'K', 'Q', 'E', 'N', 'D'},
                                          'intermediate':{'M', 'P', 'S', 'T', 'H', 'Y'}}
                 }
    
    #Make dictionaries to store all scores before appending them to the Dataframe
    all_first_occurence_scores = {} 
    all_scores_at_25 = {}
    all_scores_at_50 = {}
    all_scores_at_75 = {}
    all_scores_at_100 = {}
    
    
    for descriptor in descriptors: #Populate the dictionary with empty lists and the names of all the classes per descriptor
        for class_name in descriptors[descriptor]:
            key= f'{descriptor} {class_name}'
            all_first_occurence_scores[key]=[]
            all_scores_at_25[f'{key} 25% score']=[]
            all_scores_at_50[f'{key} 50% score']=[]
            all_scores_at_75[f'{key} 75% score']=[]
            all_scores_at_100[f'{key} 100% score']=[]
    
    
    for seq in pep_seqs:
        extracted_features={}
        first_occurence_scores={}
        scores_at_25_percent={}
        scores_at_50_percent={}
        scores_at_75_percent={}
        scores_at_100_percent={}
        
        for descriptor in descriptors.keys():
            extracted_features[descriptor]={}
            for class_name in descriptors[descriptor].keys():
                class_values=descriptors[descriptor][class_name]
                position_classes=[]
                for i,polar_aa in enumerate(seq):
                    if polar_aa in class_values:
                        position_classes.append(i+1)
                if position_classes:
                    extracted_features[descriptor][class_name]=position_classes
                else:
                    extracted_features[descriptor][class_name]=[]
                    
        
        for descriptor in extracted_features.keys():
            first_occurence_scores[descriptor]={}
            scores_at_25_percent[descriptor]={}
            scores_at_50_percent[descriptor]={}
            scores_at_75_percent[descriptor]={}
            scores_at_100_percent[descriptor]={}
            for class_name in extracted_features[descriptor].keys():
                if len(extracted_features[descriptor][class_name])!=0:
                    first=min(extracted_features[descriptor][class_name])
                    first_score=(first/len(seq))*100
                    score_at_25=((extracted_features[descriptor][class_name][math.floor(len(extracted_features[descriptor][class_name])*(25/100))-1])/len(seq))*100
                    score_at_50=((extracted_features[descriptor][class_name][math.floor(len(extracted_features[descriptor][class_name])*(50/100))-1])/len(seq))*100
                    score_at_75=((extracted_features[descriptor][class_name][math.floor(len(extracted_features[descriptor][class_name])*(75/100))-1])/len(seq))*100
                    score_at_100=((extracted_features[descriptor][class_name][math.floor(len(extracted_features[descriptor][class_name])*(100/100))-1])/len(seq))*100
                    
                    first_occurence_scores[descriptor][class_name]=first_score
                    scores_at_25_percent[descriptor][class_name]=score_at_25
                    scores_at_50_percent[descriptor][class_name]=score_at_50
                    scores_at_75_percent[descriptor][class_name]=score_at_75
                    scores_at_100_percent[descriptor][class_name]=score_at_100
                else:
                    first=0
                    first_score=(first/len(seq))*100
                    score_at_25=0
                    score_at_50=0
                    score_at_75=0
                    score_at_100=0
                    
                    first_occurence_scores[descriptor][class_name]=first_score
                    scores_at_25_percent[descriptor][class_name]=score_at_25
                    scores_at_50_percent[descriptor][class_name]=score_at_50
                    scores_at_75_percent[descriptor][class_name]=score_at_75
                    scores_at_100_percent[descriptor][class_name]=score_at_100
                
        
        for descriptor in first_occurence_scores.keys():
            for class_name in first_occurence_scores[descriptor].keys():
                first_sc=first_occurence_scores[descriptor][class_name]
                all_first_occurence_scores[f'{descriptor} {class_name}'].append(first_sc)
        for descriptor in scores_at_25_percent.keys():
            for class_name in scores_at_25_percent[descriptor].keys():
                scores25=scores_at_25_percent[descriptor][class_name]
                all_scores_at_25[f'{descriptor} {class_name} 25% score'].append(scores25)
        for descriptor in scores_at_50_percent.keys():
            for class_name in scores_at_50_percent[descriptor].keys():
                scores50=scores_at_50_percent[descriptor][class_name]
                all_scores_at_50[f'{descriptor} {class_name} 50% score'].append(scores50)
        for descriptor in scores_at_75_percent.keys():
            for class_name in scores_at_75_percent[descriptor].keys():
                scores75=scores_at_75_percent[descriptor][class_name]
                all_scores_at_75[f'{descriptor} {class_name} 75% score'].append(scores75)
        for descriptor in scores_at_100_percent.keys():
            for class_name in scores_at_100_percent[descriptor].keys():
                scores100=scores_at_100_percent[descriptor][class_name]
                all_scores_at_100[f'{descriptor} {class_name} 100% score'].append(scores100)
    
    
    for column, values in all_first_occurence_scores.items():
        Dataset_df[column] = values
    for column, values in all_scores_at_25.items():
        Dataset_df[column] = values
    for column, values in all_scores_at_50.items():
        Dataset_df[column] = values
    for column, values in all_scores_at_75.items():
        Dataset_df[column] = values
    for column, values in all_scores_at_100.items():
        Dataset_df[column] = values
    
    with pd.ExcelWriter(f"{save_path}/Extracted_AMPeP_Features.xlsx", engine='openpyxl', mode='w') as feature_writer:
        Dataset_df.to_excel(feature_writer,index=False,sheet_name='Features')


#User Model training options
def user_data_training(feature_extracted_data,model_type,oversampling,models_parent_dir,
                       user_pn_ratios=None,user_NN_layers=None, user_testing_data=None):
    
    if model_type=='all':
        os.makedirs(f'{models_parent_dir}/RF',exist_ok=True)
        os.makedirs(f'{models_parent_dir}/LogReg',exist_ok=True)
        os.makedirs(f'{models_parent_dir}/NN',exist_ok=True)
    else:
        os.makedirs(f'{models_parent_dir}/{model_type}')
    
    Amp_df = feature_extracted_data[feature_extracted_data['Type (Binary)'] == 1]#Load class 1 from user training data
    Nonamp_df = feature_extracted_data[feature_extracted_data['Type (Binary)'] == 0]#Load class 0 from user training data
    
    if user_testing_data==None:
        test_amp = Amp_df.sample(n=1250, random_state=1300)
        test_amp_2 = Amp_df.sample(n=1750, random_state=3789541)
        test_nonamp = Nonamp_df.sample(n=3000, random_state=7348908)
        test_bench = pd.concat([test_amp, test_nonamp, test_amp_2])
        X_testbench = test_bench.drop(columns=['Type (Binary)', 'SEQ', 'NAME', 'Type'])
        Y_testbench = test_bench['Type (Binary)']
    else:
        X_testbench = user_testing_data.drop(columns=['Type (Binary)', 'SEQ', 'NAME', 'Type'])
        Y_testbench = user_testing_data['Type (Binary)']
    
    all_metrics = []
         #Allows the user to check model performances of RF, LogReg and NN models with their own variable set
    if user_pn_ratios==None: #Uses default p/n ratios
        data_ratios_list=[(1,1),(1,2),(1,5),(1,10),(1,15),(1,20),(1,25),(1,30),(1,35),(1,40),(1,45),(1,50)]
    else: #Allows the user to use their own list of p/n ratios
        data_ratios_list=user_pn_ratios
    if user_NN_layers==None: #Uses default hidden layer size
        hiddenlayers_size_list=[(10,),(10,1),(10,2),(10,5),(10,10),
                                (20,),(20,1),(20,2),(20,5),(20,10),(20,15),(20,20),
                                (30,),(30,1),(30,2),(30,5),(30,10),(30,15),(30,20),(30,25),(30,30),
                                (40,),(40,1),(40,2),(40,5),(40,10),(40,15),(40,20),(40,25),(40,30),(40,35),(40,40),
                                (50,),(50,1),(50,2),(50,5),(50,10),(50,15),(50,20),(50,25),(50,30),(50,35),(50,40),(50,45),(50,50)]
    else: #Allows the user to use their own list of Hidden Layer Size
        hiddenlayers_size_list=user_NN_layers
        
    for ratio in data_ratios_list:
        if oversampling==True:
            total_amps=3000
            total_nonamps=ratio[1]*total_amps
            train_amp = Amp_df.sample(n=1250, random_state=600)  # Sample 1250 AMPs
            train_amp_2 = Amp_df.sample(n=1750, random_state=1234567)  # Sample 1750 AMPs 
            train_nonamp = Nonamp_df.sample(n=total_nonamps, random_state=54567653)# Sample {total_nonamps} Non-AMPs
            training_df = pd.concat([train_amp, train_nonamp, train_amp_2])# concatenate all the sampling
        elif oversampling==False:
            total_amps=2000
            total_nonamps=ratio[1]*total_amps
            train_amp = Amp_df.sample(n=total_amps, random_state=1234567)  # Sample 1250 AMPs
            train_nonamp = Nonamp_df.sample(n=total_nonamps, random_state=54567653)# Sample {total_nonamps} Non-AMPs
            training_df = pd.concat([train_amp, train_nonamp])
                    
        X_train = training_df.drop(columns=['Type (Binary)', 'SEQ', 'NAME', 'Type'])
        Y_train = training_df['Type (Binary)']
            
            #Model training
        if model_type=='RF':
            model = RandomForestClassifier(criterion='gini', max_depth=None, random_state=0)
            model.fit(X_train, Y_train)
            joblib.dump(model, f'{models_parent_dir}/{model_type}/RF_pn_ratio_{total_amps}_{total_nonamps}_oversampling_{oversampling}_v1.joblib')
                #Testing RF MODELS
            predictions = model.predict(X_testbench)
            report = classification_report(Y_testbench, predictions, output_dict=True)
            metrics = {f"p/n ratio {model_type}": total_amps/total_nonamps,
                        "accuracy": report["accuracy"],
                        "precision class 1": report["1"]["precision"],
                        "recall class 1": report["1"]["recall"],
                        "precision class 0": report["0"]["precision"],
                        "recall class 0": report["0"]["recall"]}
            all_metrics.append(metrics)
        elif model_type=='LogReg':
            model = LogisticRegression(penalty='l2',random_state=42, max_iter=5000, solver='liblinear')
            model.fit(X_train, Y_train)
            joblib.dump(model, f'{models_parent_dir}/{model_type}/LogReg_pn_ratio_{total_amps}_{total_nonamps}_oversampling_{oversampling}_v1.joblib')
            # Testing LogReg MODELS
            predictions = model.predict(X_testbench)
            report = classification_report(Y_testbench, predictions, output_dict=True)
            metrics = {f"p/n ratio {model_type}": total_amps/total_nonamps,
                        "accuracy": report["accuracy"],
                        "precision class 1": report["1"]["precision"],
                        "recall class 1": report["1"]["recall"],
                        "precision class 0": report["0"]["precision"],
                        "recall class 0": report["0"]["recall"]}
            all_metrics.append(metrics)
        elif model_type=='NN':
                    #We iterated over the hiddenlayersize as well to check if affects model performance
            for layer_size in hiddenlayers_size_list:
                model=MLPClassifier(hidden_layer_sizes=layer_size, solver='adam', alpha=1e-10, verbose=False, max_iter=30000)
                model.fit(X_train,Y_train)
                joblib.dump(model, f'{models_parent_dir}/{model_type}/NN_pn_ratio_{total_amps}_{total_nonamps}_layersize_{layer_size}_oversampling_{oversampling}_v1.joblib')
                # Testing NN MODELS
                predictions = model.predict(X_testbench)
                report = classification_report(Y_testbench, predictions, output_dict=True, zero_division=1)
                metrics = {'Layer Size':f'Neural Network {layer_size}',
                            f"p/n ratio {model_type}": total_amps/total_nonamps,
                            "accuracy": report["accuracy"],
                            "precision class 1": report["1"]["precision"],
                            "recall class 1": report["1"]["recall"],
                            "precision class 0": report["0"]["precision"],
                            "recall class 0": report["0"]["recall"]}
                all_metrics.append(metrics)

    scores_df = pd.DataFrame(all_metrics)
    if model_type=='RF':
        with pd.ExcelWriter(f"{models_parent_dir}/{model_type}/Metrics_RF_oversampled_{oversampling}.xlsx", engine='openpyxl', mode='w') as prediction_writer:
            scores_df.to_excel(prediction_writer,index=False,sheet_name='RF_Metrics')
    elif model_type=='LogReg':
        with pd.ExcelWriter(f"{models_parent_dir}/{model_type}/Metrics_LogReg_oversampled_{oversampling}.xlsx", engine='openpyxl', mode='w') as prediction_writer:
            scores_df.to_excel(prediction_writer,index=False,sheet_name='LogReg_Metrics')
    elif model_type=='NN':
        with pd.ExcelWriter(f"{models_parent_dir}/{model_type}/Metrics_NN_oversampled_{oversampling}.xlsx", engine='openpyxl', mode='w') as prediction_writer:
            scores_df.to_excel(prediction_writer,index=False,sheet_name='NN_Metrics')
            

def stack_helper(model_parent_path,X_data,type, feature_set): #Helper training function for the stacking training function
    ensemble_predictions=pd.DataFrame()
    for model in sorted(os.listdir(model_parent_path)):
        if model.endswith('.joblib'):
            model_path=os.path.join(model_parent_path, model)
            model=joblib.load(model_path)
            ensemble_predictions[f'{model}']=model.predict(X_data)
        ensemble_predictions[f'{type}_{feature_set}_avg']=ensemble_predictions.mean(axis=1)
    return ensemble_predictions



def stack_train(RF_path, NN_path, stack_save_path, stack_train_data, feature_set):
    os.makedirs(f'{RF_path}/RF_votes',exist_ok=True)
    #os.makedirs(f'{LogReg_path}/LogReg_votes',exist_ok=True)
    os.makedirs(f'{NN_path}/NN_votes',exist_ok=True)
    
    sMACREL_training_df=pd.DataFrame()
    
    Amp_df = stack_train_data[stack_train_data['Type (Binary)'] == 1]#Load class 1 from user training data
    Nonamp_df = stack_train_data[stack_train_data['Type (Binary)'] == 0]#Load class 0 from user training data
    
    train_amp=Amp_df.sample(n=1250,random_state=1300) # Sample class1
    train_amp_2=Amp_df.sample(n=1750,random_state=3789541)# Sample class1
    train_nonamp=Nonamp_df.sample(n=75000, random_state=7348908)# Sample class0
    training_df=pd.concat([train_amp,train_nonamp,train_amp_2])#Combine sampled data to make a training df
    X_train = training_df.drop(columns=['Type (Binary)', 'SEQ', 'NAME', 'Type']) #prep for training
    Y_train = training_df['Type (Binary)']#Prep for training
    
    rf_ensemble_predictions=stack_helper(model_parent_path=RF_path,X_data=X_train,type='RF',feature_set=feature_set)
    #logreg_ensemble_predictions=sMACREL_helper(model_parent_path=LogReg_path,X_data=X_train,type='LogReg')
    nn_ensemble_predictions=stack_helper(model_parent_path=NN_path,X_data=X_train,type='NN',feature_set=feature_set)
    
    sMACREL_training_df=pd.concat([X_train.reset_index(drop=True),rf_ensemble_predictions[f'RF_{feature_set}_avg'],nn_ensemble_predictions[f'NN_{feature_set}_avg'],Y_train.reset_index(drop=True)],axis=1)
    sMACREL_X=sMACREL_training_df.drop(columns='Type (Binary)')
    sMACREL_Y=sMACREL_training_df['Type (Binary)']
    smacrel_model= RandomForestClassifier(criterion='gini', max_depth=None, random_state=0)
    smacrel_model.fit(sMACREL_X, sMACREL_Y)
    joblib.dump(smacrel_model, f'{stack_save_path}/stack_{feature_set}.joblib')


def stack_predictions(RF_MACREL_path, RF_AMPeP_path, NN_MACREL_path, NN_AMPeP_path, sAMPeP_path, sMACREL_path, MACREL_features, AMPeP_features, test):

    stack_df_AMPeP=pd.DataFrame()#make a df to load MetaModel Features
    stack_df_MACREL=pd.DataFrame()#make a df to load MetaModel Features
    stack_out=pd.DataFrame()#make a df to store MetaModel predictions
    
    data_MACREL=excel_reader(MACREL_features)
    data_AMPeP=excel_reader(AMPeP_features)

    if test==True:
        X_MACREL = data_MACREL.drop(columns=['Names', 'SEQ', 'Class']) #prep for training
        Y_MACREL = data_MACREL['Class']#Prep for training
        X_AMPeP = data_AMPeP.drop(columns=['Names', 'SEQ', 'Class']) #prep for training
        Y_AMPeP = data_AMPeP['Class']#Prep for training

        stack_out=data_MACREL[['Names', 'SEQ', 'Class']].copy()
    elif test==False:
        X_MACREL = data_MACREL.drop(columns=['Names', 'SEQ']) #prep for training
        #Y_MACREL = data_MACREL['Class']#Prep for training
        X_AMPeP = data_AMPeP.drop(columns=['Names', 'SEQ']) #prep for training
        #Y_AMPeP = data_AMPeP['Class']#Prep for training
        stack_out=data_MACREL[['Names','SEQ']].copy()
    
    rf_ensemble_macrel=stack_helper(model_parent_path=RF_MACREL_path,X_data=X_MACREL,type='RF',feature_set='MACREL')
    nn_ensemble_macrel=stack_helper(model_parent_path=NN_MACREL_path,X_data=X_MACREL,type='NN',feature_set='MACREL')

    rf_ensemble_ampep=stack_helper(model_parent_path=RF_AMPeP_path,X_data=X_AMPeP,type='RF',feature_set='AMPeP')
    nn_ensemble_ampep=stack_helper(model_parent_path=NN_AMPeP_path,X_data=X_AMPeP,type='NN',feature_set='AMPeP')
    
    if test==True:
        stack_df_MACREL=pd.concat([X_MACREL.reset_index(drop=True), rf_ensemble_macrel['RF_MACREL_avg'], nn_ensemble_macrel['NN_MACREL_avg'],Y_MACREL.reset_index(drop=True)],axis=1)
        stack_df_AMPeP=pd.concat([X_AMPeP.reset_index(drop=True), rf_ensemble_ampep['RF_AMPeP_avg'], nn_ensemble_ampep['NN_AMPeP_avg'],Y_AMPeP.reset_index(drop=True)],axis=1)

        sMACREL_X=stack_df_MACREL.drop(columns='Class')
        sMACREL_Y=stack_df_MACREL['Class']
        sAMPeP_X=stack_df_AMPeP.drop(columns='Class')
        sAMPeP_Y=stack_df_AMPeP['Class']
    elif test==False:
        stack_df_MACREL=pd.concat([X_MACREL.reset_index(drop=True), rf_ensemble_macrel['RF_MACREL_avg'], nn_ensemble_macrel['NN_MACREL_avg']],axis=1)
        stack_df_AMPeP=pd.concat([X_AMPeP.reset_index(drop=True), rf_ensemble_ampep['RF_AMPeP_avg'], nn_ensemble_ampep['NN_AMPeP_avg']],axis=1)

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