# Change Log

---
## Changed Base Model Selection
### Date: 12/02/2026

In the previous version, the base models were manually selected based on their metrics. This manual way has now been replaced 
by `Base_model_slection.py`. The models are first filtered to only keep the best performing ones [MCC >0.9, Precision >0.8, and Recall >0.8].
Model pairs are then formed using `itertools.combinations()` and the agreement between the model pairs is measured with `skelearn.metrics.cohen_kappa_score`.
A group of models for the ensemble is then built based on the best pairs by calculating agreement using Fleiss's Kappa, following a greedy search approach.
As of now, base model selection can find the best 10 models by calculating agreement between them.

**Please note: I am still torn between selecting models based on agreement or disagreement. I do not know what is the best approach.**

---
## Major Logical Bug in helper function
### Date: 12/02/2026

The function stack_helper() in theory takes model predictions appends them to a DataFrame and then after appending all the  predictions
takes an average of the predictions.

However in practice, the averaging step takes place inside the for loop. Which makes the function take an average every single time a model
prediction is appended to the DataFrame

*Old*
```python
def stack_helper(model_parent_path,X_data,type, feature_set):
    ensemble_predictions=pd.DataFrame()
    for model in sorted(os.listdir(model_parent_path)):
        if model.endswith('.joblib'):
            model_path=os.path.join(model_parent_path, model)
            model=joblib.load(model_path)
            ensemble_predictions[f'{model}']=model.predict(X_data)
        ensemble_predictions[f'{type}_{feature_set}_avg']=ensemble_predictions.mean(axis=1) #<- notice how this line is inside the for loop
    return ensemble_predictions
```

*What was supposed to happen*
```python
def stack_helper(model_parent_path,X_data,type, feature_set):
    ensemble_predictions=pd.DataFrame()
    for model in sorted(os.listdir(model_parent_path)):
        if model.endswith('.joblib'):
            model_path=os.path.join(model_parent_path, model)
            model=joblib.load(model_path)
            ensemble_predictions[f'{model}']=model.predict(X_data)

    ensemble_predictions[f'{type}_{feature_set}_avg']=ensemble_predictions.mean(axis=1) #<- notice how this line is outside the loop
    return ensemble_predictions
```

The averaging step has now been removed from the stack_helper() function on the development branch for the time being.

---
## Neural Network Scaling
### Date: 12/02/2026

The Neural Networks on the main branch are not scaled. I had noticed that the scaled Neural Networks performed worse compared to the unscaled ones.
This was apparently due to the different p/n ratios used during training. The scale was applied to all the p/n ratios.

The scaling logic has been implemented on the development branch. The Neural Networks will be saved with their own scaler packaged in the same joblib
file. Additionally, the data_scaler() function in the `sAMPC_helpers.py` file also saves the scalers to a specified path for future work.

**Please note: I am still working on implementing this to all the codes on the development branch it is very unstable**
