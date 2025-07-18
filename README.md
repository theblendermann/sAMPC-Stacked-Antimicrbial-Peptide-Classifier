# A bit of background
## Antimicrobial Peptides
Antimicrobial Peptides (AMPs) are short peptide chains, usually 10 to 60 amino acids in length. AMPs can originate from plants, mammals, amphibians, insects, fungi, and bacteria Huan et al. (2020). Antimicrobial peptides are non-enzymatic, i.e. their mechanism of action is usually physical. Most AMPs have a positive charge on them (cationic nature) which allows them to disrupt bacterial membranes that are usually negatively charged or anionic in nature causing the cytoplasmic components to leak out of the cell Zhang et al. (2021), Benfield and Henriques (2020), Kumar et al. (2018).

Antimicrobial peptides are not highly conserved, database search tools like BLASTp often fail to find novel AMPs. Computational methods, specifically various machine learning algorithms have been used to identify new AMPs. Even though AMPs do not have a conserved amino acid sequence, they do still exhibit physicochemical properties that are characteristic to them. Various prediction tools use these properties to distinguish between AMPs and Non-AMPs. In 2009, CAMP was set up with multiple prediction tools, a Random Forest, a Support Vector Machine and a Discriminant Analysis model that were trained on 64 physicochemical features Thomas et al. (2009). In 2018, AMPeP was made public, it is a Random Forest Algorithm that was trained on 105 features the exact features will be discussed later as the AMPeP features have been used here Bhadra et al. (2018). DeepAMPeP30 is a neural network trained to predict short AMPs (<=30). In 2020, the MACREL pipeline was made available by Santos-Júnior et al. (2020), the classifier in the MACREL pipeline was trained on 22 features, two being repurposed from AMPeP with 3 being novel. While the MACREL classifier has a high precision, it did come at the cost of recall (how many AMPs it was able to detect).

### AMPeP Features
The AMPeP feature set described by Bhadra et al. (2018) can be broadly divided into 7 major descriptors, these descriptors divide amino acids into three classes, class1, class2 and class3. A score is calculated based on the position of amino acids at the first position, 25% residue, 50% residue, 75% residue, and 100% residue of a class in a peptide chain. These descriptors are
1.	**Hydrophobicity:** Divides the amino acids into three classes, class1 for polar amino acids, class2 for neutral amino acids, and class3 for hydrophobic amino acids. 
2.	**Normalized van der Waals Volume:** In this case it is the volume occupied by a molecule that is enclosed within a van der Waals surface. These are experimentally derived, but for this study they were obtained from Table 7 provided by Bhadra et al. (2018). The descriptor divides the amino acids into three classes based on the volume occupied, class1 for amino acids with a volume in the range 0 - 2.78, class2 for the range 2.95 - 4.0, and class3 for the range 4.03 - 8.08.
3.	**Polarity:** It is the distribution of charges present on an amino acid, higher polarity implies that the charges are unevenly distributed on the molecule giving it a partial positive and partial negative charge. This distributor divides the amino acids into three classes based on the polarity values, class1 has amino acids with a polarity values ranging from 4.9 - 6.2, class2 for the range 8.0 - 9.2, and class3 for the range 10.4 - 13. As one would expect these descriptors strongly correlate with the Hydrophobicity descriptors (see Figure 1)
4.	**Polarizability:** It is defined as the ability of a molecule to develop a dipole in the presence of an electric field. Higher the value, the easier it is to induce a dipole. Class1 for amino acids in the range 0 - 0.108, class2 for the range 0.128 - 0.186, and class3 for the range 0.219 - 0.409.
5.	**Charge:** Divides the amino acids based on their charge, class1 for positively charged amino acids, class2 for neutral amino acids and class3 for negatively charged amino acids
6.	**Solvent Accessibility:** Divides the amino acids based on their position in a peptide, class1 for amino acids that are buried, class2 for amino acids that are exposed, class3 for amino acids that are partially buried.
7.	**Secondary Structure:** Divides the amino acids based on their involvement in peptide secondary structure, class1 for amino acids involved in a helical structure, class2 for strands, and class3 for coils.

![Figure 1](Correlograms/AMPeP_Feature_Correlogram.svg)
Figure 1: Correlogram of AMPeP Features

### MACREL Classifier Features
The MACREL Classifier feature set described by Santos-Júnior et al. (2020) can be divided into two main groups, Global descriptors and Local descriptors. The local descriptors are calculated based on the first occurrence of a residue in a peptide sequence.
1.	**Global Descriptors**
-	**Charge and Solubility:** Charge and isoelectric point of a peptide, calculated using the peptides library in python.
-	**Indexes:**
    -	**Instability Index:** Provides an estimate on protein stability. Calculated using the peptides library.
    -	**Boman Index:** Estimates the potential of a peptide to bind to a membrane or another protein. Calculated using the peptides library.
    -	**Aliphatic Index:** It is the volume occupied by aliphatic amino acids like Alanine, Isoleucine, Leucine, and Valine. Calculated using the peptides library.
-	**Hydrophobicity and Hydrophobic Moment:** Hydrophobic Moment is an estimate of amphiphilicity of a peptide and hydrophobicity is how soluble the peptide is in water.
-	**Amino Acid Composition:** These descriptors describe the physico-chemical properties of the peptides. These properties include, charge, and type of amino acid residue.
2.	**Local Descriptors**
-	**Free Energy of Transition:** It is the estimated change in free energy when a peptide moves from a random coil in an aqueous environment to an organized structure in a lipid environment. It can help estimate the likelihood of each amino acid to promote protein folding. FET divides the amino acids into three classes, class1 low energy change, class2 for moderate energy change, and class3 for high energy change.
-	**Solvent Accessibility:** As defined by Bhadra et al. (2018)

### References
- Huan, Y., Kong, Q., Mou, H., & Yi, H. (2020). Antimicrobial Peptides: Classification, Design, Application and Research Progress in Multiple Fields. Frontiers in Microbiology, 11. https://doi.org/10.3389/fmicb.2020.582779
- Zhang, Q.-Y., Yan, Z.-B., Meng, Y.-M., Hong, X.-Y., Shao, G., Ma, J.-J., Cheng, X.-R., Liu, J., Kang, J., & Fu, C.-Y. (2021). Antimicrobial peptides: mechanism of action, activity and clinical potential. Military Medical Research, 8(48), 48. https://doi.org/10.1186/s40779-021-00343-2
- Kumar, P., Kizhakkedathu, J., & Straus, S. (2018). Antimicrobial Peptides: Diversity, Mechanism of Action and Strategies to Improve the Activity and Biocompatibility In Vivo. Biomolecules, 8(1), 4. https://doi.org/10.3390/biom8010004
- Benfield, A. H., & Henriques, S. T. (2020). Mode-of-Action of Antimicrobial Peptides: Membrane Disruption vs. Intracellular Mechanisms. Frontiers in Medical Technology, 2. https://doi.org/10.3389/fmedt.2020.610997
- Bhadra, P., Yan, J., Li, J., Fong, S., & Siu, S. W. I. (2018). AmPEP: Sequence-based prediction of antimicrobial peptides using distribution patterns of amino acid properties and random forest. Scientific Reports, 8(1). https://doi.org/10.1038/s41598-018-19752-w
- Santos-Júnior, C. D., Pan, S., Zhao, X.-M., & Coelho, L. P. (2020). Macrel: antimicrobial peptide screening in genomes and metagenomes. PeerJ, 8, e10555. https://doi.org/10.7717/peerj.10555

# The sAMPC function file and Installation
Before we proceed it must be made clear that a lot of work still needs to be done on sAMPC and the code is still very experimental.

### sAMPC: Stacked Antimicrobial Peptide Classifier
sAMPC is an ensemble model based on the feature sets described in the AMPeP and MACREL papers. The sAMPC function file contains 5 main functions, they are `macrel_feature_extract()`, `ampep_feature_extract()`, `user_data_training()`, `stack_train()`, and `stack_predictions()`.

![Figure 2](Model_images/Model_Schema.svg)
Figure 2: sAMPC Workflow

Figure 2, shown above shows the general work flow of sAMPC. First the protein sequences are stored in a `.fasta` file, for genomic dna Prodigal is first used to identify the protein sequences and then stored in a `.fasta` file. The proteins then need to be passed through the `ampep_feature_extract()` and `macrel_feature_extract()` functions. These features are used as inputs for the base models. Currently, sAMPC has 12 Neural Networks and 7 Random Forest Classifiers that predict using the MACREL feature set and 8 Neural Networks and 5 Random Forest Classifiers that predict using the AMPeP Feature set. All the base models are have different `p:n ratios` and `hidden layer sizes`. Once the base models give their predictions `predict()`, they are averaged. For example, all the Random Forest predictions for the MACREL Feature set are averaged and the results are stored in a column `RF_MACREL_avg`. The averages and the extracted features are then used as inputs for two stacks, `sMACREL` which is a random forest trained on the MACREL feature set & averages and `sAMPeP` which is a random forest classifier trained on the AMPeP feature set & averages. The predictions of these two stacks `predict_proba()` is then averaged to give the final prediction. The prediction sets are done by a single function `stack_predictions()`.

### Functions of sAMPC
**Feature Extraction**

`ampep_feature_extract()` & `macrel_feature_extract()`
-    protein_data (str): Path to your fasta file. The file extention must be `.fasta`
-    save_path (str): Path to where you want to save the extracted features. The output will be in `.xlsx`

`stack_predictions()`
-    RF_MACREL_path (str): Path to the base Random Forest Classifier trained on the MACREL Features. 
-    RF_AMPeP_path (str): Path to the base Random Forest Classifiers trained on the AMPeP Features.
-    NN_MACREL_path (str): Path to the base Neural Networks trained on MACREL features.
-    NN_AMPeP_path (str): Path to the base Neural Networks trained on AMPeP features.
-    sAMPeP_path (str): Path to the sAMPeP classifier.
-    sMACREL_path (str): Path to the sMACREL classifier.
-    MACREL_features (str): Path to MACREL extracted features. The input must be in `.xlsx` and the column names **should NOT** be tampered with, the code will break if done so.
-    AMPeP_features (str): Path to AMPeP extracted features. The input must be in `.xlsx` and the column names **should NOT** be tampered with, the code will break if done so.
-    test (bool): The input is either True or Flase. test is meant for benchmarking, when the feature set has the ground truth included. When set to `True` it will consider the True classification (AMP or Not AMP, class 1 and class 0) and will remove it from the prediction data frame. If you are benchmarking, make sure you name your ground truth column `Class`.

`user_data_training()`

This function basically allows you to input your own protein data and train multiple models (564 per feature set to be exact). It saves all the model metrics in a `.xlsx` file which you can use to select the best models for your own stack. This function does alot but works when everything is properly labeled. My advice would be to not use it yet however here are its inputs for the curious.

-    feature_extracted_data (str): Path to MACREL or AMPeP features. The input must be `.xlsx`. The column names must be ` 'Type (Binary)', 'SEQ', 'NAME', 'Type'` plus whatever column names you get from the feature extraction functions. If not in this format the code will break unless you make the necessary changes in the code itself.
-    model_type (str): Inputs must be RF, LogReg or NN. The LogReg models were the worst performing models, so they aren't trained anymore, however their code is commented. 
-    oversampling (bool): Will oversample class 1 (AMPs) when set to `True`
-    models_parent_dir (str): Path to where you want to save all your models.
-    user_pn_ratios=None (list): Allows you to test your own p:n ratios. When set to `None` it will use the default.
-    user_NN_layers=None (list): Allows you to test your own hidden layer sizes for Neural Networks. When set to `None` it will use the default. 
-    user_testing_data=None (str): Allows you to use your own testing data for base model training. When set to `None` it will use a section of the training data.

`stack_train()`

Once the best perfomring models per feature set are selected, copy them to different folders. For example Random Forest classifiers trained on MACREL features can be copied to `sMACREL/RFs` and Random Forest Classifiers trained on AMPeP features can be copied to `sAMPeP/RFs`. The point is the keep the file tree consistent.

-    RF_path (str): Path to the selected Random Forest Classifiers. AMPeP or MACREL.
-    NN_path (str): Path to the selected Nerual Networks. AMPeP or MACREL. 
-    stack_save_path (str): Path to where you want to save the stack. 
-    stack_train_data (str): Path to the training data have you have been using. The function will sample a different subset of data and will not use the exact same subset used to train the base models. 
-    feature_set (str): Input is 'AMPeP' or 'MACREL'.

### Installation and Dependencies

To use sAMPC, simply download the sAMPC_function_file.py, and the Models. As of now sAMPC can be imported like any normal python library.
`import sAMPC_function_file.py as sff` and call the functions as you would call any normal function `sff.ampep_feature_extract(protein_data='path/to/fasta/file.fasta', save_path='path/to/save/extracted/features)`.

The dependencies needed are:
-    `pandas` (https://pandas.pydata.org/docs/getting_started/install.html)
-    `openpyxl` (Because it writes `.xlsx` files, https://pypi.org/project/openpyxl/)
-    `Biopython` (https://biopython.org/wiki/Download)
-    `scikit-learn` (https://scikit-learn.org/stable/install.html)
-    `peptides` (https://pypi.org/project/peptides/)
-    `joblib` (https://pypi.org/project/joblib/#downloads , however it is usually installed)

### Performance

From our testing, using the 920/920 AMP/NAMP becnhmark by Xiao et al 2013

| Method       | Acc. | Sp.   | Sn.   | Pr.   | MCC  | Reference                         |
|--------------|------|-------|-------|-------|------|-----------------------------------|
| AmPEP        | 0.98 | -     | -     | -     | 0.92 | Bhadra et al. (2018)              |
| Macrel       | 0.95 | 0.998 | 0.9   | 0.998 | 0.9  | Santos-Junior et al (2020)        |
| **cMACREL**  | **0.965**| 0.996 | **0.935**| 0.997 | **0.933**| This study (our clone of MACREL)  |
| **sAMPC**    | 0.954| 0.996 | **0.912** | 0.995 | **0.911**| This study                        |
| **sMACREL**  | 0.946    | 0.987 | 0.907 | 0.986 | 0.896| This study(MACREL stack)          |
| MacrelX      | 0.95 | 0.97  | 0.94  | 0.97  | 0.91 | Santos-Junior et al (2020)        |
| **sAMPeP**   | 0.937| 0.965 | 0.909 | 0.963 | 0.875| This study (AMPeP stack)          |
| iAMP-2L      | 0.95 | 0.92  | 0.97  | 0.92  | 0.9  | Xiao et al. (2013)                |
| AMAP         | 0.92 | 0.86  | 0.98  | 0.88  | 0.85 | Gull, Shamim & Minhas (2019)      |
| CAMPR3-NN    | 0.8  | 0.71  | 0.89  | 0.75  | 0.61 | Waghu et al. (2016)               |
| APSv2        | 0.78 | 0.57  | 0.99  | 0.7   | 0.61 | Veltri, Kamath & Shehu (2018)     |
| CAMPR3-DA    | 0.72 | 0.49  | 0.94  | 0.65  | 0.48 | Waghu et al. (2016)               |
| CAMPR3-SVM   | 0.68 | 0.4   | 0.95  | 0.61  | 0.42 | Waghu et al. (2016)               |
| CAMPR3-RF    | 0.65 | 0.34  | 0.96  | 0.59  | 0.39 | Waghu et al. (2016)               |
| iAMPpred     | 0.64 | 0.32  | 0.96  | 0.59  | 0.37 | Meher et al. (2017)               |

### Issues that need to be fixed
-    The Neural Networks used in sAMPC are not scaled. In our testing the scaled neural networks performed worse than the unscaled Neural Networks, the scaling issue still needs to be fixed.

# Authors

## theblendermann
-    Mail: thebiotechwerido@gmail.com

## Dr. Amit Yadav
-    Site: https://nccs.res.in/research-scientist/details/41
-    Mail: amityadav@nccs.res.in

## Dr. Ranjit Kumar
-    Mail: ranjit.kumar@adypu.edu.in

**Note: Contact Dr Amit Yadav for any queries and theblendermann for any issues regarding sAMPC**
