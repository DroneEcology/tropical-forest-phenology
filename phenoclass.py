#!/usr/bin/env python
# Ross Gray
# ross.gray12@imperial.ac.uk
#
#  phenoclass.py
#
"""

Random forest model for phenology classification.

"""

###--------------------------------------------------------------------------###
###--------------------------------------------------------------------------###

##Load modules
import os, sys #For system access
import numpy as np #Array manipulation
import pandas as pd #Dataframe manipulation
from tqdm import tqdm #For status bar
from itertools import cycle

#Preprocessing and Training
from sklearn.utils import resample #For down or up sampling imbalanced data
from sklearn.model_selection import StratifiedKFold, train_test_split, validation_curve #For testing and training model
from sklearn.preprocessing import scale #For scaling data
from sklearn.decomposition import PCA #Principal Component Analysis for plotting

##Evalutaion
from sklearn.metrics import roc_auc_score, plot_roc_curve
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, precision_recall_curve, average_precision_score, PrecisionRecallDisplay
from sklearn.metrics import classification_report, confusion_matrix, plot_confusion_matrix
from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.metrics import r2_score
from scipy import linalg, stats

#Ensemble classifiers
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier

#Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import Ellipse
import seaborn as sns

###--------------------------------------------------------------------------###

### Functions

##Random Forest Model
def rf_model(vfeatures, nfeatures, vlabels, maxdepth, nest, samsplit, samleaf):

    #Convert phenology to number
    numlabels = pd.factorize(vlabels)[0] + 1 #Add one to remove zeroes

    # ##Traing test split for training the model
    X, X_test, y, y_test = train_test_split(vfeatures, numlabels, test_size=.2, random_state=42)

    ##Renaming
    # X, X_test, y = vfeatures, nfeatures, vlabels

    ##Random Forest
    rf = RandomForestClassifier(max_depth=maxdepth, n_estimators=nest, min_samples_split=samsplit, min_samples_leaf=samleaf, n_jobs=-1, class_weight="balanced").fit(X, y)
    rfscore = rf.score(X_test, y_test)
    # print('Random Forest Score:', rfscore)

    ##Use the forest's predict method on the test data
    predictions = rf.predict(X_test) #Soft voting
    # print(predictions)
    # Calculate the absolute errors
    errors = abs(predictions - y_test)
    # Calculate mean absolute percentage error (MAPE)
    mape = 100 * (errors / y_test)
    # Calculate and display accuracy
    accuracy = 100 - np.mean(mape)
    # print('Random Forest Accuracy (MAPE):', accuracy)

    return rfscore, accuracy

# def plot_rf():
#
#     ##PCA decomposition
#
#     ##Run for each pair of CC
#     for pair in ([0, 1], [0, 2], [1, 2]):

###--------------------------------------------------------------------------###

###Run the model

##Load the data and data wrangling
validdata = pd.read_csv("../Data/LongTermStudy/PixelData/PhenoWeight.csv", dtype=object).dropna(subset=["Phenology", 'RCCchange', 'GCCchange', 'BCCchange'])
alldata = pd.read_csv("../Data/LongTermStudy/PixelData/ManualITC.csv").dropna(subset=['RCCchange', 'GCCchange', 'BCCchange'])
matching = pd.read_csv("../../GroundSampling/GroundPhenology/Data/VTReeMatching.csv")

##Remove validation trees from all data
#Select validation trees
valitrees = list(map(int, list(matching['Tree_Crown_ID'].dropna())))
#Remove them from the data
nvaldata = alldata[-alldata['Tree_Crown_ID'].isin(valitrees)]

##Subsetting data
#By ITC
# validata = validdata.loc[validdata['Tree_Crown_ID'] == 45]
#By Genus
# validata = validdata.loc[validdata['TreeGenus'] == "Shorea"]
#By Species
# validata = validdata.loc[validdata['TreeTaxon'] == "S. fallax"]

##For each tree crown, genus, species
subsets = ['Tree_Crown_ID', 'TreeFamily', 'TreeGenus', 'TreeTaxon']
for sub in subsets:
    allscore = []
    allaccu = []
    allresults = []
    for i in pd.unique(validdata[sub].dropna()):

        #Subset based on species
        validata = validdata.loc[validdata[sub] == i]

        ##Creating features and labels
        #Validated
        valilabels = validata['Phenology']
        # valifeatures = validata[["RCC", "GCC", "BCC"]]
        # valifeatures = validata[["R_StDev", "G_StDev", "B_StDev"]]
        # valifeatures = validata[["R_Mean", "G_Mean", "B_Mean"]]
        # valifeatures = validata[["RCCchange", "GCCchange", "BCCchange"]]
        valifeatures = validata[["RCC", "GCC", "BCC",
                                 "R_StDev", "G_StDev", "B_StDev",
                                 "R_Mean", "G_Mean", "B_Mean",
                                 "RCCchange", "GCCchange", "BCCchange"]]
        vfeature_list = list(valifeatures.columns)

        #Non validated
        # nvalfeatures = nvaldata[["RCC", "GCC", "BCC"]]
        # nvalfeatures = nvaldata[["R_StDev", "G_StDev", "B_StDev"]]
        # nvalfeatures = nvaldata[["R_Mean", "G_Mean", "B_Mean"]]
        # nvalfeatures = nvaldata[["RCCchange", "GCCchange", "BCCchange"]]
        nvalfeatures = nvaldata[["RCC", "GCC", "BCC",
                                 "R_StDev", "G_StDev", "B_StDev",
                                 "R_Mean", "G_Mean", "B_Mean",
                                 "RCCchange", "GCCchange", "BCCchange"]]
        nfeature_list = list(nvalfeatures.columns)

        ##Reshaping
        vfeatures = np.array(valifeatures).reshape(-1,len(vfeature_list))
        vlabels = valilabels
        nfeatures = np.array(nvalfeatures).reshape(-1,len(nfeature_list))

        # print(i)
        #max_depth: None or 2-3, nest:300 or 1200, samsplit: 5, samleaf:5
        score, accuracy = rf_model(vfeatures, nfeatures, vlabels, maxdepth=3, nest=300, samsplit=5, samleaf=5)
        result = [i, score, accuracy]
        allscore.append(score)
        allaccu.append(accuracy)
        allresults.append(result)

    print(sub, 'score:', np.array(allscore).mean())
    print(sub, 'accuracy:', np.array(allaccu).mean())
# np.array(allresults)

###--------------------------------------------------------------------------###
###--------------------------------------------------------------------------###

### Apply RF model to non validated features

##Scaling features
def scale_features(data, resolution):

    ##Data Resolution
    if resolution == 'All':

        ##Resolution
        validata = data

        ##Creating features and labels
        #Validated
        valilabels = validata['Phenology']
        # valifeatures = validata[["RCC", "GCC", "BCC"]]
        # valifeatures = validata[["R_StDev", "G_StDev", "B_StDev"]]
        # valifeatures = validata[["R_Mean", "G_Mean", "B_Mean"]]
        # valifeatures = validata[["RCCchange", "GCCchange", "BCCchange"]
        valifeatures = validata[["RCC", "GCC", "BCC",
                                 "R_StDev", "G_StDev", "B_StDev",
                                 "R_Mean", "G_Mean", "B_Mean",
                                 "RCCchange", "GCCchange", "BCCchange"]]
        vfeature_list = list(valifeatures.columns)

        ##Scale features
        scalevalifeatures = scale(valifeatures)
        dfscalevalifeatures = pd.DataFrame(scalevalifeatures)
        dfscalevalifeatures.columns  = 'Scale_' + vfeature_list
        dfscalevalifeatures[sub] = list(np.repeat(i, len(dfscalevalifeatures)))

        ##Return feature
        return dfscalevalifeatures, vfeature_list, valilabels

    else:
        ##For each tree crown, genus, or species
        def for_each_res(i):

            #Subset based on species
            singlevalidata = data.loc[data[resolution] == i]

            ##Creating features and labels
            #Validated
            singlevalilabels = singlevalidata['Phenology']
            # valifeatures = validata[["RCC", "GCC", "BCC"]]
            # valifeatures = validata[["R_StDev", "G_StDev", "B_StDev"]]
            # valifeatures = validata[["R_Mean", "G_Mean", "B_Mean"]]
            # valifeatures = validata[["RCCchange", "GCCchange", "BCCchange"]
            singlevalifeatures = singlevalidata[["RCC", "GCC", "BCC",
                                                 "R_StDev", "G_StDev", "B_StDev",
                                                 "R_Mean", "G_Mean", "B_Mean",
                                                 "RCCchange", "GCCchange", "BCCchange"
                                                 ]]
            singlevfeature_list = list(singlevalifeatures.columns)

            ##Scale features
            scalevalifeatures = scale(singlevalifeatures)
            singledf = pd.DataFrame(scalevalifeatures)
            singledf.columns  = ['Scale_' + item for item in singlevfeature_list]
            singledf[sub] = list(np.repeat(i, len(singledf)))
            singledf['Phenology'] = list(singlevalilabels)
            return singledf
        allscalefeatures = pd.concat([for_each_res(i) for i in tqdm(pd.unique(data[resolution].dropna()))])
        vfeature_list = list(allscalefeatures.columns[0:12])
        valifeatures = allscalefeatures[vfeature_list]
        valilabels = allscalefeatures['Phenology']
        ##Return feature
        return valifeatures, vfeature_list, valilabels

##Features and Labels
def load_featurelabels(subset = 'None', with_none = False, with_dead = False, scale = False, resolution = 'All'):

    ##Validated
    validdata = pd.read_csv("../Data/LongTermStudy/PixelData/PhenoWeight.csv", dtype=object)

    #Subset
    if subset == 'Tree_Crown_ID':
        validata = validdata.loc[validdata['Tree_Crown_ID'] == 45]
    elif subset == 'TreeTaxon':
        validata = validdata.loc[validdata['TreeTaxon'] == "S. fallax"]
    elif subset == 'TreeGenus':
        validata = validdata.loc[validdata['TreeGenus'] == "Shorea"]
    elif subset == 'None':
        validata = validdata

    #None included
    if with_none == True:
        validata = validata.fillna(value = {'Phenology': 'None'})
    elif with_none == False:
        validata = validata.dropna(subset=['Phenology'])

    #Dead included
    if with_dead == True:
        validata = validata
    elif with_dead == False:
        validata = validata.loc[validata['Phenology'] != "Dead"]

    #Drop empty features
    validata = validata.dropna(subset=['RCCchange', 'GCCchange', 'BCCchange'])

    #Features
    if scale == True:
        valifeatures, vfeature_list, valilabels = scale_features(validata, resolution)
    else:
        valilabels = validata['Phenology']
        valifeatures = validata[[
                                 "RCC", "GCC", "BCC",
                                 "R_StDev", "G_StDev", "B_StDev",
                                 "R_Mean", "G_Mean", "B_Mean",
                                 "RCCchange", "GCCchange", "BCCchange"
                                 ]]
        vfeature_list = list(valifeatures.columns)

    ##Non validated
    alldata = pd.read_csv("../Data/LongTermStudy/PixelData/ManualITC.csv").dropna(subset=['RCCchange', 'GCCchange', 'BCCchange'])
    matching = pd.read_csv("../../GroundSampling/GroundPhenology/Data/VTReeMatching.csv")
    valitrees = list(map(int, list(matching['Tree_Crown_ID'].dropna())))
    nvaldata = alldata[-alldata['Tree_Crown_ID'].isin(valitrees)]
    nvalfeatures = nvaldata[[
                             "RCC", "GCC", "BCC",
                             "R_StDev", "G_StDev", "B_StDev",
                             "R_Mean", "G_Mean", "B_Mean",
                             "RCCchange", "GCCchange", "BCCchange"
                             ]]
    nfeature_list = list(nvalfeatures.columns)

    ##Reshaping
    vfeatures = np.array(valifeatures).reshape(-1,len(vfeature_list))
    vlabels = valilabels
    nfeatures = np.array(nvalfeatures).reshape(-1,len(nfeature_list))

    ##return
    return vfeatures, vlabels, nfeatures

##Training, Testing and Evaluating RF Model
def RF_model():

    ##Load features and labels
    vfeatures, vlabels, nfeatures = load_featurelabels(subset = 'None', with_none = True, with_dead = False, scale = True, resolution = 'Tree_Crown_ID')

    ##Convert labels to numeric
    numlabels, labels = pd.factorize(vlabels)
    numlabels = numlabels + 1 #Add one to remove zeroes

    ##Training and testing
    X, X_test, y, y_test = train_test_split(vfeatures, numlabels, test_size=.2, random_state=42) #Random stat for reproducible set
    #Random Forest
    rf = RandomForestClassifier(max_depth=3, n_estimators=300, min_samples_split=5, min_samples_leaf=5, n_jobs=-1,
                                class_weight="balanced", #Cost Sensitive Training
                                random_state=0 #For reproducible sampling
                                ).fit(X, y)

    ##Evaluating Model
    #Score
    rfscore = rf.score(X_test, y_test)
    print('Random Forest Score:', rfscore)

    #Accuracy
    predictions = rf.predict(X_test) #Soft voting
    errors = abs(predictions - y_test)
    print('RF Average absolute error:', round(np.mean(errors), 2), 'degrees.')
    mape = 100 * (errors / y_test)
    accuracy = 100 - np.mean(mape)
    print('Random Forest Accuracy:', round(accuracy, 2), '%.')

    #Area Under ROC Curve
    predprob = rf.predict_proba(X_test)
    print( roc_auc_score(y_test, predprob, multi_class="ovr", average="weighted") )
    print( roc_auc_score(y_test, predprob, multi_class="ovo", average="weighted") )
    # plot_roc_curve(rf, X_test, y_test, multi_class="ovr")

    ##Precision Recall
    # precision, recall, fscore, support = score(y_test, predictions)
    # print('precision: {}'.format(precision*100))
    # print('recall: {}'.format(recall*100))
    # print('fscore: {}'.format(fscore))
    # print('support: {}'.format(support))
    print(classification_report(y_test, predictions))

    ##Confusion Matrix
    labels_names = labels
    # plot_confusion_matrix(rf, X_test, y_test, display_labels = labels_names, normalize = 'true')
    cm = confusion_matrix(y_test, predictions, normalize = 'true')
    ax = plt.subplot()
    sns.heatmap(cm, annot=True, ax = ax, fmt='.2', cmap='Greens'); #annot=True to annotate cells
    ax.set_xlabel('Predicted Labels'); ax.set_ylabel('True Labels');
    ax.xaxis.set_ticklabels(labels_names, rotation = 45); ax.yaxis.set_ticklabels(labels_names, rotation='horizontal');
    plt.savefig("../Results/PhenoClassifier/RFConfusionMatrix.png", bbox_inches = "tight", pad_inches=0, dpi=300)
    plt.close()

    ##Using Pearson Coeefficient (r2)
    print('RF R2 Score:', r2_score(y_test, predictions))
    print('RF Pearson Coefficient', stats.pearsonr(y_test, predictions))
    np.corrcoef(y_test, predictions)

    return rf
rf = RF_model()

###--------------------------------------------------------------------------###

def apply_RFmodel():
    ##Applying to non-validated data
    predictions = rf.predict(nfeatures)

    ##Add predictions back to data frame
    #Convert numeric labels back to strings
    pred_labels = labels[predictions-1]

    #Assign empty columns
    nvaldata2 = pd.concat([nvaldata,pd.DataFrame(columns=['Phenology', 'CPWeighting', 'PhenologyPerc', 'PhenologyTS',
                                                         'PhenologyCP', 'PhenologyGMM', 'BinPheno', 'TreeFamily',
                                                         'TreeGenus', 'TreeSpecies', 'TreeTaxon', 'ClosestTrap',
                                                         'DryWeight', 'WetWeight'])], sort=False)

    #Assign new phenology to dataframe
    nvaldata2 = nvaldata2.assign(NewPhenology=pred_labels.values)

    #Load outlier validation data with other data
    validdata2 = pd.read_csv("../Data/LongTermStudy/PixelData/PhenoWeight.csv", dtype=object)
    valitrees2 = list(map(int, list(matching['Tree_Crown_ID'].dropna())))
    validdata2 = validdata2.astype({'Tree_Crown_ID': int})
    validdata3 = validdata2[validdata2['Tree_Crown_ID'].isin(valitrees2)]
    validdata3 = validdata3.assign(NewPhenology=validdata3['Phenology'])

    #Combine validation and RF classified data
    len(alldata) == len(validdata3) + len(nvaldata2)
    combdata = pd.concat([validdata3, nvaldata2])

    combdata.to_csv("../Data/LongTermStudy/PixelData/ClassPhenoWeight.csv", index=False)

###--------------------------------------------------------------------------###
###--------------------------------------------------------------------------###

###Choosing Hyperparameters
#
# ##Subset species
# validata = validdata.loc[validdata['TreeSpecies'] == "fallax"]
# valilabels = validata['Phenology']
# valifeatures = validata[["RCC", "GCC", "BCC",
#                          "R_StDev", "G_StDev", "B_StDev",
#                          "R_Mean", "G_Mean", "B_Mean",
#                          "RCCchange", "GCCchange", "BCCchange"]]
# vfeature_list = list(valifeatures.columns)
# vfeatures = np.array(valifeatures).reshape(-1,len(vfeature_list))
# vlabels = valilabels
# numlabels = pd.factorize(vlabels)[0] + 1
# X, X_test, y, y_test = train_test_split(vfeatures, numlabels, test_size=.2, random_state=42)
#
# # param_name = 'n_estimators'
# # param_range = [10,50,100,200,300,500,800,1200,1400] #Number of Estimators
# # param_name = 'max_depth'
# # param_range = [1, 2, 3, 4, 5, 8, 15, 25, 30, 50, 100] #Max Depth
# # param_name = 'min_samples_split'
# # param_range = [2, 5, 10, 15, 100] #Minimum number of samples to split a node
# param_name = 'min_samples_leaf'
# param_range = [1, 2, 5, 10] #Mininum number of samples required to be a leaf node
# train_scores, test_scores = validation_curve(
#                                 RandomForestClassifier(), X = X, y = y,
#                                 param_name = param_name, param_range = param_range,
#                                 cv=2
#                                 )
# train_scores_mean = np.mean(train_scores, axis=1)
# train_scores_std = np.std(train_scores, axis=1)
# test_scores_mean = np.mean(test_scores, axis=1)
# test_scores_std = np.std(test_scores, axis=1)
#
# plt.title("Validation Curve with RF")
# plt.xlabel(param_name)
# plt.ylabel("Score")
# plt.ylim(0.0, 1.1)
# lw = 2
# plt.plot(param_range, train_scores_mean, label="Training score",
#              color="darkorange", lw=lw)
# plt.fill_between(param_range, train_scores_mean - train_scores_std,
#                  train_scores_mean + train_scores_std, alpha=0.2,
#                  color="darkorange", lw=lw)
# plt.plot(param_range, test_scores_mean, label="Cross-validation score",
#              color="navy", lw=lw)
# plt.fill_between(param_range, test_scores_mean - test_scores_std,
#                  test_scores_mean + test_scores_std, alpha=0.2,
#                  color="navy", lw=lw)
# plt.legend(loc="best")
# plt.show()


###--------------------------------------------------------------------------###
### Gradient Boosting Classifier
# def GB_model():
#     ##Load features and labels
#     vfeatures, vlabels, nfeatures = load_featurelabels(with_none = False, with_dead = False, resolution = 'None')
#
#     ##Convert labels to numeric
#     numlabels, labels = pd.factorize(vlabels)
#     numlabels = numlabels + 1 #Add one to remove zeroes
#
#     ##Training and testing
#     X, X_test, y, y_test = train_test_split(vfeatures, numlabels, test_size=.2, random_state=42) #Random stat for reproducible set
#     #Gradient Boosting Classifier
#     gb = GradientBoostingClassifier(max_depth=3, n_estimators=300, min_samples_split=5, min_samples_leaf=5,
#                                 # min_weight_fraction_leaf="balanced", #Cost Sensitive Training
#                                 random_state=0 #For reproducible sampling
#                                 ).fit(X, y)
#
#     ##Evaluating Model
#     #Score
#     gbscore = gb.score(X_test, y_test)
#     print('Gradient Boosting Score:', gbscore)
#
#     #Accuracy
#     predictions = gb.predict(X_test) #Soft voting
#     errors = abs(predictions - y_test)
#     print('GB Average absolute error:', round(np.mean(errors), 2), 'degrees.')
#     mape = 100 * (errors / y_test)
#     accuracy = 100 - np.mean(mape)
#     print('Gradient Boosting Accuracy:', round(accuracy, 2), '%.')
#
#     #Area Under ROC Curve
#     # predprob = rf.predict_proba(X_test)
#     # print( roc_auc_score(y, predprob, multi_class="ovr") )
#     # print( roc_auc_score(y, predprob, multi_class="ovo") )
#
#     ##Precision Recall
#     precision, recall, fscore, support = score(y_test, predictions)
#     # print('precision: {}'.format(precision*100))
#     # print('recall: {}'.format(recall*100))
#     # print('fscore: {}'.format(fscore))
#     # print('support: {}'.format(support))
#     print(classification_report(y_test, predictions))
#
#     ##Confusion Matrix
#     labels_names = labels
#     plot_confusion_matrix(gb, X_test, y_test, display_labels = labels_names, normalize = 'true')
#     plt.savefig("../Results/PhenoClassifier/GBConfusionMatrix.png", bbox_inches = "tight", pad_inches=0, dpi=300)
#     plt.close()
#
#     ##Using Pearson Coeefficient (r2)
#     print('GB R2 Score:', r2_score(y_test, predictions))
#     print('GB Pearson Coefficient', stats.pearsonr(y_test, predictions))
#     np.corrcoef(y_test, predictions)
#
#     return gb
# gb = GB_Model()
