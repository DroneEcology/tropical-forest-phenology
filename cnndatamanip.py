#!/usr/bin/env python
# Ross Gray
# ross.gray12@imperial.ac.uk
#
#  cnndatamanip.py
#
"""

Tree Crown Image Masking, Validation Data and Phenology Category Splitting, Creation and Pickle of Features and Labels for Convolutional Neural Networks.

"""

###--------------------------------------------------------------------------###
### LOAD MODULES
###--------------------------------------------------------------------------###

import os, sys
import numpy as np
import pandas as pd
import re
from datetime import datetime, date, timedelta
import cv2
from tifffile import tifffile as tiff
from PIL import Image
from tqdm import tqdm
import glob
import random
from sklearn.model_selection import train_test_split


###--------------------------------------------------------------------------###

###Create dataset of individual tree crowns
##Labels
# Array
# Tree_Crown_ID - seperate training dataset with labels
# DOY
# Area in pixels

#Tree Extraction Function
def tree_extraction():

    ###Find specific tree in orthomosaics
    # figures = {}
    shapes = []
    def for_each_file(f):
    # for f in tqdm(glob.glob('../Data/LongTermStudy/Orthomosaics/Georeferenced/' + '*.tif')):

        ##Getting filenames
        fullname = os.path.split(f)[1]
        filename = os.path.splitext(fullname)[0]

        ##Getting Day of Year
        find_date = re.compile(r"_([0-9]+-[0-9]+-[0-9]+)_")
        date = find_date.search(filename).group(1)
        doy = datetime.strptime(date, '%d-%m-%y').timetuple().tm_yday

        ##Reading 16 bit ITC tiff file and bgr image
        treecrowns = tiff.imread("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif")
        img_bgr = cv2.imread(f)
        img = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB)

        ##Extract tree crowns
        print("Extracting Tree Crowns")
        def extract_trees(tree):
        # for tree in tqdm(set(rawdata["Tree_Crown_ID"])):

            ##Find tree in image
            locations = np.where(treecrowns == tree)
            xmin = locations[0].min()
            xmax = locations[0].max()
            ymin = locations[1].min()
            ymax = locations[1].max()

            ##Find boundary
            #Take tree from boundary
            boundary = treecrowns[xmin:xmax, ymin:ymax] != tree
            values = img[xmin:xmax, ymin:ymax]
            values[boundary, ] = 0

            ##Finding largest shape
            # shapes.append(boundary.shape)

            ##Expand Mask dynamically
            Max = (340, 320)
            bshape = boundary.shape
            xchange = Max[0] - bshape[0]
            ychange = Max[1] - bshape[1]
            expval = np.pad(values, ((0, xchange),(0,ychange), (0,0)), 'constant', constant_values=0)
            # expval = np.concatenate((values, np.zeros(xchange, ychange)),-1)

            ##Convert back to rgb
            Values = cv2.cvtColor(expval, cv2.COLOR_BGR2RGB)

            ##Saving file if does not exist
            newfile = "../Data/LongTermStudy/CNN/ExtractedTrees/"+format(tree, '03d')+"_"+date+"_"+str(format(doy, '03d'))+".png"
            if not os.path.isfile(newfile):
                cv2.imwrite(newfile, Values)
        [extract_trees(tree) for tree in tqdm(set(rawdata["Tree_Crown_ID"]))]
    *map(for_each_file, tqdm(glob.glob('../Data/LongTermStudy/Orthomosaics/Georeferenced/' + '*.tif'))),
    # [for_each_file(f) for f in tqdm(glob.glob('../Data/LongTermStudy/Orthomosaics/Georeferenced/' + '*.tif'))]
rawdata = pd.read_csv("../Data/LongTermStudy/PixelData/ManualITC.csv")
tree_extraction()

# ## Find largest tree crown
# #max(shapes,key=lambda item:item[0])
# #340 x 320 mask

###--------------------------------------------------------------------------###

###validation

##Splitting trees that have been validated
def validation_split():

    ##Create list of validation tree ID's
    vtreematching = pd.read_csv("../../GroundSampling/GroundPhenology/Data/VTReeMatching.csv")
    treeids = ['{:03d}'.format(i) for i in list(map(int, list(vtreematching['Tree_Crown_ID'].dropna())))]

    ##Creating folder if doesn't exist
    if not os.path.isdir("../Data/LongTermStudy/CNN/ExtractedTrees/Validation"):
        os.mkdir("../Data/LongTermStudy/CNN/ExtractedTrees/Validation")

    ##If match add to validation folder
    def add_to_folder(f):
        #Get old name
        fullname = os.path.split(f)[1]
        filename = os.path.splitext(fullname)[0]

        #Find specific tree id
        find_id = re.compile(r"([0-9]+[0-9]+[0-9]+)_")
        nameid = find_id.search(filename).group(1)

        newpath = "../Data/LongTermStudy/CNN/ExtractedTrees/Validation/"
        #Compare to list of validation trees
        for i in treeids:
            if i == nameid:
                os.rename(f, os.path.join(newpath, fullname))
    # *map(add_to_folder, tqdm(glob.glob("../Data/LongTermStudy/CNN/ExtractedTrees/"+"*.png"))), #The star and comma unzips the map object
    [add_to_folder(f) for f in tqdm(glob.glob("../Data/LongTermStudy/CNN/ExtractedTrees/"+"*.png"))]
validation_split()

###--------------------------------------------------------------------------###

###Phenophase

##Splitting images by pheno classification
def phenocat_split():

    ##Load phenoweight for phenologies
    phenoweight = pd.read_csv("../Data/LongTermStudy/PixelData/PhenoWeight.csv", dtype='object').astype({'Tree_Crown_ID':'int64'})

    ##Replace NA with "None"
    phenoweight = phenoweight.fillna(value = {'Phenology': 'None'})

    ##Add new column with converted ManID to automated Tree_Crown_ID
    # vtreematching = pd.read_csv("../../GroundSampling/GroundPhenology/Data/VTReeMatching.csv")
    # ManID = ['{:03d}'.format(i) for i in list(map(int, list(vtreematching['Manual_ID'].dropna())))]
    # TreeID = ['{:03d}'.format(i) for i in list(map(int, list(vtreematching['Tree_Crown_ID'].dropna())))]
    # combID = list(zip(ManID, TreeID))
    # AutoID = []
    # for i in tqdm(phenoweight['Tree_Crown_ID']):
    #     for j in combID:
    #         if j[0] == '{:03d}'.format(i):
    #             AutoID.append(j[1])
    # phenoweight['AutoID'] = AutoID

    ##Subset each phenology
    LF = phenoweight[phenoweight['Phenology']=='Leaf flush']
    FL = phenoweight[phenoweight['Phenology']=='Flowering']
    FR = phenoweight[phenoweight['Phenology']=='Fruiting']
    LS = phenoweight[phenoweight['Phenology']=='Leaf senescence']
    LL = phenoweight[phenoweight['Phenology']=='Leaf loss']
    NA = phenoweight[phenoweight['Phenology']=='None']
    DD = phenoweight[phenoweight['Phenology']=='Dead']
    dflist = [LF, FL, FR, LS, LL, NA, DD]
    dfnames = ['LF', 'FL', 'FR', 'LS', 'LL', 'NA', 'DD']
    for name in dfnames:
        if not os.path.isdir("../Data/LongTermStudy/CNN/ExtractedTrees/Validation/"+name):
            os.mkdir("../Data/LongTermStudy/CNN/ExtractedTrees/Validation/"+name)
    dfdict = {}
    num=0
    for df in dflist:
        dfid = df['Tree_Crown_ID']
        dfdoy = ['{:03d}'.format(i) for i in list(map(int, list(df['DOY'].dropna())))]
        comb = list(zip(dfid, dfdoy))
        dfdict[dfnames[num]] = comb
        num +=1

    ##If match add to validation folder
    # for key in tqdm(dfdict):
    def for_each_key(key):
        def add_to_folder(f):
            #Get old name
            fullname = os.path.split(f)[1]
            filename = os.path.splitext(fullname)[0]

            #Find specific tree id
            find_id = re.compile(r"([0-9]+[0-9]+[0-9]+)_")
            nameid = find_id.search(filename).group(1)

            #Find specific doy
            find_doy = re.compile(r"_([0-9]+[0-9]+[0-9]+)")
            namedoy = find_doy.search(filename).group(1)

            #Compare to list of validation trees
            newpath = "../Data/LongTermStudy/CNN/ExtractedTrees/Validation/"+key+"/"
            for i in dfdict[key]:
                if (int(i[0]) == int(nameid)) and (int(i[1]) == int(namedoy)):
                    os.rename(f, os.path.join(newpath, fullname))
        [add_to_folder(f) for f in tqdm(glob.glob("../Data/LongTermStudy/CNN/ExtractedTrees/Validation/"+"*.png"))]
    *map(for_each_key, tqdm(dfdict)), #Using list then map is fastest method
phenocat_split()

###--------------------------------------------------------------------------###
###--------------------------------------------------------------------------###

### Creating Training Data

##Directories
# DATADIR = "../Data/LongTermStudy/CNN/ExtractedTrees/Validation"
DATADIR = "../../../../../Desktop/ExtractedTrees/Validation"
CATEGORIES = ['LF', 'FL', 'FR', 'LS', 'LL', 'NA', 'DD']
training_data = []
def create_training_data():
    for category in CATEGORIES:  # do dogs and cats

        path = os.path.join(DATADIR,category)  # create path to dogs and cats
        class_num = CATEGORIES.index(category)  # get the classification

        for img in tqdm(os.listdir(path)):  # iterate over each image per dogs and cats
            try:
                img_array = cv2.imread(os.path.join(path,img), cv2.IMREAD_COLOR)  # convert to array
                new_array = cv2.resize(img_array, (320, 340)) # resize to normalize data size
                training_data.append([new_array, class_num])  # add this to our training_data
            except Exception as e:  # in the interest in keeping the output clean...
                pass
create_training_data()

###Pickle Feature and Labels for future use

##Create features and labels
def create_featurelabels():

    ##Shuffle data
    # import random
    # random.shuffle(training_data)

    ##Creating features and labels
    X = []
    y = []
    for features,label in training_data:
        X.append(features)
        y.append(label)
    #Reshaping: -1 is len of list, 1 is for greyscale / 3 is for colour
    X = np.array(X).reshape(-1, 340, 320, 3)

    ##Traing Test Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)

    ##Saving data for later
    import pickle
    pickle_out = open("../Data/LongTermStudy/CNN/CNNFeaturesTrain.pickle","wb")
    pickle.dump(X_train, pickle_out)
    pickle_out.close()
    pickle_out = open("../Data/LongTermStudy/CNN/CNNFeaturesTest.pickle","wb")
    pickle.dump(X_test, pickle_out)
    pickle_out.close()

    pickle_out = open("../Data/LongTermStudy/CNN/CNNLabelsTrain.pickle","wb")
    pickle.dump(y_train, pickle_out)
    pickle_out.close()
    pickle_out = open("../Data/LongTermStudy/CNN/CNNLabelsTest.pickle","wb")
    pickle.dump(y_test, pickle_out)
    pickle_out.close()
create_featurelabels()

###--------------------------------------------------------------------------###
