import numpy as np
import pandas as pd
import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1"
#from google.colab import files
#import model
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
import math
from keras.utils import np_utils
from tensorflow.keras.optimizers import SGD, RMSprop, Adadelta, Adam, Nadam
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.callbacks import ModelCheckpoint, EarlyStopping
from sklearn import metrics
import json
#
from sklearn.model_selection import train_test_split 
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import style
#%matplotlib inline
# Evaluations
from sklearn.metrics import classification_report,accuracy_score, confusion_matrix, precision_score, recall_score, roc_auc_score, roc_curve, f1_score
# Random Forest
from sklearn.ensemble import RandomForestClassifier
from IPython.display import clear_output

#######################################DEFINE FUNCTION #################################################################
#def function 

def getScaledData(dataMatrix):
  scaler = StandardScaler().fit(dataMatrix)
  return scaler.transform(dataMatrix)

# separating feature matrix and class label
def separateDataAndClassLabel(dataMatrix):
    featureMatrix = dataMatrix[:, :(dataMatrix.shape[1] - 1)]
    classLabelMatrix = dataMatrix[:, -1]

    return featureMatrix, classLabelMatrix


# returns the number of classes and encode it
def encodeClassLabel(classLabel):
    labelEncoder = LabelEncoder().fit(classLabel)
    labels = labelEncoder.transform(classLabel)
    classes = list(labelEncoder.classes_)

    return len(classes)


# reshaping the data to the number of bins sizes
def reshapeDataToBinSize(dataMatrix, numberOfFeatEachBin, bins):
    ReshapedData = np.zeros((len(dataMatrix), numberOfFeatEachBin, bins))
    start = 0
    end = numberOfFeatEachBin

    for i in range(1, bins + 1):
        ReshapedData[:, :, i - 1] = dataMatrix[:, start:end * i]
        start = end * i
    return ReshapedData

# returns the TP, TN, FP and FN values
def getTPTNValues(test, testPred):
    TP, TN, FP, FN = 0, 0, 0, 0
    for i in range(len(testPred)):
        if test[i] == testPred[i] == 1:
            TP += 1
        if testPred[i] == 1 and test[i] != testPred[i]:
            FP += 1
        if test[i] == testPred[i] == 0:
            TN += 1
        if testPred[i] == 0 and test[i] != testPred[i]:
            FN += 1

    return TP, TN, FP, FN
#resize data
def resizeData(table1, table2):
  matrix1 = np.array(table1)
  matrix2 = np.array(table2)
  
  matrix1Row, matrix1Col = matrix1.shape
  matrix2Row, matrix2Col = matrix2.shape

  sampleSize = matrix1Row if matrix1Row <= matrix2Row else matrix2Row

  numSampleToSelect =sampleSize # (sampleSize * 95) / 100
  reSizedMatrix1 = matrix1[np.random.choice(matrix1Row, numSampleToSelect, replace=False), :]
  reSizedMatrix2 = matrix2[np.random.choice(matrix2Row, numSampleToSelect, replace=False), :]
  
  return np.vstack((reSizedMatrix1, reSizedMatrix2))

#function to RF, and Linear Regression, Decision Tree
# Function to remove target column and create a data frame from onlydef data_prep(df):
def data_prep(df):
  feature_columns = df.columns[:-1]
  df_features = pd.DataFrame(data=df, columns=feature_columns)
  return df_features
#Standard Scaler
# Custom scaler function
def standardScaling(feature):
   scaler = StandardScaler().fit(feature)
   scaled_feature = scaler.transform(feature)
   scaled_feat = pd.DataFrame(data = scaled_feature, columns =df_features.columns)
   return scaled_feat
########################################################################################################################

#load datasets
#path = 'datasets/G20'
path = './'
#pathResult="results/G20"
pathResult="./"
#ajout du nouveau fichier du dataset
metaDataset = pd.read_csv("datasetTest.csv",sep=";")
print(metaDataset)
row,col = metaDataset.shape
count=1
for index in range(row) :
	print("############################"+metaDataset['featureGroup'][index]+"############################################")
	print(str(count)+" / "+str(row))
	#ajout de ton fichier du datastet
	#META_DATA_FILE_NAME="G20 Gene Essential Paper.xlsx"
	META_DATA_FILE_NAME="Archaeoglobus sulfaticallidus.xlsx"
	FEAT=str(metaDataset['code'][index])
	FEAT_FILE=os.path.join(path, str(metaDataset['filename'][index]))
	META_FILE=os.path.join(path, str(META_DATA_FILE_NAME))
	RESULT_FILE=os.path.join(pathResult, str(FEAT)+"_Result.json")

	print(FEAT_FILE+"\n")
	print(META_FILE+"\n")

	print(RESULT_FILE)
	#open file 
	#FILE_SAVE=open(RESULT_FILE, 'a')
	DICT_DL=dict()
	DICT_DT=dict()
	DICT_LR=dict()
	DICT_RF=dict()

	# Kegg data from G20 Paper
	feature_data = pd.read_csv(FEAT_FILE)
	print(df_full)
	genes = feature_data.index
	#
	meta_df = pd.read_excel(META_FILE, sheet_name="Sheet1",engine="openpyxl")
	#nom du champ pour identifier les gene Gene_Locus
	meta_idx = meta_df['Gene_Locus']
	meta_idx = pd.Series([x.upper() for x in meta_idx.values])
	meta_df = meta_df.set_index(keys=meta_idx)
	print(meta_df)
	#
	# get class labels for dataset
	#remplace class par Gene_essentaility
	df_full = feature_data.merge(meta_df[['Gene_essentaility']], how='inner', left_index=True, right_index=True)
	#Nucleotide feature
	if 'X' in df_full.columns:
		df_full.drop('X', axis = 1, inplace=True)

	# class mappings, 1 = Essential and 0 = Non-Essential
	#deux class NE et E 
	#mappings = {'Dispensable': 0, 'ExpectedEssential': 1, 'EdgeInsertionOnly': 0, 'Desulfovibrio-specific essential': 1, 'NotUnique': 0, 'OtherNoInsertion': 2}
	mappings = {'NE': 0, 'E': 1}
	classes = df_full.pop('Gene_essentaility')
	essential_labels = classes.map(mappings)
	df_full['essential'] = essential_labels

	# remove unknowns
	df_full = df_full[df_full['essential'] < 2]
	#copy dataset
	dataset_full=df_full.copy()
	#
	df_essential = df_full[df_full['essential'] == 1]
	df_nonEssential = df_full[df_full['essential'] == 0]
	# rebalance the classes
	df_essential_oversample = pd.concat([df_essential, df_essential], ignore_index=True)
	# sample non-essential genes
	total_essential_samples = len(df_essential_oversample)
	df_nonE_sample_RF = df_nonEssential.sample(2*total_essential_samples)
	# combine essential and non-essential sets, drop gene name column
	#balance data
	#df_full = pd.concat([df_essential_oversample.iloc[:,:],df_nonE_sample_RF.iloc[:,:]], ignore_index=True)
	#use unbalance data

	df_full = pd.concat([df_essential.iloc[:,:],df_nonEssential.iloc[:,:]], ignore_index=True)

	#df_full = pd.concat([df_essential_oversample.iloc[:,1:],df_nonE_sample_RF.iloc[:,1:]], ignore_index=True)
	# To make original data as intact, deep copy
	#df to RF
	df_RF = df_full.copy()
	#
	print(df_full )
	#

	###########RANDOM FOREST
	#calling the function prepare to get de feature columns
	df_features = data_prep(df_RF)
	#spiting the data to train and test the model
	X=df_features.copy()
	y=df_RF['essential'].copy()
	X_train,X_test, y_train,y_test=train_test_split(X,y,test_size=0.2)
	# Calling the scaler function by passing X_train and X_test to get the scaled data set
	X_train_scaled = standardScaling(X_train)
	X_test_scaled = standardScaling(X_test)
	#
	y_train = y_train.reset_index(drop=True)
	y_test = y_test.reset_index(drop=True)
	# Print X_train scaled data
	row,col = X_train_scaled.shape
	print(row)
	print(X_train_scaled.head())
	#
	#####################RANDOM FOREST BUILD AND TRAIN

	rf_model = RandomForestClassifier(n_estimators=200)
	rf_model.fit(X_train_scaled, y_train)
	#min_impurity_split=None,
	RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
	            max_depth=None, max_features='auto', max_leaf_nodes=None,
	            min_impurity_decrease=0.0, 
	            min_samples_leaf=1, min_samples_split=2,
	            min_weight_fraction_leaf=0.0, n_estimators=200, n_jobs=None,
	            oob_score=False, random_state=None, verbose=0,
	            warm_start=False)
	#
	# Prediction using Random Forest Model
	rf_prediction = rf_model.predict(X_test_scaled)# Evaluations
	train_probs = rf_model.predict_proba(X_train_scaled)[:,1]
	probs = rf_model.predict_proba(X_test_scaled)[:, 1]
	score_roc_auc=roc_auc_score(y_test, probs)
	#
	base_fpr, base_tpr, _ = roc_curve(y_test, [1 for _ in range(len(y_test))])
	model_fpr, model_tpr, _ = roc_curve(y_test, probs)
	plt.figure(figsize = (8, 6))
	plt.rcParams['font.size'] = 16
	plt.plot(base_fpr, base_tpr, 'b', label = 'baseline')
	plt.plot(model_fpr, model_tpr, 'r', label = 'model')
	plt.legend();
	plt.xlabel('False Positive Rate');
	plt.ylabel('True Positive Rate'); plt.title('ROC Curves');
	figure_save=os.path.join(pathResult, str(FEAT)+"roc_curve.png")
	plt.savefig(figure_save)
	#plt.show();
	print("ROC AUC")
	print(score_roc_auc)
	print(train_probs)
	print('Classification Report: \n')
	result=classification_report(y_test,rf_prediction,output_dict=True) #,output_dict=True
	print(result)
	print('\nConfusion Matrix: \n')
	print(confusion_matrix(y_test,rf_prediction))
	# display actual vs. predicted values
	#%load_ext google.colab.data_table
	rf_pred_table = pd.DataFrame({'Predicted': rf_prediction, 'Actual': y_test})
	DICT_RF={"FEAT":FEAT,"MODEL":"RF","METRICS":result,"AUC":score_roc_auc}

	DATA_COLUMNS = df_features.columns.values
	feature_columns = []
	print(df_features)
	for feature_name in DATA_COLUMNS:
	  feature_columns.append(tf.feature_column.numeric_column(feature_name,
	                                         dtype=tf.float32))
	print(df_features)
	#
	feature_imp = pd.Series(rf_model.feature_importances_,index=feature_columns).sort_values(ascending=False)
	FeatImp_FILE=os.path.join(pathResult, str(FEAT)+"_FeatureImportance.csv")
	feature_imp.to_csv(FeatImp_FILE, index=True, sep=";")
	print("############################Feature Importanace BEGIN#######################################")
	print(feature_imp)
	#maptplot
	sns.barplot(x=feature_imp, y=feature_imp.index)
	# Add labels to your graph
	plt.xlabel('Feature Importance Score')
	plt.ylabel('Features')
	plt.title("Visualizing Important Features")
	plt.legend()
	figure_save=os.path.join(pathResult, str(FEAT)+"_FeatureImportance.png")
	plt.savefig(figure_save)
	#plt.savefig(figure_save)
	print("############################Feature Importanace END#######################################")
	with open(RESULT_FILE, 'w') as fp:
		json.dump(DICT_RF, fp)

	print("RANDOM FOREST RESULT "+ str(rf_pred_table))
	###############################################################END RANDOM FOREST#####################################

	print("############################END "+metaDataset['featureGroup'][index]+"############################################")
