import pandas as pd
# Kegg data from G20 Paper
FEA-NUMBER=50
path = 'datasets'
pathResult = 'Results'
dataset="/content/paperwriting/EssentialFeature_collection_Sequence.csv"
metaDataset = pd.read_csv("datasetList.csv",sep=";")
print(metaDataset)
row,col = metaDataset.shape
count=1
if os.path.isdir(pathResult): 
	print("le repertoire existe")
else:
	print("Creation du repertoire "+ pathResult)
	os.mkdir(pathResult)

feature_data = pd.read_csv(dataset)
for index in range(row) :
	print("############################"+metaDataset['featureGroup'][index]+"############################################")
	print(str(count)+" / "+str(row))
	FEAT=str(metaDataset['code'][index])
	FEAT_FILE=os.path.join(path, str(metaDataset['filename'][index]))
	best_feature = pd.read_csv(FEAT_FILE,sep=";")
	RESULT_FILE=os.path.join(pathResult, "PaperG20_Fea-Best"+str(FEA-NUMBER),str(FEAT)+".csv")
	listFea=list(best_feature['Feature'][:FEA-NUMBER])
	dataset=feature_data[listFea]
	dataset.to_csv(RESULT_FILE)
print("Fin d'extraction des features")





