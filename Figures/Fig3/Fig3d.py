import xgboost
import shap
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import os
import sys
import numpy as np
import random

def read_col(fname, col=1, sep=None):
  with open(fname) as fobj:
    next(fobj)
    return [line.split(sep)[col].strip('\n') for line in fobj]

os.chdir("human/AML/DAP")

xgb_model="AML_vs_healthy_vert5.0_qval0.5.bin"
train_data="AML_vs_healthy_vert5.0_qval0.5_572T_trainSet.csv"

#reading xgb model
model_xgb_test = xgboost.Booster()
model_xgb_test.load_model(xgb_model)

#preparing data
X = pd.read_csv(train_data, index_col = 0)
y = X["binary_celltype"]
del X['binary_celltype']

#explain
explainer = shap.Explainer(model_xgb_test)
shap_values = explainer(X)

# summarized factors
gimme_annot = "summarized_gimme_annot_shorten.txt"
motif_id_0 = read_col(gimme_annot, 0, sep='\t')
tfs_0 = read_col(gimme_annot, 1, sep='\t')
# substitute "-" with "."
motif_id = [sub.replace('-', '.') for sub in motif_id_0]
# dictionary of motif IDs linking to TF summarized name
motif2factors = {}
for i in range(0,len(motif_id)):
  motif2factors.setdefault(motif_id[i],[]).append(tfs_0[i])


motifs = [s.replace("_NA", "") for s in shap_values.feature_names]
features = []
for motif in motifs:
  if motif in motif2factors:
    features.append(motif2factors[motif])
  else:
    features.append(motif)

# Transform elements into strings and remove brackets from lists
proc_features = [
    " ".join(item) if isinstance(item, list) else str(item)
    for item in features
]

# substitute IDs with TF names
shap_values.feature_names = proc_features

filename = f'SHAP_AML_vs_healthy_beswarm.pdf'
shap.summary_plot(shap_values, X, max_display=20, title=f'SHAP - AML vs healthy')
plt.savefig(filename)
plt.close()
   
filename = f'SHAP_SHAP_AML_vs_healthy.txt'
np.savetxt(filename, shap_values.values, delimiter = '\t')
