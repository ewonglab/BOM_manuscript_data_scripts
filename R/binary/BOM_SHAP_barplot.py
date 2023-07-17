import argparse
import sys
import shap

def read_col(fname, col=1, sep=None):
  with open(fname) as fobj:
    next(fobj)
    return [line.split(sep)[col].strip('\n') for line in fobj]

#print("Reading parameters...")
parser = argparse.ArgumentParser(description='Script to produce a SHAP barplot of overall motif importance in a BOM binary model.')
parser.add_argument('--xgb_model', type=str, help='Path to XGBoost model file')
parser.add_argument('--train_data', type=str, help='Path to training data file')
parser.add_argument('--out_file', type=str, help='Path to output PDF or png file')

# Additional options for shap.summary_plot
parser.add_argument('--motif_names', type=str, default='True', help='Whether to use the motif names instead of motif IDs (default: True)')
parser.add_argument('--out_SHAP', type=str, default=None, help='Path to output SHAP values')
parser.add_argument('--max_display', type=int, default=10, help='Number of motifs to display (default: 10). This determines the maximum number of motifs to include in the plot.')
parser.add_argument('--order', type=str, nargs='+', default=shap.Explanation.abs, help='Order of the motifs. It can be a list of motif names or a function that takes an shap.Explanation object and returns a list of motif names.')
# parser.add_argument('--clustering', type=str, default=None, help='Clustering object as in clustering = shap.utils.hclust(X, y).')
parser.add_argument('--hclustering', type=bool, default=False, help='Whether to perform hierarchical clustering.')
parser.add_argument('--clustering_cutoff', type=float, default=0.5, help='Cutoff threshold for clustering (default: 0.5). Motifs with a similarity below this cutoff will not be merged in the dendrogram.')
parser.add_argument('--merge_cohorts', type=bool, default=False, help='Whether to merge motifs from different cohorts (default: False). If set to True, motifs from different cohorts will be merged.')
parser.add_argument('--show_data', type=str, default='auto', help='Option to show the data on the plot. Options include "auto" (default), "data", or "None".')
parser.add_argument('--show', type=bool, default=True, help='Whether to display the plot.')
#parser.add_argument('--help', action='store_true', help='Display help message')

args = parser.parse_args()

if '--help' in sys.argv or '-h' in sys.argv:
    parser.print_help()
    exit()

#print("Loading modules...")

import xgboost
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import os
import numpy as np

# Reading xgb model

print("Reading model...")
model_xgb = xgboost.Booster()
model_xgb.load_model(args.xgb_model)

# Preparing data
print("Reading motif counts...")
X = pd.read_csv(args.train_data, index_col = 0, sep = '\t')
y = X['binary_celltype']
del X['binary_celltype']

# Define clustering if requested
if args.hclustering is True:
  clustering = shap.utils.hclust(X, y)
else:
  clustering = None


# Explain

print("Calculating SHAP values of motif contribution...")

explainer = shap.Explainer(model_xgb)
shap_values = explainer(X)

# Save SHAP values if requested
if args.out_SHAP is not None:
  shap_values_pd = pd.DataFrame.from_records(shap_values.values)
  shap_values_pd.index = X.index
  shap_values_pd.columns = X.columns
  print("Saving SHAP values...")
  shap_values_pd.to_csv(args.out_SHAP, sep='\t', index=True)


if(args.motif_names == 'True'):
  print("Reading GIMME motifs annotation...")
  gimme_annot = "/g/data/zk16/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
  motif_id_0 = read_col(gimme_annot, 0, sep='\t')
  tfs_0 = read_col(gimme_annot, 1, sep='\t')

  motif_id = [sub.replace('-', '.') for sub in motif_id_0]
  # Creating a dictionary where keys are motif IDs and values are all the TFs included in that motif
  motif2factors = {}
  for i in range(0,len(motif_id)):
    motif2factors.setdefault(motif_id[i],[]).append(tfs_0[i])

  motif2factors_selected = {}

  for key,value in motif2factors.items():
    tfs_li = value
    # random.seed(1)
    motif2factors_selected[key] = tfs_li[0]

    motifs = shap_values.feature_names
    motif2tf = dict((k, motif2factors_selected[k]) for k in motifs if k in motif2factors_selected)

  features = []
  for motif in motifs:
    if motif in motif2tf:
      features.append(motif2tf[motif])
    else:
      features.append(motif)
  # Substitute motif names in shap object
  shap_values.feature_names = features

#print(shap_values.feature_names)
print("Saving plot...")

# shap.plots.bar(shap_values, max_display=plot_n)
shap.plots.bar(
  shap_values,
  max_display = args.max_display,
  order = args.order,
  clustering = clustering,
  clustering_cutoff = args.clustering_cutoff,
  merge_cohorts = args.merge_cohorts,
  show_data = args.show_data,
  show = args.show
)

plt.savefig(args.out_file, bbox_inches='tight')
plt.close()

print("Done")

