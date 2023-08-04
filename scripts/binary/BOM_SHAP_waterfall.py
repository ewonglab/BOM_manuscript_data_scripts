import argparse
import sys

def read_col(fname, col=1, sep=None):
  with open(fname) as fobj:
    next(fobj)
    return [line.split(sep)[col].strip('\n') for line in fobj]

#print("Reading parameters...")

# Main arguments
parser = argparse.ArgumentParser(description='Script to perform operations on a model.')
parser.add_argument('--xgb_model', type=str, help='Path to XGBoost model file')
parser.add_argument('--train_data', type=str, help='Path to training data file')
parser.add_argument('--regions', type=str, help='Path to regions IDs for local explanation')

# Additional options for shap.summary_plot
parser.add_argument('--motif_names', type=str, default='True', help='Whether to use the motif names instead of motif IDs (default: True)')
parser.add_argument('--max_display', type=int, default=10, help='Number of motifs to display (default: 10)')
parser.add_argument('--show', type=bool, default=True, help='Whether to display the plot')
parser.add_argument('--out_format', type=str, default='png', help='Format for output plots. Either png or pdf.')

args = parser.parse_args()

if '--help' in sys.argv or '-h' in sys.argv:
    parser.print_help()
    exit()

# print("Loading modules...")

import xgboost
import shap
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import os
import numpy as np
import random


# Reading xgb model
print("Reading model...")
model_xgb = xgboost.Booster()
model_xgb.load_model(args.xgb_model)

# Preparing training data
print("Reading data...")
X = pd.read_csv(args.train_data, index_col = 0, sep = '\t')
del X['binary_celltype']

# Reading regions IDs to plot
regions = pd.read_csv(args.regions, header=None, squeeze=True)
regions_ids = regions.values.tolist()

# Explain

print("Calculating SHAP values...")

explainer = shap.Explainer(model_xgb)
shap_values = explainer(X)

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



# Region IDs to SHAP indexes

# Iterate over the regions and get the corresponding indices
indices = []
for region in regions_ids:
    try:
        index = X.index.get_loc(region)
        indices.append(index)
    except KeyError:
        # Handle the case when the region is not found in 'X'
        indices.append(None)
        print("Index Error: the requested sequences are not in the table of SHAP values")
        sys.exit(0)

print("Saving beeswarm plot...")

## Plot waterfall plots
for i in range(len(regions_ids)):
  tmp = regions_ids[i].replace(":", "_")
  tmp = tmp.replace("-", "_")
  if(args.out_format == 'png'):
    out_file = "_".join([tmp, "waterfall.png"])
  elif(args.out_format == 'pdf'):
    out_file = "_".join([tmp, "waterfall.pdf"])
  else:
    print("Error: Invalid output format")
    sys.exit()
  print(" ".join(["Saving", out_file, "..."]))
  shap.plots.waterfall(
  shap_values[indices[i],:],
  max_display = args.max_display,
  show = args.show)
  plt.savefig(out_file, bbox_inches='tight')
  plt.close()

print("Done")
