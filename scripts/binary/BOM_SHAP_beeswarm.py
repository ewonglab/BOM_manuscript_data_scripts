import argparse
import sys
import shap

def read_col(fname, col=1, sep=None):
  with open(fname) as fobj:
    next(fobj)
    return [line.split(sep)[col].strip('\n') for line in fobj]

#print("Reading parameters...")

# Main arguments
parser = argparse.ArgumentParser(description='Script to produce a beeswarm plot of motif SHAP values in a binary model.')
parser.add_argument('--xgb_model', type=str, help='Path to XGBoost model file')
parser.add_argument('--train_data', type=str, help='Path to training data file')
parser.add_argument('--out_file', type=str, help='Path to output PDF or png file')

# Additional options for shap.summary_plot
parser.add_argument('--motif_names', type=str, default='True', help='Whether to use the motif names instead of motif IDs (default: True)')
parser.add_argument('--out_SHAP', type=str, default=None, help='Path to output SHAP values')
parser.add_argument('--max_display', type=int, default=20, help='Number of motifs to display (default: 20)')
parser.add_argument('--features', type=str, default=None, help='Path to features file')
parser.add_argument('--feature_names', type=str, nargs='+', default=None, help='List of feature names')
parser.add_argument('--plot_type', type=str, choices=['dot', 'bar', 'violin', 'compact_dot'], default=None, help='Type of summary plot')
parser.add_argument('--color', type=str, default=None, help='Color for the plot')
parser.add_argument('--axis_color', type=str, default='#333333', help='Color for the plot axis')
parser.add_argument('--title', type=str, default=None, help='Title for the plot')
parser.add_argument('--alpha', type=float, default=1.0, help='Transparency of the plot')
parser.add_argument('--show', type=bool, default=True, help='Whether to display the plot')
parser.add_argument('--sort', type=bool, default=True, help='Whether to sort the features')
parser.add_argument('--color_bar', type=bool, default=True, help='Whether to show the color bar')
parser.add_argument('--plot_size', type=str, default='auto', help='Size of the plot. "auto" (default), float, (float, float), or None.')
parser.add_argument('--layered_violin_max_num_bins', type=int, default=20, help='Maximum number of bins for layered violin plot')
parser.add_argument('--class_names', type=str, nargs='+', default=None, help='List of class names')
parser.add_argument('--class_inds', type=int, nargs='+', default=None, help='Indices of classes to display')
parser.add_argument('--color_bar_label', type=str, default='Feature value', help='Label for the color bar')
parser.add_argument('--cmap', type=str, default=shap.plots.colors.red_blue, help='Colormap for the plot')
parser.add_argument('--auto_size_plot', type=float, nargs=2, default=None, help='Size of the plot')
parser.add_argument('--use_log_scale', type=bool, default=False, help='Whether to use a logarithmic scale')

args = parser.parse_args()

if '--help' in sys.argv or '-h' in sys.argv:
    parser.print_help()
    exit()

# print("Loading modules...")
import xgboost
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

# Explain

print("Calculating SHAP values...")

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



## Plot beeswarm plot

print("Saving beeswarm plot...")

shap.summary_plot(
    shap_values,
    features=args.features,
    feature_names=args.feature_names,
    max_display=args.max_display,
    plot_type=args.plot_type,
    color=args.color,
    axis_color=args.axis_color,
    title=args.title,
    alpha=args.alpha,
    show=args.show,
    sort=args.sort,
    color_bar=args.color_bar,
    plot_size=args.plot_size,
    layered_violin_max_num_bins=args.layered_violin_max_num_bins,
    class_names=args.class_names,
    class_inds=args.class_inds,
    color_bar_label=args.color_bar_label,
    cmap=args.cmap,
    auto_size_plot=args.auto_size_plot,
    use_log_scale=args.use_log_scale
)

plt.savefig(args.out_file)
plt.close()

print("Done")
