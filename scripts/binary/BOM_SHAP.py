import argparse
import sys
import shap

def read_col(fname, col=1, sep=None):
  with open(fname) as fobj:
    next(fobj)
    return [line.split(sep)[col].strip('\n') for line in fobj]

#print("Reading parameters...")

# Main arguments
parser = argparse.ArgumentParser(description='Script to calculate and save SHAP values for a binary model')
parser.add_argument('--xgb_model', type=str, help='Path to XGBoost model file')
parser.add_argument('--train_data', type=str, help='Path to training data file')
parser.add_argument('--out_file', type=str, help='Path to output file (.txt)')

args = parser.parse_args()

if '--help' in sys.argv or '-h' in sys.argv:
    parser.print_help()
    exit()

print("Loading modules...")

import xgboost
import pandas as pd
import os
import numpy as np
import random


# Reading xgb model
print("Reading model...")

model_xgb_test = xgboost.Booster()
model_xgb_test.load_model(args.xgb_model)

# Preparing data
print("Reading data...")
X = pd.read_csv(args.train_data, index_col = 0, sep = '\t')
del X['binary_celltype']

# Explain
print("Calculating SHAP values...")
explainer = shap.Explainer(model_xgb_test)
shap_values = explainer(X)

shap_values_pd = pd.DataFrame.from_records(shap_values.values)
shap_values_pd.index = X.index
shap_values_pd.columns = X.columns

print("Saving SHAP values...")

shap_values_pd.to_csv(args.out_file, sep='\t', index=True)

print("Done")
