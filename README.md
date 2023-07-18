# BOM
BOM, short for Bag-of-Motifs, is a powerful framework for analyzing of cis-regulatory regions. It operates on the principle that the activity of these regions relies on the binding of transcription factors (TFs) to specific TF binding motifs. Leveraging available TF binding motif profiles and the Extreme Gradient Boosting (XGBoost) algorithm, a tree-based machine learning architecture, BOM achieves exceptional performance in classifying context-specific cis-regulatory elements. Additionally, through the utilization of SHapley Additive exPlanations (SHAP), BOM enables the identification of important TF binding motifs that contribute to the classification. Notably, BOM requires smaller data sets compared to deep learning architectures, making it an efficient and accessible tool for the analysis of cis-regulatory elements. To facilitate further analysis, BOM offers several visualization options for motif counts and motif importance scores. These visualizations allow users to explore and interpret the most influential motifs learned by BOM, providing valuable insights into the regulatory landscape of the analyzed cis-regulatory regions.

## Installation 
More details

## Usage

You can used our pretrained BOM model to do prediction directly, or you can train a model for your own datasets.

### Prediction with our pre-trained model

Add how to use pretrain model for prediction

### Train model on your own data

- Model for binary classification

- Model for multi-class prediction



## Method Comparison Session

We conducted extensive comparison with the-state-of-art deeplearing methods (Basset; DeepSTARR; DeepMEl; DNAbert) and machine-learing methods(gkm-SVM). More details of comparison can be found in the *Method_Comparison* Folder.


