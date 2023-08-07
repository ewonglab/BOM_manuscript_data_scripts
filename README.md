# BOM
BOM, short for Bag-Of-Motifs, is a framework for analyzing cis-regulatory regions. 
It operates on the principle that the activity of these regions relies on the binding of transcription factors (TFs) to specific TF binding motifs. Leveraging available TF binding motif profiles and the Extreme Gradient Boosting (XGBoost) algorithm, a tree-based machine learning architecture, BOM has achieved high performance in classifying context-specific cis-regulatory elements. 
Additionally, through the use of SHapley Additive exPlanations (SHAP), BOM enables the identification of important TF binding motifs that contribute to the classification. 

To facilitate further analysis, BOM offers several visualization options for motif counts and motif importance scores. 
These visualizations allow users to explore and interpret the most influential motifs learned by BOM, providing valuable insights into the regulatory landscape of the analyzed cis-regulatory regions.

<!---
![BOM pipeline](///BOM_pipeline.png)
-->

## R dependencies

Please load the below R packages.
> rsample; xgboost; GenomicFeatures; GenomicRanges; data.table; tidyr; dplyr; cvAUC; pROC; ggplot2; yardstick

## Installation 

```
git clone https://github.com/ewonglab/BOM.git
cd BOM
chmod +x ./scripts/binary/*.sh
```

## Tutorial

For a comprehensive guide on how to fully utilize BOM, including data pre-process, model training, model prediction and model interpretation. 
Please refer to the tutorials provided below.
- [Model training and prediction](https://ewonglab.github.io/BOM/)
- [Model interpretation](https://ewonglab.github.io/BOM/binary_models_interpretation/)

## Candicate cis-regulatory elements filter

Usage:  ```filter_CREs.R  --help ```

Parameters:
> 
> --help  Display this help message
> 
> --input_bed=<file> input BED file
> 
> -- annot=<file> assembly annotation file (.gtf) 
>
> -- chr_sizes=<file> File with chromosome length values 
>
> --u=<integer>   Number of base pairs upstream of TSS for the definition of proximal regions
>
> --d=<integer>   Number of base pairs downstream of TSS for the definition of proximal regions
>
> --nbp=<integer> Number of central base pairs in adjusted CREs 
>
> --out_bed=<file> output BED file
>
> --keep_proximal=<logical> whether only proximal regions to TSS should be retained (default: FALSE)
>
> --remove_proximal=<logical> whether proximal regions to TSS should be removed (default: FALSE)
>
> --non_exonic=<logical> whether regions overlapping exons should be removed (default: FALSE)

## Motif searching

Usage: ``` run_fimo.sh  -m /gimme.vertebrate.v5.0.meme  -g /Mus_musculus_GRCm38.fa  -b /Tutorial/bed_files -o /Tutorial/motifs```

```
run_fimo.sh: 
Usage: args [-o] [-g] [-b] [-m] 
-m means path to motif database
-g means path to genome reference fasta file
-o means set a output path (Default: ./Tutorial/motifs)
-b means path to bed file folder (Default: ./Tutorial/bed_files)
```

## Motif Counting

Usage:  ```Rscript matrix_for_binary_model.R --help```

Options:

> --help                           Display this help message
>
>  --target_ct=<target_cell_type>    Name of the target cell type
>  
>  --data_path=<data_directory>     Path to the data directory
>  
>  --qval_thresh=<qvalue_threshold> Q-value threshold for filtering
>  
>  --out_filename=<output_filename> Name for the output file
  
### Model Training

- Model for binary classification

- Model for multi-class prediction

### Training a model for binary classification

Usage: ```Rscript training_binary.R --help```

Options:

> --help                    Display this help message
>
> --input_data=<file>		Path to the input matrix of motif counts (required)
>
> --data=<data>			The training data (default: dtrain)
>
> --nrounds=<n>			Number of boosting rounds (default: 10000)
>
> --eta=<value>			Learning rate (default: 0.01)
>
> --max_depth=<n>		Maximum tree depth (default: 6)
>
> --subsample=<value>		Subsample ratio of the training instances (default: 0.5)
>
> --colsample_bytree=<value>	Subsample ratio of columns when constructing each tree (default: 0.5)
>
> --objective=<name>		Objective function (default: binary:logistic)
>
> --early_stopping_rounds=<n>	Perform early stopping if no improvement for this many rounds (default: NULL)
>
> --nthread=<n>			Number of parallel threads (default: 1)
>
> --eval_metric=<name>		Evaluation metric (default: error)
>
> --maximize=<bool>		Whether to maximize the evaluation metric (default: FALSE)
>
> --save_period=<n>		Save the model for every given number of rounds (default: NULL)
>
> --save_name=<file>		Name of the saved model file (default: xgboost.model)
>
> --feval=<file>		Customized evaluation metric (default: NULL)
>
> --Verbose=<file>		How much details on the progress to print (default: 1)
>
> --print_every_n=<file>		Print evaluation messages each n-th iterations (default: 1)
>
> --save_name=<file>		Name of the saved model file (default: xgboost.model)>


### Model Prediction

Usage: ``` Rscript xgboost_predictions.R --help```

Parameters:
> --help                    Display this help message
>
> --input_data=<file>       Path to the input data file
>
> --xgb_model=<file>        Path to the xgboost model file
>
> --predictions=<file>      Path to save the predicted values
>
> --training_set=<file>     Path to save the training set (optional)


<!---
### How to cite
If you use BOM in your work, please cite:

A Bag-Of-Motif Captures Context-Specific Distal Regulatory Elements

Paola Cornejo-Paramo, Xuan Zhang, Lithin Louis, Yi-Hua Yang, Emily S. Wong
-->








