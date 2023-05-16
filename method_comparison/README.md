## Code and pre_trained-models provided for method comparison section

In this file folder, we provided all code for method comparison, including three neural network methods for 
<u>multi-classification comparison</u>; machine learning method - gkm-SVM and fine-tuned DNAbert for 
<u>binary classification comparison</u>.

### Multi-classification comparisons with neural network methods

We compared with three representative neural network methods by retraining models on our a mouse development data (E8.25), to predict 17 cell-type specific enhancers.

- Basset-2016 (see 'Basset_model.ipynb')
- DeepSTARR-2020 (see 'deepSTARR_model.ipynb')
- DeepMEL-2020 (see 'deepMel_model.ipynb')

For a fair comparison, we used 500bp sequences from the center of peaks and their reverse complement sequences to double our sample size. For each method, we re-implemented their approaches according to their neural network architecture and only modified the output layer (n=17) to predict 17 cell-type-specific enhancers. More details about these methods can be found on their GitHub pages.

In this study, we trained the model on mouse development data (E8.25) using the following parameters: Epochs = 100, Batch size = 128, Early stopping = 10. We saved the best val_loss model to perform further evaluation.  Additional details regarding model training and evaluation can be found in the scripts.

If you prefer not to retrain the models, you can directly load our trained model files (located in the ./models folder).

### Binary-classification comparison 

We conducted a comparison with gkm-SVM and a fine-tuned version of DNAbert specific to cell-types.

For the DNAbert approach, we performed fine-tuning on a pre-trained DNAbert model using data from 5 different cell-types: Cardiomyocyte, Neural Crest, Gut, Allantois, and Endothelium cells. This allowed us to predict cell-type specificity, and we selected the best-performing model based on minimal validation loss on the validation set. This final model was then used for comparison. Further information on the DNAbert fine-tuning process can be found in the [GitHub repository](https://github.com/jerryji1993/DNABERT). All the parameters utilized in the study are provided in a supplementary file.




