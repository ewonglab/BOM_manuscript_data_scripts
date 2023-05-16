## Code and pre_trained-models provided for method comparison section

In this file folder, we provided all code for method comparison, including three neural network methods for 
<u>multi-classification comparison</u>; machine learning method - gkm-SVM and fine-tuned DNAbert for 
<u>binary classification comparison</u>.

### Multi-classification comparison with neural network methods

We compared with three representative neural network methods by retraining models on our a mouse development data (E8.25), to predict 17 cell-type specific enhancers.

- Basset-2016 (see 'Basset_model.ipynb')
- DeepSTARR-2020 (see 'deepSTARR_model.ipynb')
- DeepMEL-2020 (see 'deepMel_model.ipynb')

For a fair comparison, we used 500bp sequences from the center of peaks and their reverse complement sequences to double our sample size. For each method, we re-implemented their approaches according to their neural network architecture and only modified the output layer (n=17) to predict 17 cell-type-specific enhancers. More details about these methods can be found on their GitHub pages.

In this study, we trained the model on mouse development data (E8.25) using the following parameters: Epochs = 100, Batch size = 128, Early stopping = 10. We saved the best val_loss model to perform further evaluation.  Additional details regarding model training and evaluation can be found in the scripts.

If you prefer not to retrain the models, you can directly load our trained model files (located in the ./models folder).

### Binary-classification comparison with classical machine learning method - gkm-SVM and transformer-based DNAbert method

We compared with gkm-SVM and DNAbert fine-tuned for cell-types.

For the DNAbert method, we fine-tuned pre-trained DNAbert model for 5 cell-types, including Cardiomyocyte, Nerual_crest, Gut, Allantois and Endothelium cell to predict cell-type specificity and we chose the best performance model (minimal val_loss on valid set) as the final model to compare with. More details for DNAbert finetuning can be found on the [github repo](https://github.com/jerryji1993/DNABERT). All used paramaters are provided in a supplymentary file.




