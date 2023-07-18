## Code and pre_trained-models provided for the method comparison section

In this file folder, we provided all code for method comparison, including three neural network methods for 

<u>multi-classification comparison</u>; machine learning method - gkm-SVM and fine-tuned DNAbert for 
<u>binary classification comparison</u>.

### Multi-classification comparisons with neural network methods

We compared with three representative neural network methods by retraining models on our mouse development data (E8.25) to predict 17 cell-type specific enhancers.

- Basset-2016 (see 'Basset_model.ipynb')
- DeepSTARR-2020 (see 'deepSTARR_model.ipynb')
- DeepMEL-2020 (see 'deepMel_model.ipynb')

For a fair comparison, we used 500bp sequences from the center of peaks and their reverse complement sequences to double our sample size. 
For each method, we re-implemented their approaches according to their neural network architecture and only modified the output layer (n=17) to predict 17 cell-type-specific enhancers. 
More details about these methods can be found on their GitHub pages.

In this study, we trained the model on mouse development data (E8.25) using the following parameters: Epochs = 100, Batch size = 128, Early stopping = 10. 
We saved the best val_loss model to perform the further evaluation.  Additional details regarding model training and evaluation can be found in the scripts.

If you prefer not to retrain the models, you can directly load our trained model files (located in the ./models folder).

### Binary-classification comparison 

We conducted a comparison with the gkm-SVM and a fine-tuned version of DNAbert specific to cell types.

- gkm-SVM-2016 (see 'gkm-SVM.ipynb')
- DNAbert-2021 (see 'DNAbert_finetune.ipynb')

For the DNAbert approach, we employed the 6-mer pre-trained DNABert model and fine-tuned it to perform binary prediction of cell-type specific enhancers for a total of 17 cell types. To fine-tune models for each cell type, we used the recommended parameters: epochs=5, learning_rate=0.0002, logging_steps=100, warmup_percent=0.1 and weight_decay=0.01. 
During fine-tuning, we monitored the change in loss on validation sets and selected the model with the best loss for further evaluation. However, the recommended parameters did not yield satisfactory results for most of the cell types, as we observed no significant decrease in loss or increase in accuracy during the fine-tuning process. 
Therefore, we conducted further exploration by fine-tuning with different learning rates, namely a larger learning rate of 0.0003 and a lower learning rate of 0.0001, accompanied by a larger epoch value of 20. Ultimately, we selected the model with the lowest loss for each cell type to make comparisons. Further information on the DNAbert fine-tuning process can be found in the [GitHub repository](https://github.com/jerryji1993/DNABERT). All the parameters utilized in the study are provided in a supplementary file.




