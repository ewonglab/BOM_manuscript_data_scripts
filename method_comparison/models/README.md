## Trained models for direct use

For each method, two files, namely method.json and method.hdf5, have been provided for the model's use. 
You can re-load the model by using the following code.

This is an example for re-load deepSTARR model trained on mouse development data.

```
def json_hdf5_to_model(json_filename, hdf5_filename):  
    with open(json_filename, 'r') as f:
        model = model_from_json(f.read())
    model.load_weights(hdf5_filename)
    return model

model= json_hdf5_to_model('deepSTARR.json','deepSTARR_bestloss.hdf5')
```
