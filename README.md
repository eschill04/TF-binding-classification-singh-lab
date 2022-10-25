# Transcription Factor Binding: Classification Model

## 1. Download Data ##

All training/testing data is sourced from [DeepMotif](https://github.com/QData/DeepMotif). 
To use, download and unzip data from link above, and move desired TF folders to ./data. 

## 2. Install Requirements ##

Use 
    pip3 install -r requirements.txt 
to install dependencies.

## 3. Set parameters ##

Changable parameters:
* Transcription factor
    * String that sets the savable model and graph labels
* Training/testing dataset paths
* Training/validation dataset split ratio
* Learning rate
* Batch size
* Epochs

Modify changable parameters in top section of classification_model.py

## 3. Run! ##
To run, type
    python3 classification_model.py
in command line. Your loss curves will be uploaded to the ./result_graphs, and your model will be uploaded to ./models.

## Example Results ##

For parameters: 
    
    TF = "CTCF"

    train_data_path = "./data/CTCF_K562_CTCF_-SC-15914-_Stanford/train.fa"
    test_data_path = "./data/CTCF_K562_CTCF_-SC-15914-_Stanford/test.fa"
    split = 0.7

    learning_rate = .5 
    batch_size = 300
    epochs = 100
    
We get testing accuracy of 94.1%, and the following training/validation curve:

![Loss Curves: CTCF](/result_graphs/sample_loss_curve_CTCF.png)
