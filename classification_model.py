import torch
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import random
import matplotlib.pyplot as plt
from Bio import SeqIO

###########################################################################
# MODIFY PARAMETERS:

TF = "CTCF"

train_data_path = "./data/CTCF_K562_CTCF_-SC-15914-_Stanford/train.fa"
test_data_path = "./data/CTCF_K562_CTCF_-SC-15914-_Stanford/test.fa"
split = 0.7

learning_rate = .5 
batch_size = 300
epochs = 100

###########################################################################

# Input: a sequence
# Output: an array of values for each TF

#CONSTANTS
d = 16 #number of motif detectors
m = 24 #len of motif detectors

#initialize datasets
train_sequences = []
train_data = []
test_sequences = []
test_data = []


def process_sequences(sequences):
  data = []
  key = ['A', 'C', 'G', 'T']
  for i in range(0, len(sequences)):
    label = sequences[i].id
    sequence = str(sequences[i].seq)
    char_sequence = list(sequence)
    shape = (2*m + len(char_sequence) - 2, 4)
    tensor_sequence = torch.zeros(shape)
    for i in range(0, 2*m + len(char_sequence) - 2):
      for j in range(0, 4):
        if i < (m) or i > (len(char_sequence) + m - 1):
          tensor_sequence[i][j] = 0.25
        elif char_sequence[i-m] == 'N':
          tensor_sequence[i][j] = 0.25
        elif char_sequence[i-m] == key[j]:
          tensor_sequence[i][j] = 1
        else:
          tensor_sequence[i][j] = 0
    data.append([int(label), tensor_sequence])
  return data



#PREPROCESSING METHOD
def process_data(train_data_path, test_data_path):
  train_tmp = SeqIO.parse(train_data_path, "fasta")
  test_tmp = SeqIO.parse(test_data_path, "fasta")
  for s in train_tmp:
    train_sequences.append(s)
  for s in test_tmp:
    test_sequences.append(s)

#TRAIN VALIDATIONN SPLIT METHOD
def train_val_split(data, ratio):
    train = []
    len_orig = len(data)
    while len(train) < ratio * len_orig:
        idx = random.randint(0, len(data)-1)
        train.append(data[idx])
        data.pop(idx)
    val = data
    return train, val

#RUN PREPROCESSING
print("Beginning preprocessing...")
process_data(train_data_path, test_data_path)
tv_data = process_sequences(train_sequences)
train_data, validation_data = train_val_split(tv_data, split)
test_data = process_sequences(test_sequences)


#DEFINE CUSTOM DATASET 
class CustomDataset(Dataset):
    def __init__(self, data, data_labels, transform=None, target_transform=None):
        self.data_labels = data_labels #make preprocessing fctn -- is 2d array of id, label
        self.data = data #make preprocessing fctn -- is an array of the tensors NEED PROPER SHAPE
        self.transform = transform
        self.target_transform = target_transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        datum = self.data[idx][1] #NEED DIFF FUNCTION
        label = self.data[idx][0]
        if self.transform:
            datum = self.transform(datum)
        if self.target_transform:
            label = self.target_transform(label)
        return datum, label

#INSTANTIATE DATASETS
train_dataset = CustomDataset(train_data, train_data)
validation_dataset = CustomDataset(validation_data, train_data)
test_dataset = CustomDataset(test_data, train_data)
print("Processing done! ")
print("Len training dataset: ", len(train_dataset))
print("Len validation dataset: ", len(validation_dataset))
print("Len testing dataset: ", len(test_dataset))

#WRAP DATALOADER ITERABLE
train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
validation_dataloader = DataLoader(validation_dataset, batch_size=batch_size, shuffle=True)
test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

#BUILD MODEL
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()

        # CONV LAYER
        # input 64 x n+2m-2 (n = 40, m = 24 --> 86) x 4
        self.conv1 = nn.Conv1d(4, d, m) #"one dimensional feature map using a 4 channel input"
        
        # RECTIFY LAYER
        # try leaky relu
        self.relu = nn.ReLU()

        # MAXPOOL LAYER 
        self.maxpool = nn.MaxPool1d(5) # output 64 x 16 x 12

        # NN LAYER (Note: could have hidden layer depending on which works better)
        self.linear = nn.Linear(384, 2) 

        #SOFTMAX
        #self.sigmoid = nn.Sigmoid()
        self.softmax = nn.Softmax(1)

        

    def forward(self, x):
        x = x.permute(0,2,1)
        #print(x.shape)
        c = self.conv1(x)
        #print(c.shape)
        r = self.relu(c) 
        #print(r.shape)
        m = self.maxpool(r)
        #print(m.shape)
        m = m.flatten(1)
        l = self.linear(m)
        #print(l.shape)
        out = self.softmax(l)
        #print(out.shape)
        #print(out)
        return out.double()
        
model = Net()

# TRACK LOSS AND ACCURACY DURING TRAINING
loss_train = []
loss_validation = []
accuracy_graph = []


# Initialize the loss function annd optimizer
loss_fn = nn.CrossEntropyLoss() # MSE is specific to PBM #Binary Cross Entropy
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate) #Paper specifies SGD


def train_loop(dataloader, model, loss_fn, optimizer, t):
    size = len(dataloader.dataset)
    loss_avg = 0.0
    num_batches = len(dataloader)
    for batch, (X, y) in enumerate(dataloader):
        # Compute prediction and loss
        pred = model(X)
        #print(pred)
        #print(X)
        #print(y)
        loss = loss_fn(pred, y)
        loss_avg += loss.item()
        
        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        #torch.nn.utils.clip_grad_norm_(model.parameters(), .001)
        optimizer.step()

        if batch % 100 == 0:
            loss, current = loss.item(), batch * len(X)
            if t%50 == 0:
              print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
    loss_avg /= num_batches
    loss_train.append(loss_avg)


def get_total_correct(pred, y):
  sum = 0
  for i in range(0, len(pred)):
    ans = 0
    if pred[i][1] > pred[i][0]:
      ans = 1
    if ans == y[i]:
      sum = sum + 1
  return sum

def validation_loop(dataloader, model, loss_fn, t):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    val_loss, correct = 0, 0

    with torch.no_grad():
        for X, y in dataloader:
            pred = model(X)
            val_loss += loss_fn(pred, y).item()
            #loss_test.append(loss_fn(pred, y).item())
            correct += get_total_correct(pred, y)
    val_loss /= num_batches
    correct /= size
    accuracy_graph.append(correct)
    loss_validation.append(val_loss)
    if t%50 == 0:
      print(f"Avg loss: {val_loss:>8f} \n")
      print(f"Avg accuracy: {correct:>8f} \n")

print("Beginning training...")
for t in range(epochs):
    if t%50 == 0:
      print(f"Epoch {t}\n-------------------------------")
    train_loop(train_dataloader, model, loss_fn, optimizer, t)
    validation_loop(validation_dataloader, model, loss_fn, t)
print("Training done!")

#GRAPH RESULTS

size = 100
epochs = range(0,size)
plt.plot(epochs, loss_train[0:size], 'g', label='Training loss')
plt.plot(epochs, loss_validation[0:size], 'r', label='Validation loss')
#plt.plot(epochs, accuracy_graph[0:size], 'g', label='Accuracy')
plt.title('Transcription Factor: '+ TF)
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

#get testing accuracy
with torch.no_grad():
    correct = 0
    for X, y in test_dataloader:
        pred = model(X)
        correct += get_total_correct(pred, y)
    correct /= 1000
    print("Testing accuracy:", str(100*correct) + "%")

#save model
torch.save({'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss_validation': loss_validation,'loss_train':loss_train
            }, './models/model'+ TF +'.pt')

print("Model saved to", './models/model'+ TF +'.pt')

#SAVE GRAPH
plt.savefig("./result_graphs/sample_loss_curve_"+TF)
print("Loss graph saved to", './result_graphs/sample_loss_curve_'+ TF +'.png')
