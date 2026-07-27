#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import torch.nn.functional as F
from torch import nn
from torch.autograd import Variable
from torch.optim import Adam
import torch.nn.functional as F
from torch.nn.functional import relu,tanh
#from torchvision.datasets.mnist import MNIST
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd
import random
import warnings
import csv
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
warnings.filterwarnings('ignore')
data=pd.read_csv('S1.'+'3'+'.csv',index_col=0) #'/home/ubuntu/RareCaps/adata/X.csv',header=None)#np.genfromtxt('S1.1.counts.txt',skip_header=1)# index_col=None,skiprows=1)#np.genfromtxt('/export/scratch2/sumanta/FiRE/data/preprocessedData_jurkat_two_species_1580.txt')
label=np.genfromtxt('S1.'+'3'+'.labels.txt')#'S1.1.labels.txt')#.astype('int64')
#labels = pd.read_csv('S1.1.labels.csv',index_col=0)
labels=label[1:len(label)]-1
top_k=data.shape[1]
#label=np.genfromtxt('cell_type_no_onlyrna.csv',delimiter=',').astype('int64')
data=data.to_numpy()

for i in range(9):
    print("no. of samples",len(np.where(labels==i)[0]))

  
class ConvCaps2D(nn.Module):
    def __init__(self):
        super(ConvCaps2D, self).__init__()
        # The paper suggests having 32 8D capsules
        self.capsules = nn.ModuleList([nn.Conv2d(in_channels = 1, out_channels = 8, kernel_size=(1,5), stride=2)
                                       for _ in range(32)])
        
    def squash(self, tensor, dim=-1):
        norm = (tensor**2).sum(dim=dim, keepdim = True) # norm.size() is (None, 1152, 1)
        scale = norm / (1 + norm) # scale.size()  is (None, 1152, 1)  
        return scale*tensor / torch.sqrt(norm)
        
    def forward(self, x):
        outputs = [capsule(x).view(x.size(0), 8, -1) for capsule in self.capsules] # 32 list of (None, 1, 8, 36)
        outputs = torch.cat(outputs, dim = 2).permute(0, 2, 1)  # outputs.size() is (None, 1152, 8)
        return self.squash(outputs)
    
class Caps1D(nn.Module):
    def __init__(self):
        super(Caps1D, self).__init__()
        self.num_caps = 2
        self.num_iterations = 3
        self.W = nn.Parameter(torch.randn(2, 2336, 8, 16))
        
    def softmax(self, x, dim = 1):
        transposed_input = x.transpose(dim, len(x.size()) - 1)
        softmaxed_output = F.softmax(transposed_input.contiguous().view(-1, transposed_input.size(-1)))
        return softmaxed_output.view(*transposed_input.size()).transpose(dim, len(x.size()) - 1)

    def squash(self, tensor, dim=-1):
        norm = (tensor**2).sum(dim=dim, keepdim = True) # norm.size() is (None, 1152, 1)
        scale = norm / (1 + norm)        
        return scale*tensor / torch.sqrt(norm)
   
    # Routing algorithm
    def forward(self, u):
        # u.size() is (None, 1152, 8)
        '''
        From documentation
        For example, if tensor1 is a j x 1 x n x m Tensor and tensor2 is a k x m x p Tensor, 
        out will be an j x k x n x p Tensor.
        
        We need j = None, 1, n = 1152, k = 10, m = 8, p = 16
        '''
        
        u_ji = torch.matmul(u[:, None, :, None, :], self.W) # u_ji.size() is (None, 10, 1152, 1, 16)
        
        b = Variable(torch.zeros(u_ji.size())) # b.size() is (None, 10, 1152, 1, 16)
        
        for i in range(self.num_iterations):
            c = self.softmax(b, dim=2)
            v = self.squash((c * u_ji).sum(dim=2, keepdim=True)) # v.size() is (None, 10, 1, 1, 16)

            if i != self.num_iterations - 1:
                delta_b = (u_ji * v).sum(dim=-1, keepdim=True)
                b = b + delta_b
        
        # Now we simply compute the length of the vectors and take the softmax to get probability.
        v = v.squeeze()
     #   print(v.shape)
        y=v.data.cpu().numpy()
        hook1=y[:,0:2,0:16]
    #    print(y.shape)
        y = np.reshape(y,(len(y)*2,16))
   #     print(y.shape)
        with open('test_16.csv', 'w') as outfile:
    #        for slice_2d in x:
           writer=csv.writer(outfile, delimiter='\t')
           writer.writerows(y)
        classes = (v ** 2).sum(dim=-1) ** 0.5
       # print(classes.shape)
        classes = F.softmax(classes) # This is not done in the paper, but I've done this to use CrossEntropyLoss.
      #  print(classes.shape)
        hook=c.data.cpu().numpy()
      #  print(hook.shape)
        #from numpy import savetxt
        #savetxt('test.csv', hook[0:4][0:13,0:928,-1,-1], delimiter=',')
        #with open('test.csv', 'w', newline='') as csvfile:
        #    writer=csv.writer(csvfile, delimiter='\t')
        #    writer.writerows(hook[0:4][0:13,0:928,-1,-1])
        x=hook[:,0:2,0:2336,-1,-1]
       # print(x.shape)
        # x=x.flatten()
        x = np.reshape(x,(len(x)*2,2336))
     #   print(x.shape)
        with open('test.csv', 'w') as outfile:
    #        for slice_2d in x:
           writer=csv.writer(outfile, delimiter='\t')
           writer.writerows(x)
               #np.savetxt(outfile, slice_2d)
        return classes
net = Caps1D()
    
class CapsNet(nn.Module):
    def __init__(self):
        super(CapsNet, self).__init__()
        
        #self.conv1 = nn.Conv2d(in_channels = 1, out_channels = 256, kernel_size = (1,4), stride = 1)
        self.fc1 = nn.Linear(top_k,150)
        self.dropout1 = nn.Dropout(p=0.5)
        self.primaryCaps = ConvCaps2D()
        self.digitCaps = Caps1D()
        
        
    def forward(self, x):
        x = relu(self.dropout1(self.fc1(x)))#F.relu(self.conv1(x))
        x = self.primaryCaps(x)
        x = self.digitCaps(x)
        
        return x

net = CapsNet()
    

import torch.optim as optim
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(net.parameters())
    
def evaluate(model, X, Y, batch_size = 50):
    results = []
    predicted = []
    for i in range(len(X)//batch_size):
        s = i*batch_size
        e = i*batch_size+batch_size
        
        inputs = Variable(torch.from_numpy(X[s:e]))
        pred = net(inputs)
        
        predicted += list(np.argmax(pred.data.cpu().numpy(), axis = 1))
        

    Y=Y[0:len(predicted)]
    #len(predicted)
   # acc = sum(Y == predicted)*1.0/(len(Y)) 
    f1=f1_score(Y, predicted)#,average='weighted')
    return f1 #,acc
    
    
batch_size=128
trn_acc = []
tst_acc = []
trn_loss =[]
tst_loss=[]
f1_trn=[]
f1_tst=[]
#net=CapsNet()
#data1=torch.from_numpy(data)
#fc1 = nn.Linear(top_k,350)
#dropout1 = nn.Dropout(p=0.5)

#data1=F.relu(dropout1(fc1(data1.float())))
#data2=data1.cpu().detach().numpy()
#data2=pd.DataFrame(data)
#data=pd.DataFrame(data)
skf = StratifiedShuffleSplit(n_splits=1, test_size=0.2, train_size=0.8,random_state=0)
#X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.1,stratify=1)


    #net1=net
for train, test in skf.split(data, labels):
    
    X_train=data[train,:]
    y_train=labels[train]
    X_train=pd.DataFrame(X_train)
    d=X_train.values.reshape(len(X_train),1, data.shape[1], order='F')
    d1_train=np.expand_dims(d.astype('float32'), 1)
    indices = np.random.permutation(len(d1_train))
    d1_train=d1_train[indices]
    l1_train=y_train[indices]
    
    X_test=data[test,:]
    y_test=labels[test]
    X_test=pd.DataFrame(X_test)
    d1=X_test.values.reshape(len(X_test),1, data.shape[1], order='F')
    d1_test=np.expand_dims(d1.astype('float32'), 1)
      # l1_test=y_test
    indices = np.random.permutation(len(d1_test))
    d1_test=d1_test[indices]
    l1_test=y_test[indices]
    f1_trn=[]
    f1_tst=[]
    for epoch in range(400):  # 500 epochs
        print('epoch',epoch)
        for phase in ['train', 'validation']:
            if phase == 'train':
                running_loss=0 
                for i in range(len(d1_train)//batch_size-1):    ##iteration
                  #  print(i,)
                    s = i*batch_size
                    e = i*batch_size+batch_size

                    inputs = torch.from_numpy(d1_train[s:e])
                    label = torch.LongTensor(np.array(l1_train[s:e]))

                        # wrap them in Variable
                    inputs, label = Variable(inputs), Variable(label)

                        # zero the parameter gradients
                    optimizer.zero_grad()

                        # forward + backward + optimize
                    outputs = net(inputs)

                    loss = criterion(outputs, label)
                    loss.backward()

                    optimizer.step()
                    running_loss += loss.data.item()
            #    print("Epoch, Loss - {}, {}".format(i, running_loss))
                    del inputs, label
                    #print('\n')
                   # trn_loss.append(running_loss)
            else: 
                r=random.sample(range(1, len(d1_train)), 100)
              #  trn_acc.append(evaluate(net, d1_train[r], l1_train[r], batch_size = 128)) 
               # tst_acc.append(evaluate(net, d1_test, l1_test, batch_size=128)) 
                f1_trn.append(evaluate(net, d1_train[r], l1_train[r], batch_size = 50))
                f1_tst.append(evaluate(net, d1_test, l1_test, batch_size=50))

                #out_train=net(torch.from_numpy(d1_train[r]))
                #out_test=net(torch.from_numpy(d1_test))
               # loss_trn = criterion(out_train, torch.LongTensor(np.array(l1_train[r])))
               # loss_test= criterion(out_test, torch.LongTensor(np.array(l1_test)))
                #trn_loss.append(loss_trn.data.item())
                #tst_loss.append(loss_test.data.item())
                print("f1_score train",f1_trn)
                print("f1_score_test",f1_tst)
             #   print("train_acc",trn_acc)
              #  print("test_acc",tst_acc)
                #logs['log_loss_trn'] = loss_trn.append(loss_trn)
            #logs['log_loss_tst'] = loss_test.append(loss_test)
            #logs['tr_accuracy'] = trn_acc[-1]
            #logs['tst_accuracy'] = tst_acc[-1]

        #liveloss.update(logs)
        #liveloss.draw()
        #print("Epoch, Loss - {}, {}".format(epoch, running_loss))
        #print("Train - ", trn_acc[-1])
        #print("Test - ", tst_acc[-1])
    
    

###get the coupling coefficient for the primary capsule i
print("procesing of gene specific coupling coefficient starts here:")
coup_gene_all=np.zeros(shape=(1,32))
for p in range(0,2000,1):
    print("p is ", p)
    data=pd.read_csv('S1.'+'3'+'.csv',index_col=0)#np.genfromtxt('/export/scratch2/sumanta/FiRE/data/preprocessedData_jurkat_two_species_1580.txt')
   # data=pd.DataFrame(data)
    colnames=data.columns.values
    data.loc[:, np.setdiff1d(colnames,colnames[0:p])]=0
    coup_gene=np.zeros(shape=(1,32))
    coup_gene_final=np.zeros(shape=(1,32))
    for l in range(2):
        n=np.where(labels==l)
        data=np.asarray(data)
        d=data[n[0],:]
        d=d.reshape(len(d),1, 2000)
        d=np.expand_dims(d.astype('float32'), 1)
        net(Variable(torch.from_numpy(d)))
        test_cell= pd.read_csv("test.csv", sep='\t',header=None)
        test_cell=test_cell.to_numpy()
            #test_cell=log(test_cell/(1-test_cell))
        test_cell_reformated=np.zeros(shape=(test_cell.shape[0],32))
        for i in range(np.shape(test_cell)[0]):
            k=0
            for j in range(32):
                test_cell_reformated[i,j]= np.amax(test_cell[i,k:(j+1)*73])
                k=(j+1)*73+1
        avg_s_l1=np.zeros(shape=(2,32))
        s=test_cell_reformated[0,:]
        s1=test_cell_reformated[1,:]
        for i in range(0,test_cell_reformated.shape[0],2):
            s=np.add(s,test_cell_reformated[i,:])
            s1=np.add(s1,test_cell_reformated[(i+1),:])
        avg_s_l1=np.vstack((s/n[0].shape[0],s1/n[0].shape[0]))
        coup_gene=np.vstack((coup_gene,avg_s_l1))
    coup_gene=coup_gene[1:5,:]
    #coup_gene_all.shape
    #np.savetxt('coup_gene_all.csv',coup_gene_all,delimiter=',')
    coup_gene_final=coup_gene[0:2]+coup_gene[2:4]
    coup_gene_all=np.vstack((coup_gene_all,coup_gene_final))
    
np.savetxt('coup_gene_all_jrcat.csv',coup_gene_all,delimiter=',')

