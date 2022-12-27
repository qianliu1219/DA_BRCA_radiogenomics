#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:43:43 2018

@author: root
"""

import numpy as np
from keras.models import Model
from keras.layers import Input, Dense, Conv2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras import backend as K
from keras.callbacks import TensorBoard
import matplotlib.pyplot as plt
import seaborn as sns
from keras.callbacks import Callback
import pydicom
import os
import scipy
from sklearn.model_selection import train_test_split
import pandas as pd
from scipy import stats
from keras.utils import plot_model

#path#
my_path='/home/qian/TCIA_TCGA-BRCA/DOI_filtered_AIO'
os.chdir("/home/qian/TCIA_TCGA-BRCA/DOI_filtered_AIO")
dirs = os.listdir(my_path)
numpy_int_array=np.array(dirs)
label=[]
for i in xrange(len(numpy_int_array)):
    label.append(numpy_int_array[i][0:12])
label_array=np.asarray(label)
np.savetxt("patientID.txt",label_array,fmt="%s")

#zoom#
xsize=ysize=64.0
array_all=list()
for i in numpy_int_array:
    input_im=pydicom.dcmread(i,1)
    im_array=input_im.pixel_array
    if np.max(im_array.reshape(-1)/1.0)==0:
        im_array=im_array
    else:
        im_array=im_array/(np.max(im_array.reshape(-1)/1.0))
    print im_array.shape
    xscale = xsize/im_array.shape[0]
    yscale = ysize/im_array.shape[1]
    reshaped_im_array= scipy.ndimage.interpolation.zoom(im_array,[xscale,yscale])
    print(reshaped_im_array.shape)
    array_all.append(reshaped_im_array)

#train_test_split#
x_train, x_test = train_test_split(array_all, test_size=0.3)

#add noise#
noise_factor = 0.05
x_train_noisy = x_train + noise_factor * np.random.normal(loc=0, scale=1.0, size=np.shape(x_train))
x_test_noisy = x_test + noise_factor * np.random.normal(loc=0, scale=1.0, size=np.shape(x_test))
x_train_noisy = np.clip(x_train_noisy, 0., 1.)
x_test_noisy = np.clip(x_test_noisy, 0., 1.)

os.chdir("/home/qian/TCIA_TCGA-BRCA/python")

plt.figure()
ax = plt.subplot(1, 1, 1)
plt.imshow(array_all[3] + noise_factor * np.random.normal(loc=0, scale=1.0, size=np.shape(array_all[4])))
plt.gray()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.savefig('noise.png')

#model#
x_train = np.reshape(x_train, (len(x_train),xsize,ysize,1))
x_test= np.reshape(x_test, (len(x_test),xsize,ysize,1))

x_train_noisy = np.reshape(x_train_noisy, (len(x_train_noisy),xsize,ysize,1))
x_test_noisy= np.reshape(x_test_noisy, (len(x_test_noisy),xsize,ysize,1))


input_img = Input(shape=(64,64,1),name = "Input")
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Encode_Layer1")(input_img)
x = MaxPooling2D((2, 2), padding='same',name = "Encode_Layer2")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Encode_Layer3")(x)
x = MaxPooling2D((2, 2), padding='same',name = "Encode_Layer4")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Encode_Layer5")(x)
encoded = MaxPooling2D((2, 2), padding='same',name = "Encode_Layer6")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Decode_Layer1")(encoded)
x = UpSampling2D((2, 2),name = "Decode_Layer2")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Decode_Layer3")(x)
x = UpSampling2D((2, 2),name = "Decode_Layer4")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Decode_Layer5")(x)
x = UpSampling2D((2, 2),name = "Decode_Layer6")(x)
decoded = Conv2D(1, (3, 3), activation='sigmoid', padding='same',name = "Output")(x)

autoencoder = Model(input_img, decoded)
autoencoder.compile(optimizer='adadelta', loss='cosine_proximity')
autoencoder.summary()
aa=autoencoder.fit(x_train_noisy, x_train,
                epochs=100,
                batch_size=64,
                shuffle=True,
                validation_data=(x_test_noisy, x_test),
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=0, write_graph=False)])


all_img = np.reshape(array_all, (len(array_all),xsize,ysize,1))
decoded_imgs = autoencoder.predict(all_img)
plt.figure()
ax = plt.subplot(1, 1, 1)   
plt.imshow(np.reshape(decoded_imgs[3],np.shape(array_all[3])))
plt.gray()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.savefig('denoise.png')

plot_model(autoencoder, to_file='model.png',show_shapes=True)

#extract feature
get_rd_layer_output = K.function([autoencoder.layers[0].input],
                                  [autoencoder.layers[4].output])
layer_output=list()
for i in xrange(len(all_img)):
    layer_output.append(get_rd_layer_output([np.reshape(all_img[i],(1,xsize,ysize,1))])[0])
layer_output_array=np.asarray(layer_output)
layer_output_array=np.reshape(layer_output_array,(len(all_img),16,16,16))

##quantile normalization
fea = np.reshape(layer_output_array,(len(all_img),4096))
fea = np.transpose(fea)
df=pd.DataFrame(fea)
rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
df2 = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
fea_qt = df2.values
label=[]
for lab in label_array:
    lab = lab.replace("-", ".")
    label.append(lab)
    
df2.columns=label
df2=df2.transpose()
df2.to_csv('features.csv')

fea_qt_k=np.reshape(np.transpose(fea_qt),(len(all_img),16,16,16))

##visulization
for i in xrange(16):
    feak=pd.DataFrame(fea_qt_k[4,:,:,i])
    heat = sns.heatmap(feak, square = True)
    fig = heat.get_figure()
    fig.savefig('heat{}.png'.format(i))

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(6, 5))
ax1.set_title('Before qt Normalization')
sns.kdeplot(df[1], ax=ax1)
sns.kdeplot(df[2], ax=ax1)
sns.kdeplot(df[3], ax=ax1)
ax2.set_title('After qt Normalization')
sns.kdeplot(df2[1], ax=ax2)
sns.kdeplot(df2[2], ax=ax2)
sns.kdeplot(df2[3], ax=ax2)
plt.savefig('kde.png')

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(6, 5))
ax1.set_title('Before qt Normalization')
sns.kdeplot(df.iloc[1], ax=ax1)
sns.kdeplot(df.iloc[2], ax=ax1)
sns.kdeplot(df.iloc[3], ax=ax1)
ax2.set_title('After qt Normalization')
sns.kdeplot(df2.iloc[1], ax=ax2)
sns.kdeplot(df2.iloc[2], ax=ax2)
sns.kdeplot(df2.iloc[3], ax=ax2)
plt.savefig('kde2.png')
