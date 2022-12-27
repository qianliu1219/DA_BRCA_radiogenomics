import numpy as np
from keras.models import Model
from keras.layers import Input, Dense, Conv2D, MaxPooling2D
from keras.models import Model
from keras import backend as K
from keras.callbacks import TensorBoard
from keras.callbacks import Callback
import pydicom
import os
import scipy
from sklearn.model_selection import train_test_split
import pandas as pd
from scipy import stats


#clinical label
clinic = pd.read_csv('./Data/cli60229.csv',index_col=0) 

ER_Status = clinic['ER.Status'].to_numpy()
#PR_Status = clinic['PR.Status'].to_numpy()
#HER2_Status = clinic['HER2.Final.Status'].to_numpy()
#T_Status = clinic['Tumor..T1.Coded'].to_numpy()
#N_Status = clinic['Node.Coded'].to_numpy()

label=ER_Status


#image
my_path='/home/qian/TCIA_TCGA-BRCA/DOI_filtered_AIO'
os.chdir("/home/qian/TCIA_TCGA-BRCA/DOI_filtered_AIO")
dirs = os.listdir(my_path)
numpy_int_array=np.array(dirs)

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
X_train, X_test, y_train, y_test = train_test_split(array_all, label, test_size=0.3, random_state=67)

os.chdir("/home/qian/TCIA_TCGA-BRCA/python")

#model#
x_train = np.reshape(x_train, (len(x_train),xsize,ysize,1))
x_test= np.reshape(x_test, (len(x_test),xsize,ysize,1))


input_img = Input(shape=(64,64,1),name = "Input")
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Encode_Layer1")(input_img)
x = MaxPooling2D((2, 2), padding='same',name = "Encode_Layer2")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Encode_Layer3")(x)
x = MaxPooling2D((2, 2), padding='same',name = "Encode_Layer4")(x)
x = Conv2D(16, (3, 3), activation='relu', padding='same',name = "Encode_Layer5")(x)
x = MaxPooling2D((2, 2), padding='same',name = "Encode_Layer6")(x)
output_lab = Dense(1, activation='sigmoid') (x)
classifier = Model(input_img, output_lab)
classifier.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
classifier.summary()
aa=classifier.fit(x_train, y_train,
                epochs=100,
                batch_size=64,
                shuffle=True,
                validation_data=(x_test, y_test),
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=0, write_graph=False)])


#extract feature
get_rd_layer_output = K.function([classifier.layers[0].input],
                                  [classifier.layers[4].output])
layer_output=list()
for i in xrange(len(all_img)):
    layer_output.append(get_rd_layer_output([np.reshape(all_img[i],(1,xsize,ysize,1))])[0])
layer_output_array=np.asarray(layer_output)
layer_output_array=np.reshape(layer_output_array,(len(all_img),16,16,16))

fea = np.reshape(layer_output_array,(len(all_img),4096))
fea = np.transpose(fea)
df=pd.DataFrame(fea)

#save features for LASSO prediction
df.to_csv('features.csv')