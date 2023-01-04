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
from tensorflow.keras.applications import EfficientNetB0



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

xsize=ysize=224.0
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

#model from sratch#
x_train = np.reshape(x_train, (len(x_train),xsize,ysize,3))
x_test= np.reshape(x_test, (len(x_test),xsize,ysize,3))

input_img = Input(shape=(xsize, ysize,3),name = "Input")
outputs = EfficientNetB0(include_top=True, weights=None, classes=2)(input_img)
model = tf.keras.Model(inputs, outputs)
model.compile(optimizer="adam", loss="binary_crossentropy",metrics=["accuracy"])
model.summary()

aa=model.fit(x_train, y_train,
                epochs=100,
                batch_size=64,
                shuffle=True,
                validation_data=(x_test, y_test),
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=0, write_graph=False)])

plot_hist(aa)

#pre-trained model#
def build_model(num_classes):
    inputs = layers.Input(shape=(xsize, ysize,3))
    model = EfficientNetB0(include_top=False, input_tensor=inputs, weights="imagenet")

    # Rebuild top
    x = layers.GlobalAveragePooling2D(name="avg_pool")(model.output)
    x = layers.BatchNormalization()(x)

    top_dropout_rate = 0.2
    x = layers.Dropout(top_dropout_rate, name="top_dropout")(x)
    outputs = layers.Dense(2, activation="softmax", name="pred")(x)

    # Compile
    model = tf.keras.Model(inputs, outputs, name="EfficientNet")
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-2)
    model.compile(
        optimizer=optimizer, loss="binary_crossentropy", metrics=["accuracy"]
    )
    return model


model = build_model(num_classes=2)

model.summary()

bb=model.fit(x_train, y_train,
                epochs=100,
                batch_size=64,
                shuffle=True,
                validation_data=(x_test, y_test),
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=0, write_graph=False)])


plot_hist(bb)
