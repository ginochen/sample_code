import cv2
import numpy as np
import os # to play with the dir
#import mytflearn
from random import shuffle # so to shuffle the data
from tqdm import tqdm # for professional looping with process bar

TRAIN_DIR = '/home/gino/cat-dog/competitions/dogs-vs-cats/train/' # extract the data to 
TEST_DIR = '/home/gino/cat-dog/competitions/dogs-vs-cats/test1/' # 
IMG_SIZE = 50 # 50 by 50 square images 
LR = 1e-3 # learning rate

MODEL_NAME = 'dogvscat-{}-{}.model'.format(LR, '2conv-basic') # save the model as some format

#[x, x] [1,0] [0,1]
# load in images, convert to gray scale, 
def label_img(img):
    # label image in valid syntax 
    # img (filename) is just something like dog.93.png
    word_label = img.split('.')[-3] # split by a period, and back by 3 (i.e., dog)
    if word_label == 'cat': return [1,0] # if 
    elif word_label == 'dog': return [0,1]

def create_train_data(): # create train data 
    training_data = [] # an empty list
    for img in tqdm(os.listdir(TRAIN_DIR)): # os.listdir: Return a list containing the names of the entries in the directory given by path (dog.1.png cat.2.png ...)
        label = label_img(img) # dog or cat label as [1,0] or [0,1]
        path = os.path.join(TRAIN_DIR, img) # join one or more path components intelligently
        img = cv2.resize(cv2.imread(path,cv2.IMREAD_GRAYSCALE), (IMG_SIZE,IMG_SIZE)) # read the path in as gray scale, and resize to 50by50 square
        training_data.append([np.array(img), np.array(label)]) # append numpy array of the image data
    shuffle(training_data) # shuffle the data 
    np.save('train_data.npy', training_data) # save the data in numpy file, so we don't have to rerun this file, just load it next time
    return training_data 

def process_test_data():  
    # kaggle test data, without labelling dog or cat
    testing_data = []
    for img in tqdm(os.listdir(TEST_DIR)):
        path = os.path.join(TEST_DIR, img)
        img_num = img.split('.')[0]
        img = cv2.resize(cv2.imread(path,cv2.IMREAD_GRAYSCALE), (IMG_SIZE,IMG_SIZE)) 
        testing_data.append([np.array(img), img_num]) 

    np.save('test_data.npy',testing_data)
    return testing_data

# start training data
#train_data = create_train_data()
# if you already have train data:
train_data = np.load('train_data.npy')
import tflearn
from tflearn.layers.conv import conv_2d, max_pool_2d
from tflearn.layers.core import input_data, dropout, fully_connected
from tflearn.layers.estimator import regression
#import tflearn.datasets.mnist as mnist

#X, Y, test_x, test_y = mnist.load_data(one_hot=True)

#X = X.reshape([-1, 28, 28, 1])
#test_x = test_x.reshape([-1, 28, 28, 1]) # reshape into 4d array, with 1st dimension number unclear

convnet = input_data(shape=[None, IMG_SIZE, IMG_SIZE, 1], name='input')

convnet = conv_2d(convnet, 32, 2, activation='relu') 
# 1st input: the image data
# 2nd input: the number of filters in the convolution (i.e., a filter is the
# different characteristic/feature part of a dog/cat)
# https://stackoverflow.com/questions/36243536/what-is-the-number-of-filter-in-cnn,
# 3rd input: An integer or tuple/list of 2 integers, specifying the height and width of
# the 2D convolution window. Can be a single integer to specify the same value
# for all spatial dimensions
# activation='relu': Re(ctified) L(inear) (U)nit, assigns negative values as zeros, this becomes the least favored values from the originally 1 and -1 values,
# when taking the filter and multiply by the training data (https://www.youtube.com/watch?v=umGJ30-15_A 5:57), then take an average over the pixels with multiplied values.

convnet = max_pool_2d(convnet,2) # 

convnet = fully_connected(convnet, 1024, activation='relu')
convnet = dropout(convnet, 0.8)

convnet = fully_connected(convnet, 2, activation='softmax')
convnet = regression(convnet, optimizer='adam', learning_rate=LR, loss='categorical_crossentropy',name='targets')

model = tflearn.DNN(convnet, tensorboard_dir='log')

if os.path.exists('{}.meta'.format(MODEL_NAME)):
    model.load(MODEL_NAME)
    print('model loaded!')

train = train_data[:-500] # train data is the last 500 data in the list
test = train_data[-500:] # -500 onward

X =  np.array([i[0] for i in train]).reshape(-1, IMG_SIZE, IMG_SIZE, 1) # feature set
Y =  [i[1] for i in train]

test_x =  np.array([i[0] for i in test]).reshape(-1, IMG_SIZE, IMG_SIZE, 1) # feature set
test_y =  [i[1] for i in test]

model.fit({'input': X}, {'targets': Y}, n_epoch=5, validation_set=({'input': test_x}, {'targets': test_y}), snapshot_step=500, show_metric=True, run_id=MODEL_NAME)
