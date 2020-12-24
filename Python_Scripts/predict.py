from tensorflow.keras.models import load_model
#from keras.models import load_model
import numpy as np
import cv2
from data_loader import get_image_array, get_pairs_from_paths, get_images_from_path
import os

from metrics import evaluate_metrics

def prediction(model, image_path, norm_type=None):
    img = cv2.imread(image_path, 0)
    img = get_image_array(img, norm_type)#relies on the model name
    img = np.array(img)
    img = img.reshape((1,) + img.shape)
    #print("image shape: " + str(img.shape))
    preds_test = model.predict(img, verbose=0)
    preds_test = preds_test.argmax(axis=-1) * 50

    return preds_test[0]


# method for using a model (specified by path) to output segmentations for all images in a folder to another folder
def segment_folder(model, frame_path, output_folder):
    imgs = get_images_from_path(frame_path)
    norm_type = None
    if 'divide_and_sub' in model:
        norm_type = 'divide_and_sub'
    elif 'sub_mean' in model:
        norm_type = 'sub_mean'
    elif 'divide' in model:
        norm_type = 'divide'
    model = load_model(model, compile=False)
    for img in imgs:
        #print("processing1: " + img)
        img_name = img[img.rfind('/') + 1:]#this gets the image name out of the path
        #print("processing: " + img_name)
        pred = prediction(model, img, norm_type)
        #print("saving in: " +output_folder + img_name)
        cv2.imwrite(output_folder + img_name, pred)
    return 
