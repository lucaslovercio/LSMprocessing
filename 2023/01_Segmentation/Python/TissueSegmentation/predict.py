from tensorflow.keras.models import load_model, Model
import tensorflow as tf
import numpy as np
import cv2
from TissueSegmentation.data_loader import get_image_array, get_pairs_from_paths, get_images_from_path
import os


#TODO
'''
- add documentation
- save loss curves w/ model somehow
- save text w/ learning curves and test images
'''

'''
Name: prediction
Purpose:
This function takes a preloaded, pretrained model, and uses it to create a segmentation for a 
given input image.
b) role of each formal parameter
Arguments:
- model:
 This is a preloaded model which has already been trained and is not compiled at load 
 time.
- image_path:
 The path to the raw image to segment.
- norm_type:
 One of 'divide', 'sub_mean', and 'divide_and_sub'. A string describing how to normalize the 
 images, before feeding them to the network.
Return Value:
This function returns a segmentation image as a numpy array.
'''
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
def segment_folder(model, frame_path, output_folder, filesep = '/'):
    imgs = get_images_from_path(frame_path)
    norm_type = None
    if 'divide_and_sub' in model:
        norm_type = 'divide_and_sub'
    elif 'sub_mean' in model:
        norm_type = 'sub_mean'
    elif 'divide' in model:
        norm_type = 'divide'

    #Added to clean GPU memory after using it
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)

    model = load_model(model, compile=False)

    for img in imgs:
        img_name = img[img.rfind('/') + 1:]#this gets the image name out of the path
        pred = prediction(model, img, norm_type)
        cv2.imwrite(output_folder + filesep + img_name, pred)
    return 


'''
Name: assess_model
Purpose:
This function takes a pretrained model and a set of test images and produces a set of images
showing the model's output segmentation next to the ground truth.
Arguments: 
- model:
 The path to a pretrained model (a string).
- frame_path:
 The path to the raw images for the model to segment.
- mask_path:
 The path to the ground truth segmentations for the images the model will segment.
- output_folder:
 The directory to save the output images to.
Outputs:
This function saves the segmentation/ground truth images titled the same as the raw images 
to the specified output directory.
'''
def assess_model(model, frame_path, mask_path, output_folder):
    print("enter to assess_model")
    print("model name" + model)
    img_seg_pairs = get_pairs_from_paths(frame_path, mask_path)
    #load the model
    #get normalization preprocessing type that the model was trained with
    norm_type = None
    if 'divide_and_sub' in model:
        norm_type = 'divide_and_sub'
    elif 'sub_mean' in model:
        norm_type = 'sub_mean'
    elif 'divide' in model:
        norm_type = 'divide'
    
    model = load_model(model, compile=False)
    model.summary()
    for img, mask in img_seg_pairs:
        print("entered to for in imgs")
        img_name = img[img.rfind('/') + 1:]#this gets the image name out of the path
        pred = prediction(model, img, norm_type)
        mask = cv2.imread(mask, 0) * 50
        output = np.concatenate((pred, mask), axis=1)
        print("saved in " + output_folder + img_name)
        cv2.imwrite(output_folder + img_name, output)


'''
Name: assess_models
Purpose:
Go through every pretrained model in a directory and assess it's performance on a given 
training set.
Arguments:
- models_dir:
 The directory containing all pretrained models.
- frame_path:
 The path to the raw images for the models to segment.
- mask_path:
 The path to the ground truth segmentations for the images the models will segment.
- outputs_dir:
 The top-level directory to save outputs to. The outputs for each model will be saved in a 
 subdirectory here.
Outputs:
This function outputs the assessment of every model in a directory in a subdirectory 
created in the specified output directory. The subdirectories are named after the models
and will be created if they don't already exist.
'''
def assess_models(models_dir, frame_path, mask_path, outputs_dir):
    output_folders = []
    model_paths = []
    for file in os.listdir(models_dir):
        if file.endswith('.h5'):
            print("H5 found: " + file)
            #make the subdirectory for each model
            model_subdir = outputs_dir + file[:-3] + '/'
            try:
                os.mkdir(model_subdir)
                print("Created: " + model_subdir)
            except FileExistsError:
                pass
            output_folders.append(model_subdir)
            #save the path to the model
            model_paths.append(models_dir + file)
    #assess the models
    for i in range(len(model_paths)):
        assess_model(model_paths[i], frame_path, mask_path, output_folders[i])

def assess_models_in_folders(top_dir, frame_path, mask_path, outputs_dir):
    model_dirs = []
    for file in os.listdir(top_dir):
        path = os.path.join(top_dir, file)
        if os.path.isdir(path):
            assess_models(path + '/', frame_path, mask_path, outputs_dir)
