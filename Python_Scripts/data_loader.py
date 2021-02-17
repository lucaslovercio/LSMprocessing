import itertools 
import os 
import random 
import numpy as np 
import cv2 
from augmentation import augment_pair 

#From
#https://github.com/divamgupta/image-segmentation-keras/blob/master/keras_segmentation/data_utils/data_loader.py
#
 
class DataLoaderError(Exception): 
    pass 
 
 
def get_pairs_from_paths(images_path, segs_path): 
    """ Find all the images from the images_path directory and 
        the segmentation images from the segs_path directory 
        while checking integrity of data """ 
 
    ACCEPTABLE_IMAGE_FORMATS = [".jpg", ".jpeg", ".png", ".bmp"] 
    ACCEPTABLE_SEGMENTATION_FORMATS = [".png", ".bmp"] 
 
    image_files = [] 
    segmentation_files = {} 
 
    for dir_entry in os.listdir(images_path): 
        if os.path.isfile(os.path.join(images_path, dir_entry)) and \
           os.path.splitext(dir_entry)[1] in ACCEPTABLE_IMAGE_FORMATS:
            file_name, file_extension = os.path.splitext(dir_entry) 
            image_files.append((file_name, file_extension, 
                                os.path.join(images_path, dir_entry))) 
 
    for dir_entry in os.listdir(segs_path): 
        if os.path.isfile(os.path.join(segs_path, dir_entry)) and \
           os.path.splitext(dir_entry)[1] in ACCEPTABLE_SEGMENTATION_FORMATS:
            file_name, file_extension = os.path.splitext(dir_entry) 
            full_dir_entry = os.path.join(segs_path, dir_entry) 
            if file_name in segmentation_files: 
                raise DataLoaderError("Segmentation file with filename {0}" 
                                      " already exists and is ambiguous to" 
                                      " resolve with path {1}." 
                                      " Please remove or rename the latter." 
                                      .format(file_name, full_dir_entry)) 
            segmentation_files[file_name] = (file_extension, full_dir_entry) 
 
    return_value = [] 
    # Match the images and segmentations 
    for image_file, _, image_full_path in image_files: 
        if image_file in segmentation_files: 
            return_value.append((image_full_path, 
                                segmentation_files[image_file][1])) 
        else: 
            # Error out 
            raise DataLoaderError("No corresponding segmentation " 
                                  "found for image {0}." 
                                  .format(image_full_path)) 
 
    return return_value

# this method is like get_pairs_from_paths but doesn't require any masks
def get_images_from_path(images_path):
     
    ACCEPTABLE_IMAGE_FORMATS = [".jpg", ".jpeg", ".png", ".bmp"] 
    image_files = [] 
 
    for dir_entry in os.listdir(images_path): 
        if os.path.isfile(os.path.join(images_path, dir_entry)) and \
           os.path.splitext(dir_entry)[1] in ACCEPTABLE_IMAGE_FORMATS:
            image_files.append(os.path.join(images_path, dir_entry))
            
    return image_files

 
# really as this method does is normalize and add a dimension, not load 
def get_image_array(img, norm_type): 
    """ Load image array from input """ 
    img = img.astype(np.float32) 
     
    if norm_type == 'divide': 
        img /= 255.0 
    elif norm_type == 'sub_mean': 
        img -= np.mean(img) 
    elif norm_type == 'divide_and_sub': 
        img /= 255.0 
        img -= np.mean(img) 
 
    img = img.reshape(img.shape + (1,)) 
 
    return img 
 
 
# really all this does is convert the segmentations to one-hot encodings 
def get_segmentation_array(img, n_classes, width, height): 
    """ Load segmentation array from input """ 
 
    seg_labels = np.zeros((height, width, n_classes)) 
 
    for c in range(n_classes):
        seg_labels[:, :, c] = (img == c).astype(int) 
 
    return seg_labels 
 
 
def image_segmentation_generator(images_path, segs_path, batch_size, n_classes, output_height, output_width, 
                                 norm_type=None, aug_type=None, deterministic=False): 
 
    if deterministic: 
        # if the generator is being used for visualization, fix the order it generates images in 
        random.seed(0) 
 
    img_seg_pairs = get_pairs_from_paths(images_path, segs_path) 
    random.shuffle(img_seg_pairs) 
    zipped = itertools.cycle(img_seg_pairs) 
 
    while True: 
        X = [] 
        Y = [] 
        for _ in range(batch_size): 
            img, mask = next(zipped) 
 
            img = cv2.imread(img, cv2.IMREAD_ANYDEPTH)#new changed from 0 
            mask = cv2.imread(mask, cv2.IMREAD_ANYDEPTH)#new changed from 0 
            # if the raw images are 16-bit, convert them to 8-bit for brightness augmentation and preprocessing by division 
            if img.dtype == np.uint16: 
                img = (img // 256).astype(np.uint8) 
            # apply augmentation 
            if aug_type is not None: 
                img, mask = augment_pair(img, mask, aug_type) 
             
            X.append(get_image_array(img, norm_type)) 
            Y.append(get_segmentation_array(mask, n_classes, output_width, output_height)) 
             
        yield np.array(X), np.array(Y) 
 
# feed segmentImage.prediction in batches so the normalization works 
# this function is now deprecated 
def image_generator(images_path, batch_size): 
    image_files = [] 
 
    for dir_entry in os.listdir(images_path): 
        if os.path.isfile(os.path.join(images_path, dir_entry)): 
            file_name, file_extension = os.path.splitext(dir_entry) 
            image_files.append(os.path.join(images_path, dir_entry)) 
 
    zipped = itertools.cycle(image_files) 
    while True: 
        X = [] 
        for _ in range(batch_size): 
            im = next(zipped) 
            im = cv2.imread(im, 0) 
            X.append(get_image_array(im)) 
 
        yield np.array(X)
