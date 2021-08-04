import numpy as np
import imgaug as ia
from imgaug import augmenters as iaa
from imgaug.augmentables.segmaps import SegmentationMapsOnImage
#TODO
'''

'''

# simple augmenter to perform horizontal and vertical flips on input at a rate of 50%
flip_aug = iaa.Sequential([
    iaa.Fliplr(0.5),
    iaa.Flipud(0.5)])

# simple augmenter to apply Guassian blurring
blur_aug = iaa.OneOf([
    iaa.GaussianBlur(sigma=(0,2.5)),
    iaa.GaussianBlur(sigma=(0,1.0))
    ])

# simple augmenter to vary the brightness level of images
brightness_aug = iaa.pillike.EnhanceBrightness((0.9, 1.1))

# augmenter the applies multiple affine translations
trans_aug = iaa.Affine(shear=(-16, 16),# rotate=(0, 360),
                       translate_percent={"x": (-0.10, 0.10), "y": (-0.10, 0.10)},
                       mode='reflect')
trans_aug._mode_segmentation_maps = 'reflect'

# simple augmenter to apply smooth elastic deformations
distort_aug = iaa.PiecewiseAffine(scale=(0.01, 0.075))

# simple augmenter to adjust the image contrast
contrast_aug = iaa.GammaContrast(gamma=(0.9, 1.1))

# compound augmenter that applies geometric augmenters and sporadically blurring and brightness
full_aug = iaa.Sequential([
    flip_aug, trans_aug, distort_aug,
    iaa.Sometimes(0.25, brightness_aug), iaa.Sometimes(0.25, contrast_aug), iaa.Sometimes(0.1, blur_aug)
])

full_aug_without_distortion = iaa.Sequential([
    flip_aug, trans_aug,
    iaa.Sometimes(0.25, brightness_aug), iaa.Sometimes(0.25, contrast_aug), iaa.Sometimes(0.1, blur_aug)
])

# a lookup dictionary for the different augmenters
augmenter_index = {
    'flips': flip_aug,
    'blur': blur_aug,
    'brightness': brightness_aug,
    'translation': trans_aug,
    'distortion': distort_aug,
    'contrast': contrast_aug, 
    'all': full_aug,
    'distortionless': full_aug_without_distortion
}

# takes a pair of images and applies a given augmenter to them, returning the result
def augment_pair(img, mask, augmenter):
    augmenter = augmenter_index[augmenter]
    mask = SegmentationMapsOnImage(mask, img.shape)
    img_aug, mask_aug = augmenter(image=img, segmentation_maps=mask)
    mask_aug = SegmentationMapsOnImage.get_arr(mask_aug)
    return img_aug, mask_aug
