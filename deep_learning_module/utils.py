import tensorflow as tf
import tensorlayer as tl
from tensorlayer.prepro import *
import scipy
import numpy as np


def get_imgs_fn(file_name, path):
    """ Input an image path and name, return an image array """
    # return scipy.misc.imread(path + file_name).astype(np.float)
    return scipy.misc.imread(path + file_name, mode='L')

def get_imgs_fn_lr(d, path):
    i_path = path+d
    file_list = sorted(tl.files.load_file_list(path=i_path, regx='.*.tiff', printable=False))
    # file_list = sorted(tl.files.load_file_list(path=i_path, regx='.*.TIFF', printable=False))
    file_list = list(map(lambda x: i_path+x, file_list))
    img = np.array([scipy.misc.imread(file, mode='L') for file in file_list])
    return np.moveaxis(img, 0, -1)


def crop_sub_imgs_fn(x, is_random=True):
    x = crop(x, wrg=384, hrg=384, is_random=is_random)
    x = x / (255. / 2.)
    x = x - 1.
    return x

def norm_fn(x):
    x = x / (255. / 2.)
    x = x - 1.
    return x

def norm_img(img, img_max, img_min):
    if img_max == img_min:
        return norm_fn(img)
    img = 2.0*((img-img_min)/(float(img_max)-img_min)) - 1.0
    return img

def norm_imgs(imgs):
    imgs = np.moveaxis(imgs, -1, 0)
    img_max = np.max(imgs[0])
    img_min = np.min(imgs[0])
    imgs = map(lambda x: norm_img(x, img_max, img_min),
        list(imgs))
    # imgs = map(lambda x: norm_single(x),
    #     list(imgs))    
    imgs = np.array(imgs)
    return np.moveaxis(imgs, 0, -1)

def norm_single(img):
    # print(img.shape)
    img_max = np.max(img)
    img_min = np.min(img)
    if img_max == img_min:
        return norm_fn(img)
    img = 2.0*((img-img_min)/(float(img_max)-img_min)) - 1.0
    return img

def downsample_fn_4X(x):
    # We obtained the LR images by downsampling the HR images using bicubic kernel with downsampling factor r = 2.
    x = imresize(x, size=[240, 240], interp='bicubic', mode=None)
    x = norm_single(x)
    return x

def downsample_fn_2X(x):
    # We obtained the LR images by downsampling the HR images using bicubic kernel with downsampling factor r = 4.
    x = imresize(x, size=[120, 120], interp='bicubic', mode=None)
    x = norm_single(x)
    return x