#! /usr/bin/python
# -*- coding: utf8 -*-

import os, time, pickle, random, time
from datetime import datetime
import numpy as np
from time import localtime, strftime
import logging, scipy
from sklearn.model_selection import KFold

import tensorflow as tf
import tensorlayer as tl
from model import *
from utils import *
from os.path import join

###====================== HYPER-PARAMETERS ===========================###

ni = int(np.sqrt(16))

def read_all_imgs(img_list, path='', n_threads=32):
    """ Returns all images in array by given path and name of each image file. """
    imgs = []
    for idx in range(0, len(img_list), n_threads):
        b_imgs_list = img_list[idx : idx + n_threads]
        b_imgs = tl.prepro.threading_data(b_imgs_list, fn=get_imgs_fn, path=path)
        # print(b_imgs.shape)
        imgs.extend(b_imgs)
        print('read %d from %s' % (len(imgs), path))
    return imgs

def read_all_imgs_lr(dir_list, path='', n_threads=32):
    imgs = []
    for idx in range(0, len(dir_list), n_threads):
        b_dir_list = dir_list[idx : idx + n_threads]
        b_imgs = tl.prepro.threading_data(b_dir_list, fn=get_imgs_fn_lr, path=path)
        # print(b_imgs.shape)
        imgs.extend(b_imgs)
        print('read %d from %s' % (len(imgs), path))
    return imgs


def test(i_path, out_path):
    tl.files.exists_or_mkdir(out_path)
    checkpoint_dir = "./model"

    file_list = sorted(tl.files.load_file_list(path=i_path, regx='.*', printable=False))
    file_list = list(map(lambda x: i_path+x, file_list))
    img = np.array([scipy.misc.imread(file, mode='L') for file in file_list])
    img = np.moveaxis(img, 0, -1)

    size = img.shape
    t_image = tf.placeholder('float32', [None, size[0], size[1], size[2]], name='input_image')
    net_g_6, net_g_4, net_g = SRGAN_g(t_image, is_train=False, reuse=False)

    ###========================== RESTORE G =============================###
    sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True, log_device_placement=False))
    tl.layers.initialize_global_variables(sess)
    tl.files.load_and_assign_npz(sess=sess, name=checkpoint_dir+'/g_srgan_epoch_6.npz', network=net_g)
    # tl.files.load_and_assign_npz(sess=sess, name=checkpoint_dir+'/g_srgan_init.npz', network=net_g)


    valid_lr_img = norm_imgs(img)
    for i in range(50):
        feed_dict = {t_image: [valid_lr_img]}
        feed_dict.update(net_g.all_drop)
        out, out_4, out_6 = sess.run([net_g.outputs, net_g_4.outputs, net_g_6.outputs], feed_dict=feed_dict)
        tl.vis.save_image(out[0], join(out_path,'8Xtest_gen_{}.png'.format(i)))
        # tl.vis.save_image(out_4[0], './test/4Xtest_gen_{}.png'.format(i))
        # tl.vis.save_image(out_6[0], './test/6Xtest_gen_{}.png'.format(i))
    valid_lr_img_rot = np.rot90(valid_lr_img, 1, (0,1))
    for i in range(50):
        feed_dict = {t_image: [valid_lr_img_rot]}
        feed_dict.update(net_g.all_drop)
        out, out_4, out_6 = sess.run([net_g.outputs, net_g_4.outputs, net_g_6.outputs], feed_dict=feed_dict)
        result = np.rot90(out[0], -1, (0,1))
        result_4 = np.rot90(out_4[0], -1, (0, 1))
        result_6 = np.rot90(out_6[0], -1, (0, 1))
        tl.vis.save_image(result, join(out_path,'8Xtest_gen_{}_rot.png'.format(i)))
        # tl.vis.save_image(result_4, './test/4Xtest_gen_{}_rot.png'.format(i))
        # tl.vis.save_image(result_6, './test/6Xtest_gen_{}_rot.png'.format(i))

def get_imgs_test(d, path):
    i_path = path+d
    file_list = sorted(tl.files.load_file_list(path=i_path, regx='.*', printable=False))
    
    file_list = list(map(lambda x: i_path+x, file_list))
    img = np.array([scipy.misc.imread(file, mode='L') for file in file_list])
    return np.moveaxis(img, 0, -1)

def read_all_imgs_test(dir_list, path='', n_threads=32):
    imgs = []
    for idx in range(0, len(dir_list), n_threads):
        b_dir_list = dir_list[idx : idx + n_threads]
        b_imgs = tl.prepro.threading_data(b_dir_list, fn=get_imgs_test, path=path)
        # print(b_imgs.shape)
        imgs.extend(b_imgs)
        print('read %d from %s' % (len(imgs), path))
    return imgs

def test_full(i_path, out_path):
    tl.files.exists_or_mkdir(out_path)
    checkpoint_dir = "./model"

    clip_list = sorted(os.listdir(i_path))
    clip_list_dir = map(lambda x: x+'/', clip_list)
    imgs = read_all_imgs_test(clip_list_dir, path=i_path, n_threads=32)

    size = imgs[0].shape
    t_image = tf.placeholder('float32', [None, size[0], size[1], size[2]], name='input_image')
    net_g_6, net_g_4, net_g = SRGAN_g(t_image, is_train=False, reuse=False)

    ###========================== RESTORE G =============================###
    sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True, log_device_placement=False))
    tl.layers.initialize_global_variables(sess)
    tl.files.load_and_assign_npz(sess=sess, name=checkpoint_dir+'/g_srgan_epoch_6.npz', network=net_g)
    # tl.files.load_and_assign_npz(sess=sess, name=checkpoint_dir+'/g_srgan_init.npz', network=net_g)

    for j in range(len(imgs)):
        img = imgs[j]
        valid_lr_img = norm_imgs(img)
        tl.files.exists_or_mkdir(join(out_path, '{}'.format(clip_list[j])))
        for i in range(10):
            feed_dict = {t_image: [valid_lr_img]}
            feed_dict.update(net_g.all_drop)
            out, out_4, out_6 = sess.run([net_g.outputs, net_g_4.outputs, net_g_6.outputs], feed_dict=feed_dict)
            tl.vis.save_image(out[0], join(out_path,'{}/8Xtest_gen_{}.png'.format(clip_list[j], i)))
        valid_lr_img_rot = np.rot90(valid_lr_img, 1, (0,1))
        for i in range(10):
            feed_dict = {t_image: [valid_lr_img_rot]}
            feed_dict.update(net_g.all_drop)
            out, out_4, out_6 = sess.run([net_g.outputs, net_g_4.outputs, net_g_6.outputs], feed_dict=feed_dict)
            result = np.rot90(out[0], -1, (0,1))
            result_4 = np.rot90(out_4[0], -1, (0, 1))
            result_6 = np.rot90(out_6[0], -1, (0, 1))
            tl.vis.save_image(result, join(out_path,'{}/8Xtest_gen_{}_rot.png'.format(clip_list[j], i)))

if __name__ == '__main__':
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', type=str, default='test', help='test, test_full')
    parser.add_argument('--input', type=str, required=True, help='The input path of the low resolution images.')
    parser.add_argument('--output', type=str, required=True, help='The output path of the super resolution images.')
    args = parser.parse_args()
    # print(args.mode)
    # print(args.input)

    tl.global_flag['mode'] = args.mode

    if tl.global_flag['mode'] == 'test':
        test(args.input, args.output)
    elif tl.global_flag['mode'] == 'test_full':
        test_full(args.input, args.output)
    else:
        raise Exception("Unknow --mode")
