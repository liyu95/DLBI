#! /usr/bin/python
# -*- coding: utf8 -*-

import tensorflow as tf
import tensorlayer as tl
from tensorlayer.layers import *

class RMLayer(Layer):
    def __init__(
        self,
        layer = None,
        name ='reduce_mean_layer',
    ):
        # check layer name (fixed)
        Layer.__init__(self, name=name)

        # the input of this layer is the output of previous layer (fixed)
        self.inputs = layer.outputs

        # operation (customized)
        self.outputs =  tf.reduce_mean(self.inputs, axis=-1, keep_dims=True)


        # get stuff from previous layer (fixed)
        self.all_layers = list(layer.all_layers)
        self.all_params = list(layer.all_params)
        self.all_drop = dict(layer.all_drop)

        # update layer (customized)
        self.all_layers.extend( [self.outputs] )


def SRGAN_g(t_image, is_train=False, reuse=False):
    """ Generator in Photo-Realistic Single Image Super-Resolution Using a Generative Adversarial Network
    feature maps (n) and stride (s) feature maps (n) and stride (s)
    """
    w_init = tf.random_normal_initializer(stddev=0.02)
    b_init = None # tf.constant_initializer(value=0.0)
    g_init = tf.random_normal_initializer(1., 0.02)
    with tf.variable_scope("SRGAN_g", reuse=reuse) as vs:
        tl.layers.set_name_reuse(reuse)
        n = InputLayer(t_image, name='in')
        n = Conv2d(n, 256, (7, 7), (1, 1), act=tf.nn.relu, padding='SAME', W_init=w_init, name='n128s1/c')

        # n = FlattenLayer(n, name='flatten_layer')
        n = DropoutLayer(n, keep=0.8, is_train=True, name='dropout_layer_1')
        # n = ReshapeLayer(n, shape=[-1, 60, 60, 256], name='reshape_layer')

        # The original short cut
        # temp = n

        # The denoise short cut
        temp = RMLayer(n)
        temp = TileLayer(temp, [1,1,1,256], name='tile')

        # B residual blocks
        for i in range(16):
            nn = Conv2d(n, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, b_init=b_init, name='n128s1/c1/%s' % i)
            nn = BatchNormLayer(nn, act=tf.nn.relu, is_train=is_train, gamma_init=g_init, name='n128s1/b1/%s' % i)
            nn = Conv2d(nn, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, b_init=b_init, name='n128s1/c2/%s' % i)
            nn = BatchNormLayer(nn, is_train=is_train, gamma_init=g_init, name='n128s1/b2/%s' % i)
            nn = ElementwiseLayer([n, nn], tf.add, 'b_residual_add/%s' % i)
            n = nn

        n = Conv2d(n, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, b_init=b_init, name='n128s1/c/m')
        n = BatchNormLayer(n, is_train=is_train, gamma_init=g_init, name='n128s1/b/m')
        n = ElementwiseLayer([n, temp], tf.add, 'add3')
        # B residual blocks end

        n = Conv2d(n, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, name='n256s1/1')
        n = SubpixelConv2d(n, scale=2, n_out_channel=None, act=tf.nn.relu, name='pixelshufflerx2/1')

        # n_2 = Conv2d(n, 1, (1, 1), (1, 1), act=tf.nn.tanh, padding='SAME', W_init=w_init, name='2x_out')
        n_6_temp = Conv2d(n, 270, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, name='n300s1/1')
        n_6_temp = SubpixelConv2d(n_6_temp, scale=3, n_out_channel=None, act=tf.nn.relu, name='pixelshufflerx3/1')
        n_6 = Conv2d(n_6_temp, 1, (1, 1), (1, 1), act=tf.nn.tanh, padding='SAME', W_init=w_init, name='6x_out')

        n = Conv2d(n, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, name='n256s1/2')
        n = SubpixelConv2d(n, scale=2, n_out_channel=None, act=tf.nn.relu, name='pixelshufflerx2/2')

        n_4 = Conv2d(n, 1, (1, 1), (1, 1), act=tf.nn.tanh, padding='SAME', W_init=w_init, name='4x_out')

        n = Conv2d(n, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, name='n256s1/3')
        n = SubpixelConv2d(n, scale=2, n_out_channel=None, act=tf.nn.relu, name='pixelshufflerx2/3')

        n_8 = Conv2d(n, 1, (1, 1), (1, 1), act=tf.nn.tanh, padding='SAME', W_init=w_init, name='8x_out')

        return n_6, n_4, n_8

