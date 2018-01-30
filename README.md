# Sum up
**There are three modules in our framework. Our framework contains threes parts: data simulation, deep neural network and Bayesian inference. **

Here, the simulation and deep neural network implementation is openly authorization, and the third part is partially authorization. The first part is data simulation module, whose source can be find in folder fluorescent_simulation_module. The second part is the deep learning module, who can be find in folder deep_learning_module. The third part is the Bayesian inference module, users are encouraged to contact with the authors for further authorization.

If you think the related works here are helpful, please cite "DLBI: Deep learning guided Bayesian inference for structure reconstruction of super-resolution fluorescence microscopy" and "Live cell single molecule-guided Bayesian localization super resolution microscopy".

# Prerequisites:
### For simulation module
1. cmake
2. gcc
3. g++
4. OpenCV2.4 lib

### For deep learning
1. Anaconda 2 or Miniconda 2 (https://conda.io/miniconda.html)
2. CUDA 8 and cudNN 5.1 (https://www.tensorflow.org/versions/r1.2/install/install_linux)

# Compile the simulation and Bayesian modules
```
cd ./fluorescent_simulation_module
mkdir build
cd ./build
cmake ..
make -j 2
```
A "bin" dir will be generated and the following exes will be created: calib  lapls  montage  png2pt  ranpath  sravg  srsimu  subtract


# Perform simulation:
There are three exes to do simulation (the related files can be find in the folder data_ill).

lapls: input an color image and output a gray version of the laplace filtered color image (6x), for example:
 ```
./lapls -i 0037.png -o test.pgm -s 480 -m 0
 ```
Here, test.pgm is a generated gray image (-m to say which part for translation, 0 center, 1 2 3 4 four corner).


ranpath: input a sketch image and output a size adjusted gray image.
```
./ranpath -i ./n02694662_17579-5.png -o sketch.pgm -s 480
```
Here, sketch.pgm is a rescaled sketch image (gray image).

The lapls and ranpath is used for the preparation of gray-scale images, which will be used as the fluorophore distribution and density in further simulation. 
 
srsimu: input a high resolution gray image (fluorophore distribution) and output a series of low resolution images (here the log-norm and transmitting parameters are fixed, but easy to be changed in code). srsimu uses a fully random parameter initialization.
```
./srsimu -i test.pgm -o test -n 200 -s 0.125 -p 8
```
Here, -n (number of images), -s (downsampling scale) -p (used to adjust psf width, 8 is in default).

Using these three exes, we can generate simulation data for deep learning training (in different binned size).

We totally generated about 4000 datasets from Imagenet and 8000 datasets from Sketchy. These simulated datasets are used for deep learning training.


# Get deep learning result:

### Preprocessing
The current neural network is designed for 60x60 low-resolution fluorescent images with 200 frames (these can be changed by re-training). We need to perform the following process to get the input to deep learning module.

There are several exes for the process of real data. The deep neural network need inputs of a set of image frames (If users have a image series in tif format, it can be re-saved to image sequences by ImageJ). Here, use the cell6G_ori as an example.

calib: crop a certain area and do Gaussian high pass filter. For a small patch, an example is as following:

./calib -i cell6G_ori -o cellG_60_60_60_60 -s 3 -r 60,60,60,60

The input is the cell1G_ori directory (containing 200 image frames); The output is the cellG_60_60_60_60 (containing 200 high-passed image frames with Gaussian sigma=3) and the region is 60,60,60,60  (X,Y,WIDTH,HEIGHT). Then, the cellG_60_60_60_60 can be inputted to the trained neural network to generate an initial guess of the fluorophore distribution, which can be used further refined by Bayesian technique.

If not assigned the region, it is the whole image. For the whole large field area, an example is as following:

./calib -i cell6G_ori -o cellG_GAUSS -s 3

The input is the cell1G_ori directory (containing 200 image frames); The output is the cellG_GAUSS (containing 200 high-passed image frames with Gaussian sigma=3). Here, because our network process 60x60 area each time, the large-field input should be cropped into small patches. We provide the crop.sh as an example to change a series of large-field image frames to a set of series of image frames with small patches (cellG_GAUSS as input and cellG_clip as output).

"cellG_clip" dir will contain 0-0-60-60 0-40-60-60... subdirs. Each of the subdirs will contains a 200 frames of images with 60x60 pixels. cellG_clip is the input of neural network for large-field datasets.


### Build a virtual environment (Only need to be done once)
```
conda env create -f environment.yml
```

### Usage
```
<!-- Activate the environment -->
source activate tensorflow
<!-- Run the code -->
python main.py --mode test --input path/to/lr/images
<!-- Deactivate the environment -->
source deactivate
```

### Explanation of path
There are two modes: test and test_full. In the test mode, we only need to construct one image from 200 lr images. All the 200 images should be put into the same folder. In the test_full mode, we would perform a lot of such tasks. Each 200 images would have one separate folder, which is under the folder we give the program as the input argument. 

For example, if the "cellG_clip" dir is used as the input of the test_full mode, a folder like cell6_result will be outputted. Then, with the help of previous compiled exe "montage", a large-field super-resolution image will be outputted:
```
./montage -i cell6_result -o cell6-montage.pgm
```
### Input images
Each low resolution image should be one channel, 60 by 60 grayscale image. The output will be the images with 8x super-resolution.


# Bayesian module
The Bayesian module will takes the original images (for example, cell1G_ori), the deep learning super-resolution output (for example, cell6-montage.pgm), the mask file (region of interest) and the configuration parameters as input and ouput a further refined super-resolution image. Here, the core implementation of statistic inference is developed under the kind help of our collaborators. Users are encouraged to contact with the authors for further authorization.

# Final results

For the demonstration of final results, people can access the https://drive.google.com/file/d/1DOMv2Fo9WbIyz7-biWEyVybqknJ65lid/view?usp=sharing (for people who google is blocked, please visit https://github.com/icthrm/DLBI_results) to gain the full-resolution super-resolution images, videos and the used datasets. 
• real1.avi: Movie of the first real-world dataset (Actin1).

• real2.avi: Movie of the second real-world dataset (Actin2).

• real3.avi: Movie of the third real-world dataset (ER).

• real1.png: Large-field reconstruction of the yellow area from the first real-world dataset (Actin1).

• real2.png: Large-field reconstruction of the yellow area from the second real-world dataset (Actin2).

• real3.png: Large-field reconstruction of the yellow area from the third real-world dataset (ER).

• CC.png: Overlap of the reconstruction images by DLBI (in red) and PALM (in green) (Actin1).

And real1.tif, real2.tif and real3.tif are the corresponding original high-density fluorecent images.

Please note that the real datasets real1.tif, real2.tif and real3.tif are distributed under providing of Professor Xu and Professor Zhang (http://ear.ict.ac.cn/?page_id=207). People who what to use them alone should also cite the paper "Rational design of true monomeric and bright photoactivatable fluorescent proteins". 

