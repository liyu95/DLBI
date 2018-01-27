# Prerequisites:
1. Anaconda 2 or Miniconda 2 (https://conda.io/miniconda.html)
2. CUDA 8 and cudNN 5.1 (https://www.tensorflow.org/versions/r1.2/install/install_linux)

# Build a virtual environment
```
conda env create -f environment.yml
```

# Usage
```
<!-- Activate the environment -->
source activate tensorflow
<!-- Run the code -->
python main.py --mode test --input path/to/lr/images
<!-- Deactivate the environment -->
source deactivate
```

# Explanation of path
There are two modes: test and test_full. In the test mode, we only need to construct one image from 200 lr images. All the 200 images should be put into the same folder. In the test_full mode, we would perform a lot of such tasks. Each 200 images would have one separate folder, which is under the folder we give the program as the input argument.

# Input images
Each low resolution image should be one channel, 60 by 60 grayscale image.