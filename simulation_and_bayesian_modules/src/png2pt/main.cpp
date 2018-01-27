#include <iostream>
#include <climits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <time.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <stdio.h>
#include "img_util.h"
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "opts.h"

// -i TIFF -o GAUSSPASS

using namespace std;

int main(int argc, char **argv){
    opts opt;
	opt.output[0] = '\0';
	
	GetOpts(argc, argv, &opt);
	
	IplImage* ori = cvLoadImage(opt.input, CV_LOAD_IMAGE_GRAYSCALE);
	IplImage* lawpass = cvCreateImage(cvGetSize(ori), ori->depth, ori->nChannels);
	cvSmooth(ori, lawpass, CV_GAUSSIAN, 0, 0, 3);
	cvScaleAdd(lawpass, cvScalar(-1), ori, ori);
	
	std::ofstream out(opt.output);
	
	double minimum, maximum;
	cvMinMaxLoc(ori, &minimum, &maximum);
	
	CvScalar mean, sdv;
	cvAvgSdv(ori, &mean, &sdv);
	
	for(int y = 0; y < ori->height; y++){
		uchar* ptr = (uchar*)(ori->imageData+y*ori->widthStep);
		for(int x = 0; x < ori->width; x++){
			if(*ptr < mean.val[0]){
				*ptr = 0;
			}
			else{
				*ptr = int(sqrt((*ptr-minimum)/maximum)*8);
			}
			ptr++;
		}
	}
	
	srand48(time(NULL));
	
	for(int y = 0; y < ori->height; y += 2){
		uchar* ptr = (uchar*)(ori->imageData+y*ori->widthStep);
		for(int x = 0; x < ori->width; x += 2){
			int val = int(*(ptr+x));
			if(val > 1){
				out<<x<<"\t"<<y<<std::endl;
				for(int c = val-1; c--;){
					out<<x+(8*drand48()-4)<<"\t"<<y+(8*drand48()-4)<<std::endl;
				}
			}
		}
	}
	
    return 0;
}
