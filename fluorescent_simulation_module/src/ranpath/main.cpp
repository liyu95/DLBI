#include <iostream>
#include <climits>
#include <vector>
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "opts.h"

using namespace std;

inline IplImage* GetSubImage(IplImage *image, CvRect roi)
{
    IplImage *result;
	
	if(roi.x < 0 || roi.x+roi.width >= image->width || roi.y < 0 || roi.y+roi.height >= image->height){
		return NULL;
	}
	
    cvSetImageROI(image,roi);
	CvRect roi2 = cvGetImageROI(image);
	if(roi2.width != roi.width || roi2.height != roi.height){
		cvResetImageROI(image);
		return NULL;
	}
    result = cvCreateImage( cvSize(roi.width, roi.height), image->depth, image->nChannels );
    cvCopy(image, result);
    cvResetImageROI(image);
    return result;
}

int main(int argc, char **argv) {
    opts opt;
	opt.output[0] = '\0';
	opt.size = 480;
	
    if(GetOpts(argc, argv, &opt) != 1) {
        return 0;
    }
	
	IplImage* tmp = cvLoadImage(opt.input, CV_LOAD_IMAGE_GRAYSCALE);
	IplImage* ori = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, 1);
	cvCvtScale(tmp, ori, 1, 0);
	cvReleaseImage(&tmp);
	
	double minimum, maximum;
	cvMinMaxLoc(ori, &minimum, &maximum);
// 	double dif = maximum-minimum;
	
	for(int y = 0; y < ori->height; y++){
		float* ptr = (float*)(ori->imageData+y*ori->widthStep);
		for(int x = 0; x < ori->width; x++){
			*ptr = (maximum-(*ptr-minimum)); // /dif*255;
			ptr++;
		}
	}
	
	CvRect region;
	region.x = 10000000;
	region.y = 10000000;
	region.width = -1;
	region.height = -1;
	
	for(int y = 0; y < ori->height; y++){
		float* ptr = (float*)(ori->imageData+y*ori->widthStep);
		for(int x = 0; x < ori->width; x++){
			if( *ptr > 0){
				if(region.x > x){
					region.x = x;
				}
				if(region.y > y){
					region.y = y;
				}
				
				if(region.width < x){
					region.width = x;
				}
				if(region.height < y){
					region.height = y;
				}
			}
			ptr++;
		}
	}

	region.width -= region.x;
	region.height -= region.y;
	
// 	region.width = opt.size;
// 	region.height = opt.size;
	
	if(region.width >= opt.size){
		region.width = opt.size;
	}
	if(region.height >= opt.size){
		region.height = opt.size;
	}
	
	IplImage* sub = GetSubImage(ori, region);
	cvReleaseImage(&ori);
	ori = sub;
	
	IplImage* laps = cvCreateImage(cvSize(opt.size, opt.size), ori->depth, ori->nChannels);
	
	cvResize(ori, laps, CV_INTER_CUBIC);
	
// 	cvDilate(ori, laps, NULL, 2);
// 	cvErode(laps, laps, NULL, 1);
	
	cv::RNG rnger(cv::getTickCount());
	
	float ratio = rnger.uniform(0.3, 1.0);
	
	for(int y = 0; y < laps->height; y++){
		float* ptr = (float*)(laps->imageData+y*laps->widthStep);
		for(int x = 0; x < laps->width; x++){
			*ptr = *ptr*(ratio+rnger.uniform(0.5, 1.0));
			ptr++;
		}
	}

	float* ptr = (float*)(laps->imageData);
	*ptr = 255;
	
	cvMinMaxLoc(laps, &minimum, &maximum);
	
	for(int y = 0; y < laps->height; y++){
		float* ptr = (float*)(laps->imageData+y*laps->widthStep);
		for(int x = 0; x < laps->width; x++){
			*ptr = (*ptr-minimum)/maximum*255;
			ptr++;
		}
	}

	cvSaveImage(opt.output, laps);//laps);
	
	cvReleaseImage(&ori);
	cvReleaseImage(&laps);
	
    return 0;
}
