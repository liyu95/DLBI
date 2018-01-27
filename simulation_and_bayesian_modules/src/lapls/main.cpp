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
	opt.mode = 0;
	
    if(GetOpts(argc, argv, &opt) != 1) {
        return 0;
    }
	
	IplImage* ori = cvLoadImage(opt.input, CV_LOAD_IMAGE_GRAYSCALE);
	IplImage* mag = cvCreateImage(cvSize(ori->width*6,ori->height*6), ori->depth, ori->nChannels);		//enlarge size
	cvResize(ori,mag, CV_INTER_CUBIC);
	cvReleaseImage(&ori);
	ori = mag;
		
	CvRect region;
	switch(opt.mode){
	case 0: 		//middle
		region.x = (ori->width-opt.size)/2;
		region.y = (ori->height-opt.size)/2;
		break;
		
	case 1:
		region.x = 0;
		region.y = 0;
		break;
		
	case 2:
		region.x = ori->width-opt.size-1;
		region.y = 0;
		break;
		
	case 3:
		region.x = 0;
		region.y = ori->height-opt.size-1;
		break;
	case 4:
		region.x = ori->width-opt.size-1;
		region.y = ori->height-opt.size-1;
		break;
	}
	region.width = opt.size;
	region.height = opt.size;
	
	IplImage* sub = GetSubImage(ori, region);
	cvReleaseImage(&ori);
	ori = sub;
	
	IplImage* laps = cvCreateImage(cvGetSize(ori), ori->depth, ori->nChannels);
	
	cvSmooth(ori, ori, CV_GAUSSIAN, 5);
	cvLaplace(ori, laps, 5);
	cvErode(laps, laps, NULL, 1);
	cvSmooth(ori, ori, CV_GAUSSIAN, 3);
	cvSmooth(ori, ori, CV_GAUSSIAN, 3);

	cvSaveImage(opt.output, laps);
	
	cvReleaseImage(&ori);
	cvReleaseImage(&laps);
	
    return 0;
}
