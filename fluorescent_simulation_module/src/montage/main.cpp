#include <iostream>
#include <climits>
#include <vector>
#include <algorithm>
#include <string>
#include <iterator>
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

#define SCALE 8
// -i TIFF -o GAUSSPASS

using namespace std;

void GetFilesName(const char* dirname, std::vector<std::string>& filenames)
{
	DIR* dir;
    struct dirent* ptr;
    
    dir = opendir(dirname);
	
	while((ptr = readdir(dir)) != NULL){
		if(strcmp(ptr->d_name, ".") == 0){
			continue;
		}
		if(strcmp(ptr->d_name, "..") == 0){
			continue;
		}
		
        std::string name = ptr->d_name;
		filenames.push_back(name);
	}

    closedir(dir);
}

void split(const std::string &s, char delim, int* result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = atoi(item.c_str());
    }
}

void AveragePatch(const char* dir, IplImage** output){
	std::vector<std::string> filenames;
	
	GetFilesName(dir, filenames);
	
	IplImage* tmp = cvLoadImage((std::string(dir)+"/"+filenames[0]).c_str(), CV_LOAD_IMAGE_UNCHANGED);
	IplImage* avgimg = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, 1);
	cvFillImage(avgimg, 0);
	cvReleaseImage(&tmp);

	std::cout<<(std::string(dir))<<std::endl;
	
	for(int i = 0; i < filenames.size(); i++){
// 		std::cout<<(std::string(dir)+"/"+filenames[i]).c_str()<<std::endl;
		tmp = cvLoadImage((std::string(dir)+"/"+filenames[i]).c_str(), CV_LOAD_IMAGE_UNCHANGED);
		IplImage* ori = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, 1);
		
		cvCvtScale(tmp, ori, 1, 0);
		cvReleaseImage(&tmp);

		cvScaleAdd(ori, cvScalar(1.0/filenames.size()), avgimg, avgimg);
		cvReleaseImage(&ori);
	}
	*output = avgimg;
}



int main(int argc, char **argv){
    opts opt;
	opt.output[0] = '\0';
	
	GetOpts(argc, argv, &opt);
	
	std::vector<std::string> dirnames;
	
	GetFilesName(opt.input, dirnames);
	
	int width = -1, height = -1;
	
	for(int i = 0; i < dirnames.size(); i++){
		int p[4];
		split(dirnames[i], '-', p);
		if(p[0]+p[2] > width){
			width = p[0]+p[2];
		}
		if(p[1]+p[3] > height){
			height = p[1]+p[3];
		}
	}
	width += 1; height += 1;
	width *= SCALE; height *= SCALE;
	
	IplImage* canvas = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
	IplImage* weight = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
	cvFillImage(canvas, 0);
	cvFillImage(weight, 0);
	IplImage* patch;
	
	for(int i = 0; i < dirnames.size(); i++){
		AveragePatch(std::string(opt.input+std::string("/")+dirnames[i]).c_str(), &patch);
		int p[4];
		split(dirnames[i], '-', p);
		CvRect roi; roi.width = p[2]*SCALE; roi.height = p[3]*SCALE;
		roi.x = p[0]*SCALE;
		roi.y = p[1]*SCALE;
		
		cvSetImageROI(canvas,roi);
		cvAdd(canvas, patch, canvas);
		cvResetImageROI(canvas);
		
		cvFillImage(patch, 1);
		cvSetImageROI(weight,roi);
		cvAdd(weight, patch, weight);
		cvResetImageROI(weight);
		
		cvReleaseImage(&patch);
	}
	
	cvDiv(canvas, weight, canvas);
	
	cvSaveImage(opt.output, canvas);
	
	cvReleaseImage(&canvas);
	cvReleaseImage(&weight);
	
    return 0;
}
