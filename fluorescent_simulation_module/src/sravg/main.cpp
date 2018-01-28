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

int main(int argc, char **argv){
    opts opt;
	opt.output[0] = '\0';
	
	GetOpts(argc, argv, &opt);
	
	std::vector<std::string> filenames;
	
	GetFilesName(opt.input, filenames);
	
	std::sort(filenames.begin(), filenames.end());
	
	std::cout<<filenames.size()<<std::endl;
	
	IplImage* tmp = cvLoadImage((std::string(opt.input)+"/"+filenames[0]).c_str(), CV_LOAD_IMAGE_UNCHANGED);
	IplImage* avgimg = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, 1);
	cvFillImage(avgimg, 0);
	cvReleaseImage(&tmp);
	
	for(int i = 0; i < filenames.size(); i++){
		std::cout<<(std::string(opt.input)+"/"+filenames[i]).c_str()<<std::endl;
		tmp = cvLoadImage((std::string(opt.input)+"/"+filenames[i]).c_str(), CV_LOAD_IMAGE_UNCHANGED);
		IplImage* ori = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, 1);
		
		cvCvtScale(tmp, ori, 1, 0);
		cvReleaseImage(&tmp);

		cvScaleAdd(ori, cvScalar(1.0/filenames.size()), avgimg, avgimg);
		cvReleaseImage(&ori);
	}

	cvSaveImage(opt.output, avgimg);//result);
	cvReleaseImage(&avgimg);
	
    return 0;
}
