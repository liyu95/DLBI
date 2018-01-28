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

// -i TIFF -o GAUSSPASS -s 3 -r 60,60,60,60

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
	opt.sigma = 3;
	opt.rio.x = -100;
	
	GetOpts(argc, argv, &opt);
	
	std::vector<std::string> filenames;
	
	GetFilesName(opt.input, filenames);
	
	std::sort(filenames.begin(), filenames.end());
	
	std::cout<<filenames.size()<<std::endl;
	
	if(access(opt.output,0) == -1) {		//create file folder
        mkdir(opt.output,0777);
    }
//     if(access((std::string(opt.output)+"BK").c_str(),0) == -1) {		//create file folder
//         mkdir((std::string(opt.output)+"BK").c_str(),0777);
//     }
	
	std::vector<IplImage*> images;
	std::vector<IplImage*> background;
	
	for(int i = 0; i < filenames.size(); i++){
		std::cout<<(std::string(opt.input)+"/"+filenames[i]).c_str()<<std::endl;
		IplImage* tmp = cvLoadImage((std::string(opt.input)+"/"+filenames[i]).c_str(), CV_LOAD_IMAGE_UNCHANGED);
		if(opt.rio.x >= 0){
			IplImage* tmpp = util::GetSubImage(tmp, opt.rio);
			cvReleaseImage(&tmp);
			tmp = tmpp;
		}
		IplImage* ori = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, 1);
		cvCvtScale(tmp, ori, 1, 0);
		cvReleaseImage(&tmp);
// 		util::ConvertTo1(ori, true);
		IplImage* lawpass = cvCreateImage(cvGetSize(ori), ori->depth, ori->nChannels);
		IplImage* result = cvCreateImage(cvGetSize(ori), ori->depth, ori->nChannels);
		if(opt.sigma > 0){
			cvSmooth(ori, lawpass, CV_GAUSSIAN, 0, 0, opt.sigma);
		}
		else{
			cvFillImage(lawpass, 0);
		}
		background.push_back(lawpass);
		cvScaleAdd(lawpass, cvScalar(-1), ori, result);
		images.push_back(result);
// 		util::ConvertTo1(result, true);

// 		cvReleaseImage(&ori);
// 		cvReleaseImage(&lawpass);
	}
	
	double gmin = DBL_MAX, gmax = DBL_MIN;
	
	for(int i = 0; i < images.size(); i++){
		double minimum, maximum;
// 		CvScalar mean, sdv;
// 		cvAvgSdv(images[i], &mean, &sdv);
		cvMinMaxLoc(images[i], &minimum, &maximum);
		
		if(minimum < gmin){
			gmin = minimum;
		}
		if(maximum > gmax){
			gmax = maximum;
		}
	}
	
	double dif = gmax-gmin;
	
	CvScalar mean, sdv;
	cvAvgSdv(images[0], &mean, &sdv);
	
	for(int i = 0; i < images.size(); i++){
// 		for(int y = 0; y < images[i]->height; y++){
// 			float* ptr = (float*)(images[i]->imageData+y*images[i]->widthStep);
// 			for(int x = 0; x < images[i]->width; x++){
// 				*ptr = (*ptr-gmin)/dif*255;
// 				ptr++;
// 			}
// 		}
		
		std::ostringstream oss;
		oss <<opt.output<<"/"<<filenames[i];
		cvSaveImage(oss.str().c_str(), images[i]);//result);
// 		std::ostringstream ossbk;
// 		ossbk <<opt.output<<"BK/"<<filenames[i];
// 		cvSaveImage(ossbk.str().c_str(), background[i]);
// 		cvReleaseImage(&images[i]);
	}
	
	double gdavg = 0, gdsdv = 0, gdmin = 0, gdmax = 0;
	
	cvAvgSdv(images[0], &mean, &sdv);
	
	std::cout<<"PIC AVG "<<mean.val[0]<<" SDV "<<sdv.val[0]<<std::endl;
	
	for(int i = 1; i < images.size(); i++){
		IplImage* result = cvCreateImage(cvGetSize(images[i]), images[i]->depth, images[i]->nChannels);
		cvScaleAdd(images[i-1], cvScalar(-1), images[i], result);
		
// 		std::ostringstream oss;
// 		oss <<opt.output<<"/diff"<<setfill('0') << setw(3)<<i<<".tiff";;
// 		cvSaveImage(oss.str().c_str(), result);
		
		double minimum, maximum;
		CvScalar mean, sdv;
		cvAvgSdv(result, &mean, &sdv);
		cvMinMaxLoc(result, &minimum, &maximum);
		gdavg += mean.val[0]; gdsdv += sdv.val[0]; gdmin += minimum; gdmax += maximum;
// 		std::cout<<mean.val[0]<<" "<<sdv.val[0]<<" "<<minimum<<" "<<maximum<<std::endl;
	}
	
	gdavg /= images.size()-1; gdsdv /= images.size()-1; gdmin /= images.size()-1; gdmax /= images.size()-1;
	std::cout<<"AVG  SDV  MIN  MAX"<<std::endl;
	std::cout<<gdavg<<" "<<gdsdv<<" "<<gdmin<<" "<<gdmax<<std::endl;
	
    return 0;
}
