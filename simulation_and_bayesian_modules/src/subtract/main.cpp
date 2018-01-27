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
// #include "img_util.h"
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "opts.h"

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

int main(int argc, char **argv)
{
	struct options opts;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}
	
	std::vector<std::string> filenames;
	GetFilesName(opts.input.c_str(), filenames);
	std::sort(filenames.begin(), filenames.end());
	
	std::vector<std::string> peers;
	GetFilesName(opts.peer.c_str(), peers);
	std::sort(peers.begin(), peers.end());
	
	std::vector<IplImage*> images;
	
	for(int i = 0; i < filenames.size(); i++){
		std::cout<<(std::string(opts.input)+"/"+filenames[i]).c_str()<<std::endl;
		IplImage* ori = cvLoadImage((std::string(opts.input)+"/"+filenames[i]).c_str(), CV_LOAD_IMAGE_UNCHANGED);

		IplImage* lawpass = cvLoadImage((std::string(opts.peer)+"/"+peers[i]).c_str(), CV_LOAD_IMAGE_UNCHANGED);
		IplImage* result = cvCreateImage(cvGetSize(ori), ori->depth, ori->nChannels);

		cvScaleAdd(lawpass, cvScalar(-1), ori, result);
		images.push_back(result);
// 		util::ConvertTo1(result, true);

		cvReleaseImage(&ori);
		cvReleaseImage(&lawpass);
	}
	
	if(access(opts.output.c_str(),0) == -1) {		//create file folder
        mkdir(opts.output.c_str(), 0777);
    }
	
	for(int i = 0; i < images.size(); i++){
// 		for(int y = 0; y < images[i]->height; y++){
// 			float* ptr = (float*)(images[i]->imageData+y*images[i]->widthStep);
// 			for(int x = 0; x < images[i]->width; x++){
// 				*ptr = (*ptr-gmin)/dif*255;
// 				ptr++;
// 			}
// 		}
		
		std::ostringstream oss;
		oss <<opts.output<<"/"<<filenames[i];
		cvSaveImage(oss.str().c_str(), images[i]);
	}
}

