#include <iostream>
#include <climits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include <string.h>
#include "img_util.h"
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "opts.h"


// -i test.pgm -o test -n 200 -s 0.125 -p 8

//1/(2*pi) = 0.159154943091895

double AVG = 0, SDV = 23.8*.65, MIN = -100.8, MAX = 101.36;
double SIG_MAX = 100;
double PIC_AVG = 86;



using namespace std;

#define SCALE0   8
int  MD = 12;	 //128		//(MAX-MIN)/SDV*2	//max density
#define BOUNDARY 200

#define FWHM1     3.96
#define FWHM2     4.16

double GAUSS_MAX = SIG_MAX/MD;

int SCALE = SCALE0;

class Sampling{
private:
	std::vector<int> sample_table;

private:
	float LogNormPSF(float x, float mu, float sigma){ //mu=0.242, sigma=0.1 or mu=-0.218, sigma=0.1
		float k = (log(x)-mu);
		return 1/(x*sigma*sqrt(2*M_PI))*exp(-k*k/(2*sigma*sigma));
	}
	
public:
	int InitializeRadiusSamplingTable(float mu, float sigma){
		std::vector<int> x;
		std::vector<int> y;
		for(float i = FWHM1; i < FWHM2; i+= 0.01){
			y.push_back(int(LogNormPSF(i, mu, sigma)*100+0.5));
			x.push_back(int(i*SCALE+0.5));
		}
		
		int acc = std::accumulate(y.begin(), y.end(), 0);
		sample_table.resize(acc, 0);
		
		int idx = 0;
		for(int i = 0; i < y.size(); i++){
			for(int j = y[i]; j--;){
				sample_table[idx++] = x[i];
			}
		}
	}

public:
	Sampling(){
		srand(time(NULL));
	}

	int GetRadius(){
		int idx = rand()%sample_table.size();
		return sample_table[idx];
	}
};

class Fluophor{
private:
	int val;
	bool state;		//true emitting; false not emitting
	
public:
	Fluophor(){
		val = 1;
		state = false;
	}
	
	void Set(int v, bool s = false){
		val = v;
		state = s;
	}
	
	int Transition(){
		int random;
		
		if(!val){
			return 0;
		}
		if(state){//true emitting
			random = rand()%1000;
			if(random < 210){
				return val*state;	//no change;
			}
			else{
				return val*(state = false);	//stop emitting
			}
		}
		if(!state){//false not emitting
			random = rand()%1000;
			if(random < 305){//490){//
				return val*(state = true);	//emitting
			}
			else if(random < 1000){//990){//980){//
				return val*state;	//no change;
			}
// 			else{					//for xu's data
// 				return val = 0;
// 			}
		}
		return val*state;
	}
	
	int Val(){
		return val*state;
	}
};

float gaussian2D(int x, int y, float sigma)
{
	double sig_2 = 1/(sigma*sigma);
	
	float val = sig_2*0.159154943091895*std::exp(-(x*x+y*y)/(sigma*sigma)*.5);//0.5*sig_2/M_PI*exp(-(x*x+y*y)*sig_2*.5);
	return val;
}

class GaussTable{
private:
	std::vector<IplImage*> pics;
	int start_radius;
	
public:
	GaussTable(){}
	~GaussTable(){
		if(pics.size()){
			for(int i = 0; i < pics.size(); i++){
				cvReleaseImage(&(pics[i]));
			}
		}
	}
	
	void InitializeTable(int ra, int rb){
		pics.resize(rb-ra+1);
		int size = 2*rb+1;			//*1.7 for diffison
		start_radius = ra;
		float frac = 1750.0/MD*.6; //SCALE/MD*8;
		for(int i = 0; i < pics.size(); i++){
			int r = i+ra;
			pics[i] = cvCreateImage(cvSize(size, size), IPL_DEPTH_32F, 1);
			for(int y = -rb; y <= rb; y++){
				for(int x = -rb; x <= rb; x++){
					CV_IMAGE_ELEM(pics[i], float, y+rb, x+rb) = gaussian2D(x, y, r*.65/2.355)*frac;//GAUSS_MAX;//(x/r, y/r, 1);		*.75 is a scale compensate for gaussian filter
				}
			}
		}
	}
	
	IplImage* GetDistribution(int radius){
		int idx = radius-start_radius;
		if(idx >= pics.size() || idx < 0){
			return NULL;
		}
		
		return pics[idx];
	}
};

void GenerateTimeSeries(IplImage* ori, float oscale, float pscale, int number, const char* outdir)
{
	SCALE = pscale;
	MD = rand()%16+8;
	
	double minimum, maximum;
	cvMinMaxLoc(ori, &minimum, &maximum);
	
	for(int y = 0; y < ori->height; y++){
		float* ptr = (float*)(ori->imageData+y*ori->widthStep);
		for(int x = 0; x < ori->width; x++){
			*ptr = (*ptr-minimum)/maximum*255;
			*ptr = int(*ptr/255.0*MD*4);
			ptr++;
		}
	}
	
	std::vector<Fluophor>** flv; //[ori->width][ori->height]
	flv = new std::vector<Fluophor>*[ori->height];
	for(int y = 0; y < ori->height; y++){
		flv[y] = new std::vector<Fluophor>[ori->width];
	}
	
	for(int y = 0; y < ori->height; y++){
		float* ptr = (float*)(ori->imageData+y*ori->widthStep);
		for(int x = 0; x < ori->width; x++){
			int val = int(*(ptr+x));
			if(val < 0){
				std::cout<<"!"<<std::endl;
			}
			if(val){
				flv[y][x].resize(int(val*.25));
// 				std::cout<<int(val/255.0*MD)<<std::endl;
			}
		}
	}
	
	IplImage* canvas = cvCreateImage(cvSize(ori->width+BOUNDARY*2,ori->height+BOUNDARY*2), IPL_DEPTH_32F, 1);
// 	cvFillImage(canvas, PIC_AVG-SDV);
	
	GaussTable table; 	table.InitializeTable(int(FWHM1*SCALE+0.5), int(FWHM2*SCALE+0.5));//18, 38);
	Sampling sampling; sampling.InitializeRadiusSamplingTable(1.4, 0.025);
	
	if(access(outdir,0) == -1) {		//create file folder
        mkdir(outdir,0777);
    }

    EX_TRACE("Simulate %s:", outdir);
    
	int itr = -2;
	float scale = oscale;
	
	cv::RNG rnger(cv::getTickCount());
	float nfrac = rnger.uniform(0.3, 1.);
	float nfrac2 = rnger.uniform(0.3, 1.);
	
	IplImage* bkgd = NULL;
	
	while(++itr < number){
		EX_TRACE("[%d]", itr);
// 		std::cout<<"["<<itr<<"]"<<std::endl;
		memset(canvas->imageData, 0, sizeof(char)*canvas->imageSize);
		for(int y = 0; y < ori->height; y++){
			for(int x = 0; x < ori->width; x++){
				for(int i = 0; i < flv[y][x].size(); i++){
					if(!flv[y][x][i].Transition()){
						continue;
					}
					
					int radius =sampling.GetRadius();
					IplImage* psf = table.GetDistribution(radius);
					CvRect roi; roi.width = psf->width; roi.height = psf->height;
					roi.x = BOUNDARY+x-psf->width/2;
					roi.y = BOUNDARY+y-psf->height/2;
					
					if(roi.x < 0 || roi.y < 0 || roi.x+roi.width >= canvas->width+BOUNDARY*2 || roi.y+roi.height >= canvas->height+BOUNDARY*2){
						continue;
					}
					cvSetImageROI(canvas,roi);
					cvAdd(canvas, psf, canvas);
					cvResetImageROI(canvas);
				}
			}
		}
		
		CvRect roi; roi.x = BOUNDARY; roi.y = BOUNDARY; roi.width = ori->width; roi.height = ori->height;
		cvSetImageROI(canvas,roi);
		
		std::ostringstream oss;
//         oss <<opt.output<<"/"<<opt.output<<setfill('0') << setw(3)<<itr<<".tiff";
		oss <<outdir<<"/"<<setfill('0') << setw(3)<<itr<<".tiff";
		IplImage* cpy = cvCreateImage(cvSize(ori->width*scale, ori->height*scale), canvas->depth, canvas->nChannels);
		cvResize(canvas, cpy, CV_INTER_CUBIC);
		cvResetImageROI(canvas);	
		
		if(bkgd == NULL){
			CvScalar mean, sdv;
			cvAvgSdv(cpy, &mean, &sdv);
			bkgd = cvCreateImage(cvSize(cpy->width, cpy->height), cpy->depth, cpy->nChannels);
			cvFillImage(bkgd, mean.val[0]*rnger.uniform(0.2, 0.7));
		}
		
		for(int i = 0; i < MAX/SDV*.5; i++){
			cv::Mat noise = cv::Mat(cpy->width, cpy->height, CV_32F);
			cv::RNG rnger(cv::getTickCount());
			rnger.fill(noise, cv::RNG::NORMAL, cv::Scalar::all((0.05+nfrac)*SDV), cv::Scalar::all(nfrac2*SDV));
// 			cv::randn(noise, SDV, SDV);//16*.35*.6, 16*0.7*.6);//MD*.35, MD*0.7);
			IplImage ipltemp = noise;
	// 		cvScaleAdd(&ipltemp, cvScalar(MAX/SDV), cpy, cpy);
			cvAdd(&ipltemp, cpy, cpy);
		}
		
		cvAdd(bkgd, cpy, cpy);
		
		
		if(itr >= 0 ){
			cvSaveImage(oss.str().c_str(), cpy);
		}
		
		cvReleaseImage(&cpy);
	}
	cvReleaseImage(&bkgd);
	
	EX_TRACE("\n");
// 	util::ConvertTo1(canvas, true);
// 	util::SaveImage(canvas, "gaus.pgm");
	
	for(int y = 0; y < ori->height; y++){
		delete [] flv[y];
	}
	delete [] flv;	
}

int main(int argc, char **argv){
    opts opt;
	opt.output[0] = '\0';
	opt.number = 1;
	opt.psfwidth = SCALE0;
	
    if(GetOpts(argc, argv, &opt) != 1) {
        return 0;
    }
    
    srand(time(NULL));

	IplImage* tmp = cvLoadImage(opt.input, CV_LOAD_IMAGE_GRAYSCALE);
	IplImage* ori = cvCreateImage(cvSize(tmp->width, tmp->height), IPL_DEPTH_32F, tmp->nChannels);
	double minimum, maximum;
	cvMinMaxLoc(tmp, &minimum, &maximum);
	
	for(int y = 0; y < ori->height; y++){
		float* ptr = (float*)(ori->imageData+y*ori->widthStep);
		uchar* src = (uchar*)(tmp->imageData+y*tmp->widthStep);
		for(int x = 0; x < ori->width; x++){
			*ptr = (*src-minimum)/maximum*255;
			ptr++;
			src++;
		}
	}
	
	cvReleaseImage(&tmp);
	
	if(access(opt.output,0) == -1) {		//create file folder
        mkdir(opt.output,0777);
    }
	
	std::string filename = std::string(basename(opt.input));
	filename = filename.substr(0, filename.rfind('.'));
	IplImage* cpy = cvCreateImage(cvSize(ori->width, ori->height), ori->depth, ori->nChannels);
	cvResize(ori, cpy, CV_INTER_CUBIC);
	GenerateTimeSeries(cpy, opt.sigma, opt.psfwidth, opt.number, (std::string(opt.output)+"/x8").c_str());
	cvSaveImage( (filename+"x8.TIFF").c_str(), cpy);
	cvReleaseImage(&cpy);
	
	cpy = cvCreateImage(cvSize(ori->width*6.0/8, ori->height*6.0/8), ori->depth, ori->nChannels);
	cvResize(ori, cpy, CV_INTER_CUBIC);
	GenerateTimeSeries(cpy, opt.sigma/6*8, opt.psfwidth*6.0/8, opt.number, (std::string(opt.output)+"/x6").c_str());
	cvSaveImage( (filename+"x6.TIFF").c_str(), cpy);
	cvReleaseImage(&cpy);
	
	cpy = cvCreateImage(cvSize(ori->width*4.0/8, ori->height*4.0/8), ori->depth, ori->nChannels);
	cvResize(ori, cpy, CV_INTER_CUBIC);
	GenerateTimeSeries(cpy, opt.sigma/4*8, opt.psfwidth*4.0/8, opt.number, (std::string(opt.output)+"/x4").c_str());
	cvSaveImage( (filename+"x4.TIFF").c_str(), cpy);
	cvReleaseImage(&cpy);
	
	cvReleaseImage(&ori);
	
    return 0;
}
