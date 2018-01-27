#ifndef IMG_UTIL_H__
#define IMG_UTIL_H__

#include "util/exception.h"
#include <vector>
#include <fstream>
#include "highgui.h"
#include "cxcore.h"
#include "cv.h"

namespace util {

inline void SaveImage(const IplImage* img, const char* filename) {
    EX_TRACE("Save Image %s\n", filename)
    IplImage* copy = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
	
	for(int y = 0; y < img->height; y++){
		float* src = (float*)(img->imageData+y*img->widthStep);
		char* start = (char*)(copy->imageData+y*copy->widthStep);
		for(int x = 0; x < img->width; x++){
			*start++ = (*src++)*255;
		}
	}

    cvSaveImage(filename, copy);
    cvReleaseImage(&copy);
}

// typedef struct _IplImage
// {
//     int  nSize;             /* sizeof(IplImage) */
//     int  ID;                /* version (=0)*/
//     int  nChannels;         /* Most of OpenCV functions support 1,2,3 or 4 channels */
//     int  alphaChannel;      /* Ignored by OpenCV */
//     int  depth;             /* Pixel depth in bits: IPL_DEPTH_8U, IPL_DEPTH_8S, IPL_DEPTH_16S,
//                                IPL_DEPTH_32S, IPL_DEPTH_32F and IPL_DEPTH_64F are supported.  */
//     char colorModel[4];     /* Ignored by OpenCV */
//     char channelSeq[4];     /* ditto */
//     int  dataOrder;         /* 0 - interleaved color channels, 1 - separate color channels.
//                                cvCreateImage can only create interleaved images */
//     int  origin;            /* 0 - top-left origin,
//                                1 - bottom-left origin (Windows bitmaps style).  */
//     int  align;             /* Alignment of image rows (4 or 8).
//                                OpenCV ignores it and uses widthStep instead.    */
//     int  width;             /* Image width in pixels.                           */
//     int  height;            /* Image height in pixels.                          */
//     struct _IplROI *roi;    /* Image ROI. If NULL, the whole image is selected. */
//     struct _IplImage *maskROI;      /* Must be NULL. */
//     void  *imageId;                 /* "           " */
//     struct _IplTileInfo *tileInfo;  /* "           " */
//     int  imageSize;         /* Image data size in bytes
//                                (==image->height*image->widthStep
//                                in case of interleaved data)*/
//     char *imageData;        /* Pointer to aligned image data.         */
//     int  widthStep;         /* Size of aligned image row in bytes.    */
//     int  BorderMode[4];     /* Ignored by OpenCV.                     */
//     int  BorderConst[4];    /* Ditto.                                 */
//     char *imageDataOrigin;  /* Pointer to very origin of image data
//                                (not necessarily aligned) -
//                                needed for correct deallocation */
// }
// IplImage;

inline void SeriesSaveToFile(const IplImage* img, const char* filename){
	std::ofstream ofs;
	ofs.open(filename, std::ofstream::binary);
	ofs.write(reinterpret_cast<const char*>(&(img->width)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->height)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->depth)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->nChannels)), sizeof(int));
	
	ofs.write(reinterpret_cast<const char*>(&(img->nSize)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->ID)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->alphaChannel)), sizeof(int));
	ofs.write(img->colorModel, sizeof(char)*4);
	ofs.write(img->channelSeq, sizeof(char)*4);
	ofs.write(reinterpret_cast<const char*>(&(img->dataOrder)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->origin)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->align)), sizeof(int));
	//roi   	NULL
	//maskROI   NULL
	//imageId	NULL
	//tileInfo	NULL
	ofs.write(reinterpret_cast<const char*>(&(img->imageSize)), sizeof(int));
	ofs.write(img->imageData, sizeof(char)*img->imageSize);
	ofs.write(reinterpret_cast<const char*>(&(img->widthStep)), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&(img->BorderMode)), sizeof(int)*4);
	ofs.write(reinterpret_cast<const char*>(&(img->BorderConst)), sizeof(int)*4);
	//imageDataOrigin   point to the same place of imageData
	ofs.close();
}

inline void SeriesReadFromFile(IplImage** img, const char* filename){
	std::ifstream ifs;
	ifs.open(filename, std::ofstream::binary);
	int tmp[4];//width, height, depth, nChannels;
	ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int)*4);
	
	*img = cvCreateImage(cvSize(tmp[0], tmp[1]), tmp[2], tmp[3]);
	ifs.read(reinterpret_cast<char*>(&((*img)->nSize)), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&((*img)->ID)), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&((*img)->alphaChannel)), sizeof(int));
	ifs.read((*img)->colorModel, sizeof(char)*4);
	ifs.read((*img)->channelSeq, sizeof(char)*4);
	ifs.read(reinterpret_cast<char*>(&((*img)->dataOrder)), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&((*img)->origin)), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&((*img)->align)), sizeof(int));
	//roi   	NULL
	//maskROI   NULL
	//imageId	NULL
	//tileInfo	NULL
	ifs.read(reinterpret_cast<char*>(&((*img)->imageSize)), sizeof(int));
	ifs.read((*img)->imageData, sizeof(char)*(*img)->imageSize);
	ifs.read(reinterpret_cast<char*>(&((*img)->widthStep)), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&((*img)->BorderMode)), sizeof(int)*4);
	ifs.read(reinterpret_cast<char*>(&((*img)->BorderConst)), sizeof(int)*4);
}

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

inline void ScaleImage(IplImage*& image, float scale)
{
	if(scale == 1){
		return;
	}
	
	IplImage* cpy = cvCreateImage(cvSize(image->width*scale, image->height*scale), image->depth, image->nChannels);
	cvResize(image, cpy, CV_INTER_CUBIC);
	cvReleaseImage(&image);
	image = cpy;
}

inline void Reversal(IplImage* img)
{
    double minimum, maximum;
    cvMinMaxLoc(img, &minimum, &maximum);
    double p = maximum + minimum;
    for(int y = 0; y < img->height; y++){
		float* ptr = (float*)(img->imageData+y*img->widthStep);
		for(int x = 0; x < img->width; x++){
			*ptr =p-*ptr;
			ptr++;
		}
	}
}

// 	cvErode( img, img, NULL, 5);
//      cvDilate(img, img, NULL, 5);

inline void MedianSmooth(IplImage *img)
{
    cvSmooth(img, img, CV_MEDIAN, 5);
}

#define HISTOGRAM_BIN			1024//	256
inline void HistogramStretch(IplImage *img, int binsize = 1024)
{
    cv::Mat mat(img, 0);
    double minimum, maximum;
    cvMinMaxLoc(img, &minimum, &maximum);
    minimum = minimum<0?minimum:0;
    maximum = maximum>1?maximum:1;
    float range[] = {minimum, maximum};
    const float* ranges = {range};
    cv::Mat hist;
    const int channel = 0;
	
    cv::calcHist(&mat, 1, &channel, cv::Mat(), hist, 1, &binsize, &ranges);	/*warning: if has some problem, look here */
    cv::normalize(hist, hist, 1, 0, CV_L1);
    float* vlu = (float*)hist.data;

    for(int i = 1; i < binsize; i++) {
        vlu[i] += vlu[i-1];
    }

    float delta = (maximum-minimum)/(binsize);

    float res = delta/2/(maximum-minimum);
    for(int x = 0; x < img->height; x++) {
        for(int y = 0; y < img->width; y++) {
            int loc = (int)((CV_IMAGE_ELEM(img, float, x, y)-minimum)/delta);
            CV_IMAGE_ELEM(img, float, x, y) = (maximum-minimum)*(vlu[loc]+delta)+minimum;
        }
    }
}
#undef HISTOGRAM_BIN


inline void ConvertTo1(IplImage* img, bool cutoff = false)
{
    double minimum, maximum;
    CvScalar mean, sdv;
    cvAvgSdv(img, &mean, &sdv);
    cvMinMaxLoc(img, &minimum, &maximum);

#define CUTOFF (3)

    maximum = (maximum > mean.val[0]+CUTOFF*sdv.val[0])? mean.val[0]+CUTOFF*sdv.val[0] : maximum;
    minimum = (minimum > mean.val[0]-CUTOFF*sdv.val[0])? minimum : mean.val[0]-CUTOFF*sdv.val[0];
	
	float span;
	(maximum-minimum)==0?span = 1:(span = maximum-minimum);

	for(int y = 0; y < img->height; y++){
		float* ptr = (float*)(img->imageData+y*img->widthStep);
		for(int x = 0; x < img->width; x++){
			if(cutoff){
				if(fabs(*ptr-mean.val[0]) <= CUTOFF*sdv.val[0]){
					*ptr = (*ptr-minimum)/span;
				}
				else {
					if(*ptr > mean.val[0]) {
						*ptr = 1;
					}
					else {
						*ptr = 0;
					}
				}
			}
			else{
				*ptr = (*ptr-minimum)/span;
			}
			ptr++;
		}
	}
	
#undef CUTOFF
}

inline void CopyToUChar256(const IplImage* img, IplImage** cpy, int w, int b, bool cutoff = false)
{
    *cpy = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
    double minimum, maximum;
    CvScalar mean, sdv;
    cvAvgSdv(img, &mean, &sdv);
    cvMinMaxLoc(img, &minimum, &maximum);

#define CUTOFF (3)//(4)//(5)

    maximum = (maximum > mean.val[0]+CUTOFF*sdv.val[0])? mean.val[0]+CUTOFF*sdv.val[0] : maximum;
    minimum = (minimum > mean.val[0]-CUTOFF*sdv.val[0])? minimum : mean.val[0]-CUTOFF*sdv.val[0];
    
    for(int y = 0; y < img->height; y++){
		float* ptr = (float*)(img->imageData+y*img->widthStep);
		uchar* dst = (uchar*)((*cpy)->imageData+y*(*cpy)->widthStep);
		for(int x = 0; x < img->width; x++){
			if(cutoff){
				if(fabs(*ptr-mean.val[0]) <= CUTOFF*sdv.val[0]) {
					*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
					*dst = *dst > w ? w : *dst;
				}
				else {
					if(*ptr > mean.val[0]) {
						*dst = w;
					}
					else {
						*dst = b;
					}
				}
			}
			else{
				*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
				*dst = *dst > w ? w : *dst;
			}
			ptr++;
			dst++;
		}
	}

#undef CUTOFF
}

inline void LaplaceSharpen(IplImage* img)
{
    CvMat* kernel;
    kernel = cvCreateMat(3,3,CV_32F);
    cvmSet(kernel, 0, 0, 0);
    cvmSet(kernel, 0, 1, -1);
    cvmSet(kernel, 0, 2, 0);
    cvmSet(kernel, 1, 0, -1);
    cvmSet(kernel, 1, 1, 5);
    cvmSet(kernel, 1, 2, -1);
    cvmSet(kernel, 2, 0, 0);
    cvmSet(kernel, 2, 1, -1);
    cvmSet(kernel, 2, 2, 0);

    cvFilter2D(img, img, kernel);
}

}

#endif
