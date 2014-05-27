#include "stdafx.h"
#include "nvVector.h"
#include <string>
#include <iostream>
#include <fstream>
#include <io.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
using namespace std;

#define CLAMP(v, min_v, max_v) v < min_v ? min_v : v > max_v ? max_v : v 

using namespace nv;
IplImage* convert_to_float32(IplImage* img);
IplImage *increaseImageBrightness(IplImage *img, int delta);
IplImage *decreaseImageBrightness(IplImage *img, int delta);
IplImage *IncreaseImageContrast(IplImage *img, float scopeA);
IplImage *DecreaseImageContrast(IplImage *img, float scopeA);
IplImage *GammaCorrection(IplImage *img, float gamma);
vec3f clampVector(const vec3f &v);
vec3f convertLuminanceToRGB(double luminance);
IplImage *readImagePFM(const string& fileName);
void saveImagePFM(const string& fileName, const IplImage* image);

void Compare(IplImage *img1, IplImage *img2, IplImage *ref){
	int width = ref->width, height = ref->height;
	std::vector<double> var1(width*height), var2(width*height);
	std::cout << "run1 " << std::endl;
	IplImage *rst1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	IplImage *rst2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	std::cout << "run2 " << std::endl;

	double TotalVar1 = 0, TotalVar2 = 0;
	int count = 0;
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f bgr_ref = ((vec3f*)ref->imageData)[y*width + x];
			vec3f bgr_img1 = ((vec3f*)img1->imageData)[y*width + x];
			vec3f bgr_img2 = ((vec3f*)img2->imageData)[y*width + x];
			
			var1[y*width+x] = pow((bgr_ref - bgr_img1).length(),2); 
			var2[y*width+x] = pow((bgr_ref - bgr_img2).length(),2); 
				
			int index = width * y + x;
			
			TotalVar1 += pow((bgr_img1-bgr_ref).length(),2);
			TotalVar2 += pow((bgr_img2-bgr_ref).length(),2);
		}
	}

	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f &var1C = ((vec3f*)rst1->imageData)[y*width + x];
			vec3f &var2C = ((vec3f*)rst2->imageData)[y*width + x];

			vec3f RGB = convertLuminanceToRGB((double)var1[y*width+x]*5.f);
			var1C = vec3f(255*RGB.z, 255*RGB.y, 255*RGB.x);
			RGB = convertLuminanceToRGB((double)var2[y*width+x]*5.f);
			var2C = vec3f(255*RGB.z, 255*RGB.y, 255*RGB.x);
			count++;
			
		}
	}
	
	std::cout << "run3 " << std::endl;

	cvSaveImage("rst1.jpg" , rst1);
	std::cout << "run4 " << std::endl;
	cvSaveImage("rst2.jpg" , rst2);

	double Var1 = sqrt(TotalVar1 / (count)), Var2 = sqrt(TotalVar2 / (count));
	std::cout << "Var1 = " << Var1 << " Var2 = " << Var2 << std::endl;
}

/*
int main(int argc, char* argv[])
{
	std::string file1, file2, fileRef;
	//std::cin >> file1 >> file2 >> fileRef;
	
	IplImage *img1 = convert_to_float32(readImagePFM(argv[1]));
	IplImage *img2 = convert_to_float32(readImagePFM(argv[2]));
	IplImage *imgRef = convert_to_float32(readImagePFM(argv[3]));
	
	Compare(img1, img2, imgRef);

	return 0;
}
*/
/*
int main()
{
	_finddata_t file;
	long flag;
	flag = _findfirst("D:\\Winmad\\RendererGPU\\Release\\Data\\results\\sss_vcm_5_26_gpu\\*.pfm" , &file);
	for (;;)
	{
		printf("%s\n", file.name);
		string root = "D:\\Winmad\\RendererGPU\\Release\\Data\\results\\sss_vcm_5_26_gpu\\";
		IplImage *img = convert_to_float32(readImagePFM(root + file.name));
		saveImagePFM(root + file.name , img);
		cvReleaseImage(&img);
		if (_findnext(flag , &file) == -1)
			break;
	}
	_findclose(flag);
}
*/
IplImage *readImagePFM(const string& fileName)
{
	FILE* file;
	int height, width;
	fopen_s(&file, fileName.c_str(), "rb");

	fscanf(file, "PF\n%d %d\n-1.000000\n" , &width , &height);

	IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	float* data = (float*)img->imageData;
	for (int j=img->height-1; j>=0; j--)
	{
		for (unsigned i=0; i<img->width; i++)
		{
			fread((float*)&data[3*(j*img->width+i)+2], sizeof(float), 1, file);
			fread((float*)&data[3*(j*img->width+i)+1], sizeof(float), 1, file);
			fread((float*)&data[3*(j*img->width+i)], sizeof(float), 1, file);
			if (data[3*(j*img->width+i)+2] > 50.f)
			{
				data[3*(j*img->width+i)+2] = 50.f;
				data[3*(j*img->width+i)+1] = 50.f;
				data[3*(j*img->width+i)] = 50.f;
			}
		}
	}

	fclose(file);
	return img;
}

void saveImagePFM(const string& fileName, const IplImage* image)
{
	FILE* file;
	fopen_s(&file, fileName.c_str(), "wb");

	fprintf(file, "PF\n%d %d\n-1.000000\n", image->width, image->height);

	const float* data = (float*)image->imageData;
	for(int j=image->height-1; j>=0; j--)
	{
		for(unsigned i=0; i<image->width; i++)
		{
			fwrite(&data[3*(j*image->width+i)+2], sizeof(float), 1, file);
			fwrite(&data[3*(j*image->width+i)+1], sizeof(float), 1, file);
			fwrite(&data[3*(j*image->width+i)], sizeof(float), 1, file);
		}
	}
	fclose(file);
}

vec3f convertLuminanceToRGB(double luminance){
	if(luminance < 0.25){
		float g = luminance / 0.25;
		return vec3f(0, g, 1);
	}
	else if(luminance < 0.5){
		float b = 1 - (luminance - 0.25) / 0.25;
		return vec3f(0, 1, b);
	}
	else if(luminance < 0.75){
		float r = (luminance - 0.5) / 0.25;
		return vec3f(r, 1, 0);
	}
	else{
		float g = 1 - (luminance - 0.75) / 0.25;
		return vec3f(1, g, 0);
	}

}


IplImage *increaseImageBrightness(IplImage *img, int delta){
	int width = img->width, height = img->height;
	IplImage *img_new = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f bgr = ((vec3f*)img->imageData)[y*width + x];
			vec3f &bgr_new = ((vec3f*)img_new->imageData)[y*width + x];
			bgr_new = bgr + vec3f(delta);
			bgr_new = clampVector(bgr_new);
		}
	}
	return img_new;
}

IplImage *decreaseImageBrightness(IplImage *img, int delta){
	int width = img->width, height = img->height;
	IplImage *img_new = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f bgr = ((vec3f*)img->imageData)[y*width + x];
			vec3f &bgr_new = ((vec3f*)img_new->imageData)[y*width + x];
			bgr_new = bgr - vec3f(delta);
			bgr_new = clampVector(bgr_new);
		}
	}
	return img_new;
}

IplImage *IncreaseImageContrast(IplImage *img, float scopeA){
	if(scopeA < 1.f){
		std::cerr << "assert a > 1.0 " << std::endl;
		return img;
	}
	int width = img->width, height = img->height;
	IplImage *img_new = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f bgr = ((vec3f*)img->imageData)[y*width + x];
			vec3f &bgr_new = ((vec3f*)img_new->imageData)[y*width + x];
			bgr_new = scopeA * (bgr - vec3f(127)) + vec3f(127);
			bgr_new = clampVector(bgr_new);
		}
	}
	return img_new;
}

IplImage *DecreaseImageContrast(IplImage *img, float scopeA){
	if(scopeA > 1.f){
		std::cerr << "assert a < 1.0 " << std::endl;
		return img;
	}
	int width = img->width, height = img->height;
	IplImage *img_new = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f bgr = ((vec3f*)img->imageData)[y*width + x];
			vec3f &bgr_new = ((vec3f*)img_new->imageData)[y*width + x];
			bgr_new = scopeA * (bgr - vec3f(127)) + vec3f(127);
			bgr_new = clampVector(bgr_new);
		}
	}
	return img_new;
}

IplImage *GammaCorrection(IplImage *img, float gamma){
	int width = img->width, height = img->height;
	IplImage *img_new = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			vec3f bgr = ((vec3f*)img->imageData)[y*width + x];
			vec3f &bgr_new = ((vec3f*)img_new->imageData)[y*width + x];
			for(int i = 0; i < 3; i++){
				float invGamma = 1.0 / gamma;
				bgr_new[i] = 255.0 * std::powf(bgr[i]/255.0, invGamma);
			}
			bgr_new = clampVector(bgr_new);
		}
	}
	return img_new;
}

IplImage* convert_to_float32(IplImage* img)
{
    IplImage* img32f = cvCreateImage(cvGetSize(img),IPL_DEPTH_32F,img->nChannels);

    for(int i=0; i<img->height; i++)
    {
        for(int j=0; j<img->width; j++)
        {
            cvSet2D(img32f,i,j,cvGet2D(img,i,j));
        }
    }
    return img32f;
}

vec3f clampVector(const vec3f &v){
	vec3f ans;
	for(int i = 0; i < 3; i++)
		ans[i] = CLAMP(v[i], 0, 255);
	return ans;
}