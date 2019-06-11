// flowIO.h

#include <opencv2/opencv.hpp>

// the "official" threshold - if the absolute value of either 
// flow component is greater, it's considered unknown
#define UNKNOWN_FLOW_THRESH 1e9

// value to use to represent unknown flow
#define UNKNOWN_FLOW 1e10

#define MAXWHEELCOLS 60

// return whether flow vector is unknown
bool unknown_flow(float u, float v);
bool unknown_flow(float *f);

// read a flow file into 2-band image
int ReadFlowFile(cv::Mat2f& img, const char* filename);

// write a 2-band image into flow file 
int WriteFlowFile(cv::Mat2f& img, const char* filename);


void setcols(int* colorwheel, int r, int g, int b, int k);

int createcolorwheel(int* colorwheel);

void computeColor(double fx, double fy, unsigned char *pix, int* colorwheel, int ncols);

// convert flow to color
double MotionToColor(cv::Mat2f& flow, unsigned char* fillPix, float range /*= -1*/);

// write a 2-band image into a color coded png flow file
void WriteFlowAsImage(cv::Mat2f& flow, const char* imgName, float range /*= -1*/);
