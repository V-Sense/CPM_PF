// flowIO.h
#pragma once
#ifndef FLOW_H
#define FLOW_H

#include <opencv2/opencv.hpp>

// the "official" threshold - if the absolute value of either 
// flow component is greater, it's considered unknown
#define UNKNOWN_FLOW_THRESH 1e9

// value to use to represent unknown flow
#define UNKNOWN_FLOW 1e10

#define MAXWHEELCOLS 60

const cv::Vec2i kPOSITION_INVALID = cv::Vec2i(-1, -1);
const float kMOVEMENT_UNKNOWN = 1e10;
const cv::Vec2f kFLOW_UNKNOWN = cv::Vec2f(kMOVEMENT_UNKNOWN,kMOVEMENT_UNKNOWN);

// return whether flow vector is unknown
bool unknown_flow(float u, float v);
bool unknown_flow(float *f);

/* ------------------------ OPTICAL FLOW INPUT / OUTPUT ----------------------- */
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


/* ------------------------ OPERATIONS ON OPTICAL FLOWS  ----------------------- */
// Transforms a relative flow R into an absolute flow A, checking for margins
// Ax = X + Rx , Ay = Y + Ry
// if either
//   - Ax is outside [0,w] or
//   - Ay outside [0,y]
// it returns kPOSITION_INVALID
inline cv::Vec2i getAbsoluteFlow(int x, int y, const cv::Vec2f& flow, int h, int w)
{
    cv::Vec2i result(cvRound(y + flow[1]), cvRound(x + flow[0]));
    if(result[0] >= 0 && result[0] < h && result[1] >= 0 && result[1] < w)
        return result;
    else
        return kPOSITION_INVALID;
}

inline int getAbsoluteDisp(int x, float disp, int w)
{
    int result = cvRound((float)x + disp);
    if(result >= 0 && result < w)
        return result;
    else
        return (int)kMOVEMENT_UNKNOWN;
}

cv::Mat1f getFlowConfidence(cv::Mat2f forward_flow, cv::Mat2f backward_flow);
cv::Mat1f getHorDispConfidence(cv::Mat1f forward_disp, cv::Mat1f backward_flow);
cv::Mat1f getVerDispConfidence(cv::Mat1f forward_disp, cv::Mat1f backward_flow);

#endif //!FLOW_H